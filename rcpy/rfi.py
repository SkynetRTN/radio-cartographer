"""Stage 9 — 2D RFI subtraction (Paper §3.6).

Faithful port of the C++ ProcessorRFI algorithm.

For every science sample (ra, dec) we build a local model:

  1. Gather all neighbour samples within `rfi_scale_deg` (a fraction
     of one beamwidth) — these are the "candidates" for the local
     model.
  2. Fit the cos² shape function from Paper Eq. 9 to those candidates
     analytically (Paper Eqs. B3/B4):
            flux_k = f * cos²(π·dist_k / (2·rfi_scale)) + g
     where f is the source amplitude at the anchor and g is a small
     baseline-level offset constrained by a noise-level prior.
  3. Iteratively reject the largest-positive deviation from this fit
     (one-sided rejection, biased toward removing RFI) until the RMS
     of positive residuals drops below the 2D noise level. Refit at
     each step.
  4. For every surviving candidate, emit a (LMValue, LMWeight) pair
     using paper Eq. C4 for the weight (position-dependent confidence
     of the local model).

Each science sample then ends up with a collection of LMValue/LMWeight
contributions from every anchor whose neighbourhood overlaps it.

For each sample, we run RCR `LS_MODE_DL` bulk-rejection on the
accumulated values weighted by LMWeights and take the surviving mean.
That's the RFI-subtracted flux. Samples with no surviving local model
are excised (set NaN) and the surface-model stage interpolates over the
gap.

This module skips two C++ features for now (they're additive and not
critical for Jupiter/Cyg-A):

  - Auto-centroiding around bright sources adds extra anchor centers
    on a fine grid near peaks > 75σ. Useful for high-S/N source recovery
    but doesn't change the algorithm's correctness for typical sources.
  - The correlation-map output (theta_corr, for photometry error bars)
    is set to a passable stub. Implement when photometry is wired up.

References:
    Paper §3.6, §3.6.5, Eqs. 9, B1-B4, C4.
    C++ source: src/ProcessorRFI.cpp.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.spatial import cKDTree

from rcpy.types import Survey, Composite
from rcpy.rcr2 import RCR, RejectionTech


# -----------------------------------------------------------------------------
# Analytic cos² local-model fit (Paper Eqs. B3/B4)
# -----------------------------------------------------------------------------

def _rfi_fit(
    keep: np.ndarray,
    flux: np.ndarray,
    shape: np.ndarray,
    dumps: np.ndarray,
    sigma: np.ndarray,
) -> tuple[float, float]:
    """Analytic weighted fit of  flux = f·shape + g  with a noise-level
    prior on g (pulled toward 0 with variance σ²).

    All inputs are 1D arrays of length N (candidates). `keep` selects
    which candidates are currently active. `shape` is the cos²
    pre-factor evaluated at each candidate's distance from the anchor.
    `sigma` is each candidate's 2D-noise-level estimate.

    Returns (f, g). Paper Eqs. B3/B4 (with the prior baked into the
    factor of 2 in the denominators).
    """
    sel = keep
    if not sel.any():
        return 0.0, 0.0

    f_ = flux[sel]
    s_ = shape[sel]
    d_ = dumps[sel]
    sig2 = sigma[sel] ** 2
    sig2 = np.where(sig2 > 0, sig2, 1.0)

    w = d_ / sig2
    sumW = float(w.sum())
    sumWF = float((w * f_).sum())
    sumWS = float((w * s_).sum())
    sumWFS = float((w * f_ * s_).sum())
    sumWSS = float((w * s_ * s_).sum())

    den_f = 2.0 * sumW * sumWSS - sumWS * sumWS
    if abs(den_f) < 1e-30:
        return 0.0, 0.0
    f = (2.0 * sumW * sumWFS - sumWF * sumWS) / den_f
    g = (sumWF - f * sumWS) / (2.0 * sumW)
    if not (math.isfinite(f) and math.isfinite(g)):
        return 0.0, 0.0
    return f, g


def _reject_points(
    avg_sigma: float,
    keep: np.ndarray,
    flux: np.ndarray,
    shape: np.ndarray,
    dumps: np.ndarray,
    sigma: np.ndarray,
) -> None:
    """One-sided iterative outlier rejection. For f > 0 (real source-
    like fits) we reject ABOVE-model outliers (RFI is always positive);
    for f < 0 we reject the largest |outlier| regardless of sign.

    Stops when the RMS of positive residuals (splus) drops below
    avg_sigma. Direct port of C++ ProcessorRFI::rfiRemovePoints.
    """
    last_checked: int | None = None
    prev_splus_diff = float("inf")
    splus = float("inf")

    while splus > avg_sigma:
        if int(keep.sum()) <= 1:
            return

        f, g = _rfi_fit(keep, flux, shape, dumps, sigma)
        model = f * shape + g
        delta = flux - model
        idx_kept = np.where(keep)[0]
        n_kept = idx_kept.size
        delta_k = delta[idx_kept]

        if f > 0:
            pos_mask = delta_k > 0
            n_pos = int(pos_mask.sum())
            if n_pos == 0:
                return
            splus_sq = float((delta_k[pos_mask] ** 2).sum())
            splus = math.sqrt(splus_sq / n_pos)
            if not pos_mask.any():
                return
            # Largest positive residual
            arg = int(np.argmax(delta_k))
            if delta_k[arg] <= 0:
                return
            large_local = int(idx_kept[arg])
        else:
            # f <= 0 — reject the largest absolute residual
            arg = int(np.argmax(np.abs(delta_k)))
            large_local = int(idx_kept[arg])
            n_pos = int(n_kept)
            splus_sq = float((delta_k ** 2).sum())
            splus = math.sqrt(splus_sq / n_pos)

        if splus > avg_sigma:
            keep[large_local] = False
            last_checked = large_local
            prev_splus_diff = abs(avg_sigma - splus)
        else:
            # Match C++ ProcessorRFI::rfiRemovePoints (ProcessorRFI.cpp:498-501)
            # behavior exactly. C++ does:
            #     if (splusSigmaDiff < |scatter - splus|) {
            #         checks[lastChecked] = false;  // no-op: already false
            #     }
            # This is a C++ no-op bug — `checks[lastChecked]` was already
            # set to false in the previous iteration. An earlier port of
            # this routine "fixed" the bug by RESTORING the last rejection
            # (setting to True), which made our pipeline keep MORE samples
            # than C++ kept around bright sources. For parity we match
            # C++'s observable behavior (the rejection stays in place),
            # not the "probably intended" semantics.
            return


# -----------------------------------------------------------------------------
# Per-anchor local-model assembly
# -----------------------------------------------------------------------------

def _anchor_local_model(
    anchor_xy: np.ndarray,
    anchor_flux: float,
    anchor_sigma: float,
    candidates_idx: np.ndarray,
    xy_all: np.ndarray,
    flux_all: np.ndarray,
    dumps_all: np.ndarray,
    sigma_all: np.ndarray,
    rfi_scale_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build one anchor's local model.

    Returns six parallel arrays of length = number of survivors:
       lm_idx[k]:     global index of the k-th survivor
       lm_val[k]:     local-model value at the survivor's position
       lm_w[k]:       weight per Paper Eq. C4
       lm_w_dist_sq[k]: per-survivor sum of dump * distance² for "in
                        correlation" neighbours (Paper Appendix D)
       lm_n_sum[k]:   per-survivor sum of dump for those neighbours
       lm_rfi_count[k]: cohort-member count within rfi_scale of survivor
                        (used to compute GMWeight = LMWeight/rfiCount)

    "In correlation" means a neighbour whose local-model flux falls on
    the SAME side of √2·sigma as the ANCHOR's local-model flux.

    Direct port of the C++ ProcessorRFI::rfiCollectResults (with the
    correlatedWeightMap branch enabled, which is the configuration
    used to generate the reference correlation map).
    """
    if candidates_idx.size < 2:
        empty = np.empty(0, dtype=np.float64)
        return np.empty(0, dtype=int), empty, empty, empty, empty, empty

    # Per-candidate geometry
    dx = xy_all[candidates_idx, 0] - anchor_xy[0]
    dy = xy_all[candidates_idx, 1] - anchor_xy[1]
    dist = np.hypot(dx, dy)
    shape = np.cos(math.pi * dist / (2.0 * rfi_scale_deg)) ** 2

    flux = flux_all[candidates_idx]
    dumps = dumps_all[candidates_idx]
    sigma = sigma_all[candidates_idx]

    keep = np.ones(candidates_idx.size, dtype=bool)

    avg_sigma = float(np.mean(sigma))
    if not math.isfinite(avg_sigma) or avg_sigma <= 0:
        avg_sigma = 1e-6

    _reject_points(avg_sigma, keep, flux, shape, dumps, sigma)

    if not keep.any():
        empty = np.empty(0, dtype=np.float64)
        return np.empty(0, dtype=int), empty, empty, empty, empty, empty

    f, g = _rfi_fit(keep, flux, shape, dumps, sigma)
    lm_val_full = f * shape + g  # at each candidate's position

    # Eq. C4 weight: w = sumD / (1 + 2*(pw - pwAvg)² / <pwDelta²>)
    # C++ ProcessorRFI::rfiCollectResults (lines 699-702) computes
    # <pwDelta²> as the mean squared deviation of pw across ALL
    # candidates (kept + rejected), using the KEPT-samples mean as the
    # reference. This matters: rejected samples near the anchor have
    # high pw and inflate the dispersion, which keeps LMWeights ~ sumD
    # (close to 1 in the denominator) for surviving samples. Earlier I
    # used np.var(pw[keep]) which only saw the off-source survivors'
    # narrow pw range, shrinking the dispersion and over-attenuating
    # LMWeights for any high-pw survivor.
    pw = shape  # the proximity weight = shape factor itself
    pw_avg = float(np.mean(pw[keep]))
    pw_delta_sq = float(np.mean((pw - pw_avg) ** 2))  # ALL candidates
    sumD = float(dumps[keep].sum())

    n_surviving = int(keep.sum())
    n_total_cand = int(candidates_idx.size)
    if n_total_cand < 3 or pw_delta_sq <= 0:
        weights_full = dumps.astype(np.float64)
    else:
        weights_full = sumD / (1.0 + 2.0 * (pw - pw_avg) ** 2 / pw_delta_sq)

    # === Appendix D: theta_corr accumulator ===
    # `centerLMFlux` is the LM value at the anchor itself (shape = 1).
    center_lm = f * 1.0 + g
    threshold = math.sqrt(2.0) * anchor_sigma
    count_above = center_lm > threshold

    # For each pair (survivor i, survivor j), accumulate distance²·dump
    # if j is within rfi_scale of i AND j's LM flux is on the same side
    # of √2·sigma as the anchor's LM flux.
    surv_idx = np.where(keep)[0]
    keep_xy = xy_all[candidates_idx[surv_idx]]  # positions
    keep_dumps = dumps[surv_idx]
    keep_lm = lm_val_full[surv_idx]
    flux_above = keep_lm > threshold

    # Only neighbours on the same side as the anchor are "in correlation".
    if count_above:
        cohort = flux_above
    else:
        cohort = ~flux_above

    cohort_xy = keep_xy[cohort]
    cohort_dumps = keep_dumps[cohort]

    lm_w_dist_sq = np.zeros(n_surviving, dtype=np.float64)
    lm_n_sum = np.zeros(n_surviving, dtype=np.float64)
    # rfiCount per survivor: count of cohort members (including self if
    # in cohort) within rfi_scale of survivor i. Used by C++ for
    # GMWeight = LMWeight / rfiCount and for the intersect calculation
    # in weight2. C++ ProcessorRFI.cpp:714: `rfiCountHold += 1` for
    # cohort_match within rfi_scale_deg, INCLUDING when pointDistance==0
    # (i.e., j == i). Below it then forces rfiCountHold=1 if zero.
    lm_rfi_count = np.zeros(n_surviving, dtype=np.float64)

    if cohort_xy.shape[0] > 0:
        # For each survivor i, compute its distance to every cohort point j.
        # Pairs are kept where 0 < distance < rfi_scale_deg.
        for k, i_global in enumerate(surv_idx):
            xi, yi = keep_xy[k]
            d_ij = np.hypot(cohort_xy[:, 0] - xi, cohort_xy[:, 1] - yi)
            # rfiCount: cohort members within rfi_scale (includes self at d=0)
            in_scale = d_ij < rfi_scale_deg
            lm_rfi_count[k] = float(in_scale.sum())
            valid = (d_ij > 0) & (d_ij < rfi_scale_deg)
            if not valid.any():
                continue
            lm_n_sum[k] = float(cohort_dumps[valid].sum())
            lm_w_dist_sq[k] = float((cohort_dumps[valid] * d_ij[valid] ** 2).sum())
    # C++ guard: if rfiCountHold ended up 0 (e.g. survivor not in cohort
    # at all), force it to 1 so the division doesn't blow up.
    lm_rfi_count[lm_rfi_count == 0] = 1.0

    return (
        candidates_idx[surv_idx],
        lm_val_full[surv_idx].astype(np.float64),
        weights_full[surv_idx].astype(np.float64),
        lm_w_dist_sq,
        lm_n_sum,
        lm_rfi_count,
    )


# -----------------------------------------------------------------------------
# Per-sample global-model assembly via RCR
# -----------------------------------------------------------------------------

def _rcr_mean_with_flags(vals: np.ndarray, weights: np.ndarray) -> tuple[float, np.ndarray]:
    """Run RCR LS_MODE_DL weighted bulk rejection and return BOTH the
    surviving mean AND the boolean flag array indicating which
    contributions survived. The Appendix D theta_corr accumulator
    needs the flags so it can sum only the survivors.

    Falls back to (weighted mean, all-True flags) on degenerate input.
    """
    if vals.size == 0:
        return float("nan"), np.zeros(0, dtype=bool)
    if vals.size == 1:
        return float(vals[0]), np.ones(1, dtype=bool)
    spread = float(np.max(vals) - np.min(vals))
    if spread <= 0:
        return float(vals[0]), np.ones(vals.size, dtype=bool)
    rng = np.random.default_rng(0)
    jitter = rng.normal(0, max(spread * 1e-12, 1e-15), vals.size)
    try:
        rcr = RCR(RejectionTech.LS_MODE_DL)
        rcr.perform_bulk_rejection(vals + jitter, w=weights)
        mu = float(rcr.result.mu)
        flags = rcr.result.flags
        if flags.size != vals.size or not math.isfinite(mu):
            return float(np.average(vals, weights=weights)), np.ones(vals.size, dtype=bool)
        return mu, flags.astype(bool)
    except (ZeroDivisionError, ValueError, IndexError, RuntimeError,
            OverflowError, ArithmeticError):
        return float(np.average(vals, weights=weights)), np.ones(vals.size, dtype=bool)


def _rcr_mean_safe(vals: np.ndarray, weights: np.ndarray) -> float:
    """RCR LS_MODE_DL weighted bulk-rejection mean. Falls back to a
    plain weighted mean if RCR's DL fit hits a degenerate-input
    pathology."""
    if vals.size == 0:
        return float("nan")
    if vals.size == 1:
        return float(vals[0])
    spread = float(np.max(vals) - np.min(vals))
    if spread <= 0:
        return float(vals[0])
    rng = np.random.default_rng(0)
    jitter = rng.normal(0, max(spread * 1e-12, 1e-15), vals.size)
    try:
        rcr = RCR(RejectionTech.LS_MODE_DL)
        rcr.perform_bulk_rejection(vals + jitter, w=weights)
        mu = float(rcr.result.mu)
        if not math.isfinite(mu):
            return float(np.average(vals, weights=weights))
        return mu
    except (ZeroDivisionError, ValueError, IndexError, RuntimeError,
            OverflowError, ArithmeticError):
        return float(np.average(vals, weights=weights))


# -----------------------------------------------------------------------------
# Public stage entry point
# -----------------------------------------------------------------------------

def _auto_centroid(
    xy: np.ndarray,
    flux: np.ndarray,
    sigma_1d_per_sample: np.ndarray,
    psf: float,
    has_local_model: np.ndarray,
    sigma_threshold: float = 75.0,
) -> np.ndarray:
    """Find bright-source centroids for the extra-anchor pre-pass.

    Port of ProcessorRFI::autoCentroid (Paper Footnote 21). Flags samples
    with `flux > sigma_threshold * scatter`, clusters them by 1-BW
    proximity (so a single source produces one centroid), then computes
    a weighted center-of-mass for each cluster. Used to place EXTRA
    local-model anchors on a grid around bright sources so the cos²
    shape can fit the source profile rather than rejecting it.

    Returns (n_centroids, 2) array of (x, y) centroid positions.

    Note: the C++ default of 75σ was chosen for very bright RFI spikes;
    real astronomical sources (Jupiter at 21σ in our data) sit well
    below it. The threshold is tunable from the caller.
    """
    ratios = flux / np.where(sigma_1d_per_sample > 0, sigma_1d_per_sample, 1.0)
    candidate_mask = (ratios > sigma_threshold) & has_local_model
    if not candidate_mask.any():
        return np.empty((0, 2), dtype=np.float64)

    cand_idx = np.where(candidate_mask)[0]
    # Sort by flux DESCENDING so we process the brightest source first.
    cand_sorted = cand_idx[np.argsort(-flux[cand_idx])]

    consumed = np.zeros(xy.shape[0], dtype=bool)
    centroids: list[tuple[float, float]] = []

    for idx in cand_sorted:
        if consumed[idx]:
            continue
        peak_xy = xy[idx]
        # All candidate samples within 1 BW (psf) of this peak — these
        # form ONE source cluster and get consumed together.
        dx = xy[cand_idx, 0] - peak_xy[0]
        dy = xy[cand_idx, 1] - peak_xy[1]
        within_bw = np.hypot(dx, dy) < psf
        cluster = cand_idx[within_bw]
        if cluster.size < 6:
            # C++ requires counter >= 6 in autoCentroid to accept a
            # centroid. Below that, treat the cluster as a single dot
            # (probably an isolated cal-glitch).
            consumed[cluster] = True
            continue
        # Weighted center-of-mass within the cluster (proximity-weighted
        # cos at distance from peak — same weight C++ uses inputs to
        # determineCenters).
        d = np.hypot(xy[cluster, 0] - peak_xy[0],
                     xy[cluster, 1] - peak_xy[1])
        w = np.cos(math.pi * d / (2.0 * psf)) ** 4  # peaked at d=0
        if w.sum() <= 0:
            cx, cy = float(peak_xy[0]), float(peak_xy[1])
        else:
            cx = float((xy[cluster, 0] * w).sum() / w.sum())
            cy = float((xy[cluster, 1] * w).sum() / w.sum())

        # Reject if too close to an existing centroid (avoid duplicates).
        too_close = False
        for ex, ey in centroids:
            if math.hypot(cx - ex, cy - ey) <= psf:
                too_close = True
                break
        if not too_close:
            centroids.append((cx, cy))
        consumed[cluster] = True

    return np.asarray(centroids, dtype=np.float64) if centroids \
        else np.empty((0, 2), dtype=np.float64)


def _extra_anchor_positions(
    centroid_xy: np.ndarray,
    psf: float,
    rfi_scale_bw: float,
    standard_gap_deg: float,
) -> np.ndarray:
    """Generate a grid of extra-anchor positions around one centroid.

    Port of ProcessorRFI::rfiBuildExtraModels (ProcessorRFI.cpp:531-568).
    Grid spacing is `standardGap * (1/rfi_scale² - 1)`, total extent
    bounded by `2 * sampleWidth` where `sampleWidth = psf * (1 - rfi_scale²)`.

    IMPORTANT: when extra_steps == 0 the C++ for-loop still runs ONCE
    with i=j=0 — placing a single extra anchor at the centroid itself.
    The earlier port hit early-return-empty here, which silently
    disabled auto-centroiding for any (psf, standard_gap, rfi_scale)
    config where extra_steps came out to 0 (which is nearly always the
    case for the reference cyga config — psf=0.665°, standard_gap≈0.14°,
    rfi_scale=0.35, so extra_steps = floor(0.604) = 0).
    """
    if rfi_scale_bw <= 0 or rfi_scale_bw >= 1:
        return np.empty((0, 2), dtype=np.float64)
    extra_steps = int(math.floor((psf / standard_gap_deg) * (rfi_scale_bw ** 2)))
    # NOTE: C++ does NOT early-return on extra_steps==0. The for-loops
    # iterate from -extra_steps to +extra_steps inclusive, so when
    # extra_steps==0 we still get one anchor at the centroid (i=j=0).
    if extra_steps > 0 and extra_steps % 2 != 0:
        extra_steps += 1

    # When extra_steps==0 standardGap*(1/rfi_scale²-1) is irrelevant
    # because i=j=0 anyway; guard against div-by-zero on the formula.
    step = standard_gap_deg * (1.0 / (rfi_scale_bw ** 2) - 1.0) if standard_gap_deg > 0 else 0.0
    sample_width = psf * (1.0 - rfi_scale_bw ** 2)
    max_dist = 2.0 * sample_width

    positions: list[tuple[float, float]] = []
    for i in range(-extra_steps, extra_steps + 1):
        for j in range(-extra_steps, extra_steps + 1):
            dx = i * step
            dy = j * step
            if math.hypot(dx, dy) > max_dist:
                continue
            positions.append((centroid_xy[0] + dx, centroid_xy[1] + dy))
    return np.asarray(positions, dtype=np.float64) if positions \
        else np.empty((0, 2), dtype=np.float64)


def subtract(target: Survey | Composite, scale_bw: float = 0.8,
             centroid_sigma: float = 75.0,
             photometry: bool = False) -> Survey | Composite:
    """Run 2D RFI subtraction on a Survey or Composite.

    Args:
        target: Survey or Composite. Each Scan must have `flux_bg` (or
            `flux`) populated and a 2D noise estimate (`noise_2d` or
            `noise_1d`).
        scale_bw: RFI subtraction scale in beamwidths. Default 0.8 —
            matches the C++ defaults for 20-meter L-band (paper §3.6
            Table 2).
        centroid_sigma: signal-to-noise threshold for auto-centroiding
            (Paper Footnote 21). The C++ default is 75 (Dan's value);
            real astronomical sources at 10-50σ won't trigger that, so
            the local-model RFI rejection eats them. Lower this if
            bright sources are getting subtracted as RFI.
    """
    surveys = target.surveys if isinstance(target, Composite) else [target]
    if not surveys:
        return target
    psf = surveys[0].psf_fwhm
    if psf <= 0:
        for survey in surveys:
            for scan in survey.scans:
                base = (scan.flux_bg if scan.flux_bg is not None
                        else (scan.flux if scan.flux is not None
                              else scan.flux_l))
                scan.flux_rfi = base.copy()
                scan.rfi_keep_mask = np.ones(scan.size, dtype=bool)
                scan.theta_corr = np.zeros(scan.size)
        return target

    rfi_scale_deg = scale_bw * psf

    # Flatten the scan data into one global array for kd-tree queries.
    # (For Composite, this spans all surveys.)
    if isinstance(target, Composite):
        all_scans = target.all_scans
    else:
        all_scans = target.scans

    # For each sample we need: x, y position; flux; dumps; sigma; and a
    # back-reference to (scan_idx, sample_idx) so we can write the result.
    xs, ys, fs, ds, sgs, sids, smps = [], [], [], [], [], [], []
    for sid, scan in enumerate(all_scans):
        if scan.size == 0:
            continue
        x = (scan.ra_ts if scan.ra_ts is not None else
             (scan.ra_proj if scan.ra_proj is not None else scan.ra))
        y = (scan.dec_ts if scan.dec_ts is not None else
             (scan.dec_proj if scan.dec_proj is not None else scan.dec))
        f = (scan.flux_bg if scan.flux_bg is not None else
             (scan.flux if scan.flux is not None else scan.flux_l))
        d = scan.dumps
        # Per-sample noise estimate: prefer 2D, fall back to 1D.
        if scan.noise_2d is not None:
            s = scan.noise_2d
        elif scan.noise_1d is not None:
            s = scan.noise_1d
        else:
            # Last-resort: estimate from the data itself.
            s = np.full(scan.size, max(float(np.std(f)) / 5.0, 1e-9))
        xs.append(x)
        ys.append(y)
        fs.append(f)
        ds.append(d)
        sgs.append(s)
        sids.append(np.full(scan.size, sid, dtype=np.int32))
        smps.append(np.arange(scan.size, dtype=np.int32))

    xy = np.column_stack((np.concatenate(xs), np.concatenate(ys)))
    flux_all = np.concatenate(fs)
    dumps_all = np.concatenate(ds)
    sigma_all = np.concatenate(sgs)
    scan_id_all = np.concatenate(sids)
    sample_id_all = np.concatenate(smps)

    tree = cKDTree(xy)

    # Per-sample accumulators. Each sample collects (value, weight,
    # w_dist_sq, n_sum, rfi_count) tuples from every anchor that included it.
    n_total = xy.shape[0]
    bg_vals: list[list[float]] = [[] for _ in range(n_total)]
    bg_ws: list[list[float]] = [[] for _ in range(n_total)]
    bg_dsq: list[list[float]] = [[] for _ in range(n_total)]
    bg_nsum: list[list[float]] = [[] for _ in range(n_total)]
    bg_rfic: list[list[float]] = [[] for _ in range(n_total)]

    neighbours = tree.query_ball_point(xy, r=rfi_scale_deg)

    # Pass 1: regular anchors centered at every sample position.
    has_local_model = np.zeros(n_total, dtype=bool)
    for i in range(n_total):
        cand_idx = np.asarray(neighbours[i], dtype=np.int64)
        if cand_idx.size < 2:
            continue
        out_idx, out_val, out_w, out_dsq, out_nsum, out_rfic = _anchor_local_model(
            xy[i], float(flux_all[i]), float(sigma_all[i]),
            cand_idx, xy, flux_all, dumps_all, sigma_all,
            rfi_scale_deg,
        )
        has_local_model[i] = out_idx.size > 0
        for k, idx in enumerate(out_idx):
            bg_vals[int(idx)].append(float(out_val[k]))
            bg_ws[int(idx)].append(float(out_w[k]))
            bg_dsq[int(idx)].append(float(out_dsq[k]))
            bg_nsum[int(idx)].append(float(out_nsum[k]))
            bg_rfic[int(idx)].append(float(out_rfic[k]))

    # Pass 2: auto-centroiding (Paper Footnote 21). Find bright-source
    # peaks and place EXTRA anchors on a fine grid around each one. The
    # off-center anchors can fit Jupiter's cos²-like profile with a
    # large positive f, so the source survives the rejection cycle and
    # contributes properly to the per-sample LM accumulators.
    sigma_1d_per_sample = np.concatenate([
        scan.noise_1d if scan.noise_1d is not None else
        np.full(scan.size, sigma_all[scan_id_all == sid].mean() if (scan_id_all == sid).any() else 1.0)
        for sid, scan in enumerate(all_scans) if scan.size > 0
    ])
    centroids = _auto_centroid(
        xy, flux_all, sigma_1d_per_sample, psf,
        has_local_model, sigma_threshold=centroid_sigma,
    )
    if centroids.shape[0] > 0:
        # Standard gap = median nearest-neighbour distance within the
        # data (C++'s "standardGap" — characteristic sample spacing).
        # Use kd-tree to estimate it cheaply.
        nn_dist, _ = tree.query(xy, k=2)
        standard_gap_deg = float(np.median(nn_dist[:, 1]))
        if standard_gap_deg <= 0:
            standard_gap_deg = psf * 0.1
        for ci in range(centroids.shape[0]):
            extra_pos = _extra_anchor_positions(
                centroids[ci], psf, scale_bw, standard_gap_deg,
            )
            for ap in extra_pos:
                # Neighbours of THIS extra anchor (within rfi_scale)
                cand_idx = np.asarray(tree.query_ball_point(ap, r=rfi_scale_deg),
                                      dtype=np.int64)
                if cand_idx.size < 2:
                    continue
                # Use anchor's own (interpolated) sigma — pick nearest
                # sample's noise_2d as a stand-in. C++ uses -999999 as a
                # sentinel meaning "no constraint", but for the local
                # model fit we still need a real positive sigma.
                anchor_sigma = float(np.median(sigma_all[cand_idx]))
                out_idx, out_val, out_w, out_dsq, out_nsum, out_rfic = _anchor_local_model(
                    ap, 999999.0, anchor_sigma,
                    cand_idx, xy, flux_all, dumps_all, sigma_all,
                    rfi_scale_deg,
                )
                for k, idx in enumerate(out_idx):
                    bg_vals[int(idx)].append(float(out_val[k]))
                    bg_ws[int(idx)].append(float(out_w[k]))
                    bg_dsq[int(idx)].append(float(out_dsq[k]))
                    bg_nsum[int(idx)].append(float(out_nsum[k]))
                    bg_rfic[int(idx)].append(float(out_rfic[k]))

    # Global model per sample. Also accumulate:
    #   theta_corr = 2 × 3 × sqrt(sum(w_dist_sq) / max(sum(n_sum), 1))  [Appendix D]
    #   gm_weight  = sum(LMWeight / rfiCount) over ALL contributing anchors
    #                (NOT filtered by RCR survival, per ProcessorRFI.cpp:804)
    factor_corr = 3.0
    global_flux = np.full(n_total, np.nan)
    global_theta_corr = np.zeros(n_total, dtype=np.float64)
    global_gm_weight = np.zeros(n_total, dtype=np.float64)
    for i in range(n_total):
        if not bg_vals[i]:
            continue
        v = np.array(bg_vals[i], dtype=np.float64)
        w = np.array(bg_ws[i], dtype=np.float64)
        dsq = np.array(bg_dsq[i], dtype=np.float64)
        nsum = np.array(bg_nsum[i], dtype=np.float64)
        rfic = np.array(bg_rfic[i], dtype=np.float64)
        ok = np.isfinite(v) & np.isfinite(w) & (w > 0)
        if not ok.any():
            continue
        v_ok, w_ok = v[ok], w[ok]
        dsq_ok, nsum_ok = dsq[ok], nsum[ok]
        rfic_ok = rfic[ok]

        # GMWeight = sum over ALL contributing anchors of LMWeight/rfiCount.
        # C++ ProcessorRFI.cpp:804 — NOT filtered by RCR flags.
        global_gm_weight[i] = float(
            np.where(rfic_ok > 0, w_ok / rfic_ok, w_ok).sum()
        )

        # Run RCR and use ITS flags to identify which local models
        # survived — those are the ones whose w_dist_sq / n_sum get
        # accumulated for theta_corr.
        mu, keep_flags = _rcr_mean_with_flags(v_ok, w_ok)
        global_flux[i] = mu
        # C++ ProcessorRFI.cpp:705-722: the wDistSq accumulation is
        # gated on `correlatedWeightMap` (= photometryOn). With
        # photometry off, theta_corr stays at zero and the downstream
        # correlation layer reduces to `factor_w * scale` exactly.
        if photometry and keep_flags.any():
            d_sum = float(dsq_ok[keep_flags].sum())
            n_sum_v = float(nsum_ok[keep_flags].sum())
            global_theta_corr[i] = 2.0 * factor_corr * math.sqrt(
                d_sum / max(n_sum_v, 1.0)
            )

    # Write back to each scan's flux_rfi / rfi_keep_mask / theta_corr / gm_weight.
    for sid, scan in enumerate(all_scans):
        if scan.size == 0:
            continue
        mask = scan_id_all == sid
        idx_in_scan = sample_id_all[mask]
        vals = global_flux[mask]
        tc = global_theta_corr[mask]
        gmw = global_gm_weight[mask]
        order = np.argsort(idx_in_scan)
        idx_in_scan = idx_in_scan[order]
        vals = vals[order]
        tc = tc[order]
        gmw = gmw[order]
        flux_rfi = np.full(scan.size, np.nan)
        theta_corr = np.zeros(scan.size, dtype=np.float64)
        gm_weight = np.zeros(scan.size, dtype=np.float64)
        flux_rfi[idx_in_scan] = vals
        theta_corr[idx_in_scan] = tc
        gm_weight[idx_in_scan] = gmw
        keep_mask = np.isfinite(flux_rfi)
        scan.flux_rfi = flux_rfi
        scan.rfi_keep_mask = keep_mask
        scan.theta_corr = theta_corr  # in degrees (same units as positions)
        scan.gm_weight = gm_weight

    return target
