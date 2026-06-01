"""Stage 5 — time-delay correction (Paper §3.4).

This is a faithful port of the C++ ProcessorTS algorithm:

    1. For each pair of adjacent scans, interpolate both onto a common
       dense angle grid (RA for ra-along scans, Dec otherwise) and
       cross-correlate via FFT. The argmax bin × bin-width gives the
       angular shift between the pair.
    2. Apply the sign-flip step:
           shifts[j] *= -1 * (-1)^j
       Adjacent raster scans go in opposite directions, so a real time
       delay produces alternating ± shifts. The sign flip converts the
       alternating pattern into a CONSTANT-sign series so the mean is
       meaningful.
    3. Run Footnote-28a probability rejection (shiftRejection): for
       each shift δ_i, compute the probability that |δ_i - δ_neighbor|
       could occur by chance given the scan length R. Reject shifts
       where this probability exceeds 1/(2N).
    4. Run RCR bulk rejection (LS_MODE_DL) on the survivors. The mean
       of the post-rejection set is the per-pair shift.
    5. Divide by 2 to get the per-scan offset (each scan is shifted in
       the opposite direction by half the pair-shift).
    6. Convert to a time shift by dividing by the RCR-median slew speed.

References:
    C++ source: src/ProcessorTS.cpp.
    Paper §3.4, Figure 28, Footnote 28a.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.signal import correlate

from rcpy.types import Survey, MapType
from rcpy.rcr2 import RCR, RejectionTech


# Oversampling factor for the FFT cross-correlation grid (the C++ uses
# a power-of-2 grid with a `powerOffset` of ~4 bits — i.e. 16× the
# native sampling). 30 is what my naive impl used; the C++ goes denser
# but we get sub-sample precision from the FFT anyway.
_OVERSAMPLE = 32


def _daisy_center_index(scan) -> int:
    """Find the sample index closest to projected (0, 0). Port of
    Survey::findCenters (Survey.cpp:1553-1576): for each daisy petal,
    the dish passes through the source center exactly once — that
    sample is the petal's "center index", used by daisyAngleBuilder
    to flip the sign convention so adjacent petals become like
    parallel raster scans in 1D angle space.
    """
    ra_p = scan.ra_proj if scan.ra_proj is not None else scan.ra
    dec_p = scan.dec_proj if scan.dec_proj is not None else scan.dec
    if ra_p is None or dec_p is None or ra_p.size == 0:
        return 0
    # Distance from projected origin in projected (small-angle planar) space.
    r2 = ra_p * ra_p + dec_p * dec_p
    return int(np.argmin(r2))


def _daisy_angle_builder(scan, petal_index: int) -> np.ndarray:
    """Port of Tools::daisyAngleBuilder (Tools.cpp:1166-1183).

    For each sample j, returns the signed great-circle distance from
    the petal's center sample to sample j, with a sign convention that
    makes adjacent petals look like parallel raster scans:

      negation = -(-1)^i   (initial, alternates with petal index)
      for each sample j:
          if j == centerIndex: negation *= -1   (flip at center)
          angle[j] = distance(sample_j, center_sample) * negation

    Result: petal i goes monotonically from -R to +R (or +R to -R
    depending on parity) through angle 0 at the center sample. Adjacent
    petals span the SAME angle range but in OPPOSITE directions —
    exactly the geometry that the existing pair-correlation algorithm
    expects.
    """
    ra = scan.ra_proj if scan.ra_proj is not None else scan.ra
    dec = scan.dec_proj if scan.dec_proj is not None else scan.dec
    if ra is None or dec is None or ra.size == 0:
        return np.zeros(0)
    center_idx = _daisy_center_index(scan)
    # Planar distance (degrees) from each sample to the center sample.
    # The projected coords are already locally-flat, so simple Euclidean
    # in (ra_proj, dec_proj) is the small-angle approximation of the
    # great-circle distance C++ does via getGCDistance.
    dx = ra - ra[center_idx]
    dy = dec - dec[center_idx]
    distance = np.hypot(dx, dy)
    # Signed by petal parity.
    negation = -((-1.0) ** petal_index)
    sign = np.full(ra.size, negation, dtype=np.float64)
    # Flip at and after the center sample (C++ flips when j == centerIndex,
    # and the flip persists for subsequent samples).
    sign[center_idx:] *= -1.0
    return distance * sign


def _pair_shift(s_a, s_b, i_a: int = 0, is_daisy: bool = False) -> tuple[float, float]:
    """Cross-correlate two adjacent scans along the scan axis.

    Both scans are interpolated onto a common ascending-angle grid that
    spans their overlap region with `_OVERSAMPLE`-times the native
    sample density. The grid is power-of-2 sized for FFT efficiency.

    For DAISY maps, the 1D "scan axis" coord is the per-petal signed
    angle from the center sample (see `_daisy_angle_builder`), which
    converts the daisy geometry into raster-like parallel/antiparallel
    1D scans the algorithm can cross-correlate. For RASTERS the axis
    is simply ra_proj or dec_proj depending on scan_in_ra.

    Returns (angular_shift_deg, sqrt_corr_max) where the angular shift
    is the lag at which scan_b best aligns with scan_a.
    """
    y_a = s_a.working_flux()
    y_b = s_b.working_flux()
    if y_a.size < 4 or y_b.size < 4:
        return 0.0, 0.0

    if is_daisy:
        # Port of ProcessorTS.cpp:41-44.
        x_a = _daisy_angle_builder(s_a, i_a)
        x_b = _daisy_angle_builder(s_b, i_a + 1)
    elif s_a.scan_in_ra:
        x_a = s_a.ra_proj if s_a.ra_proj is not None else s_a.ra
        x_b = s_b.ra_proj if s_b.ra_proj is not None else s_b.ra
    else:
        x_a = s_a.dec_proj if s_a.dec_proj is not None else s_a.dec
        x_b = s_b.dec_proj if s_b.dec_proj is not None else s_b.dec

    # Sort each scan ascending (alternating-direction scans need this
    # so they share a common grid orientation).
    o_a = np.argsort(x_a)
    o_b = np.argsort(x_b)
    x_a = x_a[o_a]; y_a = y_a[o_a]
    x_b = x_b[o_b]; y_b = y_b[o_b]

    # Overlap window.
    lo = max(x_a[0], x_b[0])
    hi = min(x_a[-1], x_b[-1])
    if hi <= lo:
        return 0.0, 0.0

    # Native step per scan, then chosen grid step.
    span_a = x_a[-1] - x_a[0]
    span_b = x_b[-1] - x_b[0]
    step = min(span_a / max(s_a.size - 1, 1), span_b / max(s_b.size - 1, 1)) / _OVERSAMPLE
    if step <= 0:
        return 0.0, 0.0

    # Pad up to a power-of-2 length to keep FFT happy and gain sub-bin
    # precision.
    n_bins = int((hi - lo) / step)
    if n_bins < 4:
        return 0.0, 0.0
    n_pow2 = 1 << int(math.ceil(math.log2(n_bins)))
    grid = lo + np.arange(n_pow2) * step

    y_a_i = np.interp(grid, x_a, y_a)
    y_b_i = np.interp(grid, x_b, y_b)
    y_a_i = y_a_i - y_a_i.mean()
    y_b_i = y_b_i - y_b_i.mean()

    corr = correlate(y_b_i, y_a_i, mode="full")
    norm = np.sqrt((y_a_i ** 2).sum() * (y_b_i ** 2).sum())
    if norm <= 0:
        return 0.0, 0.0
    corr = corr / norm

    peak = int(np.argmax(corr))
    lag_bins = peak - (y_a_i.size - 1)
    shift = lag_bins * step
    return float(shift), float(np.sqrt(max(corr[peak], 0.0)))


def _shift_rejection(shifts: np.ndarray, scan_span_deg: float) -> np.ndarray:
    """Footnote-28a probability rejection of pair shifts.

    For two adjacent shifts δ_i and δ_{i+1}, the probability that they
    could be coincidentally close by chance is

        p_{i,i+1} = (|δ_{i+1} - δ_i| / R) * (2 - |δ_{i+1} - δ_i| / R)

    where R is the scan length in the same units as the shifts. For a
    shift to be CONSISTENT (i.e. driven by real signal, not noise), it
    must agree with both its preceding and proceeding neighbours — the
    joint probability is 2 × p_{i-1,i} × p_{i,i+1}. We REJECT shifts
    whose joint probability exceeds 1/(2N), i.e. the Chauvenet
    threshold for this sample size.

    Note the C++ convention: probHold > 1/(2N) means REJECT (this is
    the OPPOSITE of how Chauvenet's criterion is usually stated, but
    here a "high" probability means "this shift is INconsistent with
    its neighbours, so reject"). I had to stare at the C++ for a while
    to convince myself.

    Returns boolean keep-mask of the same length as `shifts`.
    """
    n = shifts.size
    if n == 0:
        return np.zeros(0, dtype=bool)
    if scan_span_deg <= 0:
        return np.ones(n, dtype=bool)

    keep = np.ones(n, dtype=bool)
    threshold = 1.0 / (2.0 * n)

    for i in range(n):
        delta_i = shifts[i]
        if i == 0:
            j = 1 if n > 1 else 0
            delta_j = shifts[j]
            d = abs(delta_i - delta_j) / scan_span_deg
            prob = d * (2.0 - d)            # M=1 in C++, pow(2, 0) = 1
        elif i == n - 1:
            j = n - 2 if n > 1 else 0
            delta_j = shifts[j]
            d = abs(delta_i - delta_j) / scan_span_deg
            prob = d * (2.0 - d)            # M=1
        else:
            d1 = abs(delta_i - shifts[i + 1]) / scan_span_deg
            p1 = d1 * (2.0 - d1)
            d2 = abs(delta_i - shifts[i - 1]) / scan_span_deg
            p2 = d2 * (2.0 - d2)
            prob = 2.0 * p1 * p2            # M=2, pow(2, 1) = 2

        if prob > threshold:
            keep[i] = False

    return keep


def _scan_span(scan, is_daisy: bool = False, petal_index: int = 0) -> float:
    """Angular length of a scan along its scan axis.

    For daisies, the axis is the per-petal signed angle-from-center
    coord built by `_daisy_angle_builder` — a typical petal spans
    roughly [-R, +R] for R ≈ 0.5 × petal length, so the span is ~2R.
    Using the RA/Dec span instead under-counts and tightens the
    Footnote-28a rejection too much (would drop real-signal shifts).
    """
    if is_daisy:
        x = _daisy_angle_builder(scan, petal_index)
    elif scan.scan_in_ra:
        x = scan.ra_proj if scan.ra_proj is not None else scan.ra
    else:
        x = scan.dec_proj if scan.dec_proj is not None else scan.dec
    if x is None or x.size < 2:
        return 0.0
    return float(x.max() - x.min())


def _rcr_mean_safe(values: np.ndarray) -> float:
    """RCR bulk-rejection mean, robust to the degenerate-input pathology
    in RCR2's DL fit. The C++ RCR's `performBulkRejection` with
    LS_MODE_DL crashes here on inputs with many duplicate values
    (signature: ZeroDivisionError in fitDL_w). We jitter the inputs
    with ULP-scale noise that's negligible vs the data spread but
    breaks the degeneracy. The mean is robust to noise of this size.
    """
    if values.size == 0:
        return 0.0
    if values.size == 1:
        return float(values[0])
    # Use a jitter ~12 orders of magnitude below the data range so it
    # can't shift the mean meaningfully, but is large enough to break
    # exact-duplicate degeneracy.
    spread = float(np.max(values) - np.min(values))
    if spread <= 0:
        return float(values[0])
    rng = np.random.default_rng(0)
    jitter = rng.normal(0, max(spread * 1e-12, 1e-15), values.size)
    try:
        rcr = RCR(RejectionTech.LS_MODE_DL)
        rcr.perform_bulk_rejection(values + jitter)
        mu = float(rcr.result.mu)
        if not np.isfinite(mu):
            return float(np.median(values))
        return mu
    except (ZeroDivisionError, ValueError, IndexError, RuntimeError, OverflowError, ArithmeticError):
        return float(np.median(values))


def _scan_slew_speeds(survey: Survey, is_daisy: bool = False) -> np.ndarray:
    """Per-scan slew speed (deg/sec) along the scan axis.

    For rasters: RCR-median of per-sample point-to-point speeds within
    each scan along ra_proj or dec_proj. Matches the C++ `speedTemp` ->
    `performBulkRejection(speedTemp)` pattern in ProcessorTS.cpp:418-424.

    For daisies (ProcessorTS.cpp:393-403): uses the LOCAL speed at the
    petal-center crossing via cumulative arc-length (ang_dist) at
    centerIndex ± 1, i.e.

        speed = (angDist[c+1] - angDist[c-1]) / (time[c+1] - time[c-1])

    Cumulative ang_dist is monotonically increasing along time, so the
    speed is always POSITIVE regardless of petal direction — unlike the
    signed angle-from-center coord, which alternates sign per petal and
    would give a misleading mean speed.
    """
    speeds = []
    for i, scan in enumerate(survey.scans):
        if scan.size < 3:
            continue
        if is_daisy:
            ang_dist = scan.ang_dist
            if ang_dist is None or ang_dist.size != scan.size:
                continue
            c = _daisy_center_index(scan)
            if not (1 <= c <= scan.size - 2):
                continue
            dt = scan.time[c + 1] - scan.time[c - 1]
            if dt <= 0:
                continue
            speed = (ang_dist[c + 1] - ang_dist[c - 1]) / dt
            speeds.append(float(speed))
            continue
        if scan.scan_in_ra:
            x = scan.ra_proj if scan.ra_proj is not None else scan.ra
        else:
            x = scan.dec_proj if scan.dec_proj is not None else scan.dec
        if x is None or x.size < 2:
            continue
        dt = np.diff(scan.time)
        dx = np.diff(x)
        valid = dt > 0
        if not valid.any():
            continue
        sample_speeds = dx[valid] / dt[valid]
        # RCR bulk rejection on the per-sample speeds (C++ uses LS_MODE_DL).
        speeds.append(_rcr_mean_safe(sample_speeds))
    return np.array(speeds, dtype=np.float64)


def _auto_ts_daisy(survey: Survey) -> float:
    """Auto-determine the time-shift for a daisy by maximizing source
    concentration after a per-scan BG-subtracted coarse rebinning.

    The raster pair-correlation algorithm fails on daisies because
    adjacent petals are at different angles, not parallel/antiparallel.
    Here we use a sample-level metric that's robust to petal geometry:

      1. Per scan, subtract the scan's median flux (cheap BG estimate).
         Each petal passes the source briefly; the residual flux at
         on-source samples is what we're trying to localize.
      2. For each candidate TS, apply the shift to each sample's
         projected position, drop samples into a coarse 2D histogram
         weighted by BG-subbed flux, and measure the peak bin value.
      3. The TS that maximizes the peak is the one that piles the
         source-on samples into one bin (instead of smearing them
         across several due to encoder lag).

    Returns the best time-shift in seconds. Note: the daisy sign
    convention is OPPOSITE the raster convention — daisies typically
    need POSITIVE TS while rasters need NEGATIVE. See memory note
    `rcpy_daisy_handling`.
    """
    # Per-scan BG: subtract scan median flux so the source signal isn't
    # buried in the overall BG floor / elevation drift.
    bg_subbed = []  # list of (time, ra_proj, dec_proj, flux_resid) per scan
    for sc in survey.scans:
        if sc.flux is None or sc.size < 3:
            continue
        med = float(np.median(sc.flux))
        bg_subbed.append((
            sc.time, sc.ra_proj if sc.ra_proj is not None else sc.ra,
            sc.dec_proj if sc.dec_proj is not None else sc.dec,
            sc.flux - med,
        ))
    if not bg_subbed:
        return 0.0

    # Concatenate all samples and compute the rough projected map extent.
    all_ra = np.concatenate([b[1] for b in bg_subbed])
    all_dec = np.concatenate([b[2] for b in bg_subbed])
    # Use the inner 80% to set bin bounds — keeps bin width sensible
    # even when a few samples shoot out to the extreme edges.
    ra_lo, ra_hi = np.percentile(all_ra, [10, 90])
    dec_lo, dec_hi = np.percentile(all_dec, [10, 90])
    # Bins sized at ¼ of the PSF FWHM. Larger bins (full PSF) make the
    # discretization too coarse — the metric peak shifts by 0.5 s
    # depending on which side of a bin centroid the source's best-
    # alignment position lands on. ¼-PSF gives sub-PSF discretization
    # while still keeping enough samples per bin for the count≥2 cut.
    bin_size = max(0.25 * survey.psf_fwhm, 0.01)
    n_ra = max(int(np.ceil((ra_hi - ra_lo) / bin_size)), 5)
    n_dec = max(int(np.ceil((dec_hi - dec_lo) / bin_size)), 5)

    total_in_samples = sum(b[0].size for b in bg_subbed)

    # Build a HARD central mask covering bins within 1.5 PSF FWHM of
    # the projected-coord origin. For tracked moving targets the
    # source sits within ~1 beamwidth of the dish center; this hard
    # cut suppresses petal-tip artifacts whose brightness can rival
    # the source's contribution after wrong-TS pile-up. Soft Gaussian
    # weighting wasn't sharp enough — an edge spike at 1.8° from
    # center still survived the e^-3.5 ≈ 3% damping.
    ra_edges = np.linspace(ra_lo, ra_hi, n_ra + 1)
    dec_edges = np.linspace(dec_lo, dec_hi, n_dec + 1)
    ra_cen = 0.5 * (ra_edges[:-1] + ra_edges[1:])
    dec_cen = 0.5 * (dec_edges[:-1] + dec_edges[1:])
    Dgrid, Rgrid = np.meshgrid(dec_cen, ra_cen, indexing="ij")
    dist_grid = np.hypot(Rgrid, Dgrid)
    central_mask = dist_grid <= 1.5 * survey.psf_fwhm

    def _metric_at_ts(ts: float) -> float:
        # Apply TS per sample by linear interpolation in time. Mask out
        # samples whose target_t falls outside the scan's recorded
        # time range; np.interp would clamp them, piling many samples
        # onto a single endpoint position and producing a fake peak.
        peak_bins = np.zeros((n_dec, n_ra), dtype=np.float64)
        count_bins = np.zeros((n_dec, n_ra), dtype=np.int32)
        kept_total = 0
        for (t, ra_p, dec_p, resid) in bg_subbed:
            n = ra_p.size
            if n < 2:
                continue
            target_t = t - ts
            in_range = (target_t >= t[0]) & (target_t <= t[-1])
            if in_range.sum() < 3:
                continue
            tt = target_t[in_range]
            ra_shifted = np.interp(tt, t, ra_p)
            dec_shifted = np.interp(tt, t, dec_p)
            r = resid[in_range]
            ix = ((ra_shifted - ra_lo) / (ra_hi - ra_lo) * n_ra).astype(int)
            iy = ((dec_shifted - dec_lo) / (dec_hi - dec_lo) * n_dec).astype(int)
            mask = (ix >= 0) & (ix < n_ra) & (iy >= 0) & (iy < n_dec)
            ix = ix[mask]; iy = iy[mask]; r = r[mask]
            np.add.at(peak_bins, (iy, ix), r)
            np.add.at(count_bins, (iy, ix), 1)
            kept_total += ix.size
        # Require at least 2 samples per bin to count it (single-sample
        # spikes are noise, not pile-up).
        valid = count_bins >= 2
        if not valid.any():
            return 0.0
        mean_bin = np.where(valid, peak_bins / np.maximum(count_bins, 1), 0.0)
        # Max within the central PSF-sized region. Bins outside the
        # mask (petal tips etc.) cannot contribute, so an off-center
        # artifact at large radius doesn't poison the metric.
        masked = np.where(central_mask, mean_bin, -np.inf)
        if not np.any(np.isfinite(masked) & (masked > -1e30)):
            return 0.0
        frac = kept_total / max(total_in_samples, 1)
        return float(np.max(masked)) * frac

    # Coarse sweep: ±4 sec in 0.5-sec steps.
    coarse_ts = np.arange(-4.0, 4.01, 0.5)
    coarse_scores = np.array([_metric_at_ts(t) for t in coarse_ts])
    best_coarse = float(coarse_ts[np.argmax(coarse_scores)])

    # Fine sweep: ±0.6 sec around the coarse winner in 0.1-sec steps.
    fine_ts = np.arange(best_coarse - 0.6, best_coarse + 0.61, 0.1)
    fine_scores = np.array([_metric_at_ts(t) for t in fine_ts])
    best_fine = float(fine_ts[np.argmax(fine_scores)])

    return best_fine


def correct(survey: Survey, mode: str = "auto", forced_seconds: float = 0.0) -> Survey:
    """Estimate and apply the time-delay correction (Paper §3.4)."""
    if mode == "off":
        survey.time_shift = 0.0
        _apply_shift(survey)
        return survey

    if mode == "custom":
        survey.time_shift = forced_seconds
        _apply_shift(survey)
        return survey

    # AUTO mode: same C++-faithful pair-correlation algorithm for both
    # rasters and daisies. The only difference is the 1D axis used for
    # cross-correlation:
    #   - rasters: ra_proj or dec_proj (depending on scan direction)
    #   - daisies: per-petal signed angle-from-center via
    #              `_daisy_angle_builder`, which converts adjacent
    #              petals into parallel/antiparallel raster-like
    #              1D scans the algorithm can correlate. Port of
    #              Tools::daisyAngleBuilder + ProcessorTS.cpp:41-44.
    is_daisy = survey.map_type is MapType.DAISY

    n_pairs = max(len(survey.scans) - 1, 0)
    if n_pairs == 0:
        survey.time_shift = 0.0
        _apply_shift(survey)
        return survey

    # 1. Per-pair cross-correlation.
    raw_shifts = np.zeros(n_pairs, dtype=np.float64)
    weights = np.zeros(n_pairs, dtype=np.float64)
    for i in range(n_pairs):
        shift, w = _pair_shift(survey.scans[i], survey.scans[i + 1],
                               i_a=i, is_daisy=is_daisy)
        raw_shifts[i] = shift
        weights[i] = w

    # 2. Sign-flip step (C++ line 159): shifts[j] *= -1 * (-1)^j.
    # This converts the alternating ± pattern produced by raster
    # direction reversals into a consistent-sign series.
    sign = np.array([-1.0 * ((-1.0) ** j) for j in range(n_pairs)])
    signed_shifts = raw_shifts * sign

    # 3. Scan span R = mean of scan span across the data, used for
    # Footnote-28a probability calculation. The C++ uses the LAST
    # pair's overlap window (a minor bug — uses the last-computed
    # minAngle/maxAngle), but in practice all scans in a raster have
    # similar lengths.
    spans = [_scan_span(sc, is_daisy=is_daisy, petal_index=i)
             for i, sc in enumerate(survey.scans)]
    spans = [s for s in spans if s > 0]
    scan_span = float(np.mean(spans)) if spans else 1.0

    # 4. Footnote-28a probability rejection on the signed shifts.
    keep_mask = _shift_rejection(signed_shifts, scan_span)
    survivors = signed_shifts[keep_mask]

    if survivors.size == 0:
        ang_shift = 0.0
    else:
        # 5. RCR bulk-reject on survivors, take mean (C++ uses LS_MODE_DL +
        # performBulkRejection). RCR2's DL fit has a division-by-zero
        # pathology on inputs with many duplicate values; jitter with
        # tiny noise to break the degeneracy without affecting results.
        ang_shift = _rcr_mean_safe(survivors)

    # 6. Divide by 2 — pair shift = 2 × per-scan offset (C++ line 430).
    ang_shift = ang_shift / 2.0

    # 7. Convert to time via RCR-medianed slew speed.
    speeds = _scan_slew_speeds(survey, is_daisy=is_daisy)
    if speeds.size > 0:
        slew = _rcr_mean_safe(speeds)
    else:
        slew = 0.0

    if abs(slew) > 1e-12:
        survey.time_shift = ang_shift / slew
    else:
        survey.time_shift = 0.0

    _apply_shift(survey)
    return survey


def _apply_shift(survey: Survey) -> None:
    """Apply the survey's time_shift to every scan as a per-sample
    angular displacement equal to (local_slew × time_shift).

    Using a local-slew × time approach (rather than np.interp on a
    shifted time axis) avoids the endpoint-clamping pathology where
    samples near scan ends collapse onto the last sample's coordinate.
    """
    # Port of C++ Scan::cosDecTransform (Scan.cpp:455-494).
    # For each sample j, find the encoder position at time (t[j] - ts)
    # by LINEAR INTERPOLATION between the two recorded samples that
    # bracket that target time. This is critical when slew is variable
    # — many real scans have decel/accel transients near scan-changes
    # where local slew ≈ 0 but the encoder moves significantly in the
    # following few seconds. Using a local-gradient shift gets these
    # samples wrong by 50-100x, producing S-curve scatter at scan
    # endpoints. The interpolation approach correctly extrapolates
    # using the slope of the nearest two recorded samples.
    ts = survey.time_shift
    for scan in survey.scans:
        ra_src = scan.ra_proj if scan.ra_proj is not None else scan.ra
        dec_src = scan.dec_proj if scan.dec_proj is not None else scan.dec
        if scan.size < 2:
            scan.ra_ts = ra_src.copy()
            scan.dec_ts = dec_src.copy()
            continue

        t = scan.time
        n = scan.size
        ra_new = np.empty(n, dtype=np.float64)
        dec_new = np.empty(n, dtype=np.float64)
        for j in range(n):
            target_time = t[j] - ts
            # C++ uses different bracket searches for ts >= 0 vs ts < 0
            # but the net effect is the same — find adjacent samples
            # (jMin, jMin+1) whose times bracket target_time. We do
            # the simpler approach via np.searchsorted, then clamp.
            jMax = int(np.searchsorted(t, target_time, side='right'))
            jMax = max(1, min(jMax, n - 1))
            jMin = jMax - 1
            t_min = t[jMin]
            t_max = t[jMax]
            if t_max > t_min:
                frac = (target_time - t_min) / (t_max - t_min)
            else:
                frac = 0.0
            ra_new[j] = ra_src[jMin] + frac * (ra_src[jMax] - ra_src[jMin])
            dec_new[j] = dec_src[jMin] + frac * (dec_src[jMax] - dec_src[jMin])

        scan.ra_ts = ra_new
        scan.dec_ts = dec_new
