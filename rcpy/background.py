"""Stage 4 — 1D background subtraction (Paper §3.3).

Faithful port of the C++ BackgroundCUDA::calculateBG / Baseline algorithm.

For each sample (the "anchor"), in each direction (forward and backward),
we build a local quadratic model that:

  1. Passes through (x_anchor, y_anchor) exactly (the "fixed" regression).
  2. Includes all samples within `baseline` (scale_bw beamwidths) of the
     anchor on that side.
  3. Iteratively rejects the highest-residual point until the RMS of
     ABOVE-model residuals (splus) is < the 1D noise level (scatter).
  4. Drops the anchor and refits with a free intercept (unfixed
     regression), then iteratively adds back the smallest-residual
     rejected point as long as splus stays < scatter.
  5. Emits (index, value, weight) triples for every surviving sample,
     using the position-confidence weight from Paper Eq. 3.

For each sample j, all covering local models contribute. We take an RCR
LS_MODE_DL bulk-rejection weighted mean (matching C++ setBackground) and
that's the background value at j. Samples where no local model survives
are linearly interpolated between neighbors.

References:
    Paper §3.3, Eqs. 3-4.
    C++ source: src/BackgroundCUDA.cpp (Baseline class, calculateBGMulti).
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass

from rcpy.types import Survey, Scan
from rcpy.rcr2 import RCR, RejectionTech
from rcpy._pivot import _perform_pivot


# -----------------------------------------------------------------------------
# Local model fits (quadratic, with fallback to linear when quad fails)
# -----------------------------------------------------------------------------

def _fixed_quadratic_fit(
    keep: np.ndarray, w: np.ndarray, x: np.ndarray, y: np.ndarray,
    x_anchor: float, y_anchor: float,
) -> tuple[np.ndarray | None, str]:
    """Step-for-step port of C++ `Tools::fixedQuadraticRegressionPivot`
    (Tools.cpp:1030-1078) and the surrounding autoFixedWRegression
    logic (BackgroundCUDA.cpp:88-160).

    Tries weighted quadratic anchored at (x_anchor, y_anchor) via the
    real C++ performPivot algorithm. Falls back to linear (via
    fixedRegressionPivot at Tools.cpp:997) if any quad coef is NaN, and
    to no-model if linear also fails. Returns [a, b, c] in unshifted
    form (so _apply_model evaluates the polynomial directly).
    """
    if not keep.any():
        return None, "none"

    sel = keep
    xs = x[sel] - x_anchor
    ys = y[sel] - y_anchor
    ws = w[sel]

    wxx   = float(np.sum(ws * xs * xs))
    wxxx  = float(np.sum(ws * xs * xs * xs))
    wxxxx = float(np.sum(ws * xs * xs * xs * xs))
    wyx   = float(np.sum(ws * ys * xs))
    wyxx  = float(np.sum(ws * ys * xs * xs))

    # Quadratic via real C++ performPivot
    A = np.array([wxx, wxxx, wxxx, wxxxx], dtype=np.float64)
    b = np.array([wyx, wyxx], dtype=np.float64)
    try:
        coef = _perform_pivot(2, A, b)
        if np.all(np.isfinite(coef)):
            b_p, c_p = float(coef[0]), float(coef[1])
            a_u = y_anchor - b_p * x_anchor + c_p * x_anchor * x_anchor
            b_u = b_p - 2.0 * c_p * x_anchor
            c_u = c_p
            return np.array([a_u, b_u, c_u]), "quad"
    except (FloatingPointError, ZeroDivisionError):
        pass

    # Linear fallback: y - y_a = m·(x - x_a), no pivot needed (1x1)
    # C++ fixedRegressionPivot (Tools.cpp:997-1028).
    wxy_lin = float(np.sum(ws * ys * xs))
    wxx_lin = float(np.sum(ws * xs * xs))
    if wxx_lin != 0.0:
        m = wxy_lin / wxx_lin
        if np.isfinite(m):
            a_u = y_anchor - m * x_anchor
            return np.array([a_u, m, 0.0]), "lin"

    return None, "none"


def _free_quadratic_fit(
    keep: np.ndarray, w: np.ndarray, x: np.ndarray, y: np.ndarray,
) -> tuple[np.ndarray | None, str]:
    """Step-for-step port of C++ `Tools::quadraticRegressionPivot`
    (Tools.cpp:939-995) and the surrounding autoWRegression logic
    (BackgroundCUDA.cpp:161-230).

    Tries weighted quadratic via the real C++ performPivot algorithm
    (backward-Gauss + lowerTriangleSolver) — NOT np.linalg.solve, which
    uses standard LU and produces different garbage on near-singular
    matrices like the {2 distinct kept samples} case that arises after
    rejectPoints. Falls back to linear (via regressionPivot at
    Tools.cpp:887) if any quad coef is NaN. The C++ does NOT short-
    circuit based on sample count — it always attempts quadratic first.
    """
    if not keep.any():
        return None, "none"

    sel = keep
    xs_raw = x[sel]
    ys_raw = y[sel]
    ws = w[sel]

    # C++ shifts by xAxis[0]/yAxis[0] (the first element of the FULL
    # array, not the first kept sample) just for solver conditioning.
    x_shift = float(x[0]) if x.size > 0 else 0.0
    y_shift = float(y[0]) if y.size > 0 else 0.0
    xs = xs_raw - x_shift
    ys = ys_raw - y_shift

    w_sum = float(np.sum(ws))
    wx    = float(np.sum(ws * xs))
    wxx   = float(np.sum(ws * xs * xs))
    wxxx  = float(np.sum(ws * xs * xs * xs))
    wxxxx = float(np.sum(ws * xs * xs * xs * xs))
    wy    = float(np.sum(ws * ys))
    wxy   = float(np.sum(ws * xs * ys))
    wxxy  = float(np.sum(ws * xs * xs * ys))

    # Quadratic via real C++ performPivot
    A = np.array([
        w_sum, wx,   wxx,
        wx,    wxx,  wxxx,
        wxx,   wxxx, wxxxx,
    ], dtype=np.float64)
    b = np.array([wy, wxy, wxxy], dtype=np.float64)
    try:
        coef = _perform_pivot(3, A, b)
        if np.all(np.isfinite(coef)):
            a_p, b_p, c_p = float(coef[0]), float(coef[1]), float(coef[2])
            c_u = c_p
            b_u = b_p - 2.0 * c_p * x_shift
            a_u = y_shift + a_p - b_p * x_shift + c_p * x_shift * x_shift
            return np.array([a_u, b_u, c_u]), "quad"
    except (FloatingPointError, ZeroDivisionError):
        pass

    # Linear fallback via real C++ performPivot too (2x2 with intercept)
    # C++ regressionPivot (Tools.cpp:887-937).
    A_lin = np.array([w_sum, wx, wx, wxx], dtype=np.float64)
    b_lin = np.array([wy, wxy], dtype=np.float64)
    try:
        coef_lin = _perform_pivot(2, A_lin, b_lin)
        if np.all(np.isfinite(coef_lin)):
            b_p, m_p = float(coef_lin[0]), float(coef_lin[1])
            # bCoef = bPrime + (yShift - mPrime * xShift); mCoef = mPrime
            b_corr = b_p + (y_shift - m_p * x_shift)
            return np.array([b_corr, m_p, 0.0]), "lin"
    except (FloatingPointError, ZeroDivisionError):
        pass

    return None, "none"


def _apply_model(x: float | np.ndarray, coef: np.ndarray) -> float | np.ndarray:
    """Evaluate y = a + b*x + c*x²."""
    return coef[0] + coef[1] * x + coef[2] * (x ** 2)


# -----------------------------------------------------------------------------
# One local-model "Baseline" — wraps the per-anchor reject/return/setLM logic
# -----------------------------------------------------------------------------

@dataclass
class _Baseline:
    forward: bool
    scatter: float
    x_anchor: float
    y_anchor: float
    keep: np.ndarray  # bool
    flux: np.ndarray
    dumps: np.ndarray
    ang_dist: np.ndarray

    def __post_init__(self):
        self.keep = self.keep.copy()


def _reject_points(b: _Baseline) -> None:
    """Iteratively reject the largest positive residual until the RMS of
    positive residuals is <= the noise level (scatter). Anchored fits.
    Direct port of C++ Baseline::rejectPoints."""
    while True:
        if not b.keep.any():
            return
        coef, _ = _fixed_quadratic_fit(
            b.keep, b.dumps, b.ang_dist, b.flux,
            b.x_anchor, b.y_anchor,
        )
        if coef is None:
            return

        model_vals = _apply_model(b.ang_dist, coef)
        deltas = b.flux - model_vals
        idx_kept = np.where(b.keep)[0]
        if idx_kept.size == 0:
            return

        kept_deltas = deltas[idx_kept]
        # Positive residuals contribute to splus
        pos_mask = kept_deltas > 0
        n_pos = int(pos_mask.sum())
        if n_pos == 0:
            return  # No positive residuals — nothing to reject
        splus = float(np.sqrt(np.sum(kept_deltas[pos_mask] ** 2) / n_pos))

        if splus <= b.scatter:
            return

        # Find largest positive deviation
        large_local = int(np.argmax(kept_deltas))
        if kept_deltas[large_local] <= 0:
            return
        large_idx = int(idx_kept[large_local])

        # Halt if too few distinct points left.
        # C++ Baseline::sufficentPointCheck (BackgroundCUDA.cpp:37-74) caps
        # the distinct-X / distinct-Y counts at 4 and then C++ rejectPoints
        # (BackgroundCUDA.cpp:303) keeps rejecting while
        # `checkCount >= holderAVec[1]` where holderAVec[1] == 2 for the
        # quadratic model. So the C++ threshold is min(distinctX, distinctY)
        # >= 2 (i.e. stop when < 2). The earlier port used < 3, which
        # terminated the rejection loop one iteration too early on rising
        # edges of bright sources — leaving partial source-contaminated
        # contributions in bg_data that pulled the global BG model up by
        # ~4 units at the Cyg A peak.
        n_distinct_x = np.unique(b.ang_dist[b.keep]).size
        n_distinct_y = np.unique(b.flux[b.keep]).size
        if min(n_distinct_x, n_distinct_y) < 2:
            return

        b.keep[large_idx] = False


def _return_points(b: _Baseline) -> None:
    """Drop the anchor point, then iteratively add back the
    smallest-residual rejected point as long as splus stays < scatter.
    Direct port of C++ Baseline::returnPoints."""
    # Drop the anchor
    if b.forward:
        b.keep[0] = False
    else:
        b.keep[-1] = False

    last_added: int | None = None
    prev_splus_diff = float("inf")

    while True:
        if not b.keep.any():
            return
        coef, _ = _free_quadratic_fit(
            b.keep, b.dumps, b.ang_dist, b.flux,
        )
        if coef is None:
            return

        kept_idx = np.where(b.keep)[0]
        low, high = int(kept_idx.min()), int(kept_idx.max())
        lo_search = max(low - 1, 0)
        hi_search = min(high + 2, b.keep.size)

        # Compute splus on kept points in the search window.
        # C++ computes splus = sqrt(sum_pos_deltas² / counterplus). When
        # counterplus == 0, that's sqrt(0/0) = NaN, and the SUBSEQUENT
        # `while (splus < scatter)` check fails on NaN — so C++ adds
        # one more sample (because `counterplus == 0` IS the add trigger)
        # and then exits. We must match this exit-after-one-add behavior;
        # the earlier port set splus=0 and kept looping, which produced
        # extra phantom adds at the scan edges and inflated the BG model.
        deltas_window = b.flux[lo_search:hi_search] - _apply_model(b.ang_dist[lo_search:hi_search], coef)
        kept_window = b.keep[lo_search:hi_search]
        pos_kept = kept_window & (deltas_window > 0)
        n_pos = int(pos_kept.sum())
        if n_pos == 0:
            splus = float("nan")
        else:
            splus = float(np.sqrt(np.sum(deltas_window[pos_kept] ** 2) / n_pos))

        # Find smallest-deviation rejected point in window
        rej_window = ~kept_window
        rej_deltas = deltas_window.copy()
        rej_deltas[~rej_window] = np.inf
        small_idx = -1
        if np.isfinite(rej_deltas).any():
            small_local = int(np.argmin(rej_deltas))
            if np.isfinite(rej_deltas[small_local]):
                small_idx = lo_search + small_local

        # C++ enters the add branch iff `(splus < scatter || counterplus == 0)
        # && smallIndex > -1`. NaN < scatter is False, but `n_pos == 0`
        # captures the same intent.
        should_add = ((splus < b.scatter) or (n_pos == 0)) and (small_idx >= 0)
        if should_add:
            b.keep[small_idx] = True
            last_added = small_idx
            prev_splus_diff = abs(b.scatter - splus)
            # After add, C++ rechecks `while (splus < scatter)`. With splus
            # finite-and-below-scatter we continue. With splus == NaN
            # (n_pos == 0 case), the check fails → loop exits.
            if not (splus < b.scatter):
                return
        else:
            # Adding pushed us over scatter — was the previous add better?
            if last_added is not None and prev_splus_diff < abs(b.scatter - splus):
                b.keep[last_added] = False
            return


def _set_local_model(b: _Baseline) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Emit the final local model: for each surviving sample, return
    (index, value, weight) using paper Eq. 3 for the weight.

    Returns three parallel arrays. Empty if no surviving samples.
    """
    if not b.keep.any():
        return np.zeros(0, dtype=int), np.zeros(0), np.zeros(0)

    coef, _ = _free_quadratic_fit(b.keep, b.dumps, b.ang_dist, b.flux)
    if coef is None:
        return np.zeros(0, dtype=int), np.zeros(0), np.zeros(0)

    kept_idx = np.where(b.keep)[0]
    low, high = int(kept_idx.min()), int(kept_idx.max())
    counter = int(b.keep[low:high + 1].sum())

    # Compute dump-weighted moments of the kept positions
    kept_in_range = b.keep[low:high + 1]
    x_range = b.ang_dist[low:high + 1]
    d_range = b.dumps[low:high + 1]
    dump_sum = float(d_range[kept_in_range].sum())
    if dump_sum <= 0:
        return np.zeros(0, dtype=int), np.zeros(0), np.zeros(0)
    mean = float((d_range[kept_in_range] * x_range[kept_in_range]).sum() / dump_sum)
    var = float((d_range[kept_in_range] * (x_range[kept_in_range] - mean) ** 2).sum() / dump_sum)
    kurt4 = float((d_range[kept_in_range] * (x_range[kept_in_range] - mean) ** 4).sum() / dump_sum)
    stdev = np.sqrt(max(var, 0.0))
    kurtosis = max(kurt4, 0.0) ** 0.25

    # Emit model values for each sample in [low, high] (including those
    # rejected mid-range, IF the model value isn't above their flux).
    out_idx, out_val, out_w = [], [], []
    for k in range(low, high + 1):
        model_val = float(_apply_model(b.ang_dist[k], coef))
        if (model_val > b.flux[k]) and (not b.keep[k]):
            # C++ skips: model is above flux AND this sample was rejected
            continue

        if not np.isfinite(model_val):
            continue
        if stdev <= 0 or kurtosis <= 0:
            continue

        if counter == 1:
            weight = 1.0
        else:
            t1 = ((b.ang_dist[k] - mean) / stdev) ** 2
            t2 = ((b.ang_dist[k] - mean) / kurtosis) ** 4
            weight = dump_sum / (1.0 + t1 + t2)

        out_idx.append(k)
        out_val.append(model_val)
        out_w.append(float(weight))

    return (np.array(out_idx, dtype=int),
            np.array(out_val, dtype=np.float64),
            np.array(out_w, dtype=np.float64))


# -----------------------------------------------------------------------------
# Driver: find start/end indices and orchestrate per-scan
# -----------------------------------------------------------------------------

def _find_start_end_indices(ang_dist: np.ndarray, scale: float) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Find forward and backward (start, end) ranges per anchor —
    the indices of all samples within `scale` of each anchor.

    Forward: anchor at `start`, extend until ang_dist[end] - ang_dist[start] ≥ scale
    Backward: anchor at `end`, extend backward until ang_dist[end] - ang_dist[start] ≥ scale
    """
    n = ang_dist.size
    forward = []
    end = 0
    for start in range(n - 2):
        while end < n - 1 and ang_dist[end] - ang_dist[start] < scale:
            end += 1
        end -= 1
        if end == n - 2 and ang_dist[end + 1] - ang_dist[start] < scale:
            end += 1
        if end != start:
            forward.append((start, end))

    backward = []
    start = n - 1
    for end in range(n - 1, 1, -1):
        while ang_dist[end] - ang_dist[start] < scale and start > 0:
            start -= 1
        start += 1
        if start == 1 and ang_dist[end] - ang_dist[start - 1] < scale:
            start -= 1
        if start != end:
            backward.append((start, end))

    return forward, backward


def _per_pixel_rcr_mean(
    vals: np.ndarray, weights: np.ndarray,
) -> float:
    """RCR LS_MODE_DL bulk-rejection weighted mean. Falls back to
    weighted mean on degenerate inputs. Matches C++ setBackground."""
    if vals.size == 0:
        return float("nan")
    if vals.size == 1:
        return float(vals[0])
    # Jitter to avoid degenerate-input ZeroDivisionError in DL fit.
    spread = float(np.max(vals) - np.min(vals))
    if spread <= 0:
        return float(vals[0])
    rng = np.random.default_rng(0)
    jitter = rng.normal(0, max(spread * 1e-12, 1e-15), vals.size)
    try:
        rcr = RCR(RejectionTech.LS_MODE_DL)
        rcr.perform_bulk_rejection(vals + jitter, w=weights)
        mu = float(rcr.result.mu)
        if not np.isfinite(mu):
            return float(np.average(vals, weights=weights))
        return mu
    except (ZeroDivisionError, ValueError, IndexError, RuntimeError, OverflowError, ArithmeticError):
        return float(np.average(vals, weights=weights))


def _subtract_one_scan(scan: Scan, scale_deg: float) -> None:
    """Run BG subtraction on one scan in place.

    Sets `scan.background` to the modeled global background and
    `scan.flux_bg` to the original flux minus that background.
    """
    if scan.size < 5 or scan.ang_dist is None or scan.flux is None:
        scan.background = np.zeros(scan.size)
        scan.flux_bg = scan.flux.copy() if scan.flux is not None else scan.flux_l.copy()
        return

    flux = scan.flux
    ang = scan.ang_dist
    dumps = scan.dumps
    noise = (float(scan.noise_1d[0])
             if scan.noise_1d is not None and np.isfinite(scan.noise_1d[0])
             else float(np.std(flux)) / 5.0)

    # Find the per-anchor windows in forward and backward directions.
    fwd_ranges, bwd_ranges = _find_start_end_indices(ang, scale_deg)

    # Accumulate per-sample contributions from every local model.
    bg_data: list[list[float]] = [[] for _ in range(scan.size)]
    bg_weights: list[list[float]] = [[] for _ in range(scan.size)]

    for (lo, hi) in fwd_ranges:
        if hi <= lo:
            continue
        keep = np.ones(hi - lo + 1, dtype=bool)
        bl = _Baseline(
            forward=True,
            scatter=noise,
            x_anchor=float(ang[lo]),
            y_anchor=float(flux[lo]),
            keep=keep,
            flux=flux[lo:hi + 1].copy(),
            dumps=dumps[lo:hi + 1].copy(),
            ang_dist=ang[lo:hi + 1].copy(),
        )
        _reject_points(bl)
        _return_points(bl)
        idxs, vals, wts = _set_local_model(bl)
        for ii, vv, ww in zip(idxs, vals, wts):
            global_idx = int(lo + ii)
            bg_data[global_idx].append(float(vv))
            bg_weights[global_idx].append(float(ww))

    for (lo, hi) in bwd_ranges:
        if hi <= lo:
            continue
        keep = np.ones(hi - lo + 1, dtype=bool)
        bl = _Baseline(
            forward=False,
            scatter=noise,
            x_anchor=float(ang[hi]),
            y_anchor=float(flux[hi]),
            keep=keep,
            flux=flux[lo:hi + 1].copy(),
            dumps=dumps[lo:hi + 1].copy(),
            ang_dist=ang[lo:hi + 1].copy(),
        )
        _reject_points(bl)
        _return_points(bl)
        idxs, vals, wts = _set_local_model(bl)
        for ii, vv, ww in zip(idxs, vals, wts):
            global_idx = int(lo + ii)
            bg_data[global_idx].append(float(vv))
            bg_weights[global_idx].append(float(ww))

    # Global background per sample: RCR-bulk-mean of accumulated values.
    SENTINEL = 999999.0
    bg = np.full(scan.size, SENTINEL)
    for i in range(scan.size):
        if not bg_data[i]:
            continue
        v = np.array(bg_data[i], dtype=np.float64)
        w = np.array(bg_weights[i], dtype=np.float64)
        # Drop NaN/Inf
        ok = np.isfinite(v) & np.isfinite(w) & (w > 0)
        if not ok.any():
            continue
        bg[i] = _per_pixel_rcr_mean(v[ok], w[ok])

    # Linear interpolation of remaining sentinel pixels.
    valid = bg != SENTINEL
    if valid.any():
        idx_all = np.arange(scan.size)
        bg = np.interp(idx_all, idx_all[valid], bg[valid])
    else:
        bg = np.zeros(scan.size)

    scan.background = bg
    scan.flux_bg = scan.flux - bg


# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------

def subtract(survey: Survey, scale_bw: float = 6.0) -> Survey:
    """Run 1D background subtraction on every scan in the Survey.

    Args:
        survey: must have calibrated `flux` and `noise_1d` populated.
        scale_bw: background-subtraction scale in beamwidths.
            Paper §3.3 / Table 1 recommends 6-7 for 20-meter L-band.
    """
    if survey.psf_fwhm <= 0:
        # No beam size → no scale to use. Pass through unchanged.
        for scan in survey.scans:
            scan.flux_bg = scan.flux.copy() if scan.flux is not None else scan.flux_l.copy()
            scan.background = np.zeros(scan.size)
        return survey

    scale_deg = scale_bw * survey.psf_fwhm
    for scan in survey.scans:
        _subtract_one_scan(scan, scale_deg)

    return survey
