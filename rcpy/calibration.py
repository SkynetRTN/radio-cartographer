"""Stage 2 — gain calibration (Paper §3.1).

The receiver has a noise diode that can be switched on and off. By measuring
the signal level with the diode on, then off, then taking the difference,
we get a reference scale: dividing every flux sample by this delta puts
everything in "noise diode units" (1.0 = brightness of the noise diode).

Subtleties:

  - The diode-on and diode-off intervals can each contain contaminated
    samples (RFI, samples caught mid-transition, background drift). We
    fit a LINE through each interval rather than just taking the mean,
    using RCR to reject outliers (Paper Footnote 8, Fig. 8).
  - The C++ uses a 3-tech cascade: ES_MODE_DL -> LS_MODE_68 ->
    SS_MEDIAN_DL. The cascade order matters for parity.
  - Each polarization channel is calibrated independently. The composite
    flux is L + R after separate calibration (not (L+R) / combined_delta).
  - cal_method controls which calibration is applied to which samples:
      PRE          — use only the first interval's delta
      POST         — use only the last interval's delta
      INTERPOLATED — linearly interpolate delta(t) between intervals
      NONE         — no calibration (just copy raw flux)

References:
    Paper §3.1, Eq. 2, Footnote 8.
    C++ source: Survey::gainCalibration in src/Survey.cpp.
"""
from __future__ import annotations

import numpy as np

from rcpy.types import Survey, Scan, CalMethod, Channel
from rcpy.rcr2 import RCR, RejectionTech


# -----------------------------------------------------------------------------
# Calibration-interval detection
# -----------------------------------------------------------------------------

def _find_cal_intervals(cal_state: np.ndarray) -> list[tuple[int, int, int]]:
    """Return list of (start_idx, end_idx_exclusive, state) for each
    contiguous run of identical cal_state values. State 0 = off, 1 = on,
    other values may exist but are ignored for delta computation."""
    if cal_state.size == 0:
        return []
    edges = np.where(np.diff(cal_state) != 0)[0] + 1
    starts = np.concatenate(([0], edges))
    ends = np.concatenate((edges, [cal_state.size]))
    return [(int(s), int(e), int(cal_state[s])) for s, e in zip(starts, ends)]


def _fit_line_cpp_style(
    x: np.ndarray, y: np.ndarray, w: np.ndarray
) -> tuple[float, float, float]:
    """C++-parity line fit for cal lines: matches
    `LinearModel(dump, time, flux) + RCR(LS_MODE_DL).performBulkRejection`
    in Survey.cpp:1088-1163.

    Iteratively (1) fits a weighted least squares line, (2) bulk-rejects
    residual outliers using LS_MODE_DL, repeating until no further
    rejection.

    Returns (slope, intercept_at_xbar, xbar) so the caller can evaluate
    via y(t) = slope*(t - xbar) + intercept_at_xbar — matches the
    `m*(avgT - xBar) + b` form C++ uses to compute deltaStart/deltaEnd.
    """
    def _fit(xv, yv, wv):
        wsum = wv.sum()
        if wsum <= 0:
            return 0.0, 0.0, 0.0
        xbar = (wv * xv).sum() / wsum
        ybar = (wv * yv).sum() / wsum
        sxx = (wv * (xv - xbar) ** 2).sum()
        sxy = (wv * (xv - xbar) * (yv - ybar)).sum()
        slope = sxy / sxx if sxx > 0 else 0.0
        return slope, ybar, xbar

    if x.size < 2:
        return 0.0, float(np.average(y, weights=w)) if x.size else 0.0, 0.0

    slope, b, xbar = _fit(x, y, w)
    max_iters = 8
    for _ in range(max_iters):
        residuals = y - (slope * (x - xbar) + b)
        resid_spread = float(np.max(residuals) - np.min(residuals))
        ref_scale = max(abs(float(np.mean(y))), 1.0)
        if resid_spread < 1e-9 * ref_scale:
            break
        try:
            rcr = RCR(RejectionTech.LS_MODE_DL)
            rcr.perform_bulk_rejection(residuals, w=w)
            keep = rcr.result.flags if rcr.result.flags.size else np.ones(x.size, bool)
        except (ZeroDivisionError, ValueError, RuntimeError):
            break
        if keep.sum() < 3 or keep.sum() == x.size:
            # No further rejection (or too few survivors).
            break
        x = x[keep]
        y = y[keep]
        w = w[keep]
        slope, b, xbar = _fit(x, y, w)

    return slope, b, xbar


def _fit_line_with_rcr(
    x: np.ndarray, y: np.ndarray, w: np.ndarray
) -> tuple[float, float]:
    """Fit y = a*x + b via weighted least squares with RCR outlier rejection
    on residuals. Returns (slope, intercept).

    Cascade: bulk + individual rejection per Paper Footnote 8.
    """
    if x.size < 2:
        return 0.0, float(np.average(y, weights=w)) if x.size else 0.0

    # Initial weighted least-squares fit
    def _fit(xv, yv, wv):
        wsum = wv.sum()
        if wsum <= 0:
            return 0.0, 0.0
        xbar = (wv * xv).sum() / wsum
        ybar = (wv * yv).sum() / wsum
        sxx = (wv * (xv - xbar) ** 2).sum()
        sxy = (wv * (xv - xbar) * (yv - ybar)).sum()
        slope = sxy / sxx if sxx > 0 else 0.0
        intercept = ybar - slope * xbar
        return slope, intercept

    slope, intercept = _fit(x, y, w)

    # Reject residual outliers via RCR cascade.
    for tech in (RejectionTech.ES_MODE_DL,
                 RejectionTech.LS_MODE_68,
                 RejectionTech.SS_MEDIAN_DL):
        residuals = y - (slope * x + intercept)
        # If residuals are essentially constant (e.g. very clean cal
        # block where the line fit is near-perfect), RCR's `max / sigma`
        # divides by ~0 and crashes. Skip RCR in that case — the
        # initial WLS already nailed the line.
        resid_spread = float(np.max(residuals) - np.min(residuals))
        ref_scale = max(abs(float(np.mean(y))), 1.0)
        if resid_spread < 1e-9 * ref_scale:
            break
        try:
            rcr = RCR(tech)
            rcr.perform_rejection(residuals, w=w)
            keep = rcr.result.flags if rcr.result.flags.size else np.ones(x.size, bool)
        except (ZeroDivisionError, ValueError, RuntimeError):
            # RCR can fail on degenerate inputs (sigma=0). Keep current fit.
            break
        if keep.sum() < 2:
            break
        x = x[keep]
        y = y[keep]
        w = w[keep]
        slope, intercept = _fit(x, y, w)

    return slope, intercept


def _measure_delta(
    time: np.ndarray,
    flux: np.ndarray,
    dumps: np.ndarray,
    cal_state: np.ndarray,
) -> tuple[float, float] | None:
    """Measure one gain delta from a (typically short) calibration block.

    Fits separate lines through the diode-on and diode-off samples,
    evaluates each line at the dump-weighted mean time of all
    non-rejected samples in the block, and returns
    (delta, evaluation_time).

    Returns None if the block lacks both diode states.
    """
    on = cal_state == 1
    off = cal_state == 0

    if not on.any() or not off.any():
        return None

    s_on, i_on = _fit_line_with_rcr(time[on], flux[on], dumps[on])
    s_off, i_off = _fit_line_with_rcr(time[off], flux[off], dumps[off])

    # Dump-weighted mean time of all non-rejected samples in the block.
    # (Approximation: we use all samples here. C++ uses post-rejection set;
    # if bit-parity is needed, track which samples survived per fit.)
    t_eval = float(np.average(time, weights=dumps))

    y_on = s_on * t_eval + i_on
    y_off = s_off * t_eval + i_off
    delta = y_on - y_off

    return delta, t_eval


# -----------------------------------------------------------------------------
# Public stage entry point
# -----------------------------------------------------------------------------

def gain(survey: Survey) -> Survey:
    """Run gain calibration on a Survey.

    Uses the calibration samples that were filtered off in io.read_sdfits
    (SWPVALID=0 — diode-on cal blocks plus turnaround diode-off baseline
    samples). Measures one delta per cal block per channel, then applies
    the configured cal_method to populate each science scan's `flux`.

    Modifies the Survey in place; also returns it.
    """
    if survey.cal_method is CalMethod.NONE or survey.cal_time is None or survey.cal_time.size == 0:
        # No calibration: copy raw flux into the working `flux` field.
        for s in survey.scans:
            s.flux = _pick_channel(s, survey.channel).astype(np.float64).copy()
            s.flux_composite = (0.5 * (s.flux_l + s.flux_r)).astype(np.float64)
        survey.gain_delta_l = np.array([1.0])
        survey.gain_delta_r = np.array([1.0])
        survey.gain_delta_times = np.array([0.0])
        return survey

    # Diode-on samples: CALSTATE == 1
    # Diode-off baseline samples: CALSTATE == 0 (Skynet uses the
    # SWPVALID=0 turnaround samples scattered through the obs as the
    # off baseline, per pyrc/gain_calibration.py).
    cal_time = survey.cal_time
    cal_state = survey.cal_state
    cal_dumps = survey.cal_dumps
    cal_l = survey.cal_flux_l
    cal_r = survey.cal_flux_r

    on = cal_state == 1
    off = cal_state == 0

    if not on.any() or not off.any():
        for s in survey.scans:
            s.flux = _pick_channel(s, survey.channel).astype(np.float64).copy()
            s.flux_composite = (0.5 * (s.flux_l + s.flux_r)).astype(np.float64)
        survey.gain_delta_l = np.array([1.0])
        survey.gain_delta_r = np.array([1.0])
        survey.gain_delta_times = np.array([0.0])
        return survey

    # Find on-runs (each diode-on block has its own delta evaluated at
    # the block's mean time).
    on_runs = _runs(on)
    n_runs = len(on_runs)

    # For each cal-on run, the OFF samples we use to compute the delta
    # come from the SAME end of the observation as that run — NOT from
    # a global off-line fit across all off samples. This matches C++
    # Survey::gainCalibration (Survey.cpp:1082-1167) which builds
    # `lowFluxArrayStart` / `lowFluxArrayEnd` from the FIRST / LAST
    # scan's flux samples respectively, then fits separate off lines
    # `mOffStart, bOffStart` / `mOffEnd, bOffEnd` for each.
    #
    # The earlier (buggy) global off-line fit failed on cyga: the off
    # baseline rises non-linearly over the obs (RFI episodes, weather),
    # the linear fit went the WRONG direction, and predicted off at the
    # late cal block 1700 counts below the actual local median. That
    # made delta_END look 26% larger than delta_START even though the
    # true gain barely drifted, and INTERPOLATED cal then spread that
    # phantom drift across all later scans → flux baseline dropped 30%
    # by the end of the observation. With per-block local off windows
    # the two cyga deltas are 6086 vs 5948 (2% apart), matching the
    # near-constant true gain.
    # Estimate the "first scan" duration. For Skynet 20-m data, a
    # whole scan is roughly (total cal-time span) / number-of-scans.
    # OFF samples within this duration of a cal-on block are the ones
    # C++ groups into `lowFluxArrayStart` / `lowFluxArrayEnd` via the
    # `flux[0]` / `flux[last]` per-scan separation. Using a longer
    # window pulls in OFF samples from many later scans, which biases
    # the start-delta or end-delta by the cumulative off-baseline
    # drift (cyga: off rises ~1500 counts over ~9 min, biases the
    # END delta high by 1700 counts → 26% phantom gain drift).
    cal_time_span = float(cal_time.max() - cal_time.min())
    if len(survey.scans) > 0:
        scan_duration_s = cal_time_span / max(len(survey.scans), 1)
    else:
        scan_duration_s = 20.0
    # Use 1.5 scans-worth — wide enough that small per-scan timing
    # variations don't drop the OFF window entirely, narrow enough
    # that it stays "local" to the cal-on block.
    local_off_half_window_s = 1.5 * scan_duration_s

    def _local_off_mask(t_eval, t_run_lo, t_run_hi):
        """OFF samples within ~1.5 scan-durations of this cal-on block
        in cal_time. Matches C++'s `flux[0]` / `flux[last]` selection
        which only grabs OFF samples from the FIRST / LAST scan
        respectively, not from anywhere in the observation."""
        # Window: [t_run_lo - window, t_run_hi + window]
        t_min = t_run_lo - local_off_half_window_s
        t_max = t_run_hi + local_off_half_window_s
        return off & (cal_time >= t_min) & (cal_time <= t_max)

    deltas_l, deltas_r, times_eval = [], [], []
    for run_lo, run_hi in on_runs:
        # Diode-on line fit (within run) — line over time to allow drift
        s_on_l, i_on_l = _fit_line_with_rcr(
            cal_time[run_lo:run_hi], cal_l[run_lo:run_hi], cal_dumps[run_lo:run_hi]
        )
        s_on_r, i_on_r = _fit_line_with_rcr(
            cal_time[run_lo:run_hi], cal_r[run_lo:run_hi], cal_dumps[run_lo:run_hi]
        )

        # Evaluate ON line at the run's mean time.
        t_eval = float(np.average(cal_time[run_lo:run_hi], weights=cal_dumps[run_lo:run_hi]))

        # OFF line fit on samples LOCAL to this cal-on block (port of
        # C++ start/end-segregated off-line fits). If too few local
        # off samples for a robust line fit, widen the window
        # progressively until we have at least 5, falling back to all
        # off samples as a last resort.
        local_off = _local_off_mask(t_eval, cal_time[run_lo], cal_time[run_hi - 1])
        widen_factor = 1.0
        # Need at least ~20 OFF samples for a stable WLS line fit;
        # below that, the slope is dominated by the few samples present
        # and predicts wildly at the cal-block time.
        while local_off.sum() < 20 and widen_factor <= 8.0:
            widen_factor *= 2.0
            t_min = cal_time[run_lo] - widen_factor * local_off_half_window_s
            t_max = cal_time[run_hi - 1] + widen_factor * local_off_half_window_s
            local_off = off & (cal_time >= t_min) & (cal_time <= t_max)
        if local_off.sum() >= 2:
            s_off_l, i_off_l = _fit_line_with_rcr(
                cal_time[local_off], cal_l[local_off], cal_dumps[local_off])
            s_off_r, i_off_r = _fit_line_with_rcr(
                cal_time[local_off], cal_r[local_off], cal_dumps[local_off])
        else:
            # Last-resort fallback: global off fit.
            s_off_l, i_off_l = _fit_line_with_rcr(cal_time[off], cal_l[off], cal_dumps[off])
            s_off_r, i_off_r = _fit_line_with_rcr(cal_time[off], cal_r[off], cal_dumps[off])

        delta_l = (s_on_l * t_eval + i_on_l) - (s_off_l * t_eval + i_off_l)
        delta_r = (s_on_r * t_eval + i_on_r) - (s_off_r * t_eval + i_off_r)
        deltas_l.append(delta_l)
        deltas_r.append(delta_r)
        times_eval.append(t_eval)

    survey.gain_delta_l = np.array(deltas_l)
    survey.gain_delta_r = np.array(deltas_r)
    survey.gain_delta_times = np.array(times_eval)

    # Apply per cal_method.
    for s in survey.scans:
        d_l = _delta_at(s.time, deltas_l, times_eval, survey.cal_method)
        d_r = _delta_at(s.time, deltas_r, times_eval, survey.cal_method)
        l_cal = s.flux_l / d_l
        r_cal = s.flux_r / d_r
        # C++ Survey::dataProc (Survey.cpp:440): fluxComp = 0.5 * (L + R)
        # i.e. AVERAGE of the two polarizations, not sum.
        s.flux_composite = 0.5 * (l_cal + r_cal)

        if survey.channel is Channel.LEFT:
            s.flux = l_cal.copy()
        elif survey.channel is Channel.RIGHT:
            s.flux = r_cal.copy()
        else:
            s.flux = s.flux_composite.copy()

    return survey




def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Return [(start, end_exclusive), ...] for contiguous True runs."""
    if mask.size == 0:
        return []
    edges = np.where(np.diff(mask.astype(np.int8)) != 0)[0] + 1
    boundaries = np.concatenate(([0], edges, [mask.size]))
    out = []
    for i in range(boundaries.size - 1):
        s, e = int(boundaries[i]), int(boundaries[i + 1])
        if mask[s]:
            out.append((s, e))
    return out


def _pick_channel(scan: Scan, channel: Channel) -> np.ndarray:
    if channel is Channel.LEFT:
        return scan.flux_l
    if channel is Channel.RIGHT:
        return scan.flux_r
    # COMPOSITE: pre-calibration, use AVERAGE of L and R per C++
    # Survey::dataProc fluxComp = 0.5 * (L+R) (Survey.cpp:440).
    return 0.5 * (scan.flux_l + scan.flux_r)


def _identify_cal_blocks(time: np.ndarray, cal_state: np.ndarray) -> list[tuple[int, int]]:
    """Identify calibration blocks suitable for measuring a noise-diode
    delta.

    Skynet 20-m data does NOT have alternating cal_state. Instead, the
    diode is held ON for a contiguous run (typically 30 samples ~3s) at
    the start of the observation, OFF for the entire scan, then ON again
    for another run at the end. We measure delta = mean(on at start) -
    mean(off near start), and similarly for the end.

    Returns a list of (start_idx, end_idx_exclusive) ranges where each
    range contains BOTH a contiguous on-run and a temporally-adjacent
    off-stretch suitable for the difference measurement.

    Strategy:
      1. Find each contiguous on-run.
      2. Extend the block from the on-run boundary into adjacent off
         samples (up to a cap of ~10x the on-run length), to give the
         off-mean a comparable number of samples.
      3. Return one block per on-run.
    """
    if cal_state.size == 0:
        return []

    # Find on-runs (contiguous stretches of cal_state == 1).
    on = cal_state == 1
    if not on.any():
        return []

    # Identify run boundaries.
    edges = np.where(np.diff(on.astype(np.int8)) != 0)[0] + 1
    boundaries = np.concatenate(([0], edges, [on.size]))
    on_runs: list[tuple[int, int]] = []
    for i in range(boundaries.size - 1):
        s, e = int(boundaries[i]), int(boundaries[i + 1])
        if on[s]:
            on_runs.append((s, e))

    if not on_runs:
        return []

    # For each on-run, build a block that extends into adjacent off
    # samples. Cap the off extension at 10x the on-run length so we
    # don't pull in too much sky drift.
    blocks = []
    for s, e in on_runs:
        run_len = e - s
        extend = 10 * run_len
        # Find adjacent off region; prefer immediately after the on-run,
        # falling back to immediately before.
        off_after_end = min(e + extend, cal_state.size)
        off_before_start = max(s - extend, 0)

        # Pick the side with more off samples.
        n_after = int((cal_state[e:off_after_end] == 0).sum())
        n_before = int((cal_state[off_before_start:s] == 0).sum())

        if n_after >= n_before and n_after > 0:
            blocks.append((s, off_after_end))
        elif n_before > 0:
            blocks.append((off_before_start, e))
        # else: no off samples adjacent — skip this on-run.

    return blocks


def _delta_at(
    times: np.ndarray,
    deltas: list[float],
    times_eval: list[float],
    method: CalMethod,
) -> np.ndarray:
    """Evaluate the calibration delta at each sample time per cal_method."""
    arr = np.asarray(deltas)
    teval = np.asarray(times_eval)
    if method is CalMethod.PRE or arr.size == 1:
        return np.full_like(times, arr[0], dtype=np.float64)
    if method is CalMethod.POST:
        return np.full_like(times, arr[-1], dtype=np.float64)
    if method is CalMethod.CONSTANT:
        # Mean of start and end deltas. Use when calibration fits show
        # spurious drift (e.g. source contamination of cal block ON/OFF
        # line fits). For cyga, INTERPOLATED produces a ~14% phantom
        # drift; CONSTANT collapses the resulting Dec-direction gradient
        # in the raw render from ~1.1 Jy to ~0.1 Jy.
        return np.full_like(times, float(arr.mean()), dtype=np.float64)
    # INTERPOLATED — linear interpolation. Eq. 2 of the paper.
    return np.interp(times, teval, arr).astype(np.float64)
