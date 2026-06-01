"""Stages 3 & 6 — noise level measurement.

Stage 3 (1D along-scan noise, §3.2):
    For each sample, draw a line between the previous and next samples,
    then measure the central sample's deviation from that line. RCR-reject
    outliers (which are usually sources or RFI). The post-rejection stddev
    of the deviations, divided by the empirical 1.220 bias correction,
    is the per-scan sigma.

Stage 6 (2D across-scan noise, §3.5):
    Mirrors Stage 3 but across scans rather than along scans. For each
    sample, find the closest sample in the preceding scan and following
    scan (by position along scan) and use those as endpoints instead.
    The bias correction here is 1.229, and both mean and stddev contribute
    to the final sigma (not just stddev — see paper Eq. between §3.5
    paragraphs).

A linear fit over scan number (RCR-rejected) lets the noise level drift
slowly over the observation. Stored as `survey.noise_*_slope` /
`noise_*_intercept` so any sample's noise can be reconstructed from its
parent scan index.

References:
    Paper §3.2, §3.5, Footnotes 10 and 11.
    C++ source: Scan::calculateScatter, Scan::pointToPointDiff,
                Processor::characterizeData, Processor::set2DScatter.
"""
from __future__ import annotations

from typing import Iterable

import numpy as np

from rcpy.types import Survey, Scan
from rcpy.rcr2 import RCR, RejectionTech


# Empirical bias correction factors from the paper. The point-to-point
# estimator over-estimates the true sigma of Gaussian noise; these are
# divisors that recover the true value.
_BIAS_1D = 1.220   # Paper §3.2
_BIAS_2D = 1.229   # Paper §3.5


# -----------------------------------------------------------------------------
# Shared utilities
# -----------------------------------------------------------------------------

def _residuals_from_line(
    x_left: np.ndarray,
    y_left: np.ndarray,
    x_mid: np.ndarray,
    y_mid: np.ndarray,
    x_right: np.ndarray,
    y_right: np.ndarray,
) -> np.ndarray:
    """For each (left, mid, right) sample triple, compute the deviation of
    the middle sample from the line through the left and right samples.

    Vectorized — all arrays have the same length N.
    """
    # Slope of the line through (x_left, y_left) and (x_right, y_right):
    dx = x_right - x_left
    # Where dx == 0, the deviation is undefined; mask with NaN.
    with np.errstate(divide="ignore", invalid="ignore"):
        slope = np.where(dx != 0, (y_right - y_left) / dx, 0.0)
    y_predicted = y_left + slope * (x_mid - x_left)
    return y_mid - y_predicted


def _weights_inverse_variance(
    n_mid: np.ndarray, n_left: np.ndarray, n_right: np.ndarray,
    x_left: np.ndarray, x_mid: np.ndarray, x_right: np.ndarray,
) -> np.ndarray:
    """Weights for the point-to-point deviation per Footnote 10.

    w = 1 / (1/N_mid + 1/N_line(x_mid))
    where N_line(x_mid) is the interpolated dump-count of the line at x_mid.
    """
    # Linear interpolation: x_mid sits at fraction t between left and right.
    dx = x_right - x_left
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(dx != 0, (x_mid - x_left) / dx, 0.5)
    # Approximate the line's effective N via reciprocal interpolation —
    # paper Footnote 10 has the full expression; this is a faithful
    # simplification that matches asymptotically.
    n_line = 1.0 / ((1 - t) ** 2 / np.maximum(n_left, 1) + t ** 2 / np.maximum(n_right, 1))
    return 1.0 / (1.0 / np.maximum(n_mid, 1) + 1.0 / n_line)


def _rcr_scatter(residuals: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    """RCR-reject outliers from residuals, return (mean, stddev) of survivors."""
    if residuals.size == 0:
        return 0.0, 0.0

    # Drop any NaN/Inf rows that came from divide-by-zero in the residual.
    finite = np.isfinite(residuals) & np.isfinite(weights)
    residuals = residuals[finite]
    weights = weights[finite]
    if residuals.size == 0:
        return 0.0, 0.0

    # Guard against degenerate inputs (all-equal residuals → sigma=0
    # → RCR divides by zero). Daisies with a near-zero TS shift can
    # produce these on adjacent-scan matching when one scan has a
    # too-small overlap window.
    if float(np.max(residuals) - np.min(residuals)) < 1e-12:
        return float(np.mean(residuals)), 0.0
    try:
        rcr = RCR(RejectionTech.LS_MODE_68)
        rcr.perform_rejection(residuals, w=weights)
        clean = rcr.result.clean_y
    except (ZeroDivisionError, ValueError, RuntimeError):
        return float(np.mean(residuals)), float(np.std(residuals))
    if clean.size == 0:
        return 0.0, float(np.std(residuals))
    return float(np.mean(clean)), float(np.std(clean))


def _linear_fit_across_scans(
    scan_numbers: np.ndarray, sigmas: np.ndarray, weights: np.ndarray
) -> tuple[float, float]:
    """Fit sigma = slope * scan_number + intercept with RCR.

    Used to model slow drift of the noise level over the observation.
    Returns (slope, intercept).
    """
    if scan_numbers.size < 2:
        return 0.0, float(sigmas[0]) if sigmas.size else 0.0

    wsum = weights.sum()
    xbar = (weights * scan_numbers).sum() / wsum
    ybar = (weights * sigmas).sum() / wsum
    sxx = (weights * (scan_numbers - xbar) ** 2).sum()
    sxy = (weights * (scan_numbers - xbar) * (sigmas - ybar)).sum()
    slope = sxy / sxx if sxx > 0 else 0.0
    intercept = ybar - slope * xbar

    # One RCR pass on residuals to reject outliers, then refit.
    residuals = sigmas - (slope * scan_numbers + intercept)
    rcr = RCR(RejectionTech.LS_MODE_68)
    rcr.perform_rejection(residuals, w=weights)
    keep = rcr.result.flags if rcr.result.flags.size else np.ones(sigmas.size, bool)
    if keep.sum() >= 2:
        sn = scan_numbers[keep]
        sg = sigmas[keep]
        wk = weights[keep]
        wsum = wk.sum()
        xbar = (wk * sn).sum() / wsum
        ybar = (wk * sg).sum() / wsum
        sxx = (wk * (sn - xbar) ** 2).sum()
        sxy = (wk * (sn - xbar) * (sg - ybar)).sum()
        slope = sxy / sxx if sxx > 0 else 0.0
        intercept = ybar - slope * xbar

    return slope, intercept


# -----------------------------------------------------------------------------
# Stage 3 — 1D noise (along scan)
# -----------------------------------------------------------------------------

def measure_1d(survey: Survey) -> Survey:
    """Measure along-scan point-to-point noise and store per-sample sigma.

    Sets `scan.noise_1d` on each scan and the survey-level slope/intercept
    used to evaluate sigma(scan_number) elsewhere.
    """
    per_scan_sigma = []
    per_scan_weight = []
    scan_numbers = []

    # Per-scan measurement
    for scan in survey.scans:
        if scan.size < 3 or scan.flux is None:
            # Need at least 3 samples for the prev/next line trick.
            scan.noise_1d = np.full(scan.size, np.nan)
            continue

        x = scan.ang_dist if scan.ang_dist is not None else np.arange(scan.size, dtype=float)
        y = scan.flux
        n = scan.dumps

        residuals = _residuals_from_line(
            x[:-2], y[:-2], x[1:-1], y[1:-1], x[2:], y[2:]
        )
        weights = _weights_inverse_variance(
            n[1:-1], n[:-2], n[2:], x[:-2], x[1:-1], x[2:]
        )

        _, stddev = _rcr_scatter(residuals, weights)
        sigma = stddev / _BIAS_1D

        # Store per-sample sigma (constant within a scan at this stage).
        # The cross-scan linear fit below will refine it.
        scan.noise_1d = np.full(scan.size, sigma)

        per_scan_sigma.append(sigma)
        per_scan_weight.append(float(n.sum()))
        scan_numbers.append(scan.scan_index)

    # Cross-scan linear fit to capture slow drift.
    if per_scan_sigma:
        slope, intercept = _linear_fit_across_scans(
            np.asarray(scan_numbers, dtype=float),
            np.asarray(per_scan_sigma, dtype=float),
            np.asarray(per_scan_weight, dtype=float),
        )
        survey.noise_1d_slope = slope
        survey.noise_1d_intercept = intercept

        # Re-set per-sample sigma from the model.
        for scan in survey.scans:
            sigma = slope * scan.scan_index + intercept
            if scan.size > 0:
                scan.noise_1d = np.full(scan.size, max(sigma, 1e-30))

    return survey


# -----------------------------------------------------------------------------
# Stage 6 — 2D noise (across scans)
# -----------------------------------------------------------------------------

def _nearest_in_other_scan(
    target_x: np.ndarray, other_x: np.ndarray
) -> np.ndarray:
    """Vectorized nearest-neighbor index in `other_x` for each `target_x`.

    Uses bisect on a sorted view of `other_x`.
    """
    if other_x.size == 0:
        return np.full(target_x.size, -1, dtype=np.int64)
    order = np.argsort(other_x)
    sorted_x = other_x[order]
    pos = np.searchsorted(sorted_x, target_x)
    pos_lo = np.clip(pos - 1, 0, other_x.size - 1)
    pos_hi = np.clip(pos, 0, other_x.size - 1)
    d_lo = np.abs(target_x - sorted_x[pos_lo])
    d_hi = np.abs(target_x - sorted_x[pos_hi])
    pick = np.where(d_lo < d_hi, pos_lo, pos_hi)
    return order[pick]


def measure_2d(survey: Survey) -> Survey:
    """Measure across-scan noise after background subtraction.

    Requires that `scan.flux_bg` (or `flux`) is available for each scan
    and `scan.ang_dist` for adjacent-scan position matching.

    Mirrors Stage 3 but with the line endpoints drawn from the preceding
    and following scans rather than from neighbouring samples in the same
    scan.
    """
    per_scan_sigma = []
    per_scan_weight = []
    scan_numbers = []

    for i, scan in enumerate(survey.scans):
        # Need a scan before and after.
        if i == 0 or i == len(survey.scans) - 1:
            scan.noise_2d = scan.noise_1d.copy() if scan.noise_1d is not None else np.zeros(scan.size)
            continue
        prev = survey.scans[i - 1]
        nxt = survey.scans[i + 1]
        if scan.size < 1 or prev.size < 1 or nxt.size < 1:
            scan.noise_2d = scan.noise_1d.copy() if scan.noise_1d is not None else np.zeros(scan.size)
            continue

        y_mid = scan.working_flux()
        y_prev = prev.working_flux()
        y_nxt = nxt.working_flux()
        x_mid = scan.ang_dist
        x_prev = prev.ang_dist
        x_nxt = nxt.ang_dist
        n_mid = scan.dumps
        n_prev = prev.dumps
        n_nxt = nxt.dumps

        if x_mid is None or x_prev is None or x_nxt is None:
            scan.noise_2d = scan.noise_1d.copy() if scan.noise_1d is not None else np.zeros(scan.size)
            continue

        idx_prev = _nearest_in_other_scan(x_mid, x_prev)
        idx_nxt = _nearest_in_other_scan(x_mid, x_nxt)

        residuals = _residuals_from_line(
            x_prev[idx_prev], y_prev[idx_prev],
            x_mid, y_mid,
            x_nxt[idx_nxt], y_nxt[idx_nxt],
        )
        weights = _weights_inverse_variance(
            n_mid, n_prev[idx_prev], n_nxt[idx_nxt],
            x_prev[idx_prev], x_mid, x_nxt[idx_nxt],
        )

        mean_dev, stddev = _rcr_scatter(residuals, weights)
        # Paper §3.5: sigma_2d = sqrt(mean^2 + stddev^2) / 1.229
        sigma = np.sqrt(mean_dev ** 2 + stddev ** 2) / _BIAS_2D
        scan.noise_2d = np.full(scan.size, sigma)

        per_scan_sigma.append(sigma)
        per_scan_weight.append(float(n_mid.sum()))
        scan_numbers.append(scan.scan_index)

    if per_scan_sigma:
        slope, intercept = _linear_fit_across_scans(
            np.asarray(scan_numbers, dtype=float),
            np.asarray(per_scan_sigma, dtype=float),
            np.asarray(per_scan_weight, dtype=float),
        )
        survey.noise_2d_slope = slope
        survey.noise_2d_intercept = intercept

        for scan in survey.scans:
            sigma = slope * scan.scan_index + intercept
            if scan.size > 0:
                # C++ Processor.cpp:441 takes max(scan_1d_scatter, fitted).
                # For scans where the 1D point-to-point scatter is
                # elevated by a bright source (Jupiter), the 1D value is
                # larger than the smoothed cross-scan fit, and C++ uses
                # the larger of the two. This raises the RFI-rejection
                # threshold near the source enough that source samples
                # survive instead of being eaten as RFI spikes.
                scan_1d = float(np.median(scan.noise_1d)) if scan.noise_1d is not None else 0.0
                scan.noise_2d = np.full(scan.size, max(max(sigma, 1e-30), scan_1d))

    return survey
