"""Stage 1 — coordinate preparation.

The pipeline operates in a locally-flat Cartesian frame centered on the
observation. The projection used is sinusoidal / Sanson-Flamsteed (equal
area), where RA gets scaled by cos(Dec):

    x' = (RA - center_RA) * cos(Dec)
    y' = (Dec - center_Dec)

This preserves areas (so photometry doesn't bias with declination) and
keeps horizontal/vertical distances through the center distortion-free.
The center is chosen as the RCR-rejected median of all sample positions.

References:
    Paper §3.7 (sinusoidal projection block).
    C++ source: Scan::cosDecTransform / Survey coord init.
"""
from __future__ import annotations

import math
import numpy as np

from rcpy.types import Survey, Scan

_DEG2RAD = math.pi / 180.0


def _weighted_median(values: np.ndarray, weights: np.ndarray | None = None) -> float:
    """Plain weighted median. We do NOT need RCR here: the observation
    footprint (including its corners) is the thing we want centered on.
    The C++ uses RCR-median here, which can crash on degenerate inputs;
    for center-finding the difference is < 1 pixel in practice."""
    if weights is None:
        return float(np.median(values))
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cw = np.cumsum(w)
    half = cw[-1] / 2.0
    idx = int(np.searchsorted(cw, half))
    idx = min(idx, v.size - 1)
    return float(v[idx])


def project(survey: Survey, per_sample_cos: bool = True) -> Survey:
    """Apply the sinusoidal cos-Dec projection in place on every scan.

    Sets `scan.ra_proj`, `scan.dec_proj`, and `scan.ang_dist` for each scan,
    plus the Survey's `partition_sss` with center coords.

    `per_sample_cos` selects the cos-Dec variant:
      True  — Sanson-Flamsteed proper: ra_proj = (ra - center) * cos(dec)
              per sample. Equal-area, tapered shape (narrower at high
              |Dec|). Use for the MAIN flux output to match the
              reference's main HDU shape.
      False — cylindrical approximation: ra_proj = (ra - center) *
              cos(center_dec). Constant scaling, rectangular shape.
              Matches the reference's RAW HDU which appears to be
              rendered without per-sample cos-Dec — its polygon row
              widths are constant in Dec, ours (with per-sample SFL)
              tapered top to bottom.

    Returns the same Survey for chaining.
    """
    if not survey.scans:
        return survey

    # Concatenate all RA/Dec to find the center via RCR median.
    all_ra = np.concatenate([s.ra for s in survey.scans])
    all_dec = np.concatenate([s.dec for s in survey.scans])
    all_w = np.concatenate([s.dumps for s in survey.scans])

    # Handle RA wrap-around: if the spread is > 180 degrees, unwrap by
    # rotating values > 180 down by 360. (Matches C++ zeroCrossCheck.)
    if all_ra.max() - all_ra.min() > 180.0:
        all_ra = np.where(all_ra > 180.0, all_ra - 360.0, all_ra)
        for s in survey.scans:
            s.ra = np.where(s.ra > 180.0, s.ra - 360.0, s.ra)

    # Center the projection on the data — the median of the science
    # sample positions. The FITS header's RA/DEC is the COMMANDED
    # PARKING position (where the dish sits for the diode calibration),
    # NOT where the source appears during the raster. For tracked
    # moving objects, the actual raster is centered on the source's
    # apparent position at scan time, which usually differs from the
    # parking position.
    center_ra_deg = _weighted_median(all_ra, all_w)
    center_dec_deg = _weighted_median(all_dec, all_w)

    cos_dec_center = math.cos(center_dec_deg * _DEG2RAD)

    for scan in survey.scans:
        if per_sample_cos:
            cos_factor = np.cos(scan.dec * _DEG2RAD)
        else:
            cos_factor = cos_dec_center
        scan.ra_proj = (scan.ra - center_ra_deg) * cos_factor
        scan.dec_proj = scan.dec - center_dec_deg

        # Along-scan angular distance: cumulative arc length from first
        # sample. Used by every later stage that thinks in 1D (BG sub
        # windowing, etc.).
        #
        # IMPORTANT: compute from SKY coordinates via great-circle
        # distance, NOT from the projected `ra_proj`/`dec_proj`. With
        # per-sample SFL the projected `ra_proj` varies along a scan
        # even when sky-RA is constant (because cos(dec) varies sample-
        # to-sample), so `np.hypot(diff_ra_proj, diff_dec_proj)`
        # contains a spurious component proportional to (sky_ra - center)
        # times d(cos(dec))/d(sample). On scans far from the central
        # meridian that spurious term dominates and the resulting
        # ang_dist is way too large — which makes BG-sub windowing
        # collapse and the local quadratic fit produces wild outputs
        # (we saw scan 18 sample j=10 going from raw 6.9 Jy to -79.9 Jy
        # after BG sub). Matches C++ Scan::updateAngDistTemp
        # (Scan.cpp:442-454) which calls Tools::getGCDistance on sky
        # coords explicitly.
        if scan.size > 1:
            # Small-angle great-circle step: for the sub-degree distances
            # between consecutive samples, the spherical-law-of-cosines
            # form `arccos(sin·sin + cos·cos·cos)` loses precision around
            # cos≈1 (gives 0 for two close samples). Use the planar
            # small-angle form `sqrt((Δra·cos(midDec))² + Δdec²)` which
            # is accurate to <1e-9 deg per step for our typical sub-deg
            # sample spacing and never collapses to 0.
            ra_rad_mid = 0.5 * (scan.dec[:-1] + scan.dec[1:]) * _DEG2RAD
            d_ra = (scan.ra[1:] - scan.ra[:-1]) * np.cos(ra_rad_mid)
            d_dec = scan.dec[1:] - scan.dec[:-1]
            gc_steps_deg = np.hypot(d_ra, d_dec)
            scan.ang_dist = np.concatenate([[0.0], np.cumsum(gc_steps_deg)])
        elif scan.size == 1:
            scan.ang_dist = np.zeros(1)
        else:
            scan.ang_dist = np.zeros(0)

    # Stamp the partition with what we know so far. Edges come later.
    from rcpy.types import PartitionSet  # local import to avoid cycle

    ra_proj_all = np.concatenate([s.ra_proj for s in survey.scans])
    dec_proj_all = np.concatenate([s.dec_proj for s in survey.scans])

    min_ra = float(ra_proj_all.min())
    max_ra = float(ra_proj_all.max())
    min_dec = float(dec_proj_all.min())
    max_dec = float(dec_proj_all.max())

    survey.partition_sss = PartitionSet(
        map_type=survey.map_type,
        min_ra=min_ra,
        max_ra=max_ra,
        min_dec=min_dec,
        max_dec=max_dec,
        center_ra_deg=center_ra_deg,
        center_dec_deg=center_dec_deg,
        median_ra=float(np.median(ra_proj_all)),
        median_dec=float(np.median(dec_proj_all)),
        tracking=survey.tracking,
    )

    return survey


def unproject(ra_proj: np.ndarray, dec_proj: np.ndarray,
              center_ra_deg: float, center_dec_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """Invert the sinusoidal projection — for FITS WCS / output coords.

    For per-sample SFL: ra_proj = (ra - ra_c) * cos(dec). To invert we
    need cos(dec) of the original sample, which we recover from
    dec_proj first (dec_proj = dec - dec_c is exact, no projection).
    """
    dec_deg = dec_proj + center_dec_deg
    cos_dec_per_sample = np.cos(dec_deg * _DEG2RAD)
    # Guard against the pole where cos→0 (division blows up). In
    # practice we never observe within microdegrees of the pole, but
    # protect just in case.
    ra_deg = np.where(np.abs(cos_dec_per_sample) > 1e-12,
                      ra_proj / cos_dec_per_sample + center_ra_deg,
                      center_ra_deg)
    return ra_deg, dec_deg
