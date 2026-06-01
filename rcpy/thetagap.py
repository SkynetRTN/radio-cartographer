"""Stage 8 — per-sample theta_gap calculation (Paper §3.7, Fig. 49).

Faithful port of the C++ ProcessorThetaGap::circleThetaGap algorithm.

For each sample we want a robust local measure of how densely the data
is sampled around it. The C++ uses the LARGEST EMPTY CIRCLE through the
central point, computed PER QUADRANT (4 cardinal + 4 diagonal), then
takes the maximum across quadrants.

Per quadrant the algorithm:

  1. Collect all neighbour samples within 1 beamwidth of the central
     point, then keep only those that fall inside the quadrant's wedge.
  2. For each candidate neighbour, build the unique circle that
     (a) passes through both the central point and the neighbour, and
     (b) has its center on a quadrant-specific axis through the central
         point:
           TOP / BOTTOM           → center has same RA as central point
           LEFT / RIGHT           → center has same Dec as central point
           DIAG_TR / DIAG_BL      → center on the line y = -x through it
           DIAG_BR / DIAG_TL      → center on the line y = +x through it
  3. Sort circles by diameter and points by distance. The "theta_gap"
     for the quadrant is the LARGEST circle diameter that contains no
     other points inside.

Finally the per-sample theta_gap is the maximum across all 8 quadrants,
capped at 0.75 × psf_fwhm (to keep theta_w well-defined when fed
through Eq. 14: theta_w = max(theta_min, min(4/3·theta_gap, 1·psf)) ).

References:
    Paper §3.7, Fig. 49, Footnote 30.
    C++ source: src/ProcessorThetaGap.cpp.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from rcpy.types import Survey, Composite


# Per-sample theta_gap cap (C++ ProcessorThetaGap.cpp line 578-580).
# In beamwidths.
_THETA_GAP_CAP_BW = 0.75


def _quadrant_filter(x_rel: np.ndarray, y_rel: np.ndarray, quad: str) -> np.ndarray:
    """Return a boolean mask selecting points inside one of the eight
    quadrant wedges of the C++ ProcessorThetaGap::quadrantSort.

    `x_rel`, `y_rel` are the neighbour positions relative to the central
    point (so the central point sits at the origin).

    Cardinal wedges (TOP/BOTTOM/LEFT/RIGHT) are bounded by the diagonals
    y = ±x, so they're 90° wedges with their bisector along an axis.
    Diagonal wedges (DIAG_*) are bounded by the cardinal axes — they're
    the four 90° quadrants.
    """
    # Per C++ lines 179-180:
    #   raLine  = (decCheck - dec) * (-1) + ra  → x - x_self = -(y - y_self)
    #   raLine2 = (decCheck - dec) *  (1) + ra  → x - x_self =  (y - y_self)
    # In our (x_rel, y_rel) coords centered on the sample:
    #   raLine  → x_rel = -y_rel
    #   raLine2 → x_rel =  y_rel
    if quad == "TOP":            # dec > self, between the two diagonals
        return (y_rel > 0) & (x_rel < y_rel) & (x_rel > -y_rel)
    if quad == "BOTTOM":         # dec < self, between the two diagonals
        return (y_rel < 0) & (x_rel > y_rel) & (x_rel < -y_rel)
    if quad == "LEFT":           # raCheck > raLine AND raCheck > raLine2
        return (x_rel > -y_rel) & (x_rel > y_rel)
    if quad == "RIGHT":          # raCheck < raLine AND raCheck < raLine2
        return (x_rel < -y_rel) & (x_rel < y_rel)
    # The diagonal quadrants are the four 90° axis-aligned quadrants:
    if quad == "DIAG_TR":        # dec > self, ra < self  (in C++ ra<ra-self)
        return (y_rel > 0) & (x_rel < 0)
    if quad == "DIAG_BR":        # dec < self, ra < self
        return (y_rel < 0) & (x_rel < 0)
    if quad == "DIAG_TL":        # dec > self, ra > self
        return (y_rel > 0) & (x_rel > 0)
    if quad == "DIAG_BL":        # dec < self, ra > self
        return (y_rel < 0) & (x_rel > 0)
    raise ValueError(f"Unknown quadrant {quad!r}")


def _circle_center_and_diameter(
    x1: float, y1: float,
    x2: np.ndarray, y2: np.ndarray,
    quad: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute the center (x0, y0) and diameter of the circle that passes
    through (x1, y1) and each (x2[i], y2[i]) with the center constrained
    to the quadrant-specific axis through (x1, y1).

    Direct port of C++ ProcessorThetaGap::determineCircleParams.
    """
    n = x2.size
    x0 = np.empty(n)
    y0 = np.empty(n)
    if quad in ("TOP", "BOTTOM"):
        # Center stays on x = x1; solve for y0.
        y0 = (-x1 ** 2 + 2 * x1 * x2 - x2 ** 2 + y1 ** 2 - y2 ** 2) / (2 * (y1 - y2))
        x0 = np.full(n, x1)
    elif quad in ("LEFT", "RIGHT"):
        x0 = (-x1 ** 2 + x2 ** 2 + y1 ** 2 - 2 * y1 * y2 + y2 ** 2) / (-2 * x1 + 2 * x2)
        y0 = np.full(n, y1)
    elif quad in ("DIAG_TR", "DIAG_BL"):
        # Center on y - y1 = -(x - x1)  (slope -1 line, but in RA which
        # is reversed, so this is "physically" slope +1).
        denom = 2 * (x1 - x2 - y1 + y2)
        x0 = -((-x1 ** 2 + x2 ** 2 + 2 * x1 * (y1 - y2) + (y1 - y2) ** 2) / denom)
        y0 = (x1 ** 2 - 2 * x1 * x2 + x2 ** 2 + 2 * x1 * y1 - 2 * x2 * y1
              - y1 ** 2 + y2 ** 2) / denom
    elif quad in ("DIAG_BR", "DIAG_TL"):
        denom = 2 * (x1 - x2 + y1 - y2)
        x0 = -((-x1 ** 2 + x2 ** 2 - 2 * x1 * (y1 - y2) + (y1 - y2) ** 2) / denom)
        y0 = -((x1 ** 2 + x2 ** 2 + 2 * x2 * y1 - y1 ** 2
                - 2 * x1 * (x2 + y1) + y2 ** 2) / denom)
    else:
        raise ValueError(f"Unknown quadrant {quad!r}")
    radius_sq = (x1 - x0) ** 2 + (y1 - y0) ** 2
    radius_sq = np.where(np.isfinite(radius_sq) & (radius_sq >= 0), radius_sq, 0.0)
    diameter = 2.0 * np.sqrt(radius_sq)
    return x0, y0, diameter


def _max_gap_quadrant(
    x1: float, y1: float,
    x_neigh: np.ndarray, y_neigh: np.ndarray,
    quad: str,
) -> float:
    """Return the largest-empty-circle diameter for one quadrant.

    Direct port of C++ ProcessorThetaGap::maxGapQuadrant.
    """
    if x_neigh.size == 0:
        return float("inf")   # sentinel "999999" from C++

    x0, y0, diameter = _circle_center_and_diameter(x1, y1, x_neigh, y_neigh, quad)
    distance = np.sqrt((x1 - x_neigh) ** 2 + (y1 - y_neigh) ** 2)

    if diameter.size == 1:
        return float(diameter[0])

    # Sort circles by diameter ascending; sort point distances ascending.
    d_order = np.argsort(diameter)
    diameter_sorted = diameter[d_order]
    cx_sorted = x0[d_order]
    cy_sorted = y0[d_order]

    p_order = np.argsort(distance)
    distance_sorted = distance[p_order]
    px_sorted = x_neigh[p_order]
    py_sorted = y_neigh[p_order]

    # Find the smallest-diameter circle that's at least as big as the
    # nearest neighbour distance — that's where the C++ search starts.
    j = 0
    while j < diameter_sorted.size and distance_sorted[0] > diameter_sorted[j]:
        j += 1
    if j >= diameter_sorted.size:
        j = diameter_sorted.size - 1

    while True:
        # All points whose DISTANCE is less than this circle's diameter
        # are candidates for being INSIDE the circle.
        inner = np.where(distance_sorted < diameter_sorted[j])[0]
        # Check each — is the point inside the circle?
        found_inside = False
        for k in inner:
            dx = px_sorted[k] - cx_sorted[j]
            dy = py_sorted[k] - cy_sorted[j]
            # +2e-6 fudge per C++ to match exact-radius cases
            d_to_center = np.sqrt(dx * dx + dy * dy) + 2e-6
            if d_to_center < diameter_sorted[j] / 2.0:
                found_inside = True
                break
        if found_inside:
            # This circle has a point inside; the answer is the
            # PREVIOUS (smaller) diameter.
            if j > 0:
                return float(diameter_sorted[j - 1])
            return float(diameter_sorted[0])
        # No point inside this circle. Try a bigger one.
        j += 1
        if j >= diameter_sorted.size:
            return float("inf")


def _theta_gap_for_sample(
    x_self: float, y_self: float,
    x_neigh: np.ndarray, y_neigh: np.ndarray,
) -> float:
    """Compute theta_gap for one sample as the max-over-quadrants of the
    largest-empty-circle diameter in each quadrant. In the same units as
    the neighbour coordinates.

    C++ parity note: when a quadrant contains NO neighbours, C++'s
    `maxGapQuadrant` returns 999999 (effective infinity), and C++'s
    `circleThetaGap` takes the max over all 8 quadrants — so any single
    empty quadrant promotes the whole theta_gap to 999999, which then
    gets capped at the 0.75-beamwidth sentinel by `calculateThetaGapSSS`
    (ProcessorThetaGap.cpp:578-581). The earlier port silently SKIPPED
    empty quadrants and took the max over the populated ones only,
    which left edge samples (where one or more outward-facing quadrants
    are empty because there are no neighbours past the data boundary)
    with a too-small theta_gap. That under-estimate propagated into the
    Footnote-31 interpolation, lowering the per-pixel scale at edges
    and accounting for ~5% of the rcpy-vs-ref scale gap.
    """
    x_rel = x_neigh - x_self
    y_rel = y_neigh - y_self
    best = 0.0
    any_quadrant_empty = False
    for quad in ("TOP", "BOTTOM", "LEFT", "RIGHT",
                 "DIAG_TR", "DIAG_TL", "DIAG_BR", "DIAG_BL"):
        mask = _quadrant_filter(x_rel, y_rel, quad)
        if not mask.any():
            # Mirror C++ maxGapQuadrant returning 999999 for an empty
            # quadrant. We don't actually use infinity in `best` because
            # the cap at 0.75·psf is applied by the caller; instead we
            # return the cap sentinel directly. The numeric value here
            # (1.0 BW in the caller's caller frame, then capped) is what
            # C++ produces after the same cap.
            any_quadrant_empty = True
            continue
        gap = _max_gap_quadrant(x_self, y_self, x_neigh[mask], y_neigh[mask], quad)
        if gap > best and np.isfinite(gap):
            best = gap
    if any_quadrant_empty:
        # Any 999999 in C++'s max-over-quadrants forces the result to
        # 999999. Return float("inf") so the caller's eventual cap at
        # 0.75 BW kicks in.
        return float("inf")
    return best


def compute(target: Survey | Composite, psf_fwhm: float | None = None) -> Survey | Composite:
    """Populate per-sample theta_gap on every scan.

    Implements the C++ ProcessorThetaGap::circleThetaGap algorithm
    faithfully:  largest-empty-circle per quadrant, max across the 8
    quadrants, then capped at 0.75 × psf_fwhm.
    """
    if psf_fwhm is None:
        psf_fwhm = target.psf_fwhm if hasattr(target, "psf_fwhm") else 0.0
    if psf_fwhm <= 0:
        for survey in (target.surveys if isinstance(target, Composite) else [target]):
            for scan in survey.scans:
                scan.theta_gap = np.full(scan.size, 1.0)
        return target

    # Gather only RFI-KEPT samples for the kd-tree. The C++ pipeline
    # physically deletes RFI'd samples from each Scan (Scan::removeRFI,
    # Scan.cpp:655-694) and then rebuilds classificationsSSS via
    # classifySSS(composite) immediately before calling
    # calculateProcThetaGapMulti (Processor.cpp:502-519). So when
    # ProcessorThetaGap::findPossSSS walks classificationsSSS, it only
    # sees surviving samples — RFI-rejected positions show up as real
    # gaps in the spatial sampling, which is what makes the per-sample
    # theta_gap LARGER in regions where RFI has thinned the data
    # (typically around bright sources whose on-source samples were
    # falsely rejected by the cos² local model).
    #
    # An earlier port comment claimed C++ used "all science samples
    # independent of RFI" — that was wrong; the visible halo of larger
    # theta_w around Cyg A in the reference scale map is exactly the
    # signature of this rejected-sample bookkeeping.
    xs, ys, sids, smps = [], [], [], []
    surveys = target.surveys if isinstance(target, Composite) else [target]
    all_scans = target.all_scans if isinstance(target, Composite) else target.scans
    for sid, scan in enumerate(all_scans):
        if scan.ra_ts is not None:
            x, y = scan.ra_ts, scan.dec_ts
        elif scan.ra_proj is not None:
            x, y = scan.ra_proj, scan.dec_proj
        else:
            x, y = scan.ra, scan.dec
        keep = scan.rfi_keep_mask if scan.rfi_keep_mask is not None else np.ones(scan.size, bool)
        if not keep.any():
            continue
        # Track which sample index each kept-array entry maps to inside
        # the parent scan, so we can write the result back to the right
        # slot in scan.theta_gap.
        kept_indices = np.where(keep)[0].astype(np.int32)
        xs.append(x[keep])
        ys.append(y[keep])
        sids.append(np.full(kept_indices.size, sid, dtype=np.int32))
        smps.append(kept_indices)

    xy = np.column_stack((np.concatenate(xs), np.concatenate(ys)))
    scan_id_all = np.concatenate(sids)
    sample_id_all = np.concatenate(smps)

    # Per-sample min_rcr_theta_gap (in degrees, same as positions). C++
    # uses this to exclude too-close neighbours from the theta_gap
    # calculation (ProcessorThetaGap.cpp:485 — "distanceCheck >
    # scans[i_0].getRCRMinThetaGap(j_0)"). Without this filter, samples
    # at scan turnarounds (where consecutive samples are nearly on top
    # of each other) get spuriously small theta_gap values.
    # Must align with the kept-sample ordering used to build xy above.
    min_dist_per_sample = []
    for sid, scan in enumerate(all_scans):
        keep = scan.rfi_keep_mask if scan.rfi_keep_mask is not None else np.ones(scan.size, bool)
        if not keep.any():
            continue
        if scan.min_rcr_theta_gap is not None:
            min_dist_per_sample.append(scan.min_rcr_theta_gap[keep])
        else:
            min_dist_per_sample.append(np.zeros(int(keep.sum())))
    min_dist_per_sample = np.concatenate(min_dist_per_sample) if min_dist_per_sample else np.zeros(0)

    # Pre-allocate per-scan theta_gap (in BEAMWIDTHS — same as before).
    for sid, scan in enumerate(all_scans):
        scan.theta_gap = np.full(scan.size, _THETA_GAP_CAP_BW)

    tree = cKDTree(xy)
    neighbours = tree.query_ball_point(xy, r=psf_fwhm)

    for i in range(xy.shape[0]):
        idx = np.asarray(neighbours[i], dtype=np.int64)
        # Drop self
        idx = idx[idx != i]
        if idx.size == 0:
            continue
        # Filter neighbours by minimum distance threshold (RCRMinThetaGap)
        thresh = min_dist_per_sample[i]
        if thresh > 0:
            d = np.hypot(xy[idx, 0] - xy[i, 0], xy[idx, 1] - xy[i, 1])
            idx = idx[d > thresh]
            if idx.size == 0:
                continue
        gap_deg = _theta_gap_for_sample(
            xy[i, 0], xy[i, 1], xy[idx, 0], xy[idx, 1]
        )
        if not np.isfinite(gap_deg) or gap_deg <= 0:
            continue
        gap_bw = gap_deg / psf_fwhm
        gap_bw = min(gap_bw, _THETA_GAP_CAP_BW)
        sid = int(scan_id_all[i])
        smp = int(sample_id_all[i])
        all_scans[sid].theta_gap[smp] = gap_bw

    return target
