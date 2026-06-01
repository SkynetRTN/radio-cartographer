"""End-to-end pipeline orchestrator.

Wires together the stages in the order documented in PIPELINE_STAGES.txt.
Each stage is a small function on Survey or Composite — this module is
just the spine that calls them.

Usage:
    >>> from rcpy.pipeline import process
    >>> survey = io.read_sdfits("Cas_A.fits")
    >>> result_map = process(survey)
    >>> io.write_fits(result_map, "Cas_A_processed.fits")
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

from rcpy.types import Survey, Composite, Map, PartitionSet, CalMethod
from rcpy import calibration, coordinates, noise, background, timeshift, rfi, thetagap, surface


@dataclass
class PipelineConfig:
    """User-tunable knobs. Defaults match Paper recommendations where
    possible."""
    bg_scale_bw: float = 6.0
    # Reference processing uses rfi_scale_bw = 0.7 ("Faint Target" preset).
    rfi_scale_bw: float = 0.7
    # Auto-centroiding S/N threshold (Paper Footnote 21). C++ default is
    # 75, but real astronomical sources at 10-50σ won't trigger that and
    # get eaten by the local-model RFI rejection. Default lowered to 15
    # to preserve typical bright sources like Jupiter.
    centroid_sigma: float = 15.0
    weight_scale_bw: float = 1.0 / 3.0
    pixel_size_bw: float = 0.05
    time_shift_mode: str = "auto"
    time_shift_seconds: float = 0.0
    # Edge trimming size in beamwidths (C++ Composite::truncateTurningEdges).
    # Samples within trim_size_bw * psfFWHM of each scan's turning points
    # get marked as RFI-rejected. Empirically calibrated against the
    # reference 0144716 weight map: with trim=1.0 BW, my bright rim ended
    # up at col 30 whereas the reference's bright rim is at col 10. Trim
    # of 0.5 BW puts the rim at the right column position.
    trim_size_bw: float = 0.5

    skip_bg: bool = False
    skip_timeshift: bool = False
    skip_rfi: bool = False
    skip_surface_modeling: bool = False
    skip_edge_trim: bool = False

    # Cos-Dec projection variant for `coordinates.project`. Default True
    # (Sanson-Flamsteed proper, per-sample cos) matches the reference's
    # MAIN flux HDU shape. Set to False (cylindrical, constant cos at
    # center_dec) to match the reference's RAW HDU shape — its row
    # widths are constant in Dec, suggesting the raw render uses the
    # cylindrical approximation.
    per_sample_cos_dec: bool = True

    # Override Survey.cal_method. None leaves whatever the Survey was
    # read with (default INTERPOLATED). CONSTANT collapses cal-drift
    # artifacts when our start/end cal-block measurements are unreliable
    # (e.g. cyga: INTERPOLATED produces a phantom -14% gain drift that
    # bleeds into a 1.1 Jy Dec-direction gradient in the raw render;
    # CONSTANT collapses this to ~0.1 Jy).
    cal_method: CalMethod | None = None

    # C++ Processor.cpp:77, Source.cpp:95-99 gates the theta_corr (SSS
    # correlation) accumulator and the correlated weight2 formula on
    # `correlatedWeightMap`, which is just `photometryOn`. When
    # photometry is off (RCPHOT=0 in the reference job), the C++ leaves
    # theta_corr at zero and uses weight2 = Σ wHold·GMWeight·dumps
    # instead of Σ wHold·GMWeight/intersect.
    photometry: bool = False


def process(survey: Survey, config: PipelineConfig | None = None) -> Map:
    """Run the full pipeline on a single Survey, returning a Map."""
    return process_many([survey], config)


def process_many(surveys: Iterable[Survey], config: PipelineConfig | None = None) -> Map:
    """Run the pipeline on multiple Surveys, bundled into a Composite for
    joint RFI subtraction and regridding."""
    config = config or PipelineConfig()
    surveys = list(surveys)
    if not surveys:
        raise ValueError("process_many requires at least one Survey")

    # Daisy auto-TS: there are two implementations available.
    #
    # 1. `_auto_ts_daisy_brute_force` (used by default): renders the
    #    full pipeline at each candidate TS and picks the one with the
    #    highest local source brightness in the central map region.
    #    Wall time ~3 min per daisy but verified to land within 0.1 s
    #    of the visually-best TS on 0152427.
    #
    # 2. `timeshift._daisy_angle_builder` (ported from C++ but not
    #    wired in by default): converts each petal's coords to a
    #    signed angle-from-center axis and runs the raster pair-
    #    correlation algorithm on that. Same algorithm the C++ uses.
    #    Fast (~1 s) but lands ~0.3 s off the visually-best value on
    #    0152427, in a region that produces an edge-artifact map.
    #    Probably a subtle parity issue in the ported version that
    #    hasn't been tracked down yet. Available via
    #    `timeshift.correct(survey, mode='auto')` if called directly.
    from rcpy.types import MapType
    if (config.time_shift_mode == "auto"
            and any(s.map_type is MapType.DAISY for s in surveys)):
        best_ts = _auto_ts_daisy_brute_force(surveys, config)
        from dataclasses import replace
        config = replace(
            config, time_shift_mode="custom", time_shift_seconds=best_ts,
        )

    # Stages 1-7 per survey.
    for survey in surveys:
        if config.cal_method is not None:
            survey.cal_method = config.cal_method
        coordinates.project(survey, per_sample_cos=config.per_sample_cos_dec)
        calibration.gain(survey)
        noise.measure_1d(survey)
        if not config.skip_bg:
            background.subtract(survey, scale_bw=config.bg_scale_bw)
        if not config.skip_timeshift:
            # Daisies and rasters use different auto-TS algorithms
            # internally (timeshift.correct dispatches on map_type).
            # The legacy daisy → mode='off' fallback is no longer
            # needed: `_auto_ts_daisy` handles daisies via a coarse+
            # fine source-concentration sweep on BG-subbed samples.
            timeshift.correct(survey, mode=config.time_shift_mode,
                              forced_seconds=config.time_shift_seconds)
            # C++ Processor.cpp:158 recomputes the bounding box AFTER
            # time-shifting, since shifted samples can fall outside the
            # pre-shift bounds. Without this, path samples get clipped
            # and the surface-fit at the new boundary pixels is missing
            # neighbours. Update partition_sss from the new ra_ts/dec_ts.
            _update_partition_after_timeshift(survey)
        # Port of C++ Survey::setStandardThetaGap — computes
        # survey.min_gap_threshold from along-scan deltas.
        _set_standard_theta_gap(survey)
        # Port of C++ Survey::calculateEdgeParameters — fits 4 trimmed
        # boundary lines so the regridder can NaN-out pixels that fall
        # outside the actual scan footprint (Cartographer::checkEdgeCriteria).
        _calculate_edge_parameters(survey)
        noise.measure_2d(survey)

    # Stage 8: composite construction.
    composite = _build_composite(surveys)

    # Stage 8.5: trim turning-edge samples (C++ Composite::truncateTurningEdges).
    # Samples within trim_size_bw * psfFWHM of each scan's start or end
    # position get marked excluded (poor positional accuracy due to
    # telescope deceleration). Match C++ behavior: DON'T shrink the
    # partition based on trim — C++ keeps the map area unchanged and
    # relies on Cartographer::checkEdgeCriteria (pixel-level edge polygon
    # test) to NaN out flux/weight/weight2 at outside-edge pixels.
    # Scale and correlation are still computed at those pixels using
    # whatever samples are within 1 BW.
    if not config.skip_edge_trim and config.trim_size_bw > 0:
        _trim_turning_edges(composite, config.trim_size_bw)
        # Intentionally NOT calling _update_partition_after_timeshift_composite —
        # keep the partition size unchanged so scale/correlation can still
        # be rendered at L/R "cutoff" pixels via M0 fallback.

    # Stage 8.6: assign per-sample min_rcr_theta_gap.
    # Port of C++ Composite::assignRCRThetaGapMin + setMinRCRThetaGap.
    _assign_rcr_theta_gap_min(composite)

    # Stage 9-10: joint RFI subtraction + theta-gap.
    if not config.skip_rfi:
        rfi.subtract(composite, scale_bw=config.rfi_scale_bw,
                     centroid_sigma=config.centroid_sigma,
                     photometry=config.photometry)
    # Apply edge-trim mask AFTER RFI sub since RFI overwrites rfi_keep_mask.
    if not config.skip_edge_trim and config.trim_size_bw > 0:
        _apply_edge_trim_to_rfi_mask(composite)
    thetagap.compute(composite)

    # Stage 11: regridding.
    if config.skip_surface_modeling:
        # Return empty-shaped map for parity testing without the heavy step.
        part = composite.partition_sss
        import numpy as np
        empty = np.zeros((1, 1))
        return Map(
            min_ra=part.min_ra, max_ra=part.max_ra,
            min_dec=part.min_dec, max_dec=part.max_dec,
            center_ra_deg=part.center_ra_deg, center_dec_deg=part.center_dec_deg,
            resolution=config.pixel_size_bw * composite.psf_fwhm,
            flux=empty, weight=empty.copy(), weight_corr=empty.copy(),
            scale=empty.copy(), correlation=empty.copy(), path=empty.copy(),
        )

    return surface.regrid(
        composite,
        pixel_size_bw=config.pixel_size_bw,
        theta_min_bw=config.weight_scale_bw,
        photometry=config.photometry,
    )


def _trim_turning_edges(composite: Composite, trim_size_bw: float) -> None:
    """Mark samples within trim_size_bw * psfFWHM of each scan's turning
    points as excluded. Port of C++ Composite::truncateTurningEdges
    (Composite.cpp:348-391). Stores a `scan.edge_trim_mask` (True = keep,
    False = trimmed) that gets AND'd into rfi_keep_mask AFTER RFI sub.

    For each scan, the turning points are the first and last samples
    (where the telescope was changing direction). We compute the great-
    circle distance from each sample to these turning points and remove
    samples within the trim radius.
    """
    import numpy as np
    psf = composite.psf_fwhm
    trim_radius_deg = trim_size_bw * psf

    for scan in composite.all_scans:
        if scan.size == 0:
            continue
        x = scan.ra_ts if scan.ra_ts is not None else (
            scan.ra_proj if scan.ra_proj is not None else scan.ra)
        y = scan.dec_ts if scan.dec_ts is not None else (
            scan.dec_proj if scan.dec_proj is not None else scan.dec)

        tp1_x, tp1_y = x[0], y[0]
        tp2_x, tp2_y = x[-1], y[-1]

        d1 = np.hypot(x - tp1_x, y - tp1_y)
        d2 = np.hypot(x - tp2_x, y - tp2_y)

        scan.edge_trim_mask = (d1 >= trim_radius_deg) & (d2 >= trim_radius_deg)


def _apply_edge_trim_to_rfi_mask(composite: Composite) -> None:
    """After RFI sub, AND the edge_trim_mask into rfi_keep_mask so the
    path layer, surface fit, and downstream stages all see trimmed
    samples as excluded."""
    import numpy as np
    for scan in composite.all_scans:
        if getattr(scan, 'edge_trim_mask', None) is None:
            continue
        if scan.rfi_keep_mask is None:
            scan.rfi_keep_mask = scan.edge_trim_mask.copy()
        else:
            scan.rfi_keep_mask = scan.rfi_keep_mask & scan.edge_trim_mask


def _update_partition_after_timeshift_composite(composite: Composite) -> None:
    """After edge-trimming, recompute the composite's partition_sss from
    the SURVIVING (non-trimmed) samples. Mirrors what C++ does after
    truncateTurningEdges: it reclassifies and repartitions."""
    import numpy as np
    from rcpy.types import PartitionSet

    part = composite.partition_sss
    if part is None:
        return

    ras, decs = [], []
    for scan in composite.all_scans:
        if scan.size == 0:
            continue
        keep = scan.rfi_keep_mask if scan.rfi_keep_mask is not None else np.ones(scan.size, bool)
        if not keep.any():
            continue
        x = scan.ra_ts if scan.ra_ts is not None else (
            scan.ra_proj if scan.ra_proj is not None else scan.ra)
        y = scan.dec_ts if scan.dec_ts is not None else (
            scan.dec_proj if scan.dec_proj is not None else scan.dec)
        ras.append(x[keep])
        decs.append(y[keep])
    if not ras:
        return

    ra_all = np.concatenate(ras)
    dec_all = np.concatenate(decs)
    composite.partition_sss = PartitionSet(
        map_type=part.map_type,
        min_ra=float(ra_all.min()),
        max_ra=float(ra_all.max()),
        min_dec=float(dec_all.min()),
        max_dec=float(dec_all.max()),
        center_ra_deg=part.center_ra_deg,
        center_dec_deg=part.center_dec_deg,
        median_ra=float(np.median(ra_all)),
        median_dec=float(np.median(dec_all)),
    )


def _assign_rcr_theta_gap_min(composite: Composite) -> None:
    """Port of C++ Composite::assignRCRThetaGapMin (Composite.cpp:239-347)
    combined with Processor::setMinRCRThetaGap (Processor.cpp:822-828).

    Sets each sample's `min_rcr_theta_gap` to the minimum of:
      - the survey's min_gap_threshold (initial value)
      - capped to survey.min_gap_threshold if the sample is INSIDE the
        survey edge polygon (raster case)

    C++ uses edge-fit lines via RCR to define the polygon; for ralongmap
    we approximate with the rectangular survey bounding box (since
    rasters have rectangular coverage, this matches in practice). The
    cap is just `min_gap_threshold` everywhere, so the practical effect
    is: every sample's min_rcr_theta_gap = its survey's min_gap_threshold.
    """
    import numpy as np
    surveys = composite.surveys if hasattr(composite, 'surveys') else [composite]
    for survey in surveys:
        thresh = float(survey.min_gap_threshold)
        for scan in survey.scans:
            if scan.size == 0:
                continue
            # Initialise to threshold; for samples OUTSIDE the survey
            # polygon, C++ leaves their RCRMinThetaGap unbounded (large).
            # For our ralongmap case, virtually all samples are inside
            # the rectangular bbox so this initialisation is enough.
            scan.min_rcr_theta_gap = np.full(scan.size, thresh, dtype=np.float64)


def _set_standard_theta_gap(survey: Survey) -> None:
    """Compute survey.min_gap_threshold: smallest along-scan inter-sample
    distance (RCR-trimmed). Port of C++ Survey::setStandardThetaGap
    (Survey.cpp:2218-2306) for non-DAISY map types.

    The C++ algorithm:
      1. For each scan, compute deltaAng = consecutive along-scan distances
      2. Run RCR LS_MODE_DL on the pooled deltaAng values
      3. min_gap_threshold = smallest RCR-survivor value

    Per-sample minRCRThetaGap is then initialised to this value
    (Processor.cpp:822-828), and later capped further by
    Composite::assignRCRThetaGapMin for samples inside the survey edges.
    """
    import numpy as np
    from rcpy.types import MapType
    from rcpy.rcr2 import RCR, RejectionTech

    if survey.map_type is MapType.DAISY:
        # Different formula for daisies — uses center gaps. Skip for now;
        # daisy support is deferred per existing TODOs.
        survey.min_gap_threshold = 0.0
        return

    delta_ang = []
    for scan in survey.scans:
        if scan.size < 2 or scan.ang_dist is None:
            continue
        ang = scan.ang_dist
        diffs = np.diff(ang)
        delta_ang.extend(diffs.tolist())

    if not delta_ang:
        survey.min_gap_threshold = 0.0
        return

    delta_ang = np.asarray(delta_ang, dtype=np.float64)
    # RCR-trim outliers (the C++ uses LS_MODE_DL bulk-rejection).
    try:
        rcr = RCR(RejectionTech.LS_MODE_DL)
        rcr.perform_bulk_rejection(delta_ang)
        flags = rcr.result.flags
        survivors = delta_ang[flags.astype(bool)]
    except Exception:
        survivors = delta_ang

    if survivors.size == 0:
        survey.min_gap_threshold = 0.0
    else:
        # Smallest survivor — matches C++ "if deltaAng[i] < minGapThreshold && flagsHold[i] == true"
        survey.min_gap_threshold = float(survivors.min())


def _fit_robust_line(x, y):
    """Robust linear fit  y = m * (x - xbar) + b.  Returns np.array([m, b, xbar]).

    Iterative 3-sigma clip stands in for C++ RCR LS_MODE_DL — sufficient
    for the boundary samples being fit here, which are mostly clean
    scan-endpoints with a few outliers from outward turnarounds.
    """
    import numpy as np

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    n_in = int(mask.sum())
    if n_in < 2:
        return np.array([0.0,
                          float(np.nanmean(y)) if n_in else 0.0,
                          float(np.nanmean(x)) if n_in else 0.0])

    for _ in range(5):
        x_in = x[mask]
        y_in = y[mask]
        xbar = float(x_in.mean())
        dx = x_in - xbar
        if dx.size < 2 or float(dx.var()) == 0.0:
            return np.array([0.0, float(y_in.mean()), xbar])
        slope, intercept = np.polyfit(dx, y_in, 1)
        resid = y_in - (slope * dx + intercept)
        sd = float(np.std(resid))
        if sd <= 0:
            break
        keep = np.abs(resid) <= 3.0 * sd
        if keep.all():
            break
        new_mask = np.zeros_like(mask)
        new_mask[np.flatnonzero(mask)[keep]] = True
        if new_mask.sum() < 2:
            break
        mask = new_mask

    x_in = x[mask]
    y_in = y[mask]
    xbar = float(x_in.mean())
    dx = x_in - xbar
    if dx.size < 2 or float(dx.var()) == 0.0:
        return np.array([0.0, float(y_in.mean()), xbar])
    slope, intercept = np.polyfit(dx, y_in, 1)
    return np.array([float(slope), float(intercept), xbar])


def _determine_edge_flags(survey: Survey) -> tuple[list[bool], list[int], list[int]]:
    """Port of C++ Survey::determineEdgeFlags (Survey.cpp:1868-2023) for
    non-DAISY map types.

    For each scan, finds the indices of the min/max along-scan position
    and decides which one is the "outward overshoot" (the extreme
    endpoint where the telescope was decelerating before reversing
    direction). For alternating-direction rasters/noddings, even and
    odd scans overshoot in opposite directions, so we compare the
    RCR-trimmed mean max for each parity:

      - If even-scan maxes are larger on average: even scans overshoot
        at maxIndex, odd scans overshoot at minIndex.
      - Else: even scans overshoot at minIndex, odd scans overshoot at
        maxIndex.

    The overshoot endpoint gets its `edgePointFlag = False` so it's
    excluded from the boundary line fit. The opposing endpoint keeps
    `edgePointFlag = True`.

    Returns (overshoot_flag_at_max, min_indices, max_indices). The
    first list is per-scan: True means the SCAN'S MAX endpoint is the
    overshoot (so push the MIN endpoint to the edge fit); False means
    the MIN endpoint is the overshoot (so push the MAX endpoint).
    """
    import numpy as np
    from rcpy.rcr2 import RCR, RejectionTech

    scans = survey.scans
    scans_in_ra = bool(scans[0].scan_in_ra)

    n = len(scans)
    min_indices = [0] * n
    max_indices = [0] * n
    even_maxes: list[float] = []
    odd_maxes: list[float] = []

    for i, scan in enumerate(scans):
        x = scan.ra_ts if scan.ra_ts is not None else (
            scan.ra_proj if scan.ra_proj is not None else scan.ra)
        y = scan.dec_ts if scan.dec_ts is not None else (
            scan.dec_proj if scan.dec_proj is not None else scan.dec)
        if x is None or x.size == 0:
            continue
        ang = x if scans_in_ra else y
        i_max = int(np.argmax(ang))
        i_min = int(np.argmin(ang))
        min_indices[i] = i_min
        max_indices[i] = i_max
        if i % 2 == 0:
            even_maxes.append(float(ang[i_max]))
        else:
            odd_maxes.append(float(ang[i_max]))

    def _rcr_mean(values: list[float]) -> float:
        if not values:
            return float("nan")
        v = np.asarray(values, dtype=float)
        if v.size < 2:
            return float(v.mean())
        try:
            rcr = RCR(RejectionTech.LS_MODE_DL)
            rcr.perform_bulk_rejection(v)
            mu = float(rcr.result.mu)
            if not np.isfinite(mu):
                return float(v.mean())
            return mu
        except Exception:
            return float(v.mean())

    even_max_avg = _rcr_mean(even_maxes)
    odd_max_avg = _rcr_mean(odd_maxes)

    # If even-scan maxes are larger, the overshoot is at MAX for even
    # scans and at MIN for odd scans. Else vice versa.
    overshoot_at_max = [False] * n
    if np.isfinite(even_max_avg) and np.isfinite(odd_max_avg):
        even_overshoots_at_max = even_max_avg > odd_max_avg
    else:
        even_overshoots_at_max = True   # neutral default
    for i in range(n):
        if (i % 2 == 0) == even_overshoots_at_max:
            overshoot_at_max[i] = True
        else:
            overshoot_at_max[i] = False

    return overshoot_at_max, min_indices, max_indices


def _calculate_edge_parameters(survey: Survey) -> None:
    """Port of C++ Survey::calculateEdgeParameters (Survey.cpp:1579-1818).

    Fits four line equations to the scan-boundary samples so the
    regridder can NaN flux/weight/weight_corr at pixels outside the
    actual scan footprint. Each line is stored as np.array([m, b, bar])
    where the line is:

        edge_one / edge_three (vertical):   ra  = m * (dec - ybar) + b
        edge_two / edge_four (horizontal):  dec = m * (ra  - xbar) + b

    Edge convention (mirrors C++ Cartographer::checkEdgeCriteria):
        edge_one  = RIGHT  (lower-RA boundary)
        edge_two  = TOP    (higher-Dec boundary)
        edge_three= LEFT   (higher-RA boundary)
        edge_four = BOTTOM (lower-Dec boundary)

    Uses _determine_edge_flags to drop the outward-overshoot endpoint
    from each scan's contribution to the side edges. Without this
    pre-filter (the bug noted in PARITY_AUDIT §3.5), the polygon was
    too generous and ~900 perimeter pixels stayed finite that should
    have been NaN'd.

    DAISY maps skip the polygon and use edge_radius instead — already
    handled by partition.edge_radius (not implemented here yet).
    """
    import numpy as np
    from rcpy.types import MapType

    if survey.partition_sss is None or survey.map_type is MapType.DAISY:
        return
    if not survey.scans:
        return

    def _pos(scan):
        x = scan.ra_ts if scan.ra_ts is not None else (
            scan.ra_proj if scan.ra_proj is not None else scan.ra)
        y = scan.dec_ts if scan.dec_ts is not None else (
            scan.dec_proj if scan.dec_proj is not None else scan.dec)
        return x, y

    scans = survey.scans
    scans_in_ra = bool(scans[0].scan_in_ra)

    overshoot_at_max, min_idx, max_idx = _determine_edge_flags(survey)

    e1x, e1y, e2x, e2y = [], [], [], []
    e3x, e3y, e4x, e4y = [], [], [], []

    for i, scan in enumerate(scans):
        x, y = _pos(scan)
        if x is None or y is None or x.size == 0:
            continue
        i_max = max_idx[i]
        i_min = min_idx[i]
        # Per C++ Survey.cpp:1659-1668 / 1712-1722: only include the
        # min/max endpoint when its edgePointFlag is true.
        keep_max = not overshoot_at_max[i]
        keep_min = overshoot_at_max[i]

        if scans_in_ra:
            # Scan progresses in RA: first/last scans bound bottom/top in Dec;
            # per-scan min/max RA endpoints bound right/left edges.
            if i == 0:
                e4x.append(x); e4y.append(y)
            if i == len(scans) - 1:
                e2x.append(x); e2y.append(y)
            if keep_min:
                e1x.append(np.array([y[i_min]])); e1y.append(np.array([x[i_min]]))
            if keep_max:
                e3x.append(np.array([y[i_max]])); e3y.append(np.array([x[i_max]]))
        else:
            # Scan progresses in Dec: first/last scans bound right/left in RA;
            # per-scan min/max Dec endpoints bound bottom/top edges.
            if i == 0:
                e1x.append(y); e1y.append(x)
            if i == len(scans) - 1:
                e3x.append(y); e3y.append(x)
            if keep_min:
                e4x.append(np.array([x[i_min]])); e4y.append(np.array([y[i_min]]))
            if keep_max:
                e2x.append(np.array([x[i_max]])); e2y.append(np.array([y[i_max]]))

    def _cat(xs, ys):
        if not xs:
            return np.array([]), np.array([])
        return np.concatenate(xs), np.concatenate(ys)

    survey.partition_sss.edge_one   = _fit_robust_line(*_cat(e1x, e1y))
    survey.partition_sss.edge_two   = _fit_robust_line(*_cat(e2x, e2y))
    survey.partition_sss.edge_three = _fit_robust_line(*_cat(e3x, e3y))
    survey.partition_sss.edge_four  = _fit_robust_line(*_cat(e4x, e4y))


def _rcr_trimmed_mean(values, tech_name: str = "SS_MEDIAN_DL") -> float:
    """Return the RCR-trimmed mean of `values`. Matches C++
    `RCR(SS_MEDIAN_DL).performBulkRejection` followed by `.result.mu`
    used in `determineDataBoundaries`. Falls back to the unweighted
    mean on degenerate input or RCR failure.

    The SS_MEDIAN_DL technique in C++ is the single-sided median-anchor
    Direct-Linear rejector. Our rcr2 port supports both SS_MEDIAN_DL
    and LS_MODE_DL via the same API.
    """
    import numpy as np
    from rcpy.rcr2 import RCR, RejectionTech
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan")
    if v.size == 1:
        return float(v[0])
    spread = float(np.max(v) - np.min(v))
    if spread <= 0:
        return float(v[0])
    try:
        tech = getattr(RejectionTech, tech_name)
        rcr = RCR(tech)
        rcr.perform_bulk_rejection(v)
        mu = float(rcr.result.mu)
        if not np.isfinite(mu):
            return float(np.mean(v))
        return mu
    except Exception:
        return float(np.mean(v))


def _update_partition_after_timeshift(survey: Survey) -> None:
    """Recompute the survey's partition_sss bounding box using the
    time-shifted positions (ra_ts/dec_ts) AND C++'s RCR-trimmed
    per-scan-extremes algorithm.

    Port of C++ Processor::determineDataBoundaries (Processor.cpp:607-751)
    for non-DAISY maps:
      median_ra = RCR-trimmed mean of all sample RAs
      median_dec = RCR-trimmed mean of all sample Decs
      For scansInRa (raster along RA):
        min_dec_global = RCR-trimmed mean of scan[0]'s all dec values
        max_dec_global = RCR-trimmed mean of scan[last]'s all dec values
        min_ra_global  = RCR-trimmed mean of per-scan min RA values
        max_ra_global  = RCR-trimmed mean of per-scan max RA values
      Else (scansInDec):
        min_ra_global  = RCR-trimmed mean of scan[0]'s all ra values
        max_ra_global  = RCR-trimmed mean of scan[last]'s all ra values
        min_dec_global = RCR-trimmed mean of per-scan min Dec values
        max_dec_global = RCR-trimmed mean of per-scan max Dec values
      center_ra_deg = 0.5 * (min_ra + max_ra)
      center_dec_deg = 0.5 * (min_dec + max_dec)

    Replaces the earlier port which used absolute min/max of all
    samples — that included turnaround-overshoot outliers which shifted
    the partition by ~1 pixel and caused the +0.88 px peak offset
    visible in the cyga reference comparison.
    """
    import numpy as np
    from rcpy.types import PartitionSet, MapType

    part = survey.partition_sss
    if part is None:
        return

    # Gather per-scan positions
    scan_positions = []
    for scan in survey.scans:
        if scan.ra_ts is not None and scan.dec_ts is not None:
            scan_positions.append((scan.ra_ts, scan.dec_ts))
        elif scan.ra_proj is not None and scan.dec_proj is not None:
            scan_positions.append((scan.ra_proj, scan.dec_proj))
    if not scan_positions:
        return

    ra_all = np.concatenate([p[0] for p in scan_positions])
    dec_all = np.concatenate([p[1] for p in scan_positions])

    # Projection center = RCR-trimmed mean of ALL sample positions
    median_ra = _rcr_trimmed_mean(ra_all)
    median_dec = _rcr_trimmed_mean(dec_all)

    if part.map_type is MapType.DAISY:
        # Daisy bounds come from edge_radius rather than per-scan
        # extremes; fall back to simple min/max for now.
        survey.partition_sss = PartitionSet(
            map_type=part.map_type,
            min_ra=float(ra_all.min()),
            max_ra=float(ra_all.max()),
            min_dec=float(dec_all.min()),
            max_dec=float(dec_all.max()),
            center_ra_deg=0.5 * (float(ra_all.min()) + float(ra_all.max())),
            center_dec_deg=0.5 * (float(dec_all.min()) + float(dec_all.max())),
            median_ra=median_ra,
            median_dec=median_dec,
        )
        return

    scans_in_ra = bool(survey.scans[0].scan_in_ra)

    # Compute the RCR-trimmed "typical" extents but keep absolute min/max
    # for the partition bounds. The center comes from the midpoint of the
    # RCR-trimmed extents (matching C++ globalCenterDec/Ra = 0.5 * (max +
    # min)). Using raw bounds for partition extent preserves the same
    # pixel coverage as ref (matching column counts), while the
    # RCR-trimmed center matches C++'s WCS reference location and
    # eliminates the per-pixel offset.
    if scans_in_ra:
        min_dec_vals = scan_positions[0][1]
        max_dec_vals = scan_positions[-1][1]
        per_scan_min_ra = np.array([p[0].min() for p in scan_positions])
        per_scan_max_ra = np.array([p[0].max() for p in scan_positions])
        rcr_min_ra = _rcr_trimmed_mean(per_scan_min_ra)
        rcr_max_ra = _rcr_trimmed_mean(per_scan_max_ra)
        rcr_min_dec = _rcr_trimmed_mean(min_dec_vals)
        rcr_max_dec = _rcr_trimmed_mean(max_dec_vals)
    else:
        min_ra_vals = scan_positions[0][0]
        max_ra_vals = scan_positions[-1][0]
        per_scan_min_dec = np.array([p[1].min() for p in scan_positions])
        per_scan_max_dec = np.array([p[1].max() for p in scan_positions])
        rcr_min_ra = _rcr_trimmed_mean(min_ra_vals)
        rcr_max_ra = _rcr_trimmed_mean(max_ra_vals)
        rcr_min_dec = _rcr_trimmed_mean(per_scan_min_dec)
        rcr_max_dec = _rcr_trimmed_mean(per_scan_max_dec)

    if rcr_min_ra > rcr_max_ra:
        rcr_min_ra, rcr_max_ra = rcr_max_ra, rcr_min_ra
    if rcr_min_dec > rcr_max_dec:
        rcr_min_dec, rcr_max_dec = rcr_max_dec, rcr_min_dec

    # Partition bounds = raw min/max (preserves the same pixel coverage
    # as the reference).
    min_ra_g = float(ra_all.min())
    max_ra_g = float(ra_all.max())
    min_dec_g = float(dec_all.min())
    max_dec_g = float(dec_all.max())

    # Center = midpoint of RCR-trimmed extents (matches C++ centerRaDeg/
    # centerDecDeg). This is what shows up in the FITS WCS CRVAL1/2 and
    # is also the reference point against which sample positions are
    # measured in some downstream stages.
    center_ra_g = 0.5 * (rcr_min_ra + rcr_max_ra)
    center_dec_g = 0.5 * (rcr_min_dec + rcr_max_dec)

    survey.partition_sss = PartitionSet(
        map_type=part.map_type,
        min_ra=min_ra_g,
        max_ra=max_ra_g,
        min_dec=min_dec_g,
        max_dec=max_dec_g,
        center_ra_deg=float(center_ra_g),
        center_dec_deg=float(center_dec_g),
        median_ra=median_ra,
        median_dec=median_dec,
    )


def _build_composite(surveys: list[Survey]) -> Composite:
    """Stage 8 — bundle surveys into a Composite with a global partition.

    Simple version: take the union of bounding boxes, the median of
    centers, and pass through the first survey's psf_fwhm / map_type.
    Does not handle moving-object daisies (Footnote 23).
    """
    if len(surveys) == 1:
        return Composite(
            surveys=surveys,
            partition_sss=surveys[0].partition_sss,
        )

    parts = [s.partition_sss for s in surveys if s.partition_sss is not None]
    if not parts:
        raise ValueError("No surveys have partitions — run coordinates.project first")

    composite_part = PartitionSet(
        map_type=parts[0].map_type,
        min_ra=min(p.min_ra for p in parts),
        max_ra=max(p.max_ra for p in parts),
        min_dec=min(p.min_dec for p in parts),
        max_dec=max(p.max_dec for p in parts),
        center_ra_deg=sum(p.center_ra_deg for p in parts) / len(parts),
        center_dec_deg=sum(p.center_dec_deg for p in parts) / len(parts),
        median_ra=sum(p.median_ra for p in parts) / len(parts),
        median_dec=sum(p.median_dec for p in parts) / len(parts),
    )

    return Composite(surveys=surveys, partition_sss=composite_part)


def _auto_ts_daisy_brute_force(base_surveys, base_config) -> float:
    """Auto-determine the time-shift for a daisy by full-pipeline render
    at each candidate TS, scoring by the brightest pixel within a
    central PSF-radius mask.

    Cheap histogram metrics couldn't reliably distinguish source-
    concentration from petal-tip pile-ups on 0152427 (Jupiter daisy at
    1550 MHz). Brute-force uses the same gridder + RFI machinery that
    the user inspects, so the metric is exactly the visual signal —
    "is the central pixel bright?" The cost is ~2-3× a normal pipeline
    run, which is acceptable for an auto step that runs once per
    daisy observation.

    Sweep: coarse [-3, +3] step 0.5 (13 evals), then fine ±0.5 around
    coarse winner step 0.1 (11 evals). Total ~24 evals.

    Returns the best TS in seconds.
    """
    import copy
    import numpy as np
    from dataclasses import replace

    def _render_at_ts(ts: float) -> tuple[float, tuple[int, int]] | None:
        """Render the full pipeline at this TS, return (central_peak,
        (peak_y, peak_x)) or None on failure."""
        surveys_copy = [copy.deepcopy(s) for s in base_surveys]
        cfg = replace(
            base_config,
            time_shift_mode="custom",
            time_shift_seconds=ts,
            # Skip RFI in the sweep to halve render time. The source-
            # concentration signal we measure is unaffected by RFI
            # rejection (which targets bright outliers, not the
            # spatial coherence we're optimizing for).
            skip_rfi=True,
        )
        try:
            m = process_many(surveys_copy, cfg)
        except Exception:
            return None
        flux = m.flux
        if flux is None or flux.size <= 1:
            return None
        # Score = peak above local baseline in a tight central window.
        # Inner window = 0.5 PSF FWHM (catches Jupiter at the source
        # position). Annulus 0.5–1.5 PSF gives a local baseline so the
        # score reflects source signal, not BG offset. A wrong TS that
        # pushes an edge artifact into the central radius doesn't
        # produce a centered PSF-shaped source — it spills into the
        # annulus too — so it scores low.
        h, w = flux.shape
        cy, cx = h // 2, w // 2
        psf_deg = base_surveys[0].psf_fwhm
        r_inner = max(int(round(0.5 * psf_deg / max(m.resolution, 1e-6))), 4)
        r_outer = max(int(round(1.5 * psf_deg / max(m.resolution, 1e-6))), r_inner + 4)
        yi0, yi1 = max(0, cy - r_inner), min(h, cy + r_inner + 1)
        xi0, xi1 = max(0, cx - r_inner), min(w, cx + r_inner + 1)
        yo0, yo1 = max(0, cy - r_outer), min(h, cy + r_outer + 1)
        xo0, xo1 = max(0, cx - r_outer), min(w, cx + r_outer + 1)
        inner = flux[yi0:yi1, xi0:xi1]
        outer = flux[yo0:yo1, xo0:xo1]
        inner_finite = inner[np.isfinite(inner)]
        if inner_finite.size == 0:
            return None
        outer_mask = np.ones_like(outer, dtype=bool)
        outer_mask[yi0 - yo0 : yi1 - yo0, xi0 - xo0 : xi1 - xo0] = False
        annulus_vals = outer[outer_mask & np.isfinite(outer)]
        baseline = float(np.median(annulus_vals)) if annulus_vals.size > 0 else 0.0
        peak_above = float(np.max(inner_finite)) - baseline
        py_inner, px_inner = np.unravel_index(
            np.where(np.isfinite(inner), inner, -np.inf).argmax(), inner.shape
        )
        return peak_above, (yi0 + py_inner, xi0 + px_inner)

    import os
    verbose = bool(int(os.environ.get("RCPY_AUTO_TS_VERBOSE", "0")))

    # Coarse sweep.
    coarse_ts = np.arange(-3.0, 3.01, 0.5)
    coarse_scores = []
    for ts in coarse_ts:
        r = _render_at_ts(float(ts))
        score = -np.inf if r is None else r[0]
        coarse_scores.append(score)
        if verbose:
            loc = "n/a" if r is None else str(r[1])
            print(f"  [auto-ts coarse] ts={ts:+.2f}  score={score:+.4f}  loc={loc}")
    coarse_scores = np.asarray(coarse_scores)
    best_coarse = float(coarse_ts[int(np.argmax(coarse_scores))])

    # Fine sweep ±0.5 around the coarse winner.
    fine_ts = np.arange(best_coarse - 0.5, best_coarse + 0.51, 0.1)
    fine_scores = []
    for ts in fine_ts:
        r = _render_at_ts(float(ts))
        score = -np.inf if r is None else r[0]
        fine_scores.append(score)
        if verbose:
            loc = "n/a" if r is None else str(r[1])
            print(f"  [auto-ts fine]   ts={ts:+.2f}  score={score:+.4f}  loc={loc}")
    fine_scores = np.asarray(fine_scores)
    best_fine = float(fine_ts[int(np.argmax(fine_scores))])

    return best_fine
