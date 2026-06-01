"""Stage 9 — surface modeling / regridding (Paper §3.7).

STATUS: PARTIAL. Implements the per-pixel local polynomial fit at 1st and
3rd order, with the cos^alpha weighting function (Eq. 12) and a basic
sufficiency check. Skips: the noise-level prior (Footnote 27), the
2nd-order fallback, the rigorous sufficiency check (Footnote 28's per-
scan diversity requirement), and the proximity-weighted interpolation of
theta_gap to pixel centers (Footnote 31 — we use nearest-neighbour).

Good enough to produce a real-shaped pixel image from clean point data.

C++ reference: src/Cartographer.cpp.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.spatial import cKDTree

from rcpy.types import Survey, Composite, Map, PartitionSet
from rcpy._pivot import _pivot_system_const_term


def _alpha_from_fwhm(theta_w: float) -> float:
    """Power for the cos^alpha weighting so FWHM == theta_w (Eq. 13)."""
    if theta_w <= 0:
        return 1e6
    c = math.cos(math.pi * theta_w / 4.0)
    if c <= 0 or c >= 1:
        return 1.0
    return -math.log(2.0) / math.log(c)


def _cos_weight(distance_bw: np.ndarray, theta_w: float) -> np.ndarray:
    """Eq. 12: cos^alpha((pi/2) * distance), zero beyond 1 beamwidth."""
    out = np.zeros_like(distance_bw)
    mask = distance_bw < 1.0
    if not mask.any():
        return out
    alpha = _alpha_from_fwhm(theta_w)
    out[mask] = np.cos(math.pi * distance_bw[mask] / 2.0) ** alpha
    return out


def _check_scan_diversity(scan_ids: np.ndarray, min_scans: int, min_per_scan: int) -> bool:
    """Paper Footnote 28: a polynomial fit requires samples spread across
    multiple scans, with at least one scan having enough samples to
    constrain the model along the scan axis.

    Args:
        scan_ids: per-sample scan index.
        min_scans: minimum number of distinct scans required.
        min_per_scan: minimum samples in at least one scan.
    """
    unique, counts = np.unique(scan_ids, return_counts=True)
    if unique.size < min_scans:
        return False
    if counts.max() < min_per_scan:
        return False
    return True


def _solve_with_prior(basis: np.ndarray, z: np.ndarray, w: np.ndarray,
                      apply_prior: bool) -> float:
    """Weighted least-squares solve of basis · c = z with optional
    Footnote-27 noise-level prior on the constant term, returning c[0].

    Faithful port of the C++ surfaceFit10/6/3 pattern
    (Cartographer.cpp:1502): builds the explicit normal-equations
    `B^T W B · c = B^T W z`, applies the noise-level prior by doubling
    the (0,0) entry (Cartographer.cpp:1499 `A[0] = 2 * w`), then calls
    `Tools::pivotSystem` and returns `b[0] / A[0][0]` of the resulting
    lower-triangular system — which is exactly x[0] of the full solve.

    Uses `_pivot_system_const_term` (verbatim port of pivotSystem +
    forward-sub) instead of `np.linalg.solve`. The unusual backward-
    Gauss elimination order in C++ produces different garbage on near-
    singular matrices than LU/SVD; for high-M polynomial fits with
    clustered samples this routinely matters.

    C++ has NO fallback for singular cases — the divide-by-zero/NaN
    just propagates. We match that exactly (no `try/except LinAlgError`).
    """
    BtW = basis.T * w[None, :]
    A = BtW @ basis
    rhs = BtW @ z
    if apply_prior:
        A[0, 0] *= 2.0
    n = A.shape[0]
    return _pivot_system_const_term(n, A.flatten(), rhs)


def _fit_poly3(dx: np.ndarray, dy: np.ndarray, z: np.ndarray, w: np.ndarray,
               m10_plus_processing: bool = True) -> float:
    """Weighted 3rd-order 2D polynomial fit, returns constant term (a00).

    Polynomial basis: 1, dx, dy, dx^2, dx*dy, dy^2, dx^3, dx^2*dy,
    dx*dy^2, dy^3 — 10 terms (Eq. 11).

    If `m10_plus_processing` is true and the first solve yields a
    negative constant, recompute with the Footnote-27 noise-level
    prior (port of C++ Cartographer::surfaceFit10 lines 1496-1513
    recursive call with m10PlusCriteria=false).
    """
    basis = np.column_stack([
        np.ones_like(dx), dx, dy,
        dx ** 2, dx * dy, dy ** 2,
        dx ** 3, dx ** 2 * dy, dx * dy ** 2, dy ** 3,
    ])
    a00 = _solve_with_prior(basis, z, w, apply_prior=False)
    if m10_plus_processing and a00 < 0.0:
        a00 = _solve_with_prior(basis, z, w, apply_prior=True)
    return a00


def _fit_poly1(dx: np.ndarray, dy: np.ndarray, z: np.ndarray, w: np.ndarray,
               m10_plus_processing: bool = True) -> float:
    """Weighted 1st-order (plane) fit. Returns constant term.

    Mirrors C++ Cartographer::surfaceFit3 (lines 1059-1077) which also
    has the m10PlusCriteria recursive-fallback path.
    """
    basis = np.column_stack([np.ones_like(dx), dx, dy])
    a00 = _solve_with_prior(basis, z, w, apply_prior=False)
    if m10_plus_processing and a00 < 0.0:
        a00 = _solve_with_prior(basis, z, w, apply_prior=True)
    return a00


def _build_edge_mask(part, nx: int, ny: int, resolution: float) -> np.ndarray:
    """Return a bool array (ny, nx) — True for pixels inside the survey
    edge polygon. Port of C++ Cartographer::checkEdgeCriteria
    (Cartographer.cpp:686-781) for non-DAISY map types.

    Lines stored on the PartitionSet as np.array([m, b, bar]):
        edge_one  / edge_three (vertical):   ra  = m*(dec - ybar) + b
        edge_two  / edge_four  (horizontal): dec = m*(ra  - xbar) + b
    Inside polygon iff
        edge_one(ra-right) < ra < edge_three(ra-left)
        AND
        edge_four(dec-bottom) < dec < edge_two(dec-top).
    Falls back to all-True if the edges aren't populated.
    """
    if (part.edge_one is None or part.edge_two is None
            or part.edge_three is None or part.edge_four is None):
        return np.ones((ny, nx), dtype=bool)

    ra = part.min_ra + np.arange(nx) * resolution     # (nx,)
    dec = part.min_dec + np.arange(ny) * resolution   # (ny,)

    m1, b1, y1 = part.edge_one
    m2, b2, x2 = part.edge_two
    m3, b3, y3 = part.edge_three
    m4, b4, x4 = part.edge_four

    ra_right = m1 * (dec - y1) + b1   # (ny,)
    ra_left  = m3 * (dec - y3) + b3   # (ny,)
    dec_top  = m2 * (ra  - x2) + b2   # (nx,)
    dec_bot  = m4 * (ra  - x4) + b4   # (nx,)

    # Broadcast: ra (nx,) compared against ra_right/ra_left (ny,) → (ny, nx).
    in_ra = (ra_right[:, None] < ra[None, :]) & (ra[None, :] < ra_left[:, None])
    in_dec = (dec_bot[None, :] < dec[:, None]) & (dec[:, None] < dec_top[None, :])
    return in_ra & in_dec


def _fit_poly2(dx: np.ndarray, dy: np.ndarray, z: np.ndarray, w: np.ndarray,
               m10_plus_processing: bool = True) -> float:
    """Weighted 2nd-order (quadratic) 2D fit. Returns constant term (M6).

    Basis: 1, dx, dy, dx², dy², dx*dy (6 coefficients).
    Port of C++ Cartographer::surfaceFit6 (Cartographer.cpp:1091-1247),
    including the noise-level prior recursive fallback at line 1234.
    """
    basis = np.column_stack([
        np.ones_like(dx), dx, dy,
        dx ** 2, dy ** 2, dx * dy,
    ])
    a00 = _solve_with_prior(basis, z, w, apply_prior=False)
    if m10_plus_processing and a00 < 0.0:
        a00 = _solve_with_prior(basis, z, w, apply_prior=True)
    return a00


def regrid(
    target: Survey | Composite,
    pixel_size_bw: float = 0.05,
    theta_min_bw: float = 2.0 / 3.0,
    photometry: bool = False,
) -> Map:
    """Produce a Map from the RFI-subtracted, theta-gapped data.

    Args:
        pixel_size_bw: pixel size in beamwidths. Default 0.05 (= 1/20 BW).
        theta_min_bw: minimum weighting scale in beamwidths. Default 2/3.
        photometry: if True, use the correlated weight_corr formula
            (Σ wHold·GMWeight/intersect) and rely on per-sample theta_corr
            for the correlation layer. If False (the C++ default when
            photometryOn=0 / RCPHOT=0), use weight_corr = Σ wHold·GMWeight·dumps
            and skip the intersect calculation. Match RCPHOT in the FITS
            header of the reference render.
    """
    psf = target.psf_fwhm
    if psf <= 0:
        raise ValueError("psf_fwhm must be set on the Survey/Composite")

    # Resolution in projected-coordinate units
    resolution = pixel_size_bw * psf

    # Determine bounding box from partition
    if isinstance(target, Composite):
        part = target.partition_sss
    else:
        part = target.partition_sss
    if part is None:
        raise ValueError("No PartitionSet on the target — run coordinate prep first")

    nx = max(int(round((part.max_ra - part.min_ra) / resolution)) + 1, 1)
    ny = max(int(round((part.max_dec - part.min_dec) / resolution)) + 1, 1)

    # Gather all kept samples into one flat dataset.
    xs, ys, zs, ws, gaps, sids, tcs, gmws = [], [], [], [], [], [], [], []
    all_scans = target.all_scans if isinstance(target, Composite) else target.scans
    for scan in all_scans:
        keep = scan.rfi_keep_mask if scan.rfi_keep_mask is not None else np.ones(scan.size, bool)
        if not keep.any():
            continue
        x = (scan.ra_ts if scan.ra_ts is not None else
             (scan.ra_proj if scan.ra_proj is not None else scan.ra))[keep]
        y = (scan.dec_ts if scan.dec_ts is not None else
             (scan.dec_proj if scan.dec_proj is not None else scan.dec))[keep]
        z = (scan.flux_rfi if scan.flux_rfi is not None else
             (scan.flux_bg if scan.flux_bg is not None else scan.flux))[keep]
        w = scan.dumps[keep]
        g = (scan.theta_gap if scan.theta_gap is not None else
             np.full(scan.size, 1.0))[keep]
        tc = (scan.theta_corr if scan.theta_corr is not None else
              np.zeros(scan.size))[keep]
        gmw = (scan.gm_weight if scan.gm_weight is not None else
               np.ones(scan.size))[keep]
        xs.append(x)
        ys.append(y)
        zs.append(z)
        ws.append(w)
        gaps.append(g)
        sids.append(np.full(int(keep.sum()), scan.scan_index, dtype=np.int32))
        tcs.append(tc)
        gmws.append(gmw)

    if not xs:
        # Empty target: return zero map
        empty = np.zeros((ny, nx))
        return Map(
            min_ra=part.min_ra, max_ra=part.max_ra,
            min_dec=part.min_dec, max_dec=part.max_dec,
            center_ra_deg=part.center_ra_deg, center_dec_deg=part.center_dec_deg,
            resolution=resolution,
            flux=empty, weight=empty.copy(), weight_corr=empty.copy(),
            scale=empty.copy(), correlation=empty.copy(), path=empty.copy(),
        )

    xy = np.column_stack((np.concatenate(xs), np.concatenate(ys)))
    z_all = np.concatenate(zs)
    w_all = np.concatenate(ws)
    g_all = np.concatenate(gaps)
    sid_all = np.concatenate(sids)
    tc_all = np.concatenate(tcs)
    gmw_all = np.concatenate(gmws)

    tree = cKDTree(xy)

    # Initialize flux/weight/weight_corr to NaN — pixels that never get
    # a surface fit (or that fall in the M0 outside-edge region) keep
    # NaN. C++ Cartographer::checkEdgeCriteria explicitly NaN's flux,
    # weight, and weight2 at pixels outside the survey edge polygon
    # (Cartographer.cpp:263-271). Scale and correlation are computed
    # for ALL pixels, regardless of edge.
    flux = np.full((ny, nx), np.nan)
    weight = np.full((ny, nx), np.nan)
    weight_corr = np.full((ny, nx), np.nan)
    scale = np.zeros((ny, nx))
    correlation = np.zeros((ny, nx))
    # C++ initialises every pixel to -1.0 in Pixel.cpp:13; sample-painted
    # pixels are overwritten with +1.0 (even scans) or +0.5 (odd scans)
    # in Cartographer.cpp:147-151. Reference FITS has exactly 3 distinct
    # path values: {-1.0, +0.5, +1.0}.
    path = np.full((ny, nx), -1.0)

    # Pre-compute per-pixel inside/outside the survey edge polygon. If
    # the partition has the 4 edge-line parameters populated (from
    # pipeline._calculate_edge_parameters), we use them; otherwise we
    # default to "all inside" (matches old behaviour).
    edge_inside = _build_edge_mask(part, nx, ny, resolution)

    search_radius = 1.0 * psf

    # Iterate pixel centers. For real performance this should be vectorized
    # or parallelized — but for the spike, a row loop is fine.
    for iy in range(ny):
        dec_c = part.min_dec + iy * resolution
        for ix in range(nx):
            ra_c = part.min_ra + ix * resolution

            idxs = tree.query_ball_point((ra_c, dec_c), r=search_radius)
            if not idxs:
                continue
            idxs = np.asarray(idxs)
            local_xy = xy[idxs]
            # Polynomial basis offsets: keep in DEGREES so the basis values
            # match C++ Cartographer::surfaceFit10 (lines 1304-1305:
            #   xHold = inRange[6*i+2] - i_0;   // raw degrees
            #   yHold = inRange[6*i+3] - j_0;
            # ). The 3rd-order WLS solve is scale-invariant in the constant
            # term mathematically, but the conditioning of the normal
            # equations differs by orders of magnitude between BW units
            # (basis values O(1)) and degree units (basis values O(psf)),
            # which can shift the constant term by ~percent for ill-
            # conditioned point clouds.
            dx_deg = local_xy[:, 0] - ra_c
            dy_deg = local_xy[:, 1] - dec_c
            # Keep a BW-units distance for the cos^alpha radial weight, which
            # paper Eq. 12 defines on BW-scaled distance.
            dx = dx_deg / psf  # in beamwidths
            dy = dy_deg / psf
            distance_bw = np.hypot(dx, dy)

            # Interpolate theta_gap from neighbouring samples to this
            # pixel center using the log-weighted average of Paper
            # Footnote 31 (and C++ Cartographer::getPixelThetaGapSSS).
            #
            #   a    = -2.329 * ln(theta_gap/2) - 0.510
            #   w(d) = (-ln d)^a   for d in (0, 1) beamwidths
            #
            # Then theta_w = max(theta_min, min(4/3 * <theta_gap>, 1)).
            gaps_local = g_all[idxs]
            d_safe = np.clip(distance_bw, 1e-6, 0.999)
            # Avoid log(theta_gap/2) when theta_gap is tiny → blowup.
            tg_safe = np.clip(gaps_local, 1e-3, 1.0)
            exponent = -2.329 * np.log(tg_safe / 2.0) - 0.510
            log_d = -np.log(d_safe)
            log_w = np.where(log_d > 0, log_d, 0.0) ** exponent
            log_w[~np.isfinite(log_w)] = 0.0
            ws = float(log_w.sum())
            if ws > 0:
                local_gap = float((log_w * gaps_local).sum() / ws)
            else:
                local_gap = float(gaps_local.max() if gaps_local.size else 1.0)
            theta_w = max(theta_min_bw, min(4.0 / 3.0 * local_gap, 1.0))

            radial_w = _cos_weight(distance_bw, theta_w)
            # C++ Cartographer::surfaceFit3 uses pure cos^alpha (wHold) for
            # BOTH the WLS normal equations AND the weight output (toRet[1]).
            # No dumps factor — including dumps here was a parity bug that
            # made my "weight" layer differ from the reference.
            sample_w = radial_w
            if sample_w.sum() <= 0:
                continue

            # Enforce paper Footnote 28 scan-diversity requirements. Port of
            # C++ Cartographer::determineSurfaceType (Cartographer.cpp:911-925):
            #   M10 (poly3, 10 coeffs): planeApplication(5, 5, 10)
            #   M6  (poly2, 6 coeffs):  planeApplication(4, 4, 6)
            #   M3  (poly1, 3 coeffs):  planeApplication(2, 2, 3)
            #   M0  (weighted mean):    fallback when none of the above
            local_sids = sid_all[idxs][sample_w > 0]
            n_samples = int((sample_w > 0).sum())
            value = float("nan")
            model_type = "M0"
            if n_samples >= 10 and _check_scan_diversity(local_sids, min_scans=5, min_per_scan=5):
                value = _fit_poly3(dx_deg, dy_deg, z_all[idxs], sample_w)
                model_type = "M10"
            elif n_samples >= 6 and _check_scan_diversity(local_sids, min_scans=4, min_per_scan=4):
                value = _fit_poly2(dx_deg, dy_deg, z_all[idxs], sample_w)
                model_type = "M6"
            elif n_samples >= 3 and _check_scan_diversity(local_sids, min_scans=2, min_per_scan=2):
                value = _fit_poly1(dx_deg, dy_deg, z_all[idxs], sample_w)
                model_type = "M3"
            # else: value stays NaN, model_type stays "M0"

            inside_edge = bool(edge_inside[iy, ix])

            if inside_edge:
                flux[iy, ix] = value
                if model_type == "M0":
                    # C++ Cartographer.cpp:940-955: in M0 fallback, weight uses
                    # pure cos (NOT cos^alpha). This is what produces the bright
                    # top/bottom rim bands visible in C++ reference renders.
                    m0_cos = np.cos(math.pi * distance_bw / 2.0)
                    m0_cos = np.where(distance_bw < 1.0, m0_cos, 0.0)
                    weight[iy, ix] = float(m0_cos.sum())
                else:
                    weight[iy, ix] = float(sample_w.sum())
            # Else: flux, weight, weight_corr stay NaN (C++ Cartographer.cpp:263-271).
            # Scale/correlation are still computed below regardless of inside_edge.
            scale[iy, ix] = theta_w

            # weight2 / weight_corr per C++ Cartographer.cpp:1023-1026:
            # toRet[2] = wLocalW = Σ wHold · GMWeight / intersect
            # where intersect counts inRange samples j with
            # GC-distance(sample_i, sample_j) < theta_corr_i / 2.
            # Note: NOT exported to FITS in C++, but used internally for
            # photometry corrections. Compute it now that gm_weight and
            # theta_corr are available per sample.
            tc_local = tc_all[idxs]  # per-sample theta_corr in degrees
            gmw_local = gmw_all[idxs]
            dumps_local = w_all[idxs]  # per-sample DataDumps
            local_pos = xy[idxs]
            if inside_edge:
                if photometry:
                    # C++ Cartographer.cpp:1023-1024: photometryOn branch
                    # uses wHold · GMWeight / intersect, where intersect
                    # counts inRange samples j with distance from sample
                    # i < theta_corr_i / 2.
                    wlocal = 0.0
                    for i_loc in range(idxs.size):
                        if tc_local[i_loc] <= 0:
                            intersect = 1
                        else:
                            d_ij = np.hypot(
                                local_pos[:, 0] - local_pos[i_loc, 0],
                                local_pos[:, 1] - local_pos[i_loc, 1],
                            )
                            intersect = int((d_ij < tc_local[i_loc] / 2.0).sum())
                            if intersect == 0:
                                intersect = 1
                        wlocal += sample_w[i_loc] * gmw_local[i_loc] / intersect
                    weight_corr[iy, ix] = float(wlocal)
                else:
                    # C++ Cartographer.cpp:1025-1026: photometryOff branch
                    # uses wHold · GMWeight · DataDumps. Faster (no
                    # intersect loop) and matches the reference render
                    # for RCPHOT=0 jobs.
                    weight_corr[iy, ix] = float(
                        (sample_w * gmw_local * dumps_local).sum()
                    )
            # Else: weight_corr stays NaN.

            # Correlation map per Paper Appendix D / C++
            # Cartographer::getCorrMapValueSSS:
            #
            #   per_sample = sqrt((theta_corr_i/psf)² + (f_w·theta_w)²)
            #   correlation = log-weighted average of per_sample over
            #                 contributing samples (same weights as
            #                 theta_gap interpolation, Footnote 31)
            #
            # theta_corr is set per-sample by the RFI subtraction; it's
            # already in degrees so we divide by psf for BW units.
            factor_w = 3.0
            tc_bw = tc_all[idxs] / psf  # in beamwidths
            corr_per_sample = np.sqrt(tc_bw ** 2 + (factor_w * theta_w) ** 2)
            if ws > 0:
                corr_val = float((log_w * corr_per_sample).sum() / ws)
            else:
                corr_val = float(factor_w * theta_w)
            correlation[iy, ix] = corr_val

    # Path map: one pixel per sample location, alternating between 1.0
    # (even scans) and 0.5 (odd scans) per C++ Cartographer.cpp:135-152.
    # IMPORTANT: C++ iterates over scan.getSize() which is the size
    # AFTER Scan::removeRFI shrank the scan by deleting NaN-rfi_subtracted
    # samples. So path only renders samples that SURVIVED the RFI sub.
    # For very bright sources (Cygnus A) where the cos² local model
    # can't capture the full amplitude, source samples get RFI-rejected
    # and create a "missing circle" in the path layer.
    for sid, scan in enumerate(all_scans):
        if scan.size == 0:
            continue
        keep = scan.rfi_keep_mask if scan.rfi_keep_mask is not None else np.ones(scan.size, bool)
        if not keep.any():
            continue
        x = (scan.ra_ts if scan.ra_ts is not None else
             (scan.ra_proj if scan.ra_proj is not None else scan.ra))[keep]
        y = (scan.dec_ts if scan.dec_ts is not None else
             (scan.dec_proj if scan.dec_proj is not None else scan.dec))[keep]
        x_idx = np.round((x - part.min_ra) / resolution).astype(int)
        y_idx = np.round((y - part.min_dec) / resolution).astype(int)
        # Clip to map bounds
        ok = (x_idx >= 0) & (x_idx < nx) & (y_idx >= 0) & (y_idx < ny)
        # Cartographer.cpp:147-151: even scans painted +1.0, odd +0.5.
        # Background was pre-filled with -1.0 above.
        value = 1.0 if (scan.scan_index % 2 == 0) else 0.5
        path[y_idx[ok], x_idx[ok]] = value

    # Match the C++ FITS-output convention: RA decreases with column
    # index (CDELT1 < 0 per astronomy standard). C++ does this reversal
    # inside getLayerData (Map.cpp:396-398):
    #     int reversed_j = (map->getSize(1) - 1) - j;
    #     array[i * naxes0 + reversed_j] = val;
    # applied to ALL output layers (main flux, path, scale, weight,
    # correlation, raw). Our rendering loop above builds arrays where
    # column 0 corresponds to min_ra (east on sky); flipping columns
    # puts max_ra at column 0 to match the FITS metadata our
    # io.write_fits already declares (CDELT1 = -resolution). Without
    # this flip our path/scale/weight/correlation panels appeared
    # horizontally mirrored relative to the C++ reference renders.
    flux = flux[:, ::-1].copy()
    weight = weight[:, ::-1].copy()
    weight_corr = weight_corr[:, ::-1].copy()
    scale = scale[:, ::-1].copy()
    correlation = correlation[:, ::-1].copy()
    path = path[:, ::-1].copy()

    return Map(
        min_ra=part.min_ra, max_ra=part.max_ra,
        min_dec=part.min_dec, max_dec=part.max_dec,
        center_ra_deg=part.center_ra_deg, center_dec_deg=part.center_dec_deg,
        resolution=resolution,
        flux=flux, weight=weight, weight_corr=weight_corr,
        scale=scale, correlation=correlation, path=path,
    )
