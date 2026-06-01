"""SDFITS ingest + FITS output.

STATUS: PARTIAL. write_fits is functional (single-extension); read_sdfits
is wired to the column schema in src/PreProcessor.cpp::accessExtTable20m
but should be considered untested until run against a real Skynet file.

20-meter SDFITS columns expected:
    UTSECS, CRVAL2 (RA), CRVAL3 (Dec), AZIMUTH, ELEVATIO,
    NSAMPS (dump count), CALSTATE, SWPINDEX, SWPVALID, CDELT1, DATA
Primary header: OBSFREQ, MJD, TELESCOP, MAP-PATT.

C++ reference: src/PreProcessor.cpp, src/GBParser.cpp.
"""
from __future__ import annotations

from typing import Optional
import math
import numpy as np

from rcpy.types import (
    Survey, Scan, Map, Telescope, MapType, Channel,
    CalMethod, Coordinates,
)


# Beam size lookup (degrees FWHM). Approximate from Paper §2 and 20-m
# observing handbooks. Real port should pull these from the FITS header
# or a config table.
_BEAM_FWHM_DEG = {
    (Telescope.TWENTY_METER, "L"): 0.55,
    (Telescope.TWENTY_METER, "X"): 0.083,
    (Telescope.FOURTY_FOOT, "L"): 1.16,
}


def _band_from_freq(freq_ghz: float) -> str:
    if 1.0 < freq_ghz < 2.0:
        return "L"
    if 8.0 < freq_ghz < 10.0:
        return "X"
    return "?"


def read_sdfits(
    path: str,
    channel: Channel = Channel.COMPOSITE,
    inclusion_band_mhz: tuple[float, float] = (1355.0, 1455.0),
    exclusion_bands_mhz: list[tuple[float, float]] | None = None,
    include_swpvalid_0: bool = False,
) -> Survey:
    """Read a Skynet SDFITS file into a Survey.

    Skynet 20-meter SDFITS layout:
      - Two polarizations (PLNUM 0=L, 1=R) interleaved at the same
        timestamps and coordinates.
      - DATA is a (n_rows, n_channels) spectrum per row.
      - Per-row frequency mapping from CRVAL1/CRPIX1/CDELT1. CDELT1 is
        typically NEGATIVE on Skynet data (channels go high freq -> low),
        so a naive linear channel index does the wrong thing.
      - SWPINDEX groups samples into individual sweeps.
      - CALSTATE = 1 marks samples with the noise diode firing.
      - BMAJ in the primary header gives the beam FWHM directly.

    Continuum integration matches C++ PreProcessor::averageSpectra:
    SUM (not mean) all in-band channels that fall outside any exclusion
    notches. C++ comment: "Continuum should be the entire sum, not the
    average".

    Args:
        path: SDFITS file.
        channel: which polarization channel becomes the working flux.
        inclusion_band_mhz: (lo, hi) frequency range to integrate, in MHz.
            Default 1355-1455 matches the C++ FAINT-preset default.
        exclusion_bands_mhz: list of (lo, hi) tuples to NOTCH out of the
            integration (e.g. for RFI lines).
    """
    from astropy.io import fits

    with fits.open(path) as hdul:
        prim = hdul[0].header
        tbl = hdul[1].data

        # Build per-row frequency axis. CDELT1 < 0 means channels run
        # high-to-low frequency, which is the Skynet convention.
        data = np.asarray(tbl["DATA"])
        if data.ndim == 2:
            crval1 = np.asarray(tbl["CRVAL1"], dtype=np.float64)
            crpix1 = np.asarray(tbl["CRPIX1"], dtype=np.float64)
            cdelt1 = np.asarray(tbl["CDELT1"], dtype=np.float64)
            nchan = data.shape[1]
            # frequencies[i, j] in Hz
            chan_offsets = np.arange(nchan)
            freqs_hz = (
                crval1[:, None]
                + (chan_offsets[None, :] - (crpix1[:, None] - 1.0)) * cdelt1[:, None]
            )

            band_lo = inclusion_band_mhz[0] * 1e6
            band_hi = inclusion_band_mhz[1] * 1e6

            in_band = (freqs_hz >= band_lo) & (freqs_hz <= band_hi)
            for notch_lo, notch_hi in (exclusion_bands_mhz or []):
                in_notch = (freqs_hz >= notch_lo * 1e6) & (freqs_hz <= notch_hi * 1e6)
                in_band &= ~in_notch

            # Sum, not mean — matches C++ averageSpectra.
            flux_band = (data * in_band).sum(axis=1).astype(np.float64)
        else:
            flux_band = data.astype(np.float64)

        plnum = np.asarray(tbl["PLNUM"], dtype=np.int32)
        # Split by PLNUM into L (=0) and R (=1). Rows for the two
        # polarizations should be paired at the same UTSECS.
        l_mask = plnum == 0
        r_mask = plnum == 1
        if not l_mask.any() or not r_mask.any():
            # Single-polarization file: use whichever exists for both.
            l_mask = r_mask = np.ones(len(tbl), dtype=bool)

        # Use the L-pol row ordering as canonical; R-pol values get
        # matched 1:1 by index (which holds because L and R are
        # interleaved at identical timestamps in the Skynet schema).
        all_utsecs = np.asarray(tbl["UTSECS"][l_mask], dtype=np.float64)
        all_ra = np.asarray(tbl["CRVAL2"][l_mask], dtype=np.float64)
        all_dec = np.asarray(tbl["CRVAL3"][l_mask], dtype=np.float64)
        all_az = np.asarray(tbl["AZIMUTH"][l_mask], dtype=np.float64)
        all_el = np.asarray(tbl["ELEVATIO"][l_mask], dtype=np.float64)
        all_nsamps = np.asarray(tbl["NSAMPS"][l_mask], dtype=np.float64)
        all_cal = np.asarray(tbl["CALSTATE"][l_mask], dtype=np.int32)
        all_swpv = np.asarray(tbl["SWPVALID"][l_mask], dtype=np.int32)
        all_swpidx = np.asarray(tbl["SWPINDEX"][l_mask], dtype=np.int32)
        all_flux_l = flux_band[l_mask].astype(np.float64)
        all_flux_r = flux_band[r_mask].astype(np.float64)

        # Guard against L/R shape mismatch (single-pol files etc.)
        n = min(
            all_flux_l.size, all_flux_r.size, all_utsecs.size,
            all_swpv.size, all_swpidx.size, all_cal.size,
        )
        all_utsecs = all_utsecs[:n]
        all_ra = all_ra[:n]
        all_dec = all_dec[:n]
        all_az = all_az[:n]
        all_el = all_el[:n]
        all_nsamps = all_nsamps[:n]
        all_cal = all_cal[:n]
        all_swpv = all_swpv[:n]
        all_swpidx = all_swpidx[:n]
        all_flux_l = all_flux_l[:n]
        all_flux_r = all_flux_r[:n]

        # SWPVALID=1 are science samples; SWPVALID=0 are calibration
        # samples (cal-on diode samples, plus turnaround samples).
        #
        # C++ Survey::dataProc40 (Survey.cpp:506-591) actually loads ALL
        # samples regardless of SWPVALID — it only splits scans by
        # SWPINDEX. The optional `include_swpvalid_0` flag mimics this
        # behaviour. By default we exclude SWPVALID=0 (cleaner pipeline,
        # avoids cal-recovery spikes), but including them produces the
        # dense sample-pile-up regions at scan endpoints that appear as
        # the bright spots in C++ reference scale/weight maps.
        if include_swpvalid_0:
            sci_mask = np.ones_like(all_swpv, dtype=bool)
        else:
            sci_mask = all_swpv == 1
        cal_mask = all_swpv == 0

        # Science arrays — feed to scans.
        utsecs = all_utsecs[sci_mask]
        ra = all_ra[sci_mask]
        dec = all_dec[sci_mask]
        az = all_az[sci_mask]
        el = all_el[sci_mask]
        nsamps = all_nsamps[sci_mask]
        cal = all_cal[sci_mask]
        swpidx = all_swpidx[sci_mask]
        flux_l = all_flux_l[sci_mask]
        flux_r = all_flux_r[sci_mask]

        # Calibration arrays — feed to calibration stage.
        cal_time = all_utsecs[cal_mask]
        cal_state = all_cal[cal_mask]
        cal_dumps = all_nsamps[cal_mask]
        cal_flux_l = all_flux_l[cal_mask]
        cal_flux_r = all_flux_r[cal_mask]

    # Telescope + frequency from header
    telescope_name = (prim.get("TELESCOP", "") or "").upper()
    if "40" in telescope_name:
        telescope = Telescope.FOURTY_FOOT
    elif "20" in telescope_name or "GREENBANK-20" in telescope_name:
        telescope = Telescope.TWENTY_METER
    else:
        telescope = Telescope.TWENTY_METER  # best guess

    # OBSFREQ is in MHz in Skynet SDFITS — convert to GHz.
    freq_ghz = float(prim.get("OBSFREQ", 1400.0)) / 1000.0

    # Beam FWHM from the diffraction formula θ = 1.22 · λ / D, NOT from
    # the SDFITS BMAJ keyword. C++ Survey::setSdfitsParams
    # (Survey.cpp:255/294/322/353) always computes it this way and
    # IGNORES BMAJ — and BMAJ in the source files is set by the Skynet
    # data-acquisition pipeline using a slightly different recipe
    # (BMAJ=0.66538 vs diffraction=0.67599 for GBT-20 @ 1550 MHz, a 1.6%
    # discrepancy). Since `resolution = pixel_size_bw * psf_fwhm` flows
    # into every per-pixel calculation downstream, that 1.6% mismatch
    # propagates as a uniform contraction of all painted/computed
    # features relative to the reference render. Matching the C++
    # formula is required for pixel-grid parity.
    #
    # Dish diameter D varies by telescope:
    #   20-meter (GBT-20m)   → 20.0 m
    #   GBT (NRAO_GBT)       → 105.0 m
    #   40-foot              → 12.192 m
    _DISH_DIAMETER_M = {
        Telescope.TWENTY_METER: 20.0,
        Telescope.FOURTY_FOOT: 12.192,
    }
    dish_d_m = _DISH_DIAMETER_M.get(telescope, 20.0)
    freq_hz = freq_ghz * 1e9
    if freq_hz > 0 and dish_d_m > 0:
        c_m_per_s = 299792458.0
        psf_fwhm = 1.22 * c_m_per_s * 180.0 / (freq_hz * dish_d_m * math.pi)
    else:
        # Last-resort fallback if frequency is missing.
        psf_fwhm = float(prim.get("BMAJ", 0.0))
        if psf_fwhm <= 0:
            band = _band_from_freq(freq_ghz)
            psf_fwhm = _BEAM_FWHM_DEG.get((telescope, band), 0.5)

    # Map type from OBSMODE
    obsmode = (prim.get("OBSMODE", "") or "").lower()
    if "daisy" in obsmode:
        map_type = MapType.DAISY
    elif "nod" in obsmode:
        map_type = MapType.NODDING
    else:
        map_type = MapType.RASTER  # ralongmap / declongmap / raster all map here

    # Split into Scans. Rasters/noddings have one SWPINDEX value per
    # scan; daisies have all SWPINDEX=0 and need to be split by petal
    # (direction reversal through the center).
    if map_type is MapType.DAISY:
        scans = _split_scans_daisy(
            utsecs, ra, dec, az, el, nsamps, cal, flux_l, flux_r,
            center_ra_deg=_parse_header_radec(prim)[0],
            center_dec_deg=_parse_header_radec(prim)[1],
        )
    else:
        scans = _split_scans_by_sweep(
            utsecs, ra, dec, az, el, nsamps, cal, swpidx, flux_l, flux_r
        )

    # Read commanded source position from primary header. For
    # moving-object observations (Jupiter etc.), the science samples
    # may not cover the source position itself — the telescope parks
    # at the source for the cal block then scans away. The commanded
    # position from the FITS header is the only reliable way to know
    # where the source actually is.
    target_ra_deg, target_dec_deg = _parse_header_radec(prim)

    survey = Survey(
        scans=scans,
        telescope=telescope,
        frequency_ghz=freq_ghz,
        mjd=float(prim.get("MJD", 0.0)),
        map_type=map_type,
        channel=channel,
        psf_fwhm=psf_fwhm,
        source_path=path,
        cal_time=cal_time,
        cal_state=cal_state,
        cal_dumps=cal_dumps,
        cal_flux_l=cal_flux_l,
        cal_flux_r=cal_flux_r,
    )
    # Stash the commanded center on the Survey for coordinates.project
    # to use as the projection origin.
    survey.target_ra_deg = target_ra_deg
    survey.target_dec_deg = target_dec_deg
    return survey


def _parse_header_radec(header) -> tuple[float, float]:
    """Parse the OBJECT's commanded RA/Dec from a Skynet FITS primary
    header. Skynet stores RA as 'HH:MM:SS.ss' (in hours) and Dec as
    'DD:MM:SS.ss' (in degrees). Returns (ra_deg, dec_deg) or (nan, nan)
    if either is missing/unparseable."""
    ra_str = header.get("RA")
    dec_str = header.get("DEC")
    if not ra_str or not dec_str:
        return float("nan"), float("nan")
    try:
        # RA: HH:MM:SS.s -> degrees
        h, m, s = (float(x) for x in str(ra_str).split(":"))
        ra_deg = (h + m / 60.0 + s / 3600.0) * 15.0
        # Dec: DD:MM:SS.s (preserve sign on the degrees field)
        parts = str(dec_str).split(":")
        sign = -1.0 if parts[0].strip().startswith("-") else 1.0
        d, m, s = (abs(float(x)) for x in parts)
        dec_deg = sign * (d + m / 60.0 + s / 3600.0)
        return ra_deg, dec_deg
    except (ValueError, IndexError):
        return float("nan"), float("nan")


def _parse_continuum_channels(header) -> tuple[int, int]:
    """Look in FITS HISTORY cards for 'START,STOP channels' and return
    (lo, hi) channel indices. Defaults to (0, 1024) if not found."""
    history = "\n".join(str(line) for line in header.get("HISTORY", []))
    import re
    m = re.search(r"channels\s+(\d+)\s*,\s*(\d+)", history, re.IGNORECASE)
    if m:
        return int(m.group(1)), int(m.group(2))
    return 0, 1023


def _split_scans_daisy(
    utsecs, ra, dec, az, el, nsamps, cal, flux_l, flux_r,
    center_ra_deg: float, center_dec_deg: float,
) -> list[Scan]:
    """Split a daisy observation into per-petal scans.

    Daisies have all SWPINDEX=0, so we can't use the sweep column. The
    natural split point is each pass through (or close to) the daisy
    center: the radial distance from center drops to a local minimum
    at every petal-to-petal transition.

    Algorithm: compute radial distance from (center_ra, center_dec) per
    sample; find local minima with sufficient separation (each minimum
    is one petal boundary). Split between consecutive minima.
    """
    if utsecs.size == 0:
        return []

    # Radial distance from the PETAL-PATTERN center, which is the
    # median of the science sample positions — NOT the header center.
    # For tracked moving targets (Jupiter, Mars, …), the header RA/Dec
    # is the commanded parking position and can be degrees off from
    # the actual petal center (e.g. 0152427: header (118.0, 21.83) vs
    # data center (115.2, 21.86) — a 2.8° offset that wrecked find_peaks
    # with min_prominence=0.25*max_r). For stationary targets they
    # agree to within seconds of arc, so this is safe in both cases.
    cx = float(np.median(ra))
    cy = float(np.median(dec))
    cos_dec = math.cos(math.radians(cy))
    dx = (ra - cx) * cos_dec
    dy = dec - cy
    radius = np.hypot(dx, dy)

    # Find local MAXIMA of the radius (petal tips), not minima of
    # radius (center crossings). C++ daisySweepBreaker
    # (Survey.cpp:1457-1492) defines each petal as tip-to-tip, so the
    # "center" sample (closest to the daisy origin) sits in the
    # MIDDLE of each scan — which is what `daisyAngleBuilder` needs to
    # produce a monotonic 1D angle coord for cross-correlation. If you
    # split at minima instead, each scan goes center→tip→center with
    # the center sample at the BOUNDARY (or boundaries), making the
    # signed-angle coord non-monotonic and breaking the auto-TS
    # algorithm.
    from scipy.signal import find_peaks
    n = radius.size
    max_r = float(radius.max())
    min_separation = max(int(n / 50), 3)
    min_prominence = 0.25 * max_r
    tip_indices, _ = find_peaks(radius,
                                distance=min_separation,
                                prominence=min_prominence)

    # The observation starts and ends at some non-tip position (the
    # dish was already moving when recording began). The dish typically
    # starts at a TIP though, so the LEADING segment [0, tips[0]) is
    # usually about half a petal and folds into the first tip-to-tip
    # scan; same for the TRAILING segment [tips[-1], n). Drop the
    # leading and trailing partials and keep only the tip-to-tip
    # segments. Result: scans = len(tips) - 1. For a 12-petal daisy
    # we detect 12 tips (one per petal) → 11 full inter-tip scans.
    # That under-counts by one; restore by appending the leading and
    # trailing partials as their own scans only if they have enough
    # samples to be useful.
    tips = list(tip_indices.tolist())
    boundaries = list(tips) if tips else [0, n]

    scans = []
    for i in range(len(boundaries) - 1):
        lo = int(boundaries[i])
        hi = int(boundaries[i + 1])
        if hi - lo < 5:
            continue
        scans.append(Scan(
            scan_index=len(scans),
            time=utsecs[lo:hi].copy(),
            ra=ra[lo:hi].copy(),
            dec=dec[lo:hi].copy(),
            elevation=el[lo:hi].copy(),
            azimuth=az[lo:hi].copy(),
            flux_l=flux_l[lo:hi].copy(),
            flux_r=flux_r[lo:hi].copy(),
            dumps=nsamps[lo:hi].copy(),
            cal_state=cal[lo:hi].copy(),
            scan_in_ra=True,  # arbitrary for daisies
        ))

    return scans


def _split_scans_by_sweep(
    utsecs, ra, dec, az, el, nsamps, cal, swpidx, flux_l, flux_r
) -> list[Scan]:
    """Split samples into Scans using the SWPINDEX column. Each contiguous
    run of identical sweep indices becomes one Scan. Calibration samples
    (cal != 0/1 ambiguity isn't a real issue here — we just include them
    in their natural scan position; Stage 2 will find them by CALSTATE).
    """
    if utsecs.size == 0:
        return []

    # Determine working axis (RA vs Dec) by total span
    ra_span = ra.max() - ra.min()
    dec_span = dec.max() - dec.min()
    scan_in_ra = ra_span >= dec_span

    # Find boundaries where SWPINDEX changes.
    changes = np.where(np.diff(swpidx) != 0)[0] + 1
    boundaries = np.concatenate(([0], changes, [swpidx.size]))

    scans = []
    for i in range(boundaries.size - 1):
        lo, hi = int(boundaries[i]), int(boundaries[i + 1])
        if hi - lo < 3:
            continue
        scans.append(Scan(
            scan_index=len(scans),
            time=utsecs[lo:hi].copy(),
            ra=ra[lo:hi].copy(),
            dec=dec[lo:hi].copy(),
            elevation=el[lo:hi].copy(),
            azimuth=az[lo:hi].copy(),
            flux_l=flux_l[lo:hi].copy(),
            flux_r=flux_r[lo:hi].copy(),
            dumps=nsamps[lo:hi].copy(),
            cal_state=cal[lo:hi].copy(),
            scan_in_ra=bool(scan_in_ra),
        ))

    return scans


def _split_scans(
    utsecs, ra, dec, az, el, nsamps, cal, flux, map_type: MapType
) -> list[Scan]:
    """Split a flat sample list into Scans on direction reversals."""
    if utsecs.size == 0:
        return []

    # Determine which axis is the scan direction
    ra_span = ra.max() - ra.min()
    dec_span = dec.max() - dec.min()
    scan_in_ra = ra_span >= dec_span

    primary = ra if scan_in_ra else dec
    secondary = dec if scan_in_ra else ra

    # Find direction-reversal indices: zero crossings of the second derivative
    # in the primary axis. Crude but workable for rasters.
    d = np.diff(primary)
    sign = np.sign(d)
    # Locations where sign flips
    flips = np.where(np.diff(sign) != 0)[0] + 1
    boundaries = np.concatenate(([0], flips, [primary.size]))

    scans = []
    for i in range(boundaries.size - 1):
        lo, hi = int(boundaries[i]), int(boundaries[i + 1])
        if hi - lo < 3:
            # Drop sub-3-sample scans (probably noise at turning points).
            continue

        # L/R flux: for now we assume single-channel data and use it as L
        # with R = 0. Real port needs the proper channel-split logic.
        scans.append(Scan(
            scan_index=len(scans),
            time=utsecs[lo:hi].copy(),
            ra=ra[lo:hi].copy(),
            dec=dec[lo:hi].copy(),
            elevation=el[lo:hi].copy(),
            azimuth=az[lo:hi].copy(),
            flux_l=flux[lo:hi].astype(np.float64).copy(),
            flux_r=flux[lo:hi].astype(np.float64).copy(),
            dumps=nsamps[lo:hi].astype(np.float64).copy(),
            cal_state=cal[lo:hi].copy(),
            scan_in_ra=bool(scan_in_ra),
        ))

    return scans


def write_fits(rc_map: Map, path: str, header_info: Optional[dict] = None) -> None:
    """Write a Map to a multi-extension FITS file.

    Layers: PRIMARY (flux), then SCALE, WEIGHT, WEIGHT2, CORRELATION,
    PATH as image extensions.
    """
    from astropy.io import fits

    hdr = fits.Header()
    hdr["CTYPE1"] = "RA---SFL"
    hdr["CTYPE2"] = "DEC--SFL"
    hdr["CRVAL1"] = rc_map.center_ra_deg
    hdr["CRVAL2"] = rc_map.center_dec_deg
    hdr["CRPIX1"] = rc_map.flux.shape[1] // 2 + 1
    hdr["CRPIX2"] = rc_map.flux.shape[0] // 2 + 1
    hdr["CDELT1"] = -rc_map.resolution
    hdr["CDELT2"] = rc_map.resolution
    hdr["ORIGIN"] = "rcpy (Skynet 2.0)"

    if header_info:
        for k, v in header_info.items():
            try:
                hdr[k] = v
            except ValueError:
                # Skip keys that aren't valid FITS
                pass

    hdus = [
        fits.PrimaryHDU(data=rc_map.flux, header=hdr),
        fits.ImageHDU(data=rc_map.scale, name="SCALE"),
        fits.ImageHDU(data=rc_map.weight, name="WEIGHT"),
        fits.ImageHDU(data=rc_map.weight_corr, name="WEIGHT2"),
        fits.ImageHDU(data=rc_map.correlation, name="CORRELATION"),
        fits.ImageHDU(data=rc_map.path, name="PATH"),
    ]
    if rc_map.raw_flux is not None:
        hdus.append(fits.ImageHDU(data=rc_map.raw_flux, name="RAW"))

    fits.HDUList(hdus).writeto(path, overwrite=True)
