"""Synthetic Survey generators for testing.

Generates a raster scan with configurable noise, Gaussian point sources,
en-route drift, and a calibration block at the start and end. Used by
tests so we have a runnable end-to-end pipeline without needing real
SDFITS files.

The geometry roughly matches a 20-meter L-band raster: ~0.55 deg beam,
0.2-beamwidth sampling, ~24-beamwidth-wide map.
"""
from __future__ import annotations

import numpy as np

from rcpy.types import (
    Survey, Scan, Telescope, MapType, Channel, CalMethod, Coordinates,
)


def synthetic_raster(
    n_scans: int = 60,
    samples_per_scan: int = 120,
    map_width_deg: float = 13.0,   # ~24 BW at 0.55 deg
    psf_fwhm_deg: float = 0.55,
    center_ra_deg: float = 180.0,
    center_dec_deg: float = 0.0,
    noise_amplitude: float = 0.001,
    seed: int = 0,
    sources: list[tuple[float, float, float]] | None = None,
    include_cal_blocks: bool = True,
) -> Survey:
    """Construct a synthetic 20-meter raster Survey.

    Args:
        n_scans: number of horizontal scans.
        samples_per_scan: samples per scan.
        map_width_deg: full width of the mapped region in degrees.
        sources: list of (ra_offset_deg, dec_offset_deg, peak_amplitude)
                 tuples for point sources to inject.
    """
    rng = np.random.default_rng(seed)

    half = map_width_deg / 2.0
    scan_dec_centers = np.linspace(
        center_dec_deg - half, center_dec_deg + half, n_scans
    )

    # Default sources: one bright at center, one faint offset.
    if sources is None:
        sources = [(0.0, 0.0, 1.0), (3.0, 2.0, 0.1)]

    scans = []
    t_cursor = 0.0
    sample_period = 0.1  # seconds

    # Optional starting cal block (5s on, 5s off, simulated)
    if include_cal_blocks:
        t_cursor = _emit_cal_block(scans, t_cursor, sample_period, center_ra_deg,
                                   center_dec_deg, base_level=2.0, diode_jump=1.0,
                                   noise_amplitude=noise_amplitude, rng=rng)

    for i, dec in enumerate(scan_dec_centers):
        # Alternate direction
        ra_lo = center_ra_deg - half
        ra_hi = center_ra_deg + half
        if i % 2 == 1:
            ra_lo, ra_hi = ra_hi, ra_lo

        ras = np.linspace(ra_lo, ra_hi, samples_per_scan)
        decs = np.full(samples_per_scan, dec)
        times = t_cursor + np.arange(samples_per_scan) * sample_period
        t_cursor = times[-1] + sample_period

        # Flux: base + sources + noise + slow drift
        flux = np.full(samples_per_scan, 2.0)  # base level (will get calibrated out)
        flux += rng.normal(0.0, noise_amplitude, samples_per_scan)
        flux += 0.01 * np.sin(times / 30.0)  # en-route drift

        # Sources
        for src_ra_off, src_dec_off, peak in sources:
            src_ra = center_ra_deg + src_ra_off * psf_fwhm_deg
            src_dec = center_dec_deg + src_dec_off * psf_fwhm_deg
            sigma = psf_fwhm_deg / 2.355
            distance_sq = (ras - src_ra) ** 2 + (decs - src_dec) ** 2
            flux += peak * np.exp(-distance_sq / (2 * sigma ** 2))

        scans.append(Scan(
            scan_index=len(scans),
            time=times,
            ra=ras,
            dec=decs,
            elevation=np.full(samples_per_scan, 45.0),
            azimuth=np.full(samples_per_scan, 180.0),
            flux_l=flux.copy(),
            flux_r=flux.copy(),
            dumps=np.ones(samples_per_scan),
            cal_state=np.full(samples_per_scan, -1, dtype=np.int32),
            scan_in_ra=True,
        ))

    if include_cal_blocks:
        t_cursor = _emit_cal_block(scans, t_cursor, sample_period, center_ra_deg,
                                   center_dec_deg, base_level=2.0, diode_jump=1.0,
                                   noise_amplitude=noise_amplitude, rng=rng)

    # Re-index scans so they're contiguous from 0
    for i, s in enumerate(scans):
        s.scan_index = i

    return Survey(
        scans=scans,
        telescope=Telescope.TWENTY_METER,
        frequency_ghz=1.4,
        mjd=58000.0,
        map_type=MapType.RASTER,
        channel=Channel.COMPOSITE,
        cal_method=CalMethod.INTERPOLATED,
        psf_fwhm=psf_fwhm_deg,
    )


def _emit_cal_block(
    scans, t_cursor, sample_period, ra, dec, base_level, diode_jump,
    noise_amplitude, rng
):
    """Append two scans of cal-on / cal-off alternations to the scan list."""
    n_per_state = 30
    states_pattern = np.tile([1, 0], n_per_state)  # alternating
    n_total = states_pattern.size

    times = t_cursor + np.arange(n_total) * sample_period
    cal = states_pattern.copy()
    flux = np.full(n_total, base_level + 0.0)
    flux[cal == 1] += diode_jump
    flux += rng.normal(0.0, noise_amplitude, n_total)

    scans.append(Scan(
        scan_index=len(scans),
        time=times,
        ra=np.full(n_total, ra),
        dec=np.full(n_total, dec),
        elevation=np.full(n_total, 45.0),
        azimuth=np.full(n_total, 180.0),
        flux_l=flux.copy(),
        flux_r=flux.copy(),
        dumps=np.ones(n_total),
        cal_state=cal,
        scan_in_ra=True,
    ))
    return times[-1] + sample_period
