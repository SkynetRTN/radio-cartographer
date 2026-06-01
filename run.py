"""rcpy runner — RC Job parameters from the Skynet job UI.

Usage
-----
Drop one or more `.fits` files into `input/`, then::

    python run.py                           # process every .fits in input/
    python run.py input/cyg_a.fits          # process a single file
    python run.py input/a.fits input/b.fits # multiple files

Each input file produces in `output/`:
    <basename>_rcpy.fits   — processed multi-extension FITS
    <basename>_rcpy.png    — diagnostic 6-panel image
                              (flux, scale, weight, correlation, path, raw)

Job parameters
--------------
These match the "RC Job" settings tab in Skynet (see the screenshot the
lab is running against):

    Channel            : sum (= COMPOSITE, 0.5*(L+R))
    Gain Calibration   : interpolated
    Image Coordinate   : equatorial
    Time-Delay         : custom 1.0 s
    Include Raw Image  : yes
    Min Frequency      : 1355.0 MHz
    Max Frequency      : 1435.0 MHz
    BG Subtraction     : 6.0 BW
    RFI Preset         : Bright Target with Airy Rings
    RFI Scale          : 0.35 BW
    Surface Model      : 3rd-order 2D polynomial with noise prior (M10)
    Weight Scale       : 1/3 BW (measurement-quality image)

To change any setting, edit JOB below. The defaults reproduce the
"RC Job" preset on a bright continuum target (Cyg A, Cas A, Jupiter).
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt

from rcpy import io, pipeline
from rcpy.pipeline import PipelineConfig
from rcpy.types import CalMethod, Channel


# -----------------------------------------------------------------------------
# RC Job parameters (edit here to change the job)
# -----------------------------------------------------------------------------
@dataclass
class Job:
    # Display name only — written to the FITS JOBNAME header card.
    name: str = "RC Job"

    # Which polarization channel to image. Skynet's "sum" preset is
    # COMPOSITE (= 0.5*(L+R), the standard total-intensity image).
    #   Channel.LEFT       L-pol only
    #   Channel.RIGHT      R-pol only
    #   Channel.COMPOSITE  0.5*(L+R)              ← "sum" / recommended
    channel: Channel = Channel.COMPOSITE

    # How the noise-diode delta is applied to the science scans.
    #   CalMethod.PRE          use the FIRST cal block's delta everywhere
    #   CalMethod.POST         use the LAST cal block's delta everywhere
    #   CalMethod.INTERPOLATED linearly interpolate between cal blocks  ← Skynet default
    #   CalMethod.CONSTANT     mean(start, end) — use when INTERPOLATED
    #                          produces a spurious cal drift (bright
    #                          source contamination; see Cyg A notes)
    #   CalMethod.NONE         no calibration (raw counts)
    cal_method: CalMethod = CalMethod.INTERPOLATED

    # How to determine the encoder-to-data time shift.
    #   "off"     no shift (time_shift_seconds ignored)
    #   "auto"    cross-correlate adjacent scans for rasters; brute-
    #             force gridder sweep for daisies
    #   "custom"  apply the value in `time_shift_seconds` directly
    time_shift_mode: str = "custom"

    # Time shift to apply when `time_shift_mode == "custom"` (seconds).
    # Typical raster values are around -1.5 to -2 s. Daisies typically
    # want +1 s (opposite sign convention from rasters — see
    # `rcpy_daisy_handling` memory note for context).
    time_shift_seconds: float = 1.0

    # If True, runs a second BG/RFI/timeshift-disabled pass to produce
    # the "raw" image alongside the processed one. Adds about 30-60 s
    # per file. Matches Skynet's "Include Raw Image" checkbox.
    include_raw: bool = True

    # Spectral inclusion band, in MHz. Skynet reads its DATA column on
    # this channel range — used as the FITS RCMINFQ/RCMAXFQ header
    # cards. Defaults are the L-band CONTINUUM band.
    min_freq_mhz: float = 1355.0
    max_freq_mhz: float = 1435.0

    # Background-subtraction window size, in beamwidths. The fitter
    # subtracts a 4th-order polynomial fit within ±bg_scale_bw of each
    # sample. Larger values → smoother BG removal (more source flux
    # preserved at low spatial frequencies); smaller → tighter
    # baseline. 6.0 is the Skynet default for continuum maps; set to
    # 0.0 to disable BG subtraction entirely.
    bg_scale_bw: float = 6.0

    # RFI detection window size, in beamwidths. The RFI step builds
    # local models within this radius and rejects samples that deviate
    # too far. Skynet presets:
    #   0.35  "Bright Target with Airy Rings"  ← Cyg A / Cas A / Jupiter
    #   0.70  "Faint Target" (default)
    #   0.0   "No RFI"  (skip RFI rejection entirely)
    rfi_scale_bw: float = 0.35

    # Per-pixel surface-fit weighting scale, in beamwidths. Controls
    # how much the weighted least-squares fit smooths the source.
    # Skynet presets:
    #   1/3   "Measurement-quality image"  ← default; sharpest
    #   2/3   "Visualization-quality image" — smoother
    #   1.0   maximum smoothing
    weight_scale_bw: float = 1.0 / 3.0

    # Output pixel size, in beamwidths. 0.05 BW per pixel (= 20 px per
    # beam) is standard; smaller → finer grid but slower gridder;
    # larger → coarser. The Skynet UI exposes this only indirectly.
    pixel_size_bw: float = 0.05

    # Turning-edge trim, in beamwidths. Marks samples within this
    # distance of each scan's turning points (where the dish was
    # decelerating) as excluded. Skynet reference renders use 0.0
    # (rely on the edge-polygon NaN-out at the gridder stage instead).
    trim_size_bw: float = 0.0

    # Auto-centroiding S/N threshold for the RFI step (paper
    # Footnote 21). Lower → more samples flagged as potential point
    # sources (and protected from RFI rejection); higher → fewer.
    #   15    "Faint Target" (recommended for normal sources)
    #   75    "Bright Target with Airy Rings"  ← keeps the Airy
    #         rings of bright continuum sources from being eaten
    #         as RFI; matches Skynet's bright-target preset
    centroid_sigma: float = 75.0


JOB = Job()


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
ROOT = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(ROOT, "input")
OUTPUT_DIR = os.path.join(ROOT, "output")


def _ensure_dirs() -> None:
    os.makedirs(INPUT_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)


def _collect_inputs(argv: list[str]) -> list[str]:
    if argv:
        files = [os.path.abspath(p) for p in argv]
        missing = [p for p in files if not os.path.exists(p)]
        if missing:
            raise SystemExit(f"Input file(s) not found: {missing}")
        return files
    # Auto-discover from input/
    if not os.path.isdir(INPUT_DIR):
        raise SystemExit(f"No input/ directory and no files given.")
    files = sorted(
        os.path.join(INPUT_DIR, f) for f in os.listdir(INPUT_DIR)
        if f.lower().endswith(".fits")
    )
    if not files:
        raise SystemExit(
            f"No .fits files in {INPUT_DIR}. "
            f"Drop your input files there or pass paths as arguments."
        )
    return files


def _build_config(job: Job, raw_pass: bool) -> PipelineConfig:
    """Build a PipelineConfig from a Job. raw_pass=True returns the
    config for the raw-image render (BG/RFI/timeshift disabled), to
    match the "Include Raw Image" output that Skynet writes."""
    if raw_pass:
        return PipelineConfig(
            bg_scale_bw=0.0,
            rfi_scale_bw=job.rfi_scale_bw,
            weight_scale_bw=0.0,
            pixel_size_bw=job.pixel_size_bw,
            trim_size_bw=0.0,
            skip_edge_trim=True,
            skip_bg=True,
            skip_timeshift=True,
            skip_rfi=True,
            cal_method=job.cal_method,
        )
    return PipelineConfig(
        bg_scale_bw=job.bg_scale_bw,
        rfi_scale_bw=job.rfi_scale_bw,
        weight_scale_bw=job.weight_scale_bw,
        pixel_size_bw=job.pixel_size_bw,
        trim_size_bw=job.trim_size_bw,
        skip_edge_trim=(job.trim_size_bw == 0.0),
        time_shift_mode=job.time_shift_mode,
        time_shift_seconds=job.time_shift_seconds,
        centroid_sigma=job.centroid_sigma,
        cal_method=job.cal_method,
    )


def _job_header_info(job: Job) -> dict[str, str]:
    """FITS header cards recording the job settings."""
    return {
        "JOBNAME":   job.name,
        "RCCHAN":    job.channel.value,
        "RCCALI":    job.cal_method.value,
        "RCTS":      job.time_shift_mode,
        "RCTMSFT":   f"{job.time_shift_seconds}",
        "RCMINFQ":   f"{job.min_freq_mhz}",
        "RCMAXFQ":   f"{job.max_freq_mhz}",
        "RCBG":      f"{job.bg_scale_bw}",
        "RCRFI":     f"{job.rfi_scale_bw}",
        "RCWGT":     f"{job.weight_scale_bw}",
        "RCRAW":     "1" if job.include_raw else "0",
    }


def _save_diagnostic_png(rc_map, raw_flux, out_path: str, title: str) -> None:
    """6-panel diagnostic image: flux, scale, weight, correlation,
    path, raw."""
    layers = [
        ("flux",        rc_map.flux,        "Jy/beam"),
        ("scale",       rc_map.scale,       "BW"),
        ("weight",      rc_map.weight,      "rel."),
        ("correlation", rc_map.correlation, "rel."),
        ("path",        rc_map.path,        "code"),
        ("raw",         raw_flux,           "Jy/beam"),
    ]
    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    ax = ax.flatten()
    for a, (name, arr, unit) in zip(ax, layers):
        if arr is None:
            a.axis("off")
            continue
        finite = np.isfinite(arr)
        if not finite.any():
            a.set_title(f"{name} (no data)"); a.axis("off"); continue
        vmin = np.nanpercentile(arr, 2)
        vmax = np.nanpercentile(arr, 99)
        im = a.imshow(arr, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax)
        a.set_title(name)
        plt.colorbar(im, ax=a, fraction=0.046, pad=0.04, label=unit)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


def _process_one(fits_path: str, job: Job) -> None:
    name = os.path.splitext(os.path.basename(fits_path))[0]
    print(f"\n=== {name} ===")
    print(f"  reading: {fits_path}")

    survey = io.read_sdfits(fits_path)
    survey.channel = job.channel
    print(f"  map_type={survey.map_type.value}  "
          f"scans={len(survey.scans)}  "
          f"psf_fwhm={survey.psf_fwhm:.4f} deg")

    cfg = _build_config(job, raw_pass=False)
    print(f"  processing (BG={cfg.bg_scale_bw} RFI={cfg.rfi_scale_bw} "
          f"weight={cfg.weight_scale_bw:.4f} TS={job.time_shift_mode}:{job.time_shift_seconds}s)")
    rc = pipeline.process(survey, cfg)

    raw_flux = None
    if job.include_raw:
        print("  rendering raw image")
        survey_raw = io.read_sdfits(fits_path)
        survey_raw.channel = job.channel
        cfg_raw = _build_config(job, raw_pass=True)
        rc_raw = pipeline.process(survey_raw, cfg_raw)
        raw_flux = rc_raw.flux
        rc.raw_flux = raw_flux

    out_fits = os.path.join(OUTPUT_DIR, f"{name}_rcpy.fits")
    out_png = os.path.join(OUTPUT_DIR, f"{name}_rcpy.png")
    io.write_fits(rc, out_fits, header_info=_job_header_info(job))
    _save_diagnostic_png(rc, raw_flux, out_png, f"{name} — {job.name}")
    print(f"  wrote {out_fits}")
    print(f"  wrote {out_png}")


def main(argv: list[str] | None = None) -> int:
    argv = list(argv if argv is not None else sys.argv[1:])
    _ensure_dirs()
    files = _collect_inputs(argv)
    print(f"Job: {JOB.name}  ({len(files)} input file{'s' if len(files) != 1 else ''})")
    for f in files:
        _process_one(f, JOB)
    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
