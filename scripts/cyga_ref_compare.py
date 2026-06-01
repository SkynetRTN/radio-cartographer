"""Direct pixel-by-pixel comparison between the rcpy pipeline output for
Cygnus A (0144716.fits) and the reference C++ render (RC Job 7983).

Reference config (read from FITS headers):
    BGSCALE  = 6.0
    RFISCALE = 0.35
    SMSCALE  = 0.333333  (= 1/3)
    PIXLSIZE = 0.05      (BW)
    BEAM     = 0.675992  (deg)
    TIMESHIF = -1.12321  (sec)  — C++ auto time-shift result
    CENTERRA = 299.852, CENTERDE = 40.7338
    RCCHAN   = composite (= 0.5*(L+R))

Reference layers (134 x 148):
    main, path, scale, weight, correlation, raw

The rcpy pipeline currently produces (per Map):
    flux, weight, weight_corr, scale, correlation, path

Layer-pair mapping for comparison:
    rcpy.flux        <->  ref.main
    rcpy.weight      <->  ref.weight
    rcpy.scale       <->  ref.scale
    rcpy.correlation <->  ref.correlation
    rcpy.path        <->  ref.path

The two maps have slightly different bounding boxes / grids, so we
reproject the rcpy output onto the reference WCS using astropy.wcs
and bilinear interpolation (scipy.ndimage.map_coordinates).
"""
from __future__ import annotations

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import map_coordinates

# Ensure the rcpy package is importable when running from radio-cartographer/.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from rcpy import io, pipeline               # noqa: E402
from rcpy.pipeline import PipelineConfig    # noqa: E402

REF_FITS = os.path.join(ROOT, "data", "cyg a_RC Job_7983_0007673.fits")
SRC_FITS = os.path.join(ROOT, "data", "0144716.fits")
OUT_PNG  = os.path.join(ROOT, "data", "outputs", "cyga_compare.png")


def _map_to_array(rc_map, layer: str) -> np.ndarray:
    return {
        "flux":        rc_map.flux,
        "weight":      rc_map.weight,
        "weight_corr": rc_map.weight_corr,
        "scale":       rc_map.scale,
        "correlation": rc_map.correlation,
        "path":        rc_map.path,
    }[layer]


def _peak_offset(src: np.ndarray, ref: np.ndarray) -> tuple[float, float]:
    """Find the sub-pixel array-index offset (dy, dx) that aligns the
    SOURCE peak in `src` with the SOURCE peak in `ref`.

    Both pipelines' WCS metadata is inconsistent with where the source
    actually sits in the array (each anchors CRPIX at its own user-
    specified or array-center pixel), so a direct WCS lookup gives
    misleading sky coordinates. We work in raw array coordinates and
    align by the source peak instead.

    A 3x3 paraboloid fit around the integer-pixel argmax gives a
    sub-pixel refinement of the peak. The offset is src_peak - ref_peak,
    so subtracting it from src array indices maps src onto ref.
    """
    def _refine(a):
        iy, ix = np.unravel_index(np.nanargmax(a), a.shape)
        if iy == 0 or ix == 0 or iy >= a.shape[0] - 1 or ix >= a.shape[1] - 1:
            return float(iy), float(ix)
        # 1-D parabola fits in x and y
        z = a[iy, ix - 1:ix + 2]
        dx = 0.5 * (z[0] - z[2]) / (z[0] - 2 * z[1] + z[2] + 1e-30)
        z = a[iy - 1:iy + 2, ix]
        dy = 0.5 * (z[0] - z[2]) / (z[0] - 2 * z[1] + z[2] + 1e-30)
        return float(iy) + dy, float(ix) + dx

    sy, sx = _refine(src)
    ry, rx = _refine(ref)
    return sy - ry, sx - rx


def _resample_to_ref(rcpy_layer: np.ndarray,
                     peak_offset: tuple[float, float],
                     ref_shape: tuple[int, int]) -> np.ndarray:
    """Bilinear-sample rcpy_layer onto the reference array grid, after
    subtracting peak_offset so the source peaks coincide. Returns NaN
    where the reference pixel falls outside the rcpy footprint or where
    any of the 4 enclosing source pixels were NaN."""
    dy, dx = peak_offset
    ny_r, nx_r = ref_shape
    ix_r, iy_r = np.meshgrid(np.arange(nx_r), np.arange(ny_r))
    py_r = iy_r + dy
    px_r = ix_r + dx

    coords = np.stack([py_r, px_r])
    finite_in = np.isfinite(rcpy_layer)
    fill = rcpy_layer.copy()
    fill[~finite_in] = 0.0
    val = map_coordinates(fill, coords, order=1, mode="constant", cval=0.0)

    iy0 = np.clip(np.floor(py_r).astype(int), 0, finite_in.shape[0] - 1)
    ix0 = np.clip(np.floor(px_r).astype(int), 0, finite_in.shape[1] - 1)
    iy1 = np.clip(iy0 + 1, 0, finite_in.shape[0] - 1)
    ix1 = np.clip(ix0 + 1, 0, finite_in.shape[1] - 1)
    any_nan = ~(finite_in[iy0, ix0] & finite_in[iy0, ix1]
                & finite_in[iy1, ix0] & finite_in[iy1, ix1])
    out = np.where(any_nan, np.nan, val)
    oob = (px_r < 0) | (px_r > finite_in.shape[1] - 1) | \
          (py_r < 0) | (py_r > finite_in.shape[0] - 1)
    out[oob] = np.nan
    return out


def _summarise(label: str, a: np.ndarray, b: np.ndarray) -> dict:
    """Report stats on a vs b. a = rcpy resampled, b = reference."""
    both = np.isfinite(a) & np.isfinite(b)
    n = int(both.sum())
    if n == 0:
        return {"n": 0}
    resid = a[both] - b[both]
    # Relative residual where b is large enough that it dominates noise
    b_in = b[both]
    big = np.abs(b_in) > 0.05 * np.nanmax(np.abs(b_in))
    rel = np.where(big, np.abs(resid) / np.where(big, np.abs(b_in), 1.0), np.nan)
    out = {
        "n": n,
        "ref_min": float(np.nanmin(b)),
        "ref_max": float(np.nanmax(b)),
        "ref_med": float(np.nanmedian(b)),
        "src_min": float(np.nanmin(a)),
        "src_max": float(np.nanmax(a)),
        "src_med": float(np.nanmedian(a)),
        "resid_med": float(np.median(resid)),
        "resid_rms": float(np.sqrt(np.mean(resid ** 2))),
        "resid_abs_med": float(np.median(np.abs(resid))),
        "rel_abs_med": float(np.nanmedian(np.abs(rel))) if rel.size else float("nan"),
        "nan_src_unique": int(np.isfinite(b)[~np.isfinite(a)].sum()),
        "nan_ref_unique": int(np.isfinite(a)[~np.isfinite(b)].sum()),
    }
    print(f"\n--- {label} ---")
    print(f"  pixels compared : {n}")
    print(f"  ref     min/med/max : {out['ref_min']:.4g} / {out['ref_med']:.4g} / {out['ref_max']:.4g}")
    print(f"  rcpy    min/med/max : {out['src_min']:.4g} / {out['src_med']:.4g} / {out['src_max']:.4g}")
    print(f"  resid   med/rms     : {out['resid_med']:+.4g} / {out['resid_rms']:.4g}")
    print(f"  |resid| median      : {out['resid_abs_med']:.4g}")
    print(f"  relative |resid| (where ref >5% of peak): median = {out['rel_abs_med']:.4g}")
    print(f"  NaN-only-in-rcpy: {out['nan_src_unique']}   NaN-only-in-ref: {out['nan_ref_unique']}")
    return out


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", action="store_true",
                        help="Compare against ref HDU[5] ('raw') layer for flux; "
                             "run our pipeline in raw mode (no BG, no RFI, no "
                             "timeshift, theta_min=0, no noise prior). Other "
                             "layers still compare main-mode (no raw equivalent).")
    args = parser.parse_args()

    # --- Reference ---
    f_ref = fits.open(REF_FITS)
    ref_shape = f_ref[0].data.shape

    if args.raw:
        print("RAW MODE: comparing flux against ref HDU['raw']")
        ref_layers = {
            "flux":        f_ref["raw"].data,
            "path":        f_ref["path"].data,
            "scale":       f_ref["scale"].data,
            "weight":      f_ref["weight"].data,
            "correlation": f_ref["correlation"].data,
        }
    else:
        ref_layers = {
            "flux":        f_ref["main"].data,
            "path":        f_ref["path"].data,
            "scale":       f_ref["scale"].data,
            "weight":      f_ref["weight"].data,
            "correlation": f_ref["correlation"].data,
        }
    # The raw HDU is always available; include it as an extra row in the
    # diagnostic so the eye can spot pre-BG/RFI distortions independently
    # of the processed-map distortions.
    ref_raw = f_ref["raw"].data

    # --- rcpy run ---
    # Reference config from the FITS header:
    #   BGSCALE  = 6.0     SMSCALE = 0.333333    RFISCALE = 0.35
    #   PIXLSIZE = 0.05    RCTRIM  = 0.0         RCRAW = 1
    #   RCCHAN   = composite  RCTS = auto  TIMESHIF = -1.12321
    survey = io.read_sdfits(SRC_FITS)
    if args.raw:
        # C++ raw render disables BG, RFI, timeshift; uses theta_min=0.
        cfg = PipelineConfig(
            bg_scale_bw=0.0,
            rfi_scale_bw=0.35,
            weight_scale_bw=0.0,
            pixel_size_bw=0.05,
            trim_size_bw=0.0,
            skip_edge_trim=True,
            skip_bg=True,
            skip_timeshift=True,
            skip_rfi=True,
        )
    else:
        cfg = PipelineConfig(
            bg_scale_bw=6.0,
            rfi_scale_bw=0.35,
            weight_scale_bw=1.0 / 3.0,
            pixel_size_bw=0.05,
            trim_size_bw=0.0,
            skip_edge_trim=True,
        )
    rc = pipeline.process(survey, cfg)
    print("rcpy map shape    :", rc.flux.shape)
    print("reference shape   :", ref_shape)

    # Second pipeline pass for the "raw" panel — match C++'s raw render
    # mode (skip BG/RFI/timeshift). We always run this so the raw row
    # in the diagnostic shows ref vs rcpy directly even when args.raw is
    # False. The processed-comparison rows come from `rc`; the raw row
    # uses `rc_raw`.
    survey_raw = io.read_sdfits(SRC_FITS)
    cfg_raw = PipelineConfig(
        bg_scale_bw=0.0, rfi_scale_bw=0.35, weight_scale_bw=0.0,
        pixel_size_bw=0.05, trim_size_bw=0.0,
        skip_edge_trim=True, skip_bg=True, skip_timeshift=True,
        skip_rfi=True,
    )
    rc_raw = pipeline.process(survey_raw, cfg_raw)
    print("rcpy raw shape    :", rc_raw.flux.shape)

    # Align by source-peak position. Both pipelines write a CRPIX that
    # doesn't reliably point at the source, so we ignore WCS metadata
    # entirely and align maps in raw array index space using the flux
    # peak as the registration point.
    peak_offset = _peak_offset(rc.flux, ref_layers["flux"])
    print(f"peak offset (dy, dx) src - ref: ({peak_offset[0]:+.2f}, {peak_offset[1]:+.2f}) pixels")

    # --- Per-layer compare ---
    # All layers use the same flux peak-offset alignment. Only the flux
    # layer has a strong registration signal; using its offset for the
    # other layers is fine in main mode (they all share the same pixel
    # grid). For raw mode the raw HDU has a different shape (134x142 vs
    # the main HDUs' 134x148) — in that case we use the raw layer's own
    # shape for the flux compare and the main shape for the others.
    summaries = {}
    src_resampled = {}
    for layer in ("flux", "weight", "scale", "correlation", "path"):
        src_arr = _map_to_array(rc, layer)
        layer_ref_shape = ref_layers[layer].shape
        if args.raw and layer == "flux":
            # The raw HDU shape differs from main HDUs and has its own
            # source peak position in its own pixel grid.
            try:
                layer_peak_offset = _peak_offset(src_arr, ref_layers[layer])
            except Exception:
                layer_peak_offset = peak_offset
        else:
            layer_peak_offset = peak_offset
        re = _resample_to_ref(src_arr, layer_peak_offset, layer_ref_shape)
        src_resampled[layer] = re
        summaries[layer] = _summarise(layer, re, ref_layers[layer])

    # Raw row: peak-align rc_raw.flux to ref_raw, resample onto ref_raw's grid.
    try:
        raw_peak_offset = _peak_offset(rc_raw.flux, ref_raw)
    except Exception:
        raw_peak_offset = peak_offset
    raw_resampled = _resample_to_ref(rc_raw.flux, raw_peak_offset, ref_raw.shape)
    summaries["raw"] = _summarise("raw", raw_resampled, ref_raw)

    # --- Render side-by-side diagnostic ---
    panel_layers = ("flux", "weight", "scale", "correlation", "path", "raw")
    fig, axes = plt.subplots(len(panel_layers), 3, figsize=(13, 21))
    for i, layer in enumerate(panel_layers):
        if layer == "raw":
            ref = ref_raw
            src = raw_resampled
        else:
            ref = ref_layers[layer]
            src = src_resampled[layer]
        resid = src - ref
        cmap = "viridis" if layer != "path" else "gray"
        # Symmetric limits for residual
        finite_resid = resid[np.isfinite(resid)]
        if finite_resid.size:
            r_lim = float(np.nanpercentile(np.abs(finite_resid), 99))
        else:
            r_lim = 1.0
        # Match colour scales on ref/src
        finite_both = np.concatenate([ref[np.isfinite(ref)].ravel(),
                                       src[np.isfinite(src)].ravel()])
        vmin = float(np.nanpercentile(finite_both, 1)) if finite_both.size else 0
        vmax = float(np.nanpercentile(finite_both, 99)) if finite_both.size else 1
        axes[i, 0].imshow(ref, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
        axes[i, 0].set_title(f"reference {layer}")
        axes[i, 1].imshow(src, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
        axes[i, 1].set_title(f"rcpy {layer}")
        axes[i, 2].imshow(resid, origin="lower", cmap="RdBu_r",
                           vmin=-r_lim, vmax=r_lim)
        axes[i, 2].set_title(f"rcpy - ref  ({layer})")
        for ax in axes[i]:
            ax.set_xticks([]); ax.set_yticks([])

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=120)
    print(f"\nSaved {OUT_PNG}")


if __name__ == "__main__":
    main()
