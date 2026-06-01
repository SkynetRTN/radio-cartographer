"""Compare our pipeline against the Jupiter reference render.

Reference RC Job 7978 settings (from FITS header + UI screenshot):
    RCBG     = 6.0          (1D background subtraction scale, BW)
    RCRFI    = 0.7          (2D RFI subtraction scale, "Faint Target" preset)
    RCWGT    = 0.333333     (surface model weighting scale, 1/3 BW)
    PIXLSIZE = 0.05         (pixel scale, BW)
    BEAM     = 0.675992     (beam FWHM, deg)
    RCPHOT   = 0            (photometry off)
    RCTS     = custom       (time-shift mode)
    TIMESHIF = -1.807       (time-shift, seconds)
    RCSFM    = brightest    (centroid method)
    RCCHAN   = composite    (channel = 0.5*(L+R))

Surface model: "3rd-Order 2D Polynomial" — NO noise prior, unlike the
cyga reference. This means m10_plus_processing should be set to False.
"""
from __future__ import annotations

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import map_coordinates

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from rcpy import io, pipeline                                      # noqa: E402
from rcpy.pipeline import PipelineConfig                           # noqa: E402

REF_FITS = os.path.join(ROOT, "data", "jupiter_RC Job - Chelsea_7978_0007668.fits")
SRC_FITS = os.path.join(ROOT, "data", "0144954.fits")
OUT_PNG  = os.path.join(ROOT, "data", "outputs", "jupiter_compare.png")


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
    def _refine(a):
        flat = a.copy()
        flat[~np.isfinite(flat)] = -np.inf
        iy, ix = np.unravel_index(np.argmax(flat), flat.shape)
        if iy == 0 or ix == 0 or iy >= a.shape[0] - 1 or ix >= a.shape[1] - 1:
            return float(iy), float(ix)
        z = flat[iy, ix - 1:ix + 2]
        dx = 0.5 * (z[0] - z[2]) / (z[0] - 2 * z[1] + z[2] + 1e-30)
        z = flat[iy - 1:iy + 2, ix]
        dy = 0.5 * (z[0] - z[2]) / (z[0] - 2 * z[1] + z[2] + 1e-30)
        return float(iy) + dy, float(ix) + dx
    sy, sx = _refine(src)
    ry, rx = _refine(ref)
    return sy - ry, sx - rx


def _resample_to_ref(rcpy_layer, peak_offset, ref_shape):
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


def _summarise(label, a, b):
    both = np.isfinite(a) & np.isfinite(b)
    n = int(both.sum())
    if n == 0:
        return {"n": 0}
    resid = a[both] - b[both]
    b_in = b[both]
    big = np.abs(b_in) > 0.05 * np.nanmax(np.abs(b_in))
    rel = np.where(big, np.abs(resid) / np.where(big, np.abs(b_in), 1.0), np.nan)
    out = {
        "ref_min": float(np.nanmin(b)), "ref_max": float(np.nanmax(b)),
        "ref_med": float(np.nanmedian(b)),
        "src_min": float(np.nanmin(a)), "src_max": float(np.nanmax(a)),
        "src_med": float(np.nanmedian(a)),
        "resid_med": float(np.median(resid)),
        "resid_rms": float(np.sqrt(np.mean(resid ** 2))),
        "resid_abs_med": float(np.median(np.abs(resid))),
        "rel_abs_med": float(np.nanmedian(np.abs(rel))) if rel.size else float("nan"),
        "nan_src_unique": int(np.isfinite(b)[~np.isfinite(a)].sum()),
        "nan_ref_unique": int(np.isfinite(a)[~np.isfinite(b)].sum()),
        "n": n,
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
    f_ref = fits.open(REF_FITS)
    ref_shape = f_ref[0].data.shape

    ref_layers = {
        "flux":        f_ref["main"].data,
        "path":        f_ref["path"].data,
        "scale":       f_ref["scale"].data,
        "weight":      f_ref["weight"].data,
        "correlation": f_ref["correlation"].data,
    }
    ref_raw = f_ref["raw"].data  # extra "raw" panel for the eye check

    # Jupiter reference config: RCBG=6.0, RCRFI=0.7 ("Faint Target"),
    # RCWGT=1/3, PIXLSIZE=0.05, custom time-shift = -1.807 sec, no
    # noise prior. We currently lack a config knob for the prior; cyga
    # parity uses prior=True (the cyga ref also had prior=True per the
    # UI). For jupiter the UI shows "3rd-Order 2D Polynomial" (no prior).
    # We override surface._fit_poly* to disable the prior just for this
    # render so it matches the reference.
    survey = io.read_sdfits(SRC_FITS)
    cfg = PipelineConfig(
        bg_scale_bw=6.0,
        rfi_scale_bw=0.7,
        weight_scale_bw=1.0 / 3.0,
        pixel_size_bw=0.05,
        trim_size_bw=0.0,
        skip_edge_trim=True,
        time_shift_mode="custom",
        time_shift_seconds=-1.807,
    )

    # Disable the noise prior for jupiter parity (UI screenshot says
    # "3rd-Order 2D Polynomial" without "with Noise Prior").
    from rcpy import surface
    orig_poly3 = surface._fit_poly3
    orig_poly2 = surface._fit_poly2
    orig_poly1 = surface._fit_poly1
    surface._fit_poly3 = lambda dx, dy, z, w, m10_plus_processing=True: \
        orig_poly3(dx, dy, z, w, m10_plus_processing=False)
    surface._fit_poly2 = lambda dx, dy, z, w, m10_plus_processing=True: \
        orig_poly2(dx, dy, z, w, m10_plus_processing=False)
    surface._fit_poly1 = lambda dx, dy, z, w, m10_plus_processing=True: \
        orig_poly1(dx, dy, z, w, m10_plus_processing=False)
    try:
        rc = pipeline.process(survey, cfg)
    finally:
        surface._fit_poly3 = orig_poly3
        surface._fit_poly2 = orig_poly2
        surface._fit_poly1 = orig_poly1

    print("rcpy map shape    :", rc.flux.shape)
    print("reference shape   :", ref_shape)

    # Second pipeline pass for the "raw" panel — match C++'s raw render
    # (skip BG/RFI/timeshift). Time-shift is left off so the raw render
    # uses the raw encoder positions as the C++ raw HDU does.
    survey_raw = io.read_sdfits(SRC_FITS)
    cfg_raw = PipelineConfig(
        bg_scale_bw=0.0, rfi_scale_bw=0.7, weight_scale_bw=0.0,
        pixel_size_bw=0.05, trim_size_bw=0.0,
        skip_edge_trim=True, skip_bg=True, skip_timeshift=True,
        skip_rfi=True,
    )
    rc_raw = pipeline.process(survey_raw, cfg_raw)
    print("rcpy raw shape    :", rc_raw.flux.shape)

    peak_offset = _peak_offset(rc.flux, ref_layers["flux"])
    print(f"peak offset (dy, dx) src - ref: ({peak_offset[0]:+.2f}, {peak_offset[1]:+.2f}) pixels")

    summaries = {}
    src_resampled = {}
    for layer in ("flux", "weight", "scale", "correlation", "path"):
        src_arr = _map_to_array(rc, layer)
        layer_ref_shape = ref_layers[layer].shape
        re = _resample_to_ref(src_arr, peak_offset, layer_ref_shape)
        src_resampled[layer] = re
        summaries[layer] = _summarise(layer, re, ref_layers[layer])

    # Jupiter raw is faint (peak ~6.3 vs background ~5.7), so the
    # parabola peak-refinement can land on a NaN border and produce
    # NaN offsets. Fall back to a NaN-free centroid in that case.
    raw_peak_offset = _peak_offset(rc_raw.flux, ref_raw)
    if not all(np.isfinite(raw_peak_offset)):
        def _centroid(a):
            mask = np.isfinite(a)
            if not mask.any():
                return float(a.shape[0]) / 2.0, float(a.shape[1]) / 2.0
            ys, xs = np.where(mask)
            return float(ys.mean()), float(xs.mean())
        sy, sx = _centroid(rc_raw.flux)
        ry, rx = _centroid(ref_raw)
        raw_peak_offset = (sy - ry, sx - rx)
    raw_resampled = _resample_to_ref(rc_raw.flux, raw_peak_offset, ref_raw.shape)
    summaries["raw"] = _summarise("raw", raw_resampled, ref_raw)

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
        finite_resid = resid[np.isfinite(resid)]
        if finite_resid.size:
            r_lim = float(np.nanpercentile(np.abs(finite_resid), 99))
        else:
            r_lim = 1.0
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
