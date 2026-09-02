#!/usr/bin/env python3
"""Shared utilities for CRD_DRP post-stacking validation."""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import math
import numpy as np
from astropy import units as u
from astropy.io import fits
from scipy import ndimage

MAD_TO_SIGMA = 1.482602218505602


@dataclass
class StackPaths:
    icube: Path
    vcube: Path | None
    mcube: Path | None
    ecube: Path | None


@dataclass
class Cube:
    flux: np.ndarray
    variance: np.ndarray | None
    mask: np.ndarray | None
    exposure: np.ndarray | None
    wavelength: np.ndarray
    header: fits.Header


def resolve_path(path):
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(p)
    return p


def infer_companion(icube: Path, letter: str):
    suffix = "_icubes.fits"
    if not icube.name.endswith(suffix):
        return None
    p = icube.with_name(icube.name[:-len(suffix)] + f"_{letter}cubes.fits")
    return p if p.is_file() else None


def resolve_stack_paths(icube, vcube=None, mcube=None, ecube=None):
    i = resolve_path(icube)
    return StackPaths(
        icube=i,
        vcube=resolve_path(vcube) if vcube else infer_companion(i, "v"),
        mcube=resolve_path(mcube) if mcube else infer_companion(i, "m"),
        ecube=resolve_path(ecube) if ecube else infer_companion(i, "e"),
    )


def _find_image_hdu(path, min_ndim=1):
    with fits.open(path, memmap=True) as hdul:
        for hdu in hdul:
            data = getattr(hdu, "data", None)
            if data is None:
                continue
            arr = np.asarray(data)
            if arr.ndim >= min_ndim and np.issubdtype(arr.dtype, np.number):
                return arr.copy(), hdu.header.copy()
    raise ValueError(f"No numeric image HDU with ndim>={min_ndim}: {path}")


def _spectral_numpy_axis(header, ndim):
    tokens = ("WAVE", "AWAV", "LAMB", "SPEC", "FREQ", "VELO", "VRAD")
    for fits_axis in range(1, ndim + 1):
        ctype = str(header.get(f"CTYPE{fits_axis}", "")).upper()
        if any(t in ctype for t in tokens):
            return ndim - fits_axis, fits_axis
    if ndim == 3 and "CRVAL3" in header:
        return 0, 3
    if ndim == 1 and "CRVAL1" in header:
        return 0, 1
    raise ValueError("Could not identify spectral axis.")


def _wave_axis(header, n, fits_axis):
    crval = header.get(f"CRVAL{fits_axis}")
    crpix = header.get(f"CRPIX{fits_axis}", 1.0)
    step = header.get(f"CDELT{fits_axis}", header.get(f"CD{fits_axis}_{fits_axis}"))
    if crval is None or step is None:
        raise ValueError(f"Missing wavelength WCS on FITS axis {fits_axis}.")
    pix = np.arange(n, dtype=float) + 1
    wave = float(crval) + (pix - float(crpix)) * float(step)
    unit_text = str(header.get(f"CUNIT{fits_axis}", "Angstrom")).strip()
    try:
        unit = u.Unit(unit_text) if unit_text else u.AA
        wave = (wave * unit).to_value(u.AA)
    except Exception:
        pass
    return np.asarray(wave, float)


def _read_cube(path):
    raw, hdr = _find_image_hdu(path, 3)
    saxis, faxis = _spectral_numpy_axis(hdr, raw.ndim)
    canon = raw if saxis == 0 else np.moveaxis(raw, saxis, 0)
    return canon, hdr, saxis, faxis


def load_stack(paths: StackPaths):
    flux, hdr, saxis, faxis = _read_cube(paths.icube)
    companions = {}
    for key, p in (("variance", paths.vcube), ("mask", paths.mcube), ("exposure", paths.ecube)):
        if p is None:
            companions[key] = None
            continue
        arr, _, csaxis, _ = _read_cube(p)
        if arr.shape != flux.shape:
            raise ValueError(f"{key} cube shape {arr.shape} != science shape {flux.shape}")
        companions[key] = arr
    wave = _wave_axis(hdr, flux.shape[0], faxis)
    if wave[0] > wave[-1]:
        wave = wave[::-1]
        flux = flux[::-1]
        for k in companions:
            if companions[k] is not None:
                companions[k] = companions[k][::-1]
    return Cube(flux, companions["variance"], companions["mask"], companions["exposure"], wave, hdr)


def robust_center_sigma(values):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if not v.size:
        return np.nan, np.nan
    med = float(np.median(v))
    sig = float(MAD_TO_SIGMA * np.median(np.abs(v - med)))
    return med, sig


def local_negative_metrics(j, spec, valid, inner=3, outer=24):
    lo = max(0, j - outer)
    hi = min(spec.size, j + outer + 1)
    idx = np.arange(lo, hi)
    side = idx[(np.abs(idx - j) >= inner) & valid[idx] & np.isfinite(spec[idx])]
    if side.size < 8:
        return None
    vals = np.asarray(spec[side], float)
    med, sig = robust_center_sigma(vals)
    level = max(float(np.median(np.abs(vals))), np.finfo(float).tiny)
    sig = max(sig if np.isfinite(sig) else 0.0, level * 1e-3, np.finfo(float).tiny)
    f = float(spec[j])
    return {
        "local_median": med,
        "local_sigma": sig,
        "local_z": (f - med) / sig,
        "amplitude_factor": (med - f) / level,
    }


def group_contiguous(indices, gap=1):
    vals = sorted(set(int(v) for v in indices))
    if not vals:
        return []
    out = [[vals[0]]]
    for v in vals[1:]:
        if v - out[-1][-1] <= gap:
            out[-1].append(v)
        else:
            out.append([v])
    return out


def scan_negative_candidates(
    cube,
    warning_local_z=15.0,
    catastrophic_local_z=50.0,
    warning_amplitude_factor=8.0,
    catastrophic_amplitude_factor=25.0,
    rough_z_threshold=10.0,
    sideband_inner_channels=3,
    sideband_outer_channels=24,
):
    flux, exp, var = cube.flux, cube.exposure, cube.variance
    nw, ny, nx = flux.shape
    rows = []
    min_map = np.full((ny, nx), np.nan)
    n_nonfinite_covered = 0
    n_negvar_covered = 0
    n_spax = 0

    for y in range(ny):
        for x in range(nx):
            spec = np.asarray(flux[:, y, x], float)
            covered = np.ones(nw, bool) if exp is None else (np.isfinite(exp[:, y, x]) & (exp[:, y, x] > 0))
            finite = np.isfinite(spec)
            valid = finite & covered
            n_nonfinite_covered += int(np.sum(covered & ~finite))
            if var is not None:
                vv = np.asarray(var[:, y, x], float)
                n_negvar_covered += int(np.sum(covered & np.isfinite(vv) & (vv < 0)))
            if np.count_nonzero(valid) < max(32, 2*sideband_outer_channels+1):
                continue
            n_spax += 1
            min_map[y, x] = float(np.nanmin(spec[valid]))
            med, sig = robust_center_sigma(spec[valid])
            abs_med = max(float(np.median(np.abs(spec[valid]))), np.finfo(float).tiny)
            sig = max(sig if np.isfinite(sig) else 0.0, abs_med*1e-3, np.finfo(float).tiny)
            rz = (spec - med) / sig
            rough = valid & (spec < 0) & (rz <= -abs(float(rough_z_threshold)))
            rough |= valid & (spec < 0) & (np.abs(spec) >= float(warning_amplitude_factor)*abs_med)

            local = []
            for j in np.flatnonzero(rough):
                m = local_negative_metrics(j, spec, valid, sideband_inner_channels, sideband_outer_channels)
                if m is None:
                    continue
                warn = m["local_z"] <= -abs(float(warning_local_z)) and m["amplitude_factor"] >= float(warning_amplitude_factor)
                fail = m["local_z"] <= -abs(float(catastrophic_local_z)) and m["amplitude_factor"] >= float(catastrophic_amplitude_factor)
                if not warn and not fail:
                    continue
                local.append({
                    "WAVE_INDEX": int(j),
                    "WAVELENGTH_ANGSTROM": float(cube.wavelength[j]),
                    "Y": int(y), "X": int(x),
                    "FLUX": float(spec[j]),
                    "LOCAL_Z_ROBUST": float(m["local_z"]),
                    "AMPLITUDE_FACTOR": float(m["amplitude_factor"]),
                    "RAW_SEVERITY": "FAIL" if fail else "WARN",
                })
            if local:
                byj = {r["WAVE_INDEX"]: r for r in local}
                for g in group_contiguous(byj):
                    rs = [byj[j] for j in g]
                    rep = dict(min(rs, key=lambda r: r["LOCAL_Z_ROBUST"]))
                    rep["GROUP_START_INDEX"] = min(g)
                    rep["GROUP_END_INDEX"] = max(g)
                    rep["GROUP_N_CHANNELS"] = len(g)
                    rows.append(rep)

    return rows, {
        "n_scanned_spaxels": n_spax,
        "n_positive_exposure_nonfinite_science_voxels": n_nonfinite_covered,
        "n_positive_exposure_negative_variance_voxels": n_negvar_covered,
        "minimum_flux_map": min_map,
    }


def recurrent_wavelength_classification(rows, n_scanned_spaxels, recurrence_gap_channels=2, recurrent_spatial_fraction=0.03, recurrent_min_spaxels=25):
    if not rows:
        return [], []
    wave_groups = group_contiguous([r["WAVE_INDEX"] for r in rows], gap=recurrence_gap_channels)
    summaries, lookup = [], {}
    for gid, g in enumerate(wave_groups, 1):
        lo, hi = min(g), max(g)
        members = [r for r in rows if lo <= r["WAVE_INDEX"] <= hi]
        spax = {(r["X"], r["Y"]) for r in members}
        frac = len(spax) / max(n_scanned_spaxels, 1)
        recurrent = len(spax) >= recurrent_min_spaxels and frac >= recurrent_spatial_fraction
        summaries.append({"RECURRENCE_GROUP": gid, "INDEX_LO": lo, "INDEX_HI": hi, "N_SPAXELS": len(spax), "SPATIAL_FRACTION": frac, "RECURRENT_COMMON_MODE": recurrent})
        for j in range(lo, hi+1):
            lookup[j] = (gid, recurrent, len(spax), frac)
    out = []
    for row in rows:
        r = dict(row)
        gid, recurrent, nsp, frac = lookup[r["WAVE_INDEX"]]
        r["RECURRENCE_GROUP"] = gid
        r["RECURRENT_COMMON_MODE"] = recurrent
        r["RECURRENT_N_SPAXELS"] = nsp
        r["RECURRENT_SPATIAL_FRACTION"] = frac
        if r["RAW_SEVERITY"] == "FAIL" and not recurrent:
            r["INTEGRITY_CLASS"] = "LOCALIZED_STACK_CORRUPTION_FAIL"
        elif recurrent:
            r["INTEGRITY_CLASS"] = "COMMON_MODE_SPECTRAL_WARN"
        else:
            r["INTEGRITY_CLASS"] = "LOCALIZED_NEGATIVE_WARN"
        out.append(r)
    return out, summaries


def infer_wavelength_medium(header, fallback="vacuum"):
    for key in ("WAVEMED", "AIRORVAC", "WAVEFRAM"):
        if key in header:
            val = str(header[key]).lower()
            if "vac" in val:
                return "vacuum"
            if "air" in val:
                return "air"
    return fallback


def vacuum_to_air(vac):
    vac = np.asarray(vac, float)
    s2 = (1e4 / vac)**2
    n = 1 + 0.0000834254 + 0.02406147/(130-s2) + 0.00015998/(38.9-s2)
    return vac/n


def air_to_vacuum(air):
    air = np.asarray(air, float)
    vac = air*1.00028
    for _ in range(5):
        vac *= air / vacuum_to_air(vac)
    return vac


def convert_medium(wave, src, dst):
    if src == dst:
        return np.asarray(wave, float)
    if src == "air" and dst == "vacuum":
        return air_to_vacuum(wave)
    if src == "vacuum" and dst == "air":
        return vacuum_to_air(wave)
    raise ValueError(f"Unsupported medium conversion {src}->{dst}")


def extract_spatial_median_spectrum(path, object_mask=None):
    cube, hdr, _, faxis = _read_cube(path)
    wave = _wave_axis(hdr, cube.shape[0], faxis)
    spatial = cube.shape[1:]
    sky_sel = np.ones(spatial, bool)
    if object_mask is not None:
        arr, _ = _find_image_hdu(resolve_path(object_mask), 2)
        arr = np.squeeze(arr)
        if arr.shape != spatial:
            raise ValueError(f"Object mask shape {arr.shape} != cube spatial shape {spatial}")
        sky_sel = ~(arr > 0)
    flat = np.asarray(cube, float).reshape(cube.shape[0], -1)
    spec = np.nanmedian(flat[:, sky_sel.ravel()], axis=1)
    if wave[0] > wave[-1]:
        wave, spec = wave[::-1], spec[::-1]
    return wave, spec, hdr, int(np.count_nonzero(sky_sel))


def robust_highpass_z(wave, spectrum, continuum_window_angstrom=75.0):
    w, s = np.asarray(wave, float), np.asarray(spectrum, float)
    finite = np.isfinite(w) & np.isfinite(s)
    dw = float(np.nanmedian(np.diff(w[finite])))
    size = max(5, int(round(continuum_window_angstrom/dw)))
    if size % 2 == 0:
        size += 1
    filled = s.copy()
    filled[~finite] = np.nanmedian(s[finite])
    cont = ndimage.median_filter(filled, size=size, mode="nearest")
    hp = filled - cont
    med, sig = robust_center_sigma(hp[finite])
    if not np.isfinite(sig) or sig <= 0:
        sig = float(np.nanstd(hp[finite]))
    return hp, (hp-med)/sig, sig


def bridge_small_false_gaps(mask, max_gap=1):
    out = np.asarray(mask, bool).copy()
    false_groups = group_contiguous(np.flatnonzero(~out))
    for g in false_groups:
        if min(g) == 0 or max(g) == len(out)-1:
            continue
        if len(g) <= max_gap and out[min(g)-1] and out[max(g)+1]:
            out[min(g):max(g)+1] = True
    return out


def dilate_1d(mask, pixels=0):
    return np.asarray(mask, bool) if pixels <= 0 else ndimage.binary_dilation(mask, iterations=pixels)


def wavelength_edges(center):
    c = np.asarray(center, float)
    e = np.empty(c.size+1)
    e[1:-1] = 0.5*(c[:-1]+c[1:])
    e[0] = c[0]-0.5*(c[1]-c[0])
    e[-1] = c[-1]+0.5*(c[-1]-c[-2])
    return e


def save_1d_mask_fits(path, mask, wave, arm, medium, source_icube):
    m = np.asarray(mask, np.uint8)
    hdr = fits.Header()
    hdr["PIPELINE"] = "CRD_DRP"
    hdr["ARM"] = arm.upper()
    hdr["WAVEMED"] = medium
    hdr["MASKTRUE"] = 1
    hdr["CRPIX1"] = 1.0
    hdr["CRVAL1"] = float(wave[0])
    hdr["CDELT1"] = float(np.nanmedian(np.diff(wave)))
    hdr["CTYPE1"] = "WAVE"
    hdr["CUNIT1"] = "Angstrom"
    hdr["SRCICUBE"] = Path(source_icube).name[:68]
    fits.PrimaryHDU(m, hdr).writeto(path, overwrite=True)


def load_1d_mask_fits(path, expected_wave=None):
    p = resolve_path(path)
    with fits.open(p, memmap=True) as h:
        mask = np.asarray(h[0].data, bool)
        hdr = h[0].header.copy()
    wave = _wave_axis(hdr, mask.size, 1)
    if expected_wave is not None:
        ew = np.asarray(expected_wave, float)
        if ew.shape != wave.shape or not np.allclose(ew, wave, rtol=0, atol=1e-6):
            raise ValueError(f"Atmospheric mask wavelength grid does not match stack: {p}")
    return mask, wave, hdr
