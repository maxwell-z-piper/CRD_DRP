#!/usr/bin/env python3
"""
Validate final KcwiKit BLUE/RED stacks for catastrophic finite negative artifacts.

Run this immediately after KcwiKit stacking and before any CRD_DAP science
pipeline scripts.

This guardrail targets the specific failure mode found during CRD_DAP
development: non-finite trimmed samples could enter KcwiKit's integer
scaling/reprojection path and reappear in the final science stack as sparse,
enormously negative *finite* voxels.

The test is deliberately independent of pPXF and target-specific atmospheric
masks. It asks whether the final stacked cubes contain negative samples that are
orders of magnitude more extreme than their own local spectra.

Typical usage
-------------
python validate_kcwikit_stacks.py \
    --blue-icube /path/blue_icubes.fits \
    --red-icube  /path/red_icubes.fits \
    --output-dir stack_validation

Matching vcube/mcube/ecube files are inferred automatically when standard
KcwiKit filenames are used.

Exit codes
----------
0 : PASS, or WARN unless --fail-on-warning is set
1 : FAIL
2 : invocation/runtime error
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.io import fits
from astropy.table import Table


SCRIPT_VERSION = "1.1.0"
MAD_TO_SIGMA = 1.482602218505602


def parser():
    p = argparse.ArgumentParser(
        description="Validate final KcwiKit BLUE/RED stacks for catastrophic negative artifacts."
    )
    p.add_argument("--blue-icube", required=True)
    p.add_argument("--red-icube", required=True)

    for arm in ("blue", "red"):
        p.add_argument(f"--{arm}-vcube", default=None)
        p.add_argument(f"--{arm}-mcube", default=None)
        p.add_argument(f"--{arm}-ecube", default=None)

    p.add_argument("--output-dir", default="stack_validation")
    p.add_argument("--rough-z-threshold", type=float, default=10.0)
    p.add_argument("--warning-local-z", type=float, default=15.0)
    p.add_argument("--catastrophic-local-z", type=float, default=50.0)
    p.add_argument("--warning-amplitude-factor", type=float, default=8.0)
    p.add_argument("--catastrophic-amplitude-factor", type=float, default=25.0)
    p.add_argument("--sideband-inner-channels", type=int, default=3)
    p.add_argument("--sideband-outer-channels", type=int, default=24)
    p.add_argument("--exposure-deficit-fraction", type=float, default=0.99)
    p.add_argument("--max-candidates-per-arm", type=int, default=5000)
    p.add_argument("--plot-worst", type=int, default=8)
    p.add_argument("--plot-half-width-channels", type=int, default=35)
    p.add_argument("--fail-on-warning", action="store_true")
    return p


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


def infer_companion(icube: Path, letter: str) -> Path | None:
    suffix = "_icubes.fits"
    if not icube.name.endswith(suffix):
        return None
    p = icube.with_name(icube.name[:-len(suffix)] + f"_{letter}cubes.fits")
    return p if p.is_file() else None


def resolve_paths(args, arm: str) -> StackPaths:
    icube = Path(getattr(args, f"{arm}_icube")).expanduser().resolve()
    if not icube.is_file():
        raise FileNotFoundError(f"{arm.upper()} icube not found: {icube}")

    result = {"icube": icube}
    for key, letter in (("vcube", "v"), ("mcube", "m"), ("ecube", "e")):
        explicit = getattr(args, f"{arm}_{key}")
        if explicit:
            p = Path(explicit).expanduser().resolve()
            if not p.is_file():
                raise FileNotFoundError(f"{arm.upper()} {key} not found: {p}")
        else:
            p = infer_companion(icube, letter)
        result[key] = p
    return StackPaths(**result)


def spectral_axis(header, ndim=3):
    tokens = ("WAVE", "AWAV", "FREQ", "VELO", "VRAD")
    for fits_axis in range(1, ndim + 1):
        ctype = str(header.get(f"CTYPE{fits_axis}", "")).upper()
        if any(t in ctype for t in tokens):
            return ndim - fits_axis
    if ndim == 3 and "CRVAL3" in header and ("CDELT3" in header or "CD3_3" in header):
        return 0
    raise ValueError("Could not determine spectral axis from FITS WCS.")


def wavelength_axis(header, nwave, numpy_spectral_axis):
    fits_axis = 3 - numpy_spectral_axis
    crval = header.get(f"CRVAL{fits_axis}")
    crpix = header.get(f"CRPIX{fits_axis}")
    step = header.get(f"CDELT{fits_axis}", header.get(f"CD{fits_axis}_{fits_axis}"))
    if crval is None or crpix is None or step is None:
        return np.arange(nwave, dtype=float)
    pix = np.arange(nwave, dtype=float) + 1.0
    wave = float(crval) + (pix - float(crpix)) * float(step)
    unit_text = str(header.get(f"CUNIT{fits_axis}", "Angstrom")).strip()
    try:
        unit = u.Unit(unit_text) if unit_text else u.AA
        wave = (wave * unit).to_value(u.AA)
    except Exception:
        pass
    return np.asarray(wave, dtype=float)


def read_primary(path: Path):
    with fits.open(path, memmap=True) as h:
        data = np.asarray(h[0].data)
        hdr = h[0].header.copy()
    if data.ndim != 3:
        raise ValueError(f"Expected 3-D cube, got {data.shape}: {path}")
    return data, hdr


def canonical(data, saxis):
    return data if saxis == 0 else np.moveaxis(data, saxis, 0)


def load_stack(paths: StackPaths) -> Cube:
    raw, hdr = read_primary(paths.icube)
    saxis = spectral_axis(hdr)
    flux = canonical(raw, saxis)

    comps = {}
    for key, p in (("variance", paths.vcube), ("mask", paths.mcube), ("exposure", paths.ecube)):
        if p is None:
            comps[key] = None
            continue
        a, _ = read_primary(p)
        if a.shape != raw.shape:
            raise ValueError(f"{key} cube shape {a.shape} != science shape {raw.shape}: {p}")
        comps[key] = canonical(a, saxis)

    wave = wavelength_axis(hdr, flux.shape[0], saxis)
    return Cube(flux, comps["variance"], comps["mask"], comps["exposure"], wave)


def robust_center_sigma(values):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if not v.size:
        return np.nan, np.nan
    med = float(np.median(v))
    sig = float(MAD_TO_SIGMA * np.median(np.abs(v - med)))
    return med, sig


def cube_amplitude_floor(flux):
    nw, ny, nx = flux.shape
    ws = np.unique(np.linspace(0, nw - 1, min(nw, 128)).astype(int))
    ys = np.unique(np.linspace(0, ny - 1, min(ny, 12)).astype(int))
    xs = np.unique(np.linspace(0, nx - 1, min(nx, 12)).astype(int))
    sample = np.abs(np.asarray(flux[np.ix_(ws, ys, xs)], float).ravel())
    sample = sample[np.isfinite(sample)]
    if not sample.size:
        return np.finfo(float).tiny
    return max(float(np.median(sample)) * 0.05, np.finfo(float).tiny)


def local_metrics(j, spec, valid, inner, outer, amplitude_floor):
    lo = max(0, j - outer)
    hi = min(spec.size, j + outer + 1)
    idx = np.arange(lo, hi)
    sep = np.abs(idx - j)
    side_idx = idx[(sep >= inner) & (sep <= outer)]
    use = valid[side_idx] & np.isfinite(spec[side_idx])
    side = np.asarray(spec[side_idx[use]], float)
    if side.size < 8:
        return None

    med, sig = robust_center_sigma(side)
    if not np.isfinite(med):
        return None

    level = max(float(np.median(np.abs(side))), amplitude_floor)
    sigma_floor = max(level * 1e-3, amplitude_floor)
    sig_eff = max(sig if np.isfinite(sig) else 0.0, sigma_floor)

    f = float(spec[j])
    depth = med - f
    return {
        "local_median": med,
        "local_sigma": sig_eff,
        "local_z": (f - med) / sig_eff,
        "amplitude_factor": depth / level,
    }


def groups(indices):
    vals = sorted(set(int(v) for v in indices))
    if not vals:
        return []
    out = [[vals[0]]]
    for v in vals[1:]:
        if v - out[-1][-1] <= 1:
            out[-1].append(v)
        else:
            out.append([v])
    return out


def scan_arm(arm, cube, args):
    flux, var, mask, exp, wave = (
        cube.flux, cube.variance, cube.mask, cube.exposure, cube.wavelength
    )
    nw, ny, nx = flux.shape
    amp_floor = cube_amplitude_floor(flux)

    if exp is not None:
        pe = np.asarray(exp[np.isfinite(exp) & (exp > 0)], float)
        nominal_exp = float(np.percentile(pe, 99.5)) if pe.size else np.nan
    else:
        nominal_exp = np.nan

    min_map = np.full((ny, nx), np.nan)
    fail_map = np.zeros((ny, nx), dtype=int)
    candidates = []

    # Keep complete counts even if a very bad/old reduction generates thousands
    # of candidates.  The ECSV is capped only to keep the output manageable.
    n_warning_total = 0
    n_catastrophic_total = 0
    candidate_table_truncated = False

    nonfinite_with_coverage = 0
    negative_variance_with_coverage = 0
    scanned_spaxels = 0

    for y in range(ny):
        for x in range(nx):
            spec = np.asarray(flux[:, y, x], float)
            if exp is not None:
                es = np.asarray(exp[:, y, x], float)
                covered = np.isfinite(es) & (es > 0)
            else:
                es = None
                covered = np.ones(nw, bool)

            finite = np.isfinite(spec)
            valid = finite & covered
            nonfinite_with_coverage += int(np.sum(covered & ~finite))

            if var is not None:
                vs = np.asarray(var[:, y, x], float)
                negative_variance_with_coverage += int(
                    np.sum(covered & np.isfinite(vs) & (vs < 0))
                )

            if np.count_nonzero(valid) < max(32, 2 * args.sideband_outer_channels + 1):
                continue

            scanned_spaxels += 1
            min_map[y, x] = float(np.nanmin(spec[valid]))

            center, gsig = robust_center_sigma(spec[valid])
            abs_med = max(float(np.median(np.abs(spec[valid]))), amp_floor)
            gsig = max(gsig if np.isfinite(gsig) else 0.0, max(abs_med * 1e-3, amp_floor))

            rz = (spec - center) / gsig
            rough = valid & (spec < 0) & (rz <= -abs(args.rough_z_threshold))
            rough |= valid & (spec < 0) & (
                np.abs(spec) >= args.warning_amplitude_factor * abs_med
            )

            rows_here = []
            for j in np.flatnonzero(rough):
                met = local_metrics(
                    int(j), spec, valid,
                    args.sideband_inner_channels,
                    args.sideband_outer_channels,
                    amp_floor,
                )
                if met is None:
                    continue

                warn = (
                    met["local_z"] <= -abs(args.warning_local_z)
                    and met["amplitude_factor"] >= args.warning_amplitude_factor
                )
                fail = (
                    met["local_z"] <= -abs(args.catastrophic_local_z)
                    and met["amplitude_factor"] >= args.catastrophic_amplitude_factor
                )
                if not warn and not fail:
                    continue

                if exp is not None and np.isfinite(es[j]):
                    e_sec = float(es[j])
                    e_frac = e_sec / nominal_exp if np.isfinite(nominal_exp) and nominal_exp > 0 else np.nan
                else:
                    e_sec = e_frac = np.nan

                partial = bool(
                    np.isfinite(e_frac) and e_frac < args.exposure_deficit_fraction
                )

                if mask is not None and np.isfinite(mask[j, y, x]):
                    mraw = float(mask[j, y, x])
                    mflag = bool(mraw != 0)
                else:
                    mraw = np.nan
                    mflag = False

                if var is not None and np.isfinite(var[j, y, x]):
                    vv = float(var[j, y, x])
                    ss = math.sqrt(vv) if vv >= 0 else np.nan
                    formal_z = (
                        (float(spec[j]) - met["local_median"]) / ss
                        if np.isfinite(ss) and ss > 0 else np.nan
                    )
                else:
                    vv = formal_z = np.nan

                rows_here.append({
                    "ARM": arm.upper(),
                    "SEVERITY": "FAIL" if fail else "WARN",
                    "WAVE_INDEX": int(j),
                    "WAVELENGTH_ANGSTROM": float(wave[j]),
                    "Y": int(y),
                    "X": int(x),
                    "FLUX": float(spec[j]),
                    "LOCAL_MEDIAN": float(met["local_median"]),
                    "LOCAL_SIGMA_ROBUST": float(met["local_sigma"]),
                    "LOCAL_Z_ROBUST": float(met["local_z"]),
                    "AMPLITUDE_FACTOR": float(met["amplitude_factor"]),
                    "FORMAL_Z_IF_VARIANCE": float(formal_z) if np.isfinite(formal_z) else np.nan,
                    "VARIANCE": vv,
                    "MASK_RAW": mraw,
                    "MASK_FLAGGED": mflag,
                    "EXPOSURE_SECONDS": e_sec,
                    "EXPOSURE_FRACTION": float(e_frac) if np.isfinite(e_frac) else np.nan,
                    "PARTIAL_COVERAGE": partial,
                })

            if rows_here:
                byj = {r["WAVE_INDEX"]: r for r in rows_here}
                for g in groups(byj):
                    rs = [byj[j] for j in g]
                    rep = dict(min(rs, key=lambda r: r["LOCAL_Z_ROBUST"]))
                    rep["GROUP_START_INDEX"] = min(g)
                    rep["GROUP_END_INDEX"] = max(g)
                    rep["GROUP_N_CHANNELS"] = len(g)
                    if rep["SEVERITY"] == "FAIL":
                        n_catastrophic_total += 1
                        fail_map[y, x] += 1
                    else:
                        n_warning_total += 1

                    # Do NOT abort an arm when the detailed-candidate table gets
                    # large.  This is expected for an old/unfixed reduction.
                    # Continue scanning the entire arm so RED can still be tested
                    # even if BLUE is severely corrupted, and so the summary
                    # counts/maps remain complete.
                    if len(candidates) < args.max_candidates_per_arm:
                        candidates.append(rep)
                    else:
                        candidate_table_truncated = True

    finite_flux = np.asarray(flux[np.isfinite(flux)], float)
    summary = {
        "arm": arm.upper(),
        "shape_wave_y_x": [int(v) for v in flux.shape],
        "wavelength_min_angstrom": float(np.nanmin(wave)),
        "wavelength_max_angstrom": float(np.nanmax(wave)),
        "finite_flux_fraction": float(np.mean(np.isfinite(flux))),
        "minimum_finite_flux": float(np.min(finite_flux)) if finite_flux.size else np.nan,
        "median_finite_flux": float(np.median(finite_flux)) if finite_flux.size else np.nan,
        "nominal_effective_exposure_seconds": nominal_exp,
        "n_scanned_spaxels": scanned_spaxels,
        "n_warning_candidate_groups": int(n_warning_total),
        "n_catastrophic_candidate_groups": int(n_catastrophic_total),
        "n_spaxels_with_catastrophic_candidate": int(np.sum(fail_map > 0)),
        "candidate_table_rows_stored": int(len(candidates)),
        "candidate_table_truncated": bool(candidate_table_truncated),
        "candidate_table_row_limit": int(args.max_candidates_per_arm),
        "n_positive_exposure_nonfinite_science_voxels": nonfinite_with_coverage,
        "n_positive_exposure_negative_variance_voxels": negative_variance_with_coverage,
    }
    return candidates, summary, min_map, fail_map


COLS = [
    "ARM","SEVERITY","WAVE_INDEX","WAVELENGTH_ANGSTROM","Y","X","FLUX",
    "LOCAL_MEDIAN","LOCAL_SIGMA_ROBUST","LOCAL_Z_ROBUST","AMPLITUDE_FACTOR",
    "FORMAL_Z_IF_VARIANCE","VARIANCE","MASK_RAW","MASK_FLAGGED",
    "EXPOSURE_SECONDS","EXPOSURE_FRACTION","PARTIAL_COVERAGE",
    "GROUP_START_INDEX","GROUP_END_INDEX","GROUP_N_CHANNELS",
]


def write_candidate_table(rows, path):
    if rows:
        tab = Table(rows=[[r.get(c) for c in COLS] for r in rows], names=COLS)
    else:
        tab = Table(names=COLS, dtype=[
            "U8","U8","i8","f8","i8","i8","f8","f8","f8","f8","f8",
            "f8","f8","f8","bool","f8","f8","bool","i8","i8","i8"
        ])
    tab.meta["SCRIPT_VERSION"] = SCRIPT_VERSION
    tab.write(path, format="ascii.ecsv", overwrite=True)


def plot_map(data, title, label, path):
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(data, origin="lower", aspect="auto")
    fig.colorbar(im, ax=ax, label=label)
    ax.set_xlabel("x pixel")
    ax.set_ylabel("y pixel")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_worst(arm, cube, candidates, nplot, half, path):
    if not candidates or nplot <= 0:
        return
    rows = sorted(candidates, key=lambda r: r["LOCAL_Z_ROBUST"])[:nplot]
    fig, axes = plt.subplots(len(rows), 1, figsize=(12, max(3, 2.5*len(rows))), squeeze=False)
    for i, r in enumerate(rows):
        ax = axes[i, 0]
        j, y, x = r["WAVE_INDEX"], r["Y"], r["X"]
        lo, hi = max(0, j-half), min(cube.flux.shape[0], j+half+1)
        ax.plot(cube.wavelength[lo:hi], cube.flux[lo:hi, y, x], lw=1)
        ax.axvline(cube.wavelength[j], ls="--", lw=0.8)
        ax.axhline(r["LOCAL_MEDIAN"], ls=":", lw=0.8)
        ax.set_ylabel("flux")
        ax.set_title(
            f"{r['SEVERITY']} | (x,y)=({x},{y}) | wave={r['WAVELENGTH_ANGSTROM']:.2f} A | "
            f"flux={r['FLUX']:.5g} | z={r['LOCAL_Z_ROBUST']:.1f} | factor={r['AMPLITUDE_FACTOR']:.1f}"
        )
    axes[-1, 0].set_xlabel("Observed wavelength [A]")
    fig.suptitle(f"{arm.upper()} stack: worst negative validation candidates")
    fig.tight_layout(rect=(0,0,1,0.985))
    fig.savefig(path, dpi=160)
    plt.close(fig)


def status(summary, fail_on_warning):
    if summary.get("arm_error"):
        return "ERROR"

    structural = (
        summary["n_positive_exposure_nonfinite_science_voxels"] > 0
        or summary["n_positive_exposure_negative_variance_voxels"] > 0
    )
    if structural or summary["n_catastrophic_candidate_groups"] > 0:
        return "FAIL"
    if summary["n_warning_candidate_groups"] > 0:
        return "FAIL" if fail_on_warning else "WARN"
    return "PASS"


def main():
    args = parser().parse_args()
    if args.sideband_inner_channels < 1:
        raise ValueError("--sideband-inner-channels must be >= 1.")
    if args.sideband_outer_channels <= args.sideband_inner_channels + 4:
        raise ValueError("--sideband-outer-channels is too small.")

    out = Path(args.output_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    all_rows = []
    summaries = {}

    for arm in ("blue", "red"):
        try:
            paths = resolve_paths(args, arm)
            print(f"\nLoading {arm.upper()}")
            print(f"  icube: {paths.icube}")
            print(f"  vcube: {paths.vcube or 'not supplied/found'}")
            print(f"  mcube: {paths.mcube or 'not supplied/found'}")
            print(f"  ecube: {paths.ecube or 'not supplied/found'}")

            cube = load_stack(paths)
            rows, summary, min_map, fail_map = scan_arm(arm, cube, args)
            arm_status = status(summary, args.fail_on_warning)
            summary["status"] = arm_status
            summary["arm_error"] = None
            summary["paths"] = {
                "icube": str(paths.icube),
                "vcube": str(paths.vcube) if paths.vcube else None,
                "mcube": str(paths.mcube) if paths.mcube else None,
                "ecube": str(paths.ecube) if paths.ecube else None,
            }
            summaries[arm] = summary
            all_rows.extend(rows)

            plot_map(
                min_map,
                f"{arm.upper()} final stack: minimum flux per spatial pixel",
                "minimum finite flux",
                out / f"{arm}_minimum_flux_map.png",
            )
            plot_map(
                fail_map,
                f"{arm.upper()} final stack: catastrophic negative groups",
                "number of catastrophic groups",
                out / f"{arm}_catastrophic_count_map.png",
            )
            plot_worst(
                arm, cube, rows, args.plot_worst, args.plot_half_width_channels,
                out / f"{arm}_worst_negative_candidates.png"
            )

            print(f"\n{arm.upper()} STACK VALIDATION: {arm_status}")
            print(f"  minimum finite flux:              {summary['minimum_finite_flux']:.6g}")
            print(f"  warning candidate groups:         {summary['n_warning_candidate_groups']}")
            print(f"  catastrophic candidate groups:    {summary['n_catastrophic_candidate_groups']}")
            print(f"  spaxels with catastrophic groups: {summary['n_spaxels_with_catastrophic_candidate']}")
            print(f"  candidate rows stored:            {summary['candidate_table_rows_stored']}")
            if summary["candidate_table_truncated"]:
                print(
                    f"  NOTE: candidate ECSV truncated at {summary['candidate_table_row_limit']} rows; "
                    "full-arm counts/maps are still complete."
                )
            print(f"  nonfinite science with e>0:       {summary['n_positive_exposure_nonfinite_science_voxels']}")
            print(f"  negative variance with e>0:       {summary['n_positive_exposure_negative_variance_voxels']}")

        except Exception as exc:
            # One bad/old arm must not prevent validation of the other arm.
            # Record the failure in the final JSON and continue.
            print(
                f"\n{arm.upper()} STACK VALIDATION: ERROR\n"
                f"  {type(exc).__name__}: {exc}\n"
                f"  Continuing to the other arm..."
            )
            summaries[arm] = {
                "arm": arm.upper(),
                "status": "ERROR",
                "arm_error": {
                    "type": type(exc).__name__,
                    "message": str(exc),
                },
                "paths": {
                    "icube": str(getattr(args, f"{arm}_icube")),
                    "vcube": getattr(args, f"{arm}_vcube"),
                    "mcube": getattr(args, f"{arm}_mcube"),
                    "ecube": getattr(args, f"{arm}_ecube"),
                },
            }

    write_candidate_table(all_rows, out / "stack_validation_candidates.ecsv")

    states = [summaries[a]["status"] for a in ("blue", "red")]
    overall = (
        "FAIL"
        if ("FAIL" in states or "ERROR" in states)
        else ("WARN" if "WARN" in states else "PASS")
    )

    payload = {
        "script": "validate_kcwikit_stacks.py",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "overall_status": overall,
        "thresholds": {
            "rough_z_threshold": args.rough_z_threshold,
            "warning_local_z": args.warning_local_z,
            "catastrophic_local_z": args.catastrophic_local_z,
            "warning_amplitude_factor": args.warning_amplitude_factor,
            "catastrophic_amplitude_factor": args.catastrophic_amplitude_factor,
            "sideband_inner_channels": args.sideband_inner_channels,
            "sideband_outer_channels": args.sideband_outer_channels,
            "exposure_deficit_fraction": args.exposure_deficit_fraction,
        },
        "arms": summaries,
    }
    (out / "stack_validation_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )

    print("\n" + "="*72)
    print(f"OVERALL KcwiKit STACK VALIDATION: {overall}")
    print("="*72)
    print(f"BLUE: {summaries['blue']['status']}")
    print(f"RED:  {summaries['red']['status']}")
    print(f"Candidates: {out/'stack_validation_candidates.ecsv'}")
    print(f"Summary:    {out/'stack_validation_summary.json'}")

    return 1 if overall == "FAIL" else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
