#!/usr/bin/env python3
"""CRD_DRP post-KcwiKit stack-integrity validation."""

from __future__ import annotations
import argparse, json, sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table
from crd_drp_utils import resolve_stack_paths, load_stack, scan_negative_candidates, recurrent_wavelength_classification

VERSION = "2.0.0"

def parser():
    p = argparse.ArgumentParser()
    for arm in ("blue","red"):
        p.add_argument(f"--{arm}-icube", required=True)
        p.add_argument(f"--{arm}-vcube")
        p.add_argument(f"--{arm}-mcube")
        p.add_argument(f"--{arm}-ecube")
    p.add_argument("--output-dir", default="stack_validation")
    p.add_argument("--recurrent-spatial-fraction", type=float, default=0.03)
    p.add_argument("--recurrent-min-spaxels", type=int, default=25)
    p.add_argument("--recurrence-gap-channels", type=int, default=2)
    return p

def write_table(rows, path):
    if rows:
        Table(rows=rows).write(path, format="ascii.ecsv", overwrite=True)
    else:
        Table().write(path, format="ascii.ecsv", overwrite=True)

def plot_map(data, title, label, path):
    fig, ax = plt.subplots(figsize=(8,6))
    im = ax.imshow(data, origin="lower", aspect="auto")
    fig.colorbar(im, ax=ax, label=label)
    ax.set_xlabel("x pixel"); ax.set_ylabel("y pixel"); ax.set_title(title)
    fig.tight_layout(); fig.savefig(path, dpi=160); plt.close(fig)

def run_arm(arm, args, out):
    paths = resolve_stack_paths(getattr(args,f"{arm}_icube"), getattr(args,f"{arm}_vcube"), getattr(args,f"{arm}_mcube"), getattr(args,f"{arm}_ecube"))
    cube = load_stack(paths)
    raw, structural = scan_negative_candidates(cube)
    classified, recurrence = recurrent_wavelength_classification(
        raw, structural["n_scanned_spaxels"],
        recurrence_gap_channels=args.recurrence_gap_channels,
        recurrent_spatial_fraction=args.recurrent_spatial_fraction,
        recurrent_min_spaxels=args.recurrent_min_spaxels,
    )

    fails = [r for r in classified if r["INTEGRITY_CLASS"]=="LOCALIZED_STACK_CORRUPTION_FAIL"]
    common = [r for r in classified if r["INTEGRITY_CLASS"]=="COMMON_MODE_SPECTRAL_WARN"]
    localw = [r for r in classified if r["INTEGRITY_CLASS"]=="LOCALIZED_NEGATIVE_WARN"]
    structural_fail = structural["n_positive_exposure_nonfinite_science_voxels"]>0 or structural["n_positive_exposure_negative_variance_voxels"]>0
    status = "FAIL" if structural_fail or fails else ("WARN" if common or localw else "PASS")

    fail_map = np.zeros(cube.flux.shape[1:], int)
    for r in fails:
        fail_map[r["Y"], r["X"]] += 1

    plot_map(structural["minimum_flux_map"], f"{arm.upper()} final stack: minimum finite flux", "minimum finite flux", out/f"{arm}_minimum_flux_map.png")
    plot_map(fail_map, f"{arm.upper()} final stack: localized integrity failures", "failure groups", out/f"{arm}_localized_integrity_failures.png")

    for r in classified:
        r["ARM"] = arm.upper()
    write_table(classified, out/f"{arm}_stack_integrity_candidates.ecsv")

    rr = []
    for rec in recurrence:
        row = dict(rec)
        row["WAVE_LO_ANGSTROM"] = float(cube.wavelength[row["INDEX_LO"]])
        row["WAVE_HI_ANGSTROM"] = float(cube.wavelength[row["INDEX_HI"]])
        rr.append(row)
    write_table(rr, out/f"{arm}_recurrent_spectral_features.ecsv")

    finite = cube.flux[np.isfinite(cube.flux)]
    return {
        "arm": arm.upper(), "status": status,
        "paths": {k:str(getattr(paths,k)) if getattr(paths,k) else None for k in ("icube","vcube","mcube","ecube")},
        "shape_wave_y_x": list(map(int,cube.flux.shape)),
        "wavelength_min_angstrom": float(cube.wavelength[0]),
        "wavelength_max_angstrom": float(cube.wavelength[-1]),
        "minimum_finite_flux": float(np.min(finite)) if finite.size else None,
        "n_localized_stack_corruption_failures": len(fails),
        "n_common_mode_spectral_warnings": len(common),
        "n_localized_negative_warnings": len(localw),
        "n_positive_exposure_nonfinite_science_voxels": structural["n_positive_exposure_nonfinite_science_voxels"],
        "n_positive_exposure_negative_variance_voxels": structural["n_positive_exposure_negative_variance_voxels"],
    }

def main():
    args = parser().parse_args()
    out = Path(args.output_dir).expanduser().resolve(); out.mkdir(parents=True, exist_ok=True)
    arms = {}
    for arm in ("blue","red"):
        try:
            arms[arm] = run_arm(arm,args,out)
        except Exception as exc:
            arms[arm] = {"arm":arm.upper(),"status":"ERROR","error":{"type":type(exc).__name__,"message":str(exc)}}
            print(f"{arm.upper()} ERROR: {exc}", file=sys.stderr)
    states = [arms[a]["status"] for a in ("blue","red")]
    overall = "FAIL" if ("FAIL" in states or "ERROR" in states) else ("WARN" if "WARN" in states else "PASS")
    payload = {"pipeline":"CRD_DRP","script":Path(__file__).name,"version":VERSION,"created_utc":datetime.now(timezone.utc).isoformat(),"overall_status":overall,"arms":arms}
    (out/"stack_integrity_summary.json").write_text(json.dumps(payload,indent=2)+"\n")
    print(f"\nCRD_DRP STACK INTEGRITY: {overall}\nBLUE: {arms['blue']['status']}\nRED:  {arms['red']['status']}")
    return 1 if overall=="FAIL" else 0

if __name__=="__main__":
    raise SystemExit(main())
