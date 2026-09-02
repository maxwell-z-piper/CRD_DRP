#!/usr/bin/env python3
"""
Validate the final CRD_DRP package and write CRD_DRP_reduction_manifest.json.

A reduction package fails if:
  * stack-integrity Step 01 failed/errored;
  * cube/mask grids are inconsistent;
  * structural science/variance problems remain;
  * a catastrophic negative group remains outside the finalized atmospheric mask.

WARN-level negative groups outside the atmospheric mask are recorded explicitly
for provenance but do not fail the package.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from astropy.table import Table

from crd_drp_utils import (
    resolve_stack_paths,
    load_stack,
    load_1d_mask_fits,
    scan_negative_candidates,
)

VERSION = "1.1.0"


def parser():
    p = argparse.ArgumentParser()

    for arm in ("blue", "red"):
        p.add_argument(f"--{arm}-icube", required=True)
        p.add_argument(f"--{arm}-vcube")
        p.add_argument(f"--{arm}-mcube")
        p.add_argument(f"--{arm}-ecube")
        p.add_argument(
            f"--{arm}-atmospheric-mask",
            required=True,
        )
        p.add_argument(
            f"--{arm}-atmospheric-intervals",
            required=True,
        )

    p.add_argument("--stack-integrity-summary")
    p.add_argument("--atmospheric-mask-summary")
    p.add_argument("--max-mask-fraction", type=float, default=0.45)
    p.add_argument(
        "--output-dir",
        default="reduction_package_validation",
    )
    return p


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_arm(arm, args, out):
    paths = resolve_stack_paths(
        getattr(args, f"{arm}_icube"),
        getattr(args, f"{arm}_vcube"),
        getattr(args, f"{arm}_mcube"),
        getattr(args, f"{arm}_ecube"),
    )
    cube = load_stack(paths)

    mask_path = Path(
        getattr(args, f"{arm}_atmospheric_mask")
    ).expanduser().resolve()
    interval_path = Path(
        getattr(args, f"{arm}_atmospheric_intervals")
    ).expanduser().resolve()

    if not interval_path.is_file():
        raise FileNotFoundError(interval_path)

    mask, _, _ = load_1d_mask_fits(
        mask_path,
        cube.wavelength,
    )

    mask_fraction = float(mask.mean())
    if mask_fraction > float(args.max_mask_fraction):
        raise ValueError(
            f"{arm.upper()} mask fraction {mask_fraction:.3f} "
            f"exceeds {args.max_mask_fraction:.3f}"
        )

    rows, structural = scan_negative_candidates(cube)

    annotated = []
    fail_inside = []
    fail_outside = []
    warn_inside = []
    warn_outside = []

    for row in rows:
        r = dict(row)

        i0 = int(r["GROUP_START_INDEX"])
        i1 = int(r["GROUP_END_INDEX"])
        inside = bool(mask[i0 : i1 + 1].any())

        r["INSIDE_ATMOSPHERIC_MASK"] = inside
        r["ARM"] = arm.upper()
        annotated.append(r)

        if r["RAW_SEVERITY"] == "FAIL":
            (fail_inside if inside else fail_outside).append(r)
        else:
            (warn_inside if inside else warn_outside).append(r)

    out_table = out / f"{arm}_postmask_negative_candidates.ecsv"
    if annotated:
        Table(rows=annotated).write(
            out_table,
            format="ascii.ecsv",
            overwrite=True,
        )
    else:
        Table().write(
            out_table,
            format="ascii.ecsv",
            overwrite=True,
        )

    structural_fail = (
        structural[
            "n_positive_exposure_nonfinite_science_voxels"
        ] > 0
        or structural[
            "n_positive_exposure_negative_variance_voxels"
        ] > 0
    )

    status = (
        "FAIL"
        if structural_fail or fail_outside
        else "PASS"
    )

    files = {
        "icube": paths.icube,
        "vcube": paths.vcube,
        "mcube": paths.mcube,
        "ecube": paths.ecube,
        "atmospheric_mask_fits": mask_path,
        "atmospheric_intervals_ecsv": interval_path,
    }

    records = {}
    for key, path in files.items():
        if path is None:
            records[key] = None
            continue

        p = Path(path).resolve()
        records[key] = {
            "path": str(p),
            "sha256": sha256(p),
            "size_bytes": int(p.stat().st_size),
        }

    return {
        "arm": arm.upper(),
        "status": status,
        "shape_wave_y_x": list(map(int, cube.flux.shape)),
        "wavelength_min_angstrom": float(cube.wavelength[0]),
        "wavelength_max_angstrom": float(cube.wavelength[-1]),
        "atmospheric_mask_fraction": mask_fraction,
        "atmospheric_masked_native_pixels": int(mask.sum()),

        "n_extreme_fail_groups_inside_atmospheric_mask":
            len(fail_inside),
        "n_unexplained_extreme_fail_groups_outside_atmospheric_mask":
            len(fail_outside),

        "n_warning_groups_inside_atmospheric_mask":
            len(warn_inside),
        "n_unexplained_warning_groups_outside_atmospheric_mask":
            len(warn_outside),

        "n_positive_exposure_nonfinite_science_voxels":
            structural[
                "n_positive_exposure_nonfinite_science_voxels"
            ],
        "n_positive_exposure_negative_variance_voxels":
            structural[
                "n_positive_exposure_negative_variance_voxels"
            ],

        "postmask_candidate_table": str(out_table.resolve()),
        "files": records,
    }


def main():
    args = parser().parse_args()

    out = Path(args.output_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    upstream = {}
    upstream_fail = False

    if args.stack_integrity_summary:
        p = Path(
            args.stack_integrity_summary
        ).expanduser().resolve()
        d = json.loads(p.read_text())

        upstream["stack_integrity_summary"] = {
            "path": str(p),
            "overall_status": d.get("overall_status"),
        }
        upstream_fail = (
            d.get("overall_status")
            in {"FAIL", "ERROR"}
        )

    if args.atmospheric_mask_summary:
        p = Path(
            args.atmospheric_mask_summary
        ).expanduser().resolve()
        d = json.loads(p.read_text())

        upstream["atmospheric_mask_summary"] = {
            "path": str(p),
            "script_version": d.get("version"),
            "selection_principle": d.get(
                "selection_principle"
            ),
            "blue_telluric_reference": d.get(
                "arms", {}
            ).get("blue", {}).get("telluric_reference"),
            "red_telluric_reference": d.get(
                "arms", {}
            ).get("red", {}).get("telluric_reference"),
        }

    arms = {}
    for arm in ("blue", "red"):
        try:
            arms[arm] = validate_arm(
                arm,
                args,
                out,
            )
        except Exception as exc:
            arms[arm] = {
                "arm": arm.upper(),
                "status": "ERROR",
                "error": {
                    "type": type(exc).__name__,
                    "message": str(exc),
                },
            }

    states = [
        arms[a]["status"]
        for a in ("blue", "red")
    ]

    overall = (
        "FAIL"
        if (
            upstream_fail
            or "FAIL" in states
            or "ERROR" in states
        )
        else "PASS"
    )

    summary = {
        "pipeline": "CRD_DRP",
        "script": Path(__file__).name,
        "version": VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "overall_status": overall,
        "upstream_validation": upstream,
        "arms": arms,
    }

    summary_path = out / "reduction_package_validation.json"
    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    manifest = {
        "schema": "CRD_DRP_reduction_manifest_v1",
        "pipeline": "CRD_DRP",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "validation_status": overall,
        "validation_summary": str(summary_path.resolve()),
        "wavelength_frame_contract": (
            "Atmospheric masks are 1-D boolean masks sampled on each "
            "final stack's native observed wavelength grid."
        ),
        "mask_semantics": {
            "true": "exclude from stellar analysis",
            "false": "not excluded by atmospheric mask",
        },
        "blue": arms["blue"],
        "red": arms["red"],
    }

    manifest_path = (
        out / "CRD_DRP_reduction_manifest.json"
    )
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n"
    )

    print(
        f"\nCRD_DRP REDUCTION PACKAGE: {overall}\n"
        f"BLUE: {arms['blue']['status']}\n"
        f"RED:  {arms['red']['status']}\n"
        f"Manifest: {manifest_path}"
    )

    for arm in ("blue", "red"):
        if arms[arm].get("status") in {"PASS", "FAIL"}:
            print(
                f"{arm.upper()} post-mask QC: "
                f"FAIL inside/outside="
                f"{arms[arm]['n_extreme_fail_groups_inside_atmospheric_mask']}/"
                f"{arms[arm]['n_unexplained_extreme_fail_groups_outside_atmospheric_mask']}; "
                f"WARN inside/outside="
                f"{arms[arm]['n_warning_groups_inside_atmospheric_mask']}/"
                f"{arms[arm]['n_unexplained_warning_groups_outside_atmospheric_mask']}"
            )

    return 1 if overall == "FAIL" else 0


if __name__ == "__main__":
    raise SystemExit(main())
