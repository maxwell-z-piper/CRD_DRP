#!/usr/bin/env python3
"""
Optional CRD_DRP audit of Ca II triplet (CaT) wavelength retention.

Purpose
-------
This script does NOT modify the atmospheric mask and does NOT decide whether a
pPXF fit is scientifically valid.  It answers a narrower reduction-QC question:

    "Given the finalized RED atmospheric mask, how much of each CaT feature and
     its surrounding continuum remains available to CRD_DAP?"

The audit is intentionally downstream of CRD_DRP Step 03 and can read the
CRD_DRP reduction manifest directly.

It reports:
  * expected observed-frame CaT centers for the target redshift;
  * retained/masked native-pixel fractions in a narrow line-core window;
  * retained/masked fractions in a wider context window;
  * whether the nearest native pixel to the nominal line center is masked;
  * longest contiguous retained wavelength span in each window;
  * if the Step-02 atmospheric-reference ECSV is available:
      - empirical-sky, known-line, and telluric mask fractions separately;
      - minimum telluric transmission;
      - maximum empirical sky significance;
  * retained fraction across a configurable rest-frame CaT analysis window.

Default rest wavelengths are the standard AIR Ca II triplet wavelengths:
8498.02, 8542.09, and 8662.14 Angstrom.

The default windows are:
  core    = +/- 5 A observed
  context = +/- 25 A observed

These are diagnostic defaults, not universal scientific thresholds.

Outputs
-------
CaT_mask_audit.ecsv
CaT_mask_audit_summary.json
CaT_mask_overview.png
CaT_8498_mask_zoom.png
CaT_8542_mask_zoom.png
CaT_8662_mask_zoom.png

Exit code
---------
0 : audit completed
2 : invocation/runtime error

Example
-------
python validation/audit_CaT_mask_coverage.py \
    --manifest validation_products/03_reduction_package/CRD_DRP_reduction_manifest.json \
    --redshift 0.04138 \
    --output-dir validation_products/CaT_audit
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy.table import Table

from crd_drp_utils import (
    resolve_stack_paths,
    load_stack,
    load_1d_mask_fits,
    convert_medium,
)

VERSION = "1.0.0"

CAT_LINES_AIR_ANGSTROM = (
    ("CaT_8498", 8498.02),
    ("CaT_8542", 8542.09),
    ("CaT_8662", 8662.14),
)


def parser():
    p = argparse.ArgumentParser()

    p.add_argument(
        "--manifest",
        default=None,
        help=(
            "CRD_DRP_reduction_manifest.json. Preferred input mode. "
            "If omitted, --red-icube and --red-atmospheric-mask are required."
        ),
    )

    p.add_argument("--red-icube", default=None)
    p.add_argument("--red-vcube", default=None)
    p.add_argument("--red-mcube", default=None)
    p.add_argument("--red-ecube", default=None)
    p.add_argument("--red-atmospheric-mask", default=None)
    p.add_argument("--red-atmospheric-intervals", default=None)
    p.add_argument(
        "--atmospheric-reference",
        default=None,
        help=(
            "Optional Step-02 red_atmospheric_reference.ecsv. "
            "If omitted, it is inferred next to the atmospheric mask when present."
        ),
    )

    p.add_argument(
        "--redshift",
        type=float,
        required=True,
        help="Target systemic redshift.",
    )

    p.add_argument(
        "--rest-medium",
        choices=("air", "vacuum"),
        default="air",
        help="Wavelength medium of the supplied/rest CaT wavelengths (default air).",
    )
    p.add_argument(
        "--science-medium",
        choices=("auto", "air", "vacuum"),
        default="auto",
        help=(
            "Final RED science wavelength medium. With auto, use WAVEMED from "
            "the atmospheric-mask FITS header, falling back to vacuum."
        ),
    )

    p.add_argument(
        "--core-half-width-angstrom",
        type=float,
        default=5.0,
    )
    p.add_argument(
        "--context-half-width-angstrom",
        type=float,
        default=25.0,
    )
    p.add_argument(
        "--cat-rest-window",
        nargs=2,
        type=float,
        metavar=("LO", "HI"),
        default=(8400.0, 8750.0),
        help=(
            "Broad rest-frame CaT audit window in --rest-medium Angstrom "
            "(default 8400 8750)."
        ),
    )

    # These thresholds only create human-readable flags. They do not cause a
    # nonzero exit status.
    p.add_argument(
        "--review-core-retained-fraction",
        type=float,
        default=0.50,
        help=(
            "Flag a line for REVIEW if less than this fraction of its core "
            "window remains unmasked (default 0.50)."
        ),
    )
    p.add_argument(
        "--review-context-retained-fraction",
        type=float,
        default=0.50,
        help=(
            "Flag a line for REVIEW if less than this fraction of its context "
            "window remains unmasked (default 0.50)."
        ),
    )

    p.add_argument(
        "--output-dir",
        default="validation_products/CaT_audit",
    )

    return p


def _manifest_file_path(record):
    if record is None:
        return None
    if isinstance(record, str):
        return record
    if isinstance(record, dict):
        return record.get("path")
    return None


def resolve_inputs(args):
    manifest = None

    if args.manifest is not None:
        mp = Path(args.manifest).expanduser().resolve()
        if not mp.is_file():
            raise FileNotFoundError(mp)

        manifest = json.loads(mp.read_text())
        if manifest.get("pipeline") != "CRD_DRP":
            raise ValueError(
                f"{mp} does not identify itself as a CRD_DRP manifest."
            )

        red = manifest.get("red", {})
        files = red.get("files", {})

        red_icube = _manifest_file_path(files.get("icube"))
        red_vcube = _manifest_file_path(files.get("vcube"))
        red_mcube = _manifest_file_path(files.get("mcube"))
        red_ecube = _manifest_file_path(files.get("ecube"))
        mask_path = _manifest_file_path(
            files.get("atmospheric_mask_fits")
        )
        intervals_path = _manifest_file_path(
            files.get("atmospheric_intervals_ecsv")
        )

        if red_icube is None or mask_path is None:
            raise ValueError(
                "Manifest does not contain RED icube + atmospheric-mask paths."
            )
    else:
        red_icube = args.red_icube
        red_vcube = args.red_vcube
        red_mcube = args.red_mcube
        red_ecube = args.red_ecube
        mask_path = args.red_atmospheric_mask
        intervals_path = args.red_atmospheric_intervals

        if red_icube is None or mask_path is None:
            raise ValueError(
                "Without --manifest, provide --red-icube and "
                "--red-atmospheric-mask."
            )

    # Explicit CLI values override manifest values if supplied.
    red_icube = args.red_icube or red_icube
    red_vcube = args.red_vcube or red_vcube
    red_mcube = args.red_mcube or red_mcube
    red_ecube = args.red_ecube or red_ecube
    mask_path = args.red_atmospheric_mask or mask_path
    intervals_path = (
        args.red_atmospheric_intervals or intervals_path
    )

    mask_path = Path(mask_path).expanduser().resolve()

    if args.atmospheric_reference is not None:
        reference_path = Path(
            args.atmospheric_reference
        ).expanduser().resolve()
    else:
        candidate = mask_path.with_name(
            "red_atmospheric_reference.ecsv"
        )
        reference_path = candidate if candidate.is_file() else None

    return {
        "manifest": manifest,
        "red_icube": red_icube,
        "red_vcube": red_vcube,
        "red_mcube": red_mcube,
        "red_ecube": red_ecube,
        "mask_path": mask_path,
        "intervals_path": (
            Path(intervals_path).expanduser().resolve()
            if intervals_path is not None
            else None
        ),
        "reference_path": reference_path,
    }


def longest_true_run(mask, wave):
    mask = np.asarray(mask, bool)
    if not np.any(mask):
        return 0, 0.0

    idx = np.flatnonzero(mask)
    split = np.flatnonzero(np.diff(idx) > 1) + 1
    groups = np.split(idx, split)
    best = max(groups, key=len)

    n = int(best.size)

    if n == 1:
        if wave.size > 1:
            dw = float(np.nanmedian(np.diff(wave)))
            span = abs(dw)
        else:
            span = 0.0
    else:
        # Include approximately one native pixel width across the run.
        dw = float(np.nanmedian(np.diff(wave)))
        span = float(abs(wave[best[-1]] - wave[best[0]]) + abs(dw))

    return n, span


def window_metrics(
    wave,
    atmospheric_mask,
    center,
    half_width,
    reference=None,
):
    sel = (
        (wave >= center - half_width)
        & (wave <= center + half_width)
    )
    n_total = int(np.count_nonzero(sel))

    if n_total == 0:
        return {
            "n_total": 0,
            "n_masked": 0,
            "n_retained": 0,
            "masked_fraction": np.nan,
            "retained_fraction": np.nan,
            "longest_retained_run_pixels": 0,
            "longest_retained_run_angstrom": 0.0,
            "empirical_sky_mask_fraction": np.nan,
            "known_line_mask_fraction": np.nan,
            "telluric_mask_fraction": np.nan,
            "min_telluric_transmission": np.nan,
            "max_empirical_sky_z": np.nan,
        }

    masked = atmospheric_mask[sel]
    retained = ~masked

    local_wave = wave[sel]
    run_n, run_ang = longest_true_run(retained, local_wave)

    out = {
        "n_total": n_total,
        "n_masked": int(np.count_nonzero(masked)),
        "n_retained": int(np.count_nonzero(retained)),
        "masked_fraction": float(np.mean(masked)),
        "retained_fraction": float(np.mean(retained)),
        "longest_retained_run_pixels": int(run_n),
        "longest_retained_run_angstrom": float(run_ang),
        "empirical_sky_mask_fraction": np.nan,
        "known_line_mask_fraction": np.nan,
        "telluric_mask_fraction": np.nan,
        "min_telluric_transmission": np.nan,
        "max_empirical_sky_z": np.nan,
    }

    if reference is not None:
        if "EMPIRICAL_SKY_MASK" in reference.colnames:
            out["empirical_sky_mask_fraction"] = float(
                np.mean(
                    np.asarray(
                        reference["EMPIRICAL_SKY_MASK"][sel],
                        bool,
                    )
                )
            )

        if "KNOWN_LINE_MASK" in reference.colnames:
            out["known_line_mask_fraction"] = float(
                np.mean(
                    np.asarray(
                        reference["KNOWN_LINE_MASK"][sel],
                        bool,
                    )
                )
            )

        if "TELLURIC_MASK" in reference.colnames:
            out["telluric_mask_fraction"] = float(
                np.mean(
                    np.asarray(
                        reference["TELLURIC_MASK"][sel],
                        bool,
                    )
                )
            )

        if "TELLURIC_TRANSMISSION" in reference.colnames:
            vals = np.asarray(
                reference["TELLURIC_TRANSMISSION"][sel],
                float,
            )
            if np.any(np.isfinite(vals)):
                out["min_telluric_transmission"] = float(
                    np.nanmin(vals)
                )

        if "EMPIRICAL_SKY_Z" in reference.colnames:
            vals = np.asarray(
                reference["EMPIRICAL_SKY_Z"][sel],
                float,
            )
            if np.any(np.isfinite(vals)):
                out["max_empirical_sky_z"] = float(
                    np.nanmax(vals)
                )

    return out


def overlapping_interval_reasons(intervals, lo, hi):
    if intervals is None or len(intervals) == 0:
        return ""

    required = {
        "WAVE_LO_ANGSTROM",
        "WAVE_HI_ANGSTROM",
        "REASON",
    }
    if not required.issubset(set(intervals.colnames)):
        return ""

    reasons = []
    for row in intervals:
        rlo = float(row["WAVE_LO_ANGSTROM"])
        rhi = float(row["WAVE_HI_ANGSTROM"])

        if rhi >= lo and rlo <= hi:
            reasons.append(str(row["REASON"]))

    return "; ".join(sorted(set(reasons)))


def plot_overview(
    wave,
    atmospheric_mask,
    centers,
    broad_lo,
    broad_hi,
    path,
):
    sel = (wave >= broad_lo) & (wave <= broad_hi)

    fig, ax = plt.subplots(figsize=(13, 4))

    ax.step(
        wave[sel],
        (~atmospheric_mask[sel]).astype(float),
        where="mid",
        lw=1.0,
        label="retained native pixel",
    )

    for label, center in centers:
        ax.axvline(center, ls="--", lw=1.0)
        ax.text(
            center,
            1.04,
            label.replace("CaT_", ""),
            ha="center",
            va="bottom",
            transform=ax.get_xaxis_transform(),
        )

    ax.set_ylim(-0.08, 1.18)
    ax.set_ylabel("retained (1=yes)")
    ax.set_xlabel("Observed wavelength [A]")
    ax.set_title(
        "CRD_DRP CaT atmospheric-mask retention"
    )
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def plot_line_zoom(
    label,
    wave,
    atmospheric_mask,
    center,
    context_half_width,
    reference,
    path,
):
    sel = (
        (wave >= center - context_half_width)
        & (wave <= center + context_half_width)
    )

    fig, ax = plt.subplots(figsize=(10, 4))

    ax.step(
        wave[sel],
        (~atmospheric_mask[sel]).astype(float),
        where="mid",
        lw=1.1,
        label="retained",
    )
    ax.axvline(
        center,
        ls="--",
        lw=1.0,
        label="expected CaT center",
    )

    if (
        reference is not None
        and "TELLURIC_TRANSMISSION" in reference.colnames
    ):
        t = np.asarray(
            reference["TELLURIC_TRANSMISSION"][sel],
            float,
        )
        if np.any(np.isfinite(t)):
            ax.plot(
                wave[sel],
                t,
                lw=0.8,
                alpha=0.7,
                label="telluric transmission",
            )

    ax.set_ylim(-0.08, 1.18)
    ax.set_xlabel("Observed wavelength [A]")
    ax.set_ylabel("retained / transmission")
    ax.set_title(
        f"{label.replace('CaT_', 'CaT ')} mask audit"
    )
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main():
    args = parser().parse_args()

    if args.redshift < 0:
        raise ValueError("--redshift must be non-negative.")

    if args.core_half_width_angstrom <= 0:
        raise ValueError("--core-half-width-angstrom must be positive.")

    if (
        args.context_half_width_angstrom
        <= args.core_half_width_angstrom
    ):
        raise ValueError(
            "--context-half-width-angstrom must exceed the core half-width."
        )

    inputs = resolve_inputs(args)

    paths = resolve_stack_paths(
        inputs["red_icube"],
        inputs["red_vcube"],
        inputs["red_mcube"],
        inputs["red_ecube"],
    )
    cube = load_stack(paths)

    atmospheric_mask, mask_wave, mask_hdr = load_1d_mask_fits(
        inputs["mask_path"],
        expected_wave=cube.wavelength,
    )
    wave = cube.wavelength

    if args.science_medium == "auto":
        science_medium = str(
            mask_hdr.get("WAVEMED", "vacuum")
        ).strip().lower()
        if science_medium not in {"air", "vacuum"}:
            science_medium = "vacuum"
    else:
        science_medium = args.science_medium

    reference = None
    if (
        inputs["reference_path"] is not None
        and inputs["reference_path"].is_file()
    ):
        reference = Table.read(
            inputs["reference_path"],
            format="ascii.ecsv",
        )

        if "WAVELENGTH_ANGSTROM" not in reference.colnames:
            raise ValueError(
                "Atmospheric reference table lacks WAVELENGTH_ANGSTROM."
            )

        rw = np.asarray(
            reference["WAVELENGTH_ANGSTROM"],
            float,
        )
        if (
            rw.shape != wave.shape
            or not np.allclose(
                rw,
                wave,
                rtol=0,
                atol=1e-6,
            )
        ):
            raise ValueError(
                "Atmospheric reference ECSV wavelength grid does not "
                "match the final RED stack."
            )

    intervals = None
    if (
        inputs["intervals_path"] is not None
        and inputs["intervals_path"].is_file()
    ):
        intervals = Table.read(
            inputs["intervals_path"],
            format="ascii.ecsv",
        )

    # Convert rest CaT centers into the final science wavelength medium first,
    # then apply the systemic redshift.
    rest_air = np.asarray(
        [w for _, w in CAT_LINES_AIR_ANGSTROM],
        float,
    )

    if args.rest_medium == "air":
        rest_input = rest_air
    else:
        # The canonical constants are AIR. Convert once if the user requests
        # that the audit's rest coordinate be treated as vacuum.
        rest_input = convert_medium(
            rest_air,
            "air",
            "vacuum",
        )

    rest_science_medium = convert_medium(
        rest_input,
        args.rest_medium,
        science_medium,
    )
    observed_centers = (
        rest_science_medium
        * (1.0 + float(args.redshift))
    )

    centers = [
        (label, float(center))
        for (label, _), center
        in zip(
            CAT_LINES_AIR_ANGSTROM,
            observed_centers,
        )
    ]

    rows = []

    for (label, canonical_air), center in zip(
        CAT_LINES_AIR_ANGSTROM,
        observed_centers,
    ):
        center = float(center)

        nearest = int(
            np.nanargmin(
                np.abs(wave - center)
            )
        )
        center_masked = bool(
            atmospheric_mask[nearest]
        )
        center_offset = float(
            wave[nearest] - center
        )

        core = window_metrics(
            wave,
            atmospheric_mask,
            center,
            float(args.core_half_width_angstrom),
            reference,
        )
        context = window_metrics(
            wave,
            atmospheric_mask,
            center,
            float(args.context_half_width_angstrom),
            reference,
        )

        review_reasons = []

        if (
            np.isfinite(core["retained_fraction"])
            and core["retained_fraction"]
            < float(args.review_core_retained_fraction)
        ):
            review_reasons.append(
                "LOW_CORE_RETENTION"
            )

        if (
            np.isfinite(context["retained_fraction"])
            and context["retained_fraction"]
            < float(args.review_context_retained_fraction)
        ):
            review_reasons.append(
                "LOW_CONTEXT_RETENTION"
            )

        if center_masked:
            review_reasons.append(
                "NOMINAL_CENTER_MASKED"
            )

        assessment = (
            "REVIEW"
            if review_reasons
            else "GOOD_RETENTION"
        )

        reasons = overlapping_interval_reasons(
            intervals,
            center - float(args.context_half_width_angstrom),
            center + float(args.context_half_width_angstrom),
        )

        row = {
            "LINE": label,
            "REST_AIR_ANGSTROM": float(canonical_air),
            "OBSERVED_CENTER_ANGSTROM": center,
            "NEAREST_NATIVE_WAVELENGTH_ANGSTROM":
                float(wave[nearest]),
            "CENTER_NATIVE_OFFSET_ANGSTROM":
                center_offset,
            "CENTER_PIXEL_MASKED":
                center_masked,

            "CORE_HALF_WIDTH_ANGSTROM":
                float(args.core_half_width_angstrom),
            "CORE_N_PIXELS":
                core["n_total"],
            "CORE_N_RETAINED":
                core["n_retained"],
            "CORE_RETAINED_FRACTION":
                core["retained_fraction"],
            "CORE_LONGEST_RETAINED_RUN_PIXELS":
                core["longest_retained_run_pixels"],
            "CORE_LONGEST_RETAINED_RUN_ANGSTROM":
                core["longest_retained_run_angstrom"],

            "CONTEXT_HALF_WIDTH_ANGSTROM":
                float(args.context_half_width_angstrom),
            "CONTEXT_N_PIXELS":
                context["n_total"],
            "CONTEXT_N_RETAINED":
                context["n_retained"],
            "CONTEXT_RETAINED_FRACTION":
                context["retained_fraction"],
            "CONTEXT_LONGEST_RETAINED_RUN_PIXELS":
                context["longest_retained_run_pixels"],
            "CONTEXT_LONGEST_RETAINED_RUN_ANGSTROM":
                context["longest_retained_run_angstrom"],

            "CONTEXT_EMPIRICAL_SKY_MASK_FRACTION":
                context["empirical_sky_mask_fraction"],
            "CONTEXT_KNOWN_LINE_MASK_FRACTION":
                context["known_line_mask_fraction"],
            "CONTEXT_TELLURIC_MASK_FRACTION":
                context["telluric_mask_fraction"],
            "CONTEXT_MIN_TELLURIC_TRANSMISSION":
                context["min_telluric_transmission"],
            "CONTEXT_MAX_EMPIRICAL_SKY_Z":
                context["max_empirical_sky_z"],

            "OVERLAPPING_MASK_REASONS":
                reasons,
            "ASSESSMENT":
                assessment,
            "REVIEW_REASONS":
                ";".join(review_reasons),
        }
        rows.append(row)

    # Broad CaT rest-frame window.
    broad_rest = np.asarray(
        args.cat_rest_window,
        float,
    )
    broad_science_rest = convert_medium(
        broad_rest,
        args.rest_medium,
        science_medium,
    )
    broad_obs = broad_science_rest * (
        1.0 + float(args.redshift)
    )
    broad_lo = float(np.min(broad_obs))
    broad_hi = float(np.max(broad_obs))

    broad = window_metrics(
        wave,
        atmospheric_mask,
        0.5 * (broad_lo + broad_hi),
        0.5 * (broad_hi - broad_lo),
        reference,
    )

    out = Path(
        args.output_dir
    ).expanduser().resolve()
    out.mkdir(
        parents=True,
        exist_ok=True,
    )

    table = Table(rows=rows)
    table.meta["PIPELINE"] = "CRD_DRP"
    table.meta["SCRIPT"] = Path(__file__).name
    table.meta["VERSION"] = VERSION
    table.meta["REDSHIFT"] = float(args.redshift)
    table.meta["SCIENCE_WAVELENGTH_MEDIUM"] = science_medium
    table.meta["REST_WAVELENGTH_MEDIUM"] = args.rest_medium
    table.write(
        out / "CaT_mask_audit.ecsv",
        format="ascii.ecsv",
        overwrite=True,
    )

    n_good = sum(
        row["ASSESSMENT"] == "GOOD_RETENTION"
        for row in rows
    )
    n_review = len(rows) - n_good

    summary = {
        "pipeline": "CRD_DRP",
        "script": Path(__file__).name,
        "version": VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": (
            "Optional reduction-QC audit of CaT wavelength retention; "
            "does not replace pPXF convergence/model validation."
        ),
        "redshift": float(args.redshift),
        "science_wavelength_medium": science_medium,
        "rest_wavelength_medium": args.rest_medium,
        "red_icube": str(paths.icube),
        "atmospheric_mask": str(
            inputs["mask_path"]
        ),
        "atmospheric_intervals": (
            str(inputs["intervals_path"])
            if inputs["intervals_path"] is not None
            else None
        ),
        "atmospheric_reference": (
            str(inputs["reference_path"])
            if inputs["reference_path"] is not None
            else None
        ),
        "core_half_width_angstrom":
            float(args.core_half_width_angstrom),
        "context_half_width_angstrom":
            float(args.context_half_width_angstrom),
        "review_core_retained_fraction":
            float(args.review_core_retained_fraction),
        "review_context_retained_fraction":
            float(args.review_context_retained_fraction),

        "cat_rest_window_angstrom":
            [float(v) for v in args.cat_rest_window],
        "cat_observed_window_angstrom":
            [broad_lo, broad_hi],
        "cat_window_n_native_pixels":
            int(broad["n_total"]),
        "cat_window_n_retained_pixels":
            int(broad["n_retained"]),
        "cat_window_retained_fraction":
            float(broad["retained_fraction"]),
        "cat_window_masked_fraction":
            float(broad["masked_fraction"]),
        "cat_window_longest_retained_run_angstrom":
            float(broad["longest_retained_run_angstrom"]),

        "n_lines_good_retention": int(n_good),
        "n_lines_flagged_for_review": int(n_review),
        "lines": rows,
        "interpretation_note": (
            "GOOD_RETENTION/REVIEW are mask-geometry flags only. "
            "Scientific use of the CaT for one- or two-component pPXF "
            "must still be established by the CRD_DAP fit diagnostics."
        ),
    }

    (out / "CaT_mask_audit_summary.json").write_text(
        json.dumps(
            summary,
            indent=2,
        )
        + "\n"
    )

    # Plot a little wider than the broad audit window so boundaries are visible.
    pad = 20.0
    plot_overview(
        wave,
        atmospheric_mask,
        centers,
        max(float(wave[0]), broad_lo - pad),
        min(float(wave[-1]), broad_hi + pad),
        out / "CaT_mask_overview.png",
    )

    for label, center in centers:
        plot_line_zoom(
            label,
            wave,
            atmospheric_mask,
            center,
            float(args.context_half_width_angstrom),
            reference,
            out / f"{label}_mask_zoom.png",
        )

    print("\nCRD_DRP OPTIONAL CaT MASK AUDIT")
    print("=" * 72)
    print(
        f"science medium: {science_medium}; redshift={args.redshift:.6f}"
    )
    print(
        f"broad CaT window retained: "
        f"{broad['n_retained']}/{broad['n_total']} "
        f"({100*broad['retained_fraction']:.1f}%)"
    )

    for row in rows:
        print(
            f"{row['LINE']}: center={row['OBSERVED_CENTER_ANGSTROM']:.2f} A | "
            f"core retained={100*row['CORE_RETAINED_FRACTION']:.1f}% | "
            f"context retained={100*row['CONTEXT_RETAINED_FRACTION']:.1f}% | "
            f"center masked={row['CENTER_PIXEL_MASKED']} | "
            f"{row['ASSESSMENT']}"
        )

    print("=" * 72)
    print(
        f"Outputs: {out}"
    )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as exc:
        print(
            f"ERROR: {type(exc).__name__}: {exc}",
            file=sys.stderr,
        )
        raise SystemExit(2)
