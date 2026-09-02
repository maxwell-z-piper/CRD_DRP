#!/usr/bin/env python3
"""
Build fixed observed-frame BLUE/RED atmospheric masks for CRD_DRP.

Mask-selection inputs are independent of CRD_DAP/pPXF residuals:

  1. empirical sky-emission structure from reduction-stage reference cubes;
  2. known bright sky-emission lines;
  3. optional modeled telluric transmission, normally the RED
     KSkyWizard/PypeIt *_tellmodel.fits product.

The telluric handling is ported directly from the old CRD_DAP Script 03d:
the model is converted into the final stack's wavelength medium, interpolated
onto the exact native KcwiKit wavelength grid, and transmission below a fixed
threshold is independently masked.

The final science cubes are never modified.
"""

from __future__ import annotations

import argparse
import json
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
    extract_spatial_median_spectrum,
    robust_highpass_z,
    infer_wavelength_medium,
    convert_medium,
    bridge_small_false_gaps,
    dilate_1d,
    wavelength_edges,
    save_1d_mask_fits,
    resolve_path,
    load_telluric_spectrum,
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
            f"--{arm}-reference-cube",
            action="append",
            default=[],
            help=(
                "Independent reduction-stage cube used to identify sky-emission "
                "structure. Repeat for multiple exposures."
            ),
        )
        p.add_argument(
            f"--{arm}-reference-object-mask",
            action="append",
            default=[],
            help=(
                "Optional 2-D object mask(s), object>0 and sky=0. Supply none, "
                "one shared mask, or one per reference cube."
            ),
        )
        p.add_argument(
            f"--{arm}-science-medium",
            choices=("auto", "air", "vacuum"),
            default="auto",
        )
        p.add_argument(
            f"--{arm}-reference-medium",
            choices=("auto", "air", "vacuum"),
            default="auto",
        )

        p.add_argument(
            f"--{arm}-telluric-spectrum",
            default=None,
            help=(
                "Optional 1-D telluric transmission model. For RED, use the "
                "KSkyWizard/PypeIt *_tellmodel.fits produced from the standard "
                "star. BLUE normally leaves this unset."
            ),
        )
        p.add_argument(
            f"--{arm}-telluric-medium",
            choices=("air", "vacuum"),
            default="air",
            help=(
                "Fallback wavelength medium for the telluric product when its "
                "header/table does not declare one (default: air, matching old 03d)."
            ),
        )

    p.add_argument(
        "--known-line-list",
        default=None,
        help="Optional ECSV/ASCII table adding known sky-emission lines.",
    )
    p.add_argument("--sky-line-sigma", type=float, default=8.0)
    p.add_argument("--sky-reference-min-fraction", type=float, default=0.50)
    p.add_argument("--continuum-window-angstrom", type=float, default=75.0)
    p.add_argument("--sky-padding-pixels", type=int, default=2)
    p.add_argument("--bridge-pixels", type=int, default=1)
    p.add_argument("--known-line-half-width-angstrom", type=float, default=3.5)

    # Ported from old CRD_DAP 03d.
    p.add_argument(
        "--telluric-threshold",
        type=float,
        default=0.97,
        help=(
            "Mask native wavelengths with modeled atmospheric transmission below "
            "this value (default 0.97, matching old CRD_DAP 03d)."
        ),
    )
    p.add_argument(
        "--telluric-padding-pixels",
        type=int,
        default=2,
        help="Dilate telluric-mask regions by this many native pixels (default 2).",
    )

    p.add_argument("--output-dir", default="atmospheric_masks")
    return p


def built_in_lines():
    # Canonical bright optical night-sky [O I] lines, AIR wavelengths.
    return [
        (5577.338, "[O I] 5577", "air"),
        (6300.304, "[O I] 6300", "air"),
        (6363.776, "[O I] 6364", "air"),
    ]


def load_extra_lines(path):
    if path is None:
        return []

    t = Table.read(path)
    low = {n.lower(): n for n in t.colnames}

    wcol = next(
        (v for k, v in low.items() if "wave" in k or "lambda" in k),
        None,
    )
    lcol = next(
        (v for k, v in low.items() if "label" in k or "name" in k),
        None,
    )
    mcol = next(
        (v for k, v in low.items() if "medium" in k),
        None,
    )

    if wcol is None:
        raise ValueError("Known-line table needs a wavelength column.")

    return [
        (
            float(r[wcol]),
            str(r[lcol]) if lcol else "user sky line",
            str(r[mcol]).lower() if mcol else "air",
        )
        for r in t
    ]


def resolve_masks(refs, masks):
    if not masks:
        return [None] * len(refs)
    if len(masks) == 1:
        return [resolve_path(masks[0])] * len(refs)
    if len(masks) != len(refs):
        raise ValueError(
            "Supply zero, one shared object mask, or one mask per reference cube."
        )
    return [resolve_path(m) for m in masks]


def make_interval_table(
    wave,
    final_mask,
    empirical_mask,
    known_mask,
    telluric_mask,
    known_labels,
    empirical_z,
    empirical_fraction,
    reference_coverage,
    telluric_transmission,
):
    edges = wavelength_edges(wave)
    idx = np.flatnonzero(final_mask)
    rows = []

    if idx.size:
        starts = np.r_[idx[0], idx[1:][np.diff(idx) > 1]]
        ends = np.r_[idx[:-1][np.diff(idx) > 1], idx[-1]]

        for iid, (i0, i1) in enumerate(zip(starts, ends), 1):
            sl = slice(i0, i1 + 1)

            labels = sorted(
                {
                    label
                    for labels_here in known_labels[sl]
                    for label in labels_here
                }
            )

            reasons = []
            if np.any(empirical_mask[sl]):
                reasons.append("EMPIRICAL_SKY_EMISSION")
            if np.any(known_mask[sl]):
                reasons.append("KNOWN_SKY_LINE")
            if np.any(telluric_mask[sl]):
                reasons.append("TELLURIC_REFERENCE")

            rows.append(
                {
                    "INTERVAL_ID": iid,
                    "INDEX_LO": int(i0),
                    "INDEX_HI": int(i1),
                    "N_PIXELS": int(i1 - i0 + 1),
                    "WAVE_LO_ANGSTROM": float(edges[i0]),
                    "WAVE_HI_ANGSTROM": float(edges[i1 + 1]),
                    "REASON": "+".join(reasons) if reasons else "ATMOSPHERIC_REFERENCE",
                    "KNOWN_LINES": "; ".join(labels),
                    "PEAK_EMPIRICAL_SKY_Z": (
                        float(np.nanmax(empirical_z[sl]))
                        if np.any(np.isfinite(empirical_z[sl]))
                        else np.nan
                    ),
                    "PEAK_SKY_REFERENCE_DETECTION_FRACTION": (
                        float(np.nanmax(empirical_fraction[sl]))
                        if np.any(np.isfinite(empirical_fraction[sl]))
                        else np.nan
                    ),
                    "MAX_REFERENCE_COVERAGE_COUNT": (
                        int(np.nanmax(reference_coverage[sl]))
                        if reference_coverage.size
                        else 0
                    ),
                    "MIN_TELLURIC_TRANSMISSION": (
                        float(np.nanmin(telluric_transmission[sl]))
                        if telluric_transmission is not None
                        and np.any(np.isfinite(telluric_transmission[sl]))
                        else np.nan
                    ),
                }
            )

    tab = Table(rows=rows)
    tab.meta["PIPELINE"] = "CRD_DRP"
    tab.meta["MASK_SELECTION_BASIS"] = (
        "REDUCTION_STAGE_SKY_PLUS_KNOWN_LINES_PLUS_OPTIONAL_TELLURIC_MODEL"
    )
    return tab


def plot_qc(
    arm,
    wave,
    empirical_z,
    empirical_fraction,
    coverage_fraction,
    telluric_transmission,
    final_mask,
    telluric_threshold,
    path,
):
    nrows = 4 if telluric_transmission is not None else 3
    fig, axes = plt.subplots(
        nrows,
        1,
        figsize=(13, 10 if nrows == 4 else 8),
        sharex=True,
    )

    axes[0].plot(wave, empirical_z, lw=0.8)
    axes[0].axhline(0.0, ls=":", lw=0.7)
    axes[0].set_ylabel("empirical sky z")

    axes[1].plot(wave, empirical_fraction, lw=0.8)
    axes[1].set_ylim(-0.05, 1.05)
    axes[1].set_ylabel("detection fraction")

    axes[2].plot(wave, coverage_fraction, lw=0.8)
    axes[2].set_ylim(-0.05, 1.05)
    axes[2].set_ylabel("reference coverage")

    if telluric_transmission is not None:
        axes[3].plot(wave, telluric_transmission, lw=0.8)
        axes[3].axhline(telluric_threshold, ls="--", lw=0.8)
        axes[3].set_ylim(-0.05, 1.05)
        axes[3].set_ylabel("telluric T")
        bottom = axes[3]
    else:
        bottom = axes[2]

    for ax in axes:
        y0, y1 = ax.get_ylim()
        ax.fill_between(
            wave,
            y0,
            y1,
            where=final_mask,
            alpha=0.18,
        )

    bottom.set_xlabel("Observed wavelength [A]")
    fig.suptitle(
        f"{arm.upper()} atmospheric mask | "
        f"{np.count_nonzero(final_mask)}/{final_mask.size} "
        f"({100*np.mean(final_mask):.2f}%) native pixels masked"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=170)
    plt.close(fig)


def build_arm(arm, args, lines, out):
    paths = resolve_stack_paths(
        getattr(args, f"{arm}_icube"),
        getattr(args, f"{arm}_vcube"),
        getattr(args, f"{arm}_mcube"),
        getattr(args, f"{arm}_ecube"),
    )
    cube = load_stack(paths)
    wave = cube.wavelength

    requested_science_medium = getattr(args, f"{arm}_science_medium")
    science_medium = (
        infer_wavelength_medium(cube.header, "vacuum")
        if requested_science_medium == "auto"
        else requested_science_medium
    )

    refs = [
        resolve_path(p)
        for p in getattr(args, f"{arm}_reference_cube")
    ]
    object_masks = resolve_masks(
        refs,
        getattr(args, f"{arm}_reference_object_mask"),
    )

    zrefs = []
    ref_meta = []

    for p, omask in zip(refs, object_masks):
        rw, rs, hdr, nsky = extract_spatial_median_spectrum(p, omask)

        requested_reference_medium = getattr(args, f"{arm}_reference_medium")
        reference_medium = (
            infer_wavelength_medium(hdr, science_medium)
            if requested_reference_medium == "auto"
            else requested_reference_medium
        )

        rw = convert_medium(rw, reference_medium, science_medium)

        try:
            _, z, sigma = robust_highpass_z(
                rw,
                rs,
                args.continuum_window_angstrom,
            )
        except ValueError as exc:
            ref_meta.append(
                {
                    "path": str(p),
                    "object_mask": str(omask) if omask else None,
                    "n_sky_spaxels": int(nsky),
                    "medium": reference_medium,
                    "status": "skipped",
                    "reason": str(exc),
                }
            )
            continue

        z_on_science = np.interp(
            wave,
            rw,
            z,
            left=np.nan,
            right=np.nan,
        )
        z_on_science[(wave < np.nanmin(rw)) | (wave > np.nanmax(rw))] = np.nan

        zrefs.append(z_on_science)
        ref_meta.append(
            {
                "path": str(p),
                "object_mask": str(omask) if omask else None,
                "n_sky_spaxels": int(nsky),
                "medium": reference_medium,
                "status": "used",
                "robust_sigma": float(sigma),
            }
        )

    n_reference_requested = len(refs)
    n_reference_used = len(zrefs)

    if zrefs:
        zarr = np.vstack(zrefs)
        finite_ref = np.isfinite(zarr)

        coverage_count = np.sum(finite_ref, axis=0).astype(int)
        detection_count = np.sum(
            finite_ref & (zarr >= float(args.sky_line_sigma)),
            axis=0,
        ).astype(int)

        detection_fraction = np.full(wave.size, np.nan, dtype=float)
        covered = coverage_count > 0
        detection_fraction[covered] = (
            detection_count[covered] / coverage_count[covered]
        )

        empirical_z = np.full(wave.size, np.nan, dtype=float)
        for j in np.flatnonzero(covered):
            empirical_z[j] = float(np.nanmedian(zarr[:, j]))

        empirical_mask = (
            covered
            & np.isfinite(detection_fraction)
            & (
                detection_fraction
                >= float(args.sky_reference_min_fraction)
            )
        )
    else:
        coverage_count = np.zeros(wave.size, dtype=int)
        detection_count = np.zeros(wave.size, dtype=int)
        detection_fraction = np.full(wave.size, np.nan, dtype=float)
        empirical_z = np.full(wave.size, np.nan, dtype=float)
        empirical_mask = np.zeros(wave.size, dtype=bool)

    coverage_fraction = (
        coverage_count / max(n_reference_used, 1)
        if n_reference_used > 0
        else np.zeros(wave.size, dtype=float)
    )

    empirical_mask = bridge_small_false_gaps(
        empirical_mask,
        args.bridge_pixels,
    )
    empirical_mask = dilate_1d(
        empirical_mask,
        args.sky_padding_pixels,
    )

    known_mask = np.zeros_like(wave, dtype=bool)
    known_labels = np.empty(wave.size, dtype=object)
    for i in range(wave.size):
        known_labels[i] = []

    for line_wave, label, line_medium in lines:
        center = float(
            convert_medium(
                np.asarray([line_wave]),
                line_medium,
                science_medium,
            )[0]
        )
        sel = (
            np.abs(wave - center)
            <= float(args.known_line_half_width_angstrom)
        )
        known_mask |= sel
        for i in np.flatnonzero(sel):
            known_labels[i].append(label)

    # ------------------------------------------------------------------
    # Telluric absorption: direct port of the old CRD_DAP 03d logic.
    # ------------------------------------------------------------------
    telluric_path_arg = getattr(args, f"{arm}_telluric_spectrum")
    telluric_transmission = None
    telluric_mask = np.zeros(wave.size, dtype=bool)
    telluric_meta = None

    if telluric_path_arg is not None:
        telluric_path = resolve_path(telluric_path_arg)
        fallback_medium = getattr(args, f"{arm}_telluric_medium")

        tw, tt, telluric_medium = load_telluric_spectrum(
            telluric_path,
            fallback_medium=fallback_medium,
        )
        tw_science = convert_medium(
            tw,
            telluric_medium,
            science_medium,
        )

        telluric_transmission = np.interp(
            wave,
            tw_science,
            tt,
            left=np.nan,
            right=np.nan,
        )
        telluric_transmission[
            (wave < np.nanmin(tw_science))
            | (wave > np.nanmax(tw_science))
        ] = np.nan

        telluric_mask = (
            np.isfinite(telluric_transmission)
            & (
                telluric_transmission
                < float(args.telluric_threshold)
            )
        )
        telluric_mask = bridge_small_false_gaps(
            telluric_mask,
            args.bridge_pixels,
        )
        telluric_mask = dilate_1d(
            telluric_mask,
            args.telluric_padding_pixels,
        )

        telluric_meta = {
            "path": str(telluric_path),
            "input_medium": telluric_medium,
            "science_medium": science_medium,
            "threshold": float(args.telluric_threshold),
            "padding_native_pixels": int(args.telluric_padding_pixels),
            "n_masked_native_pixels": int(np.count_nonzero(telluric_mask)),
        }

    final_mask = empirical_mask | known_mask | telluric_mask

    mask_path = out / f"{arm}_atmospheric_mask.fits"
    interval_path = out / f"{arm}_atmospheric_mask.ecsv"
    reference_path = out / f"{arm}_atmospheric_reference.ecsv"
    plot_path = out / f"{arm}_atmospheric_mask_qc.png"

    save_1d_mask_fits(
        mask_path,
        final_mask,
        wave,
        arm,
        science_medium,
        paths.icube,
    )

    tab = make_interval_table(
        wave=wave,
        final_mask=final_mask,
        empirical_mask=empirical_mask,
        known_mask=known_mask,
        telluric_mask=telluric_mask,
        known_labels=known_labels,
        empirical_z=empirical_z,
        empirical_fraction=detection_fraction,
        reference_coverage=coverage_count,
        telluric_transmission=telluric_transmission,
    )
    tab.meta["ARM"] = arm.upper()
    tab.meta["WAVELENGTH_MEDIUM"] = science_medium
    tab.meta["SKY_LINE_SIGMA_THRESHOLD"] = float(args.sky_line_sigma)
    tab.meta["SKY_REFERENCE_MIN_FRACTION"] = float(
        args.sky_reference_min_fraction
    )
    tab.meta["TELLURIC_THRESHOLD"] = float(args.telluric_threshold)
    tab.write(
        interval_path,
        format="ascii.ecsv",
        overwrite=True,
    )

    ref = Table()
    ref["WAVELENGTH_ANGSTROM"] = wave
    ref["EMPIRICAL_SKY_Z"] = empirical_z
    ref["N_REFERENCE_COVERAGE"] = coverage_count
    ref["REFERENCE_COVERAGE_FRACTION"] = coverage_fraction
    ref["N_REFERENCE_DETECTIONS"] = detection_count
    ref["REFERENCE_DETECTION_FRACTION"] = detection_fraction
    ref["EMPIRICAL_SKY_MASK"] = empirical_mask
    ref["KNOWN_LINE_MASK"] = known_mask
    ref["TELLURIC_MASK"] = telluric_mask
    ref["FINAL_MASK"] = final_mask

    if telluric_transmission is not None:
        ref["TELLURIC_TRANSMISSION"] = telluric_transmission

    ref.meta["PIPELINE"] = "CRD_DRP"
    ref.meta["ARM"] = arm.upper()
    ref.meta["WAVELENGTH_MEDIUM"] = science_medium
    ref.meta["N_REFERENCE_REQUESTED"] = int(n_reference_requested)
    ref.meta["N_REFERENCE_USED"] = int(n_reference_used)
    ref.write(
        reference_path,
        format="ascii.ecsv",
        overwrite=True,
    )

    plot_qc(
        arm=arm,
        wave=wave,
        empirical_z=empirical_z,
        empirical_fraction=detection_fraction,
        coverage_fraction=coverage_fraction,
        telluric_transmission=telluric_transmission,
        final_mask=final_mask,
        telluric_threshold=float(args.telluric_threshold),
        path=plot_path,
    )

    return {
        "arm": arm.upper(),
        "science_medium": science_medium,
        "source_stack": {
            k: str(getattr(paths, k)) if getattr(paths, k) else None
            for k in ("icube", "vcube", "mcube", "ecube")
        },
        "reference_cubes": ref_meta,
        "n_reference_requested": int(n_reference_requested),
        "n_reference_used": int(n_reference_used),
        "telluric_reference": telluric_meta,
        "mask_fits": str(mask_path.resolve()),
        "mask_intervals_ecsv": str(interval_path.resolve()),
        "mask_reference_ecsv": str(reference_path.resolve()),
        "mask_qc_png": str(plot_path.resolve()),
        "n_native_pixels": int(final_mask.size),
        "n_empirical_sky_masked_pixels": int(
            np.count_nonzero(empirical_mask)
        ),
        "n_known_line_masked_pixels": int(
            np.count_nonzero(known_mask)
        ),
        "n_telluric_masked_pixels": int(
            np.count_nonzero(telluric_mask)
        ),
        "n_masked_pixels": int(np.count_nonzero(final_mask)),
        "masked_fraction": float(np.mean(final_mask)),
        "n_intervals": int(len(tab)),
    }


def main():
    args = parser().parse_args()
    out = Path(args.output_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    lines = built_in_lines() + load_extra_lines(
        args.known_line_list
    )

    arms = {
        arm: build_arm(arm, args, lines, out)
        for arm in ("blue", "red")
    }

    payload = {
        "pipeline": "CRD_DRP",
        "script": Path(__file__).name,
        "version": VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selection_principle": (
            "Fixed observed-frame masks constructed independently of "
            "CRD_DAP/pPXF residuals from reduction-stage sky emission, "
            "known sky lines, and optional modeled telluric transmission."
        ),
        "arms": arms,
        "parameters": vars(args),
    }
    (out / "atmospheric_mask_summary.json").write_text(
        json.dumps(payload, indent=2) + "\n"
    )

    for arm in ("blue", "red"):
        s = arms[arm]
        tell = s["n_telluric_masked_pixels"]
        print(
            f"{arm.upper()}: "
            f"{s['n_masked_pixels']}/{s['n_native_pixels']} "
            f"({100*s['masked_fraction']:.2f}%) masked; "
            f"telluric component={tell} native pixels"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
