#!/usr/bin/env python3
"""
compare_standard_invsens.py

Compare multiple KCWI/KCRM *_invsens.fits files produced from repeated
exposures of the same spectrophotometric standard.

For each exposure this script measures:

  1. DRP calibration residuals reconstructed from:
         100 * (fitted_invsens / raw_invsens - 1)

  2. Robust scatter of those residuals.

  3. Fraction of wavelengths with large residuals.

  4. Agreement of the fitted sensitivity-curve SHAPE with the ensemble
     median of all standard exposures.

  5. Relative normalization / throughput of each standard exposure.

Outputs:
  standard_invsens_QC.csv
  standard_invsens_curves.png
  standard_invsens_residuals.png
  standard_invsens_shape_comparison.png

Usage examples
--------------
Explicit files:

python compare_standard_invsens.py \
    redux/kr241226_00150_invsens.fits \
    redux/kr241226_00151_invsens.fits \
    redux/kr241226_00152_invsens.fits

or with a glob:

python compare_standard_invsens.py \
    "redux/kr241226_0015[0-2]_invsens.fits"
"""

from pathlib import Path
import argparse
import glob
import csv

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


# ============================================================
# Utility functions
# ============================================================

def robust_sigma(x):
    """1.4826 * MAD about the median."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]

    if x.size == 0:
        return np.nan

    med = np.nanmedian(x)
    return 1.4826 * np.nanmedian(np.abs(x - med))


def wavelength_from_header(header, npix):
    """
    Construct the wavelength axis using FITS 1-indexed CRPIX convention.
    """
    crval = float(header["CRVAL1"])
    crpix = float(header["CRPIX1"])
    cdelt = float(header["CDELT1"])

    pix_fits = np.arange(npix, dtype=float) + 1.0

    wave = crval + (pix_fits - crpix) * cdelt

    return wave


def load_invsens(filename):
    """
    Read a KCWI DRP inverse-sensitivity FITS file.

    Expected primary array:
        [0] raw inverse sensitivity
        [1] fitted inverse sensitivity
        [2] observed standard spectrum
    """

    with fits.open(filename) as hdul:
        data = np.asarray(hdul[0].data, dtype=float)
        hdr = hdul[0].header.copy()

    if data.ndim != 2:
        raise RuntimeError(
            "{}: expected a 2-D invsens array, got shape {}".format(
                filename, data.shape
            )
        )

    # Normally shape = (3, Nwave).
    # Protect against an unexpected transpose.
    if data.shape[0] == 3:
        pass
    elif data.shape[1] == 3:
        data = data.T
    else:
        raise RuntimeError(
            "{}: expected one dimension of length 3, got {}".format(
                filename, data.shape
            )
        )

    raw_invsens = data[0]
    fit_invsens = data[1]
    obs_spec = data[2]

    wave = wavelength_from_header(hdr, raw_invsens.size)

    # DRP fitting limits
    fit_w0 = float(hdr.get("INVFW0", np.nanmin(wave)))
    fit_w1 = float(hdr.get("INVFW1", np.nanmax(wave)))

    # Also retain the nominal valid inverse-sensitivity interval.
    good_w0 = float(hdr.get("INVSW0", fit_w0))
    good_w1 = float(hdr.get("INVSW1", fit_w1))

    return {
        "filename": str(filename),
        "name": Path(filename).name,
        "header": hdr,
        "wave": wave,
        "raw": raw_invsens,
        "fit": fit_invsens,
        "obs": obs_spec,
        "fit_w0": fit_w0,
        "fit_w1": fit_w1,
        "good_w0": good_w0,
        "good_w1": good_w1,
    }


# ============================================================
# File handling
# ============================================================

def expand_files(arguments):
    """
    Expand literal filenames and shell-style glob patterns.
    """

    files = []

    for item in arguments:
        matches = sorted(glob.glob(item))

        if matches:
            files.extend(matches)
        elif Path(item).exists():
            files.append(item)

    # remove duplicates while preserving order
    unique = []

    for f in files:
        if f not in unique:
            unique.append(f)

    return unique


# ============================================================
# Main analysis
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "files",
        nargs="+",
        help="KCWI *_invsens.fits files or glob patterns",
    )

    parser.add_argument(
        "--out-prefix",
        default="standard_invsens",
        help="Prefix for output files",
    )

    args = parser.parse_args()

    filenames = expand_files(args.files)

    if len(filenames) < 2:
        raise RuntimeError(
            "Need at least two inverse-sensitivity files."
        )

    print()
    print("Loading {} standard-star invsens files".format(len(filenames)))
    print()

    standards = [load_invsens(f) for f in filenames]

    # --------------------------------------------------------
    # Determine common wavelength interval
    # --------------------------------------------------------

    common_lo = max(s["fit_w0"] for s in standards)
    common_hi = min(s["fit_w1"] for s in standards)

    if common_hi <= common_lo:
        raise RuntimeError(
            "The inverse-sensitivity fitting intervals do not overlap."
        )

    print(
        "Common DRP fitting interval: {:.2f} - {:.2f} Angstrom".format(
            common_lo, common_hi
        )
    )
    print()

    # --------------------------------------------------------
    # Define common wavelength grid
    # --------------------------------------------------------

    # Use the first standard's wavelength sampling within common interval.
    wref = standards[0]["wave"]

    common = (
        np.isfinite(wref)
        & (wref >= common_lo)
        & (wref <= common_hi)
    )

    common_wave = wref[common]

    if common_wave.size < 10:
        raise RuntimeError("Too few common wavelength pixels.")

    # Interpolate all fitted curves onto common wavelength grid.
    interp_fits = []
    interp_raws = []
    interp_obs = []

    for s in standards:

        good_fit = (
            np.isfinite(s["wave"])
            & np.isfinite(s["fit"])
        )

        good_raw = (
            np.isfinite(s["wave"])
            & np.isfinite(s["raw"])
        )

        good_obs = (
            np.isfinite(s["wave"])
            & np.isfinite(s["obs"])
        )

        fit_i = np.interp(
            common_wave,
            s["wave"][good_fit],
            s["fit"][good_fit],
            left=np.nan,
            right=np.nan,
        )

        raw_i = np.interp(
            common_wave,
            s["wave"][good_raw],
            s["raw"][good_raw],
            left=np.nan,
            right=np.nan,
        )

        obs_i = np.interp(
            common_wave,
            s["wave"][good_obs],
            s["obs"][good_obs],
            left=np.nan,
            right=np.nan,
        )

        interp_fits.append(fit_i)
        interp_raws.append(raw_i)
        interp_obs.append(obs_i)

    interp_fits = np.asarray(interp_fits)
    interp_raws = np.asarray(interp_raws)
    interp_obs = np.asarray(interp_obs)

    # --------------------------------------------------------
    # Normalize fitted curves to compare SHAPE separately
    # from absolute throughput normalization.
    # --------------------------------------------------------

    fit_norm = np.full_like(interp_fits, np.nan)

    fit_medians = []

    for i in range(len(standards)):

        good = (
            np.isfinite(interp_fits[i])
            & (interp_fits[i] > 0)
        )

        med = np.nanmedian(interp_fits[i][good])
        fit_medians.append(med)

        fit_norm[i] = interp_fits[i] / med

    fit_medians = np.asarray(fit_medians)

    ensemble_shape = np.nanmedian(fit_norm, axis=0)

    # Median normalization among exposures.
    ensemble_scale = np.nanmedian(fit_medians)

    # --------------------------------------------------------
    # Calculate QC metrics
    # --------------------------------------------------------

    rows = []

    for i, s in enumerate(standards):

        raw = interp_raws[i]
        fit = interp_fits[i]
        obs = interp_obs[i]

        valid = (
            np.isfinite(raw)
            & np.isfinite(fit)
            & np.isfinite(obs)
            & (raw != 0)
            & (fit > 0)
        )

        # Exact calibration residual implied by DRP invsens product:
        #
        # reference = observed * raw_invsens
        # calibrated = observed * fitted_invsens
        #
        residual = np.full_like(raw, np.nan)

        residual[valid] = (
            100.0 * (fit[valid] / raw[valid] - 1.0)
        )

        r = residual[valid]

        resid_mean = np.nanmean(r)
        resid_std = np.nanstd(r)
        resid_median = np.nanmedian(r)
        resid_robust = robust_sigma(r)

        abs_r = np.abs(r)

        p95_abs = np.nanpercentile(abs_r, 95.0)
        p99_abs = np.nanpercentile(abs_r, 99.0)

        frac_gt_5 = np.mean(abs_r > 5.0)
        frac_gt_10 = np.mean(abs_r > 10.0)
        frac_gt_20 = np.mean(abs_r > 20.0)

        # Shape difference from ensemble median.
        shape_ratio = fit_norm[i] / ensemble_shape

        shape_good = (
            np.isfinite(shape_ratio)
            & (ensemble_shape != 0)
        )

        shape_delta_pct = (
            100.0 * (shape_ratio[shape_good] - 1.0)
        )

        shape_med_offset = np.nanmedian(shape_delta_pct)
        shape_robust = robust_sigma(shape_delta_pct)
        shape_p95 = np.nanpercentile(
            np.abs(shape_delta_pct), 95.0
        )

        # Absolute sensitivity normalization relative to the median
        # of the three standards.
        #
        # >1 means larger inverse sensitivity = lower measured throughput.
        scale_ratio = fit_medians[i] / ensemble_scale

        # Observed standard counts are useful as an independent
        # relative-throughput diagnostic because all three exposures
        # have the same exposure time.
        obs_good = (
            np.isfinite(obs)
            & (obs > 0)
        )

        obs_median = np.nanmedian(obs[obs_good])

        # Header information
        hdr = s["header"]

        frame = hdr.get("FRAMENO", "")
        obj = hdr.get("OBJECT", hdr.get("TARGNAME", ""))
        airmass = hdr.get("AIRMASS", np.nan)
        dateobs = hdr.get("DATE-OBS", "")
        utc = hdr.get("UTC", hdr.get("UT", ""))

        row = {
            "filename": s["name"],
            "frame": frame,
            "object": obj,
            "date_obs": dateobs,
            "utc": utc,
            "airmass": airmass,
            "fit_wave_min": common_lo,
            "fit_wave_max": common_hi,
            "n_valid": int(np.sum(valid)),
            "resid_mean_pct": resid_mean,
            "resid_std_pct": resid_std,
            "resid_median_pct": resid_median,
            "resid_robust_sigma_pct": resid_robust,
            "resid_p95_abs_pct": p95_abs,
            "resid_p99_abs_pct": p99_abs,
            "frac_abs_resid_gt_5pct": frac_gt_5,
            "frac_abs_resid_gt_10pct": frac_gt_10,
            "frac_abs_resid_gt_20pct": frac_gt_20,
            "shape_median_offset_pct": shape_med_offset,
            "shape_robust_sigma_pct": shape_robust,
            "shape_p95_abs_pct": shape_p95,
            "invsens_scale_ratio_to_median": scale_ratio,
            "median_observed_standard_counts": obs_median,
        }

        rows.append(row)

        # Save residual vector for plots.
        s["common_residual"] = residual
        s["common_fit"] = fit
        s["common_raw"] = raw
        s["common_obs"] = obs
        s["shape_delta_pct"] = np.full_like(
            common_wave, np.nan
        )
        s["shape_delta_pct"][shape_good] = shape_delta_pct

    # --------------------------------------------------------
    # Construct a simple ranking
    # --------------------------------------------------------
    #
    # We deliberately do NOT rank by absolute sensitivity
    # normalization. A different normalization can reflect a real
    # throughput change and should be inspected, not automatically
    # penalized.
    #
    # Ranking uses:
    #   robust calibration residual
    #   95th percentile calibration residual
    #   curve-shape agreement with other standards
    #

    metric_names = [
        "resid_robust_sigma_pct",
        "resid_p95_abs_pct",
        "shape_robust_sigma_pct",
        "shape_p95_abs_pct",
    ]

    for row in rows:
        row["rank_score"] = 0

    for metric in metric_names:

        values = np.array(
            [row[metric] for row in rows],
            dtype=float
        )

        order = np.argsort(values)

        for rank, index in enumerate(order):
            rows[index]["rank_score"] += rank

    rows_sorted = sorted(
        rows,
        key=lambda r: r["rank_score"]
    )

    # --------------------------------------------------------
    # Print summary
    # --------------------------------------------------------

    print("=" * 100)
    print("STANDARD STAR INVERSE-SENSITIVITY QC")
    print("=" * 100)

    for row in rows_sorted:

        print()
        print(row["filename"])

        print(
            "  frame                     : {}".format(
                row["frame"]
            )
        )

        print(
            "  airmass                   : {}".format(
                row["airmass"]
            )
        )

        print(
            "  residual mean/std         : "
            "{:+.3f} / {:.3f} %".format(
                row["resid_mean_pct"],
                row["resid_std_pct"],
            )
        )

        print(
            "  residual robust sigma     : {:.3f} %".format(
                row["resid_robust_sigma_pct"]
            )
        )

        print(
            "  |residual| 95th percentile: {:.3f} %".format(
                row["resid_p95_abs_pct"]
            )
        )

        print(
            "  fraction |resid| > 5%     : {:.3%}".format(
                row["frac_abs_resid_gt_5pct"]
            )
        )

        print(
            "  fraction |resid| > 10%    : {:.3%}".format(
                row["frac_abs_resid_gt_10pct"]
            )
        )

        print(
            "  shape robust sigma        : {:.3f} %".format(
                row["shape_robust_sigma_pct"]
            )
        )

        print(
            "  shape |delta| p95         : {:.3f} %".format(
                row["shape_p95_abs_pct"]
            )
        )

        print(
            "  invsens scale / median    : {:.5f}".format(
                row["invsens_scale_ratio_to_median"]
            )
        )

        print(
            "  observed median counts    : {:.6g}".format(
                row["median_observed_standard_counts"]
            )
        )

        print(
            "  QC rank score             : {}".format(
                row["rank_score"]
            )
        )

    print()
    print("=" * 100)
    print(
        "Lowest numerical QC score: {}".format(
            rows_sorted[0]["filename"]
        )
    )
    print("=" * 100)
    print()

    print(
        "IMPORTANT: the lowest score is a QC recommendation, "
        "not an automatic scientific decision."
    )

    print(
        "Check the plots, especially for telluric-region artifacts "
        "or an exposure with clearly different throughput."
    )

    # --------------------------------------------------------
    # Write CSV
    # --------------------------------------------------------

    csv_name = args.out_prefix + "_QC.csv"

    fieldnames = list(rows[0].keys())

    with open(csv_name, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames
        )
        writer.writeheader()
        writer.writerows(rows)

    print()
    print("Saved:", csv_name)

    # --------------------------------------------------------
    # Plot 1: raw + fitted inverse sensitivity
    # --------------------------------------------------------

    plt.figure(figsize=(11, 7))

    for s in standards:

        good = (
            np.isfinite(s["common_raw"])
            & (s["common_raw"] > 0)
        )

        plt.plot(
            common_wave[good],
            s["common_raw"][good],
            alpha=0.30,
            linewidth=0.8,
        )

        good = (
            np.isfinite(s["common_fit"])
            & (s["common_fit"] > 0)
        )

        plt.plot(
            common_wave[good],
            s["common_fit"][good],
            linewidth=1.8,
            label=s["name"],
        )

    plt.xlabel("Wavelength (Angstrom)")
    plt.ylabel("Inverse sensitivity")
    plt.title(
        "KCWI standard-star inverse sensitivity\n"
        "faint = raw; solid = DRP fitted curve"
    )
    plt.legend(fontsize=8)
    plt.tight_layout()

    fn = args.out_prefix + "_curves.png"
    plt.savefig(fn, dpi=180)
    plt.close()

    print("Saved:", fn)

    # --------------------------------------------------------
    # Plot 2: reconstructed calibration residuals
    # --------------------------------------------------------

    plt.figure(figsize=(11, 7))

    for s in standards:

        plt.plot(
            common_wave,
            s["common_residual"],
            linewidth=0.9,
            label=s["name"],
        )

    plt.axhline(0.0, linewidth=1.0)

    plt.xlabel("Wavelength (Angstrom)")
    plt.ylabel("(calibrated - reference) / reference (%)")
    plt.title(
        "Reconstructed DRP standard-star calibration residuals"
    )
    plt.legend(fontsize=8)
    plt.tight_layout()

    fn = args.out_prefix + "_residuals.png"
    plt.savefig(fn, dpi=180)
    plt.close()

    print("Saved:", fn)

    # --------------------------------------------------------
    # Plot 3: fitted-curve shape comparison
    # --------------------------------------------------------

    plt.figure(figsize=(11, 7))

    for s in standards:

        plt.plot(
            common_wave,
            s["shape_delta_pct"],
            linewidth=1.2,
            label=s["name"],
        )

    plt.axhline(0.0, linewidth=1.0)

    plt.xlabel("Wavelength (Angstrom)")
    plt.ylabel(
        "Normalized sensitivity shape / ensemble median - 1 (%)"
    )

    plt.title(
        "Standard-star sensitivity-curve shape agreement"
    )

    plt.legend(fontsize=8)
    plt.tight_layout()

    fn = args.out_prefix + "_shape_comparison.png"
    plt.savefig(fn, dpi=180)
    plt.close()

    print("Saved:", fn)

    print()
    print("Done.")


if __name__ == "__main__":
    main()
