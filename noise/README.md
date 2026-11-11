# Stack ZAP science cubes and replay the same stack with Gaussian noise

This records the completed 8143-1902 setup, so its two named JSON files are fixed snapshots. This is mainly for my own notes, but could be helpful for troubleshooting For a new reduction, generate both recipes from the final blue and red `.list` files with `make_recipe.py` as shown in `INSTRUCTIONS_COMPLETE_NOISE_REDUCTION.md`; the number of exposures then comes from the lists rather than these examples.

The current red `_zap_icubes.fits` exposures are the intended science inputs.
KSkyWizard first ran `crop_cube` on the earlier DRP cubes, selecting the header's
`WAVGOOD0`–`WAVGOOD1` range, then ran ZAP and wrote the flux-calibrated ZAP
products. For exposure 172 the earlier cube spans 5131–9356 Å; the ZAP product
spans 5567–8929 Å. This crop is why the later product is shorter. ZAP itself
does not recover the excluded ends. The red files are labeled `RL` in their
headers; they remain development data for a planned RH3 analysis.

The recorded old red stack list used four earlier `_icubes.fits` exposures, so
it is not the correct reference for these ZAP outputs. The recorded blue list
used only `kb241226_00119_icubes.fits`. A direct blue science replay matched
its saved variance, mask, and exposure cubes exactly, but differed in roughly
14,000 flux samples. The provenance of those flux-only changes is unresolved.
A fresh blue science stack therefore becomes the reference for blue mocks. The
recipes preserve exposure selection, `3 3` spatial trim, measured shift lists,
stack parameter files, and drizzle value 0.9; they do not recompute shifts.

Use the `kcwikit` environment on this Mac. All output paths below are new;
the existing stacks and reduced exposures are left intact. First audit the
inputs:

```bash
/Users/maxpiper/miniforge3/envs/kcwikit/bin/python /Users/maxpiper/CRD_DRP/noise/stack_noise.py \
  --recipe /Users/maxpiper/CRD_DRP/noise/8143-1902_red_zap.json \
  --output-dir /Users/maxpiper/KCWI_Data/noise_calibration/red_zap_science \
  --science-only --audit-only
```

Build both reproducible science stacks:

```bash
/Users/maxpiper/miniforge3/envs/kcwikit/bin/python /Users/maxpiper/CRD_DRP/noise/stack_noise.py \
  --recipe /Users/maxpiper/CRD_DRP/noise/8143-1902_red_zap.json \
  --output-dir /Users/maxpiper/KCWI_Data/noise_calibration/red_zap_science \
  --science-only

/Users/maxpiper/miniforge3/envs/kcwikit/bin/python /Users/maxpiper/CRD_DRP/noise/stack_noise.py \
  --recipe /Users/maxpiper/CRD_DRP/noise/8143-1902_blue.json \
  --output-dir /Users/maxpiper/KCWI_Data/noise_calibration/blue_science \
  --science-only
```

Each run writes `science_stack/stack_{i,v,m,e}cubes.fits` and a
`recipe_for_noise.json` pointing to its **new** science stack. In QFitsView,
collapse clean continuum planes of `stack_icubes.fits` to check the white-light
shape and centering. Inspect a central plane of `stack_ecubes.fits` for exposure
coverage, and compare the same plane of `stack_vcubes.fits` and
`stack_mcubes.fits` for unexpected variance or masking. Browse a few central
and outer spectra for spikes and verify the FITS wavelength limits. The
alignment notebook supplies alignment plots; `validation/01` supplies
minimum-flux and corruption maps. There is no dedicated automatic white-light
or coverage QA report yet. The DAP target config points to the new blue and
red stack filenames.

Run the stack-integrity check on these exact reference stacks before the
noise ensemble. It writes minimum-flux and localized-failure maps plus a
machine-readable summary:

```bash
cd /Users/maxpiper/CRD_DRP
/Users/maxpiper/miniforge3/envs/crd_dap/bin/python validation/01_validate_kcwikit_stacks.py \
  --blue-icube /Users/maxpiper/KCWI_Data/noise_calibration/blue_science/science_stack/stack_icubes.fits \
  --red-icube /Users/maxpiper/KCWI_Data/noise_calibration/red_zap_science/science_stack/stack_icubes.fits \
  --output-dir validation_products/01_stack_integrity
```

Generate independent zero-mean Gaussian draws from the `UNCERT` extension of
each reduced exposure, keeping headers, quality flags, masks and all other
extensions fixed. Each mock is stacked by KcwiKit with the same fixed science
shifts and drizzle configuration. The tool first replays the true science
inputs and requires them to match the designated science reference; it then
checks that every mock keeps the same variance, exposure and mask products.

```bash
/Users/maxpiper/miniforge3/envs/kcwikit/bin/python /Users/maxpiper/CRD_DRP/noise/stack_noise.py \
  --recipe /Users/maxpiper/KCWI_Data/noise_calibration/red_zap_science/recipe_for_noise.json \
  --output-dir /Users/maxpiper/KCWI_Data/noise_calibration/red_zap_noise \
  --realizations 32 --seed 20261021 --workers 2

/Users/maxpiper/miniforge3/envs/kcwikit/bin/python /Users/maxpiper/CRD_DRP/noise/stack_noise.py \
  --recipe /Users/maxpiper/KCWI_Data/noise_calibration/blue_science/recipe_for_noise.json \
  --output-dir /Users/maxpiper/KCWI_Data/noise_calibration/blue_noise \
  --realizations 32 --seed 20261022 --workers 2
```

Each `noise_manifest.json` records exact source, KcwiKit, recipe, science and
mock hashes. The mock flux stacks are intentionally retained. A failed replay
is a provenance problem to investigate, not a value to tune until it passes.
If a run fails, use a fresh output directory after resolving its cause.
`--workers 2` runs two independent stacks concurrently in separate processes;
use `--workers 1` when RAM or temporary disk is tight. Four draws check the
workflow and meet Script 2's minimum; 32 were chosen as a practical way to
reduce Monte Carlo noise in each bin's variance estimate. Check free disk
space before increasing the count, since every flux cube is retained.
Script 2 takes a median of per-channel sample variances. In an idealized
Gaussian constant-variance case, four draws make that scale about 11% low;
32 draws reduce the analogous finite-draw bias to about 1%. Four is useful
for testing the workflow, but the resulting calibration is preliminary.

After changing the science stacks, regenerate the CRD_DRP atmospheric mask
and reduction-package manifest using the new red stack. The old red mask has
the wrong wavelength grid. From `/Users/maxpiper/CRD_DRP`, use the same sky
references and telluric model recorded in `pipeline_calls.txt`, but supply the
new red stack and a new output directory:

```bash
cd /Users/maxpiper/CRD_DRP
/Users/maxpiper/miniforge3/envs/crd_dap/bin/python validation/02_build_atmospheric_masks.py \
  --blue-icube /Users/maxpiper/KCWI_Data/noise_calibration/blue_science/science_stack/stack_icubes.fits \
  --red-icube /Users/maxpiper/KCWI_Data/noise_calibration/red_zap_science/science_stack/stack_icubes.fits \
  --blue-reference-cube /Users/maxpiper/KCWI_Data/redux/kb241226_00118_icubes.fits \
  --blue-reference-cube /Users/maxpiper/KCWI_Data/redux/kb241226_00119_icubes.fits \
  --red-reference-cube /Users/maxpiper/KCWI_Data/redux/kr241226_00172_icubed.fits \
  --red-reference-cube /Users/maxpiper/KCWI_Data/redux/kr241226_00173_icubed.fits \
  --red-reference-cube /Users/maxpiper/KCWI_Data/redux/kr241226_00174_icubed.fits \
  --red-reference-cube /Users/maxpiper/KCWI_Data/redux/kr241226_00175_icubed.fits \
  --red-telluric-spectrum /Users/maxpiper/KCWI_Data/redux/kskywizard/kr241226_00110_standard_spec1d_tellmodel.fits \
  --red-telluric-medium air \
  --known-line-list data/known_sky_lines.ecsv \
  --output-dir validation_products/02_atmospheric_masks_red_zap

/Users/maxpiper/miniforge3/envs/crd_dap/bin/python validation/03_validate_reduction_package.py \
  --blue-icube /Users/maxpiper/KCWI_Data/noise_calibration/blue_science/science_stack/stack_icubes.fits \
  --red-icube /Users/maxpiper/KCWI_Data/noise_calibration/red_zap_science/science_stack/stack_icubes.fits \
  --blue-atmospheric-mask validation_products/02_atmospheric_masks_red_zap/blue_atmospheric_mask.fits \
  --blue-atmospheric-intervals validation_products/02_atmospheric_masks_red_zap/blue_atmospheric_mask.ecsv \
  --red-atmospheric-mask validation_products/02_atmospheric_masks_red_zap/red_atmospheric_mask.fits \
  --red-atmospheric-intervals validation_products/02_atmospheric_masks_red_zap/red_atmospheric_mask.ecsv \
  --atmospheric-mask-summary validation_products/02_atmospheric_masks_red_zap/atmospheric_mask_summary.json \
  --output-dir validation_products/03_reduction_package_red_zap
```

Run fresh DAP stages in this order:

```bash
cd /Users/maxpiper/CRD_DAP
/Users/maxpiper/miniforge3/envs/crd_dap/bin/python scripts/01_prepare_and_register_cubes.py \
  --config config/8143-1902.py \
  --reduction-manifest /Users/maxpiper/CRD_DRP/validation_products/03_reduction_package_red_zap/CRD_DRP_reduction_manifest.json \
  --run-name 8143-1902_scalar_S01

/Users/maxpiper/miniforge3/envs/crd_dap/bin/python scripts/02_make_master_BL_bins.py \
  --config config/8143-1902.py \
  --script1-run /Users/maxpiper/CRD_DAP/runs/8143-1902_scalar_S01 \
  --run-name 8143-1902_scalar_S02

/Users/maxpiper/miniforge3/envs/crd_dap/bin/python scripts/03_build_RH3_likelihood_cubes.py \
  --config config/8143-1902.py \
  --script1-run /Users/maxpiper/CRD_DAP/runs/8143-1902_scalar_S01 \
  --script2-run /Users/maxpiper/CRD_DAP/runs/8143-1902_scalar_S02 \
  --run-name 8143-1902_scalar_S03 \
  --workers 3
```

Use new run names if these directories already exist. Do not resume
checkpoints from the earlier stacks. Script 2 reads both noise manifests from
`/Users/maxpiper/CRD_DAP/config/8143-1902.py` and calibrates every actual
PowerBin. Script 3 applies an additional residual multiplier per bin. The
PowerBin tessellation target S/N remains based on the formal BL input noise;
its reported achieved S/N uses the calibrated uncertainty.

For a later observing run, use the final red and blue KcwiKit `.list` files
to generate new, exposure-count-independent recipes with `noise/make_recipe.py`.
The complete sequence and matching directory names are in
`INSTRUCTIONS_COMPLETE_NOISE_REDUCTION.md`. A list root ending in `_zap`
selects the corresponding `_zap_icubes.fits`; leave it off for final DRP
`_icubes.fits` products. The current list files must be updated before using
them for recipe generation: they still point to older non-ZAP or external
locations.

These realizations cover independent uncertainty at the already reduced
exposure level propagated through the fixed stack. They do not simulate the
uncertainty of the detector reduction, fitted ZAP sky model, telluric correction,
or alignment fit. A diagonal residual scalar also does not establish coverage
of parameter confidence intervals; injection and recovery checks are still
needed before interpreting widths as final.
