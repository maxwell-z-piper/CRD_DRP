
# CRD_DRP

`CRD_DRP` is the reduction and reduction-validation pipeline that prepares
KCWI/KCRM observations for the `CRD_DAP` science-analysis pipeline.

The intended boundary is:

- **CRD_DRP**: raw-data reduction, sky/telluric treatment, KcwiKit stacking,
  stack-integrity validation, atmospheric wavelength masks, and final reduction
  package validation.
- **CRD_DAP**: spatial binning, stellar kinematics, two-component decomposition,
  populations, uncertainties, and final science products.

The final KcwiKit science/variance/mask/effective-exposure cubes are never
modified by the atmospheric-mask stage. Instead, CRD_DRP creates a 1-D boolean
mask on each arm's exact native observed wavelength grid.

## Workflow

```text
Raw KCWI/KCRM data
        |
        v
KCWI DRP
        |
        +--> KCRM multi-exposure cosmic-ray masking where appropriate
        |
        v
KSkyWizard
        |
        v
KcwiKit alignment + stacking
        |
        v
validation/01_validate_kcwikit_stacks.py
        |
        v
validation/02_build_atmospheric_masks.py
        |
        v
validation/03_validate_reduction_package.py
        |
        v
CRD_DRP_reduction_manifest.json
        |
        v
CRD_DAP
```

The existing reduction instructions remain the detailed guide through KcwiKit
stacking. The scripts below begin after final BLUE and RED KcwiKit stacks exist.

## `validation/01_validate_kcwikit_stacks.py`

Checks the final KcwiKit stacks for the historical finite-negative
stack/reprojection corruption.

It distinguishes:

- `LOCALIZED_STACK_CORRUPTION_FAIL`: sparse/localized catastrophic values;
- `COMMON_MODE_SPECTRAL_WARN`: extreme features recurring at the same
  wavelength over many spatial pixels, such as a sky-line subtraction residual;
- `LOCALIZED_NEGATIVE_WARN`: unusual but non-catastrophic local structure.

A common-mode feature is not declared harmless. It is passed to the atmospheric
mask stage rather than being mislabeled as KcwiKit corruption.

Example:

```bash
python validation/01_validate_kcwikit_stacks.py \
    --blue-icube /path/blue_kskywizard_icubes.fits \
    --red-icube /path/red_kskywizard_icubes.fits \
    --output-dir validation_products/01_stack_integrity
```

## `validation/02_build_atmospheric_masks.py`

Creates fixed observed-frame BLUE and RED atmospheric masks independently of
CRD_DAP/pPXF residuals.

The production mask combines:

1. empirical sky-emission structure from one or more reduction-stage reference
   cubes;
2. an editable known-sky-line list in `data/known_sky_lines.ecsv`.

The BLUE mask therefore handles features such as [O I] 5577 using the same
general machinery as the RED mask, while allowing the actual masked intervals
to differ between arms.

Example:

```bash
python validation/02_build_atmospheric_masks.py \
    --blue-icube /path/blue_kskywizard_icubes.fits \
    --red-icube /path/red_kskywizard_icubes.fits \
    --blue-reference-cube /path/blue_reference_1_icubes.fits \
    --red-reference-cube /path/red_reference_1_icubed.fits \
    --known-line-list data/known_sky_lines.ecsv \
    --output-dir validation_products/02_atmospheric_masks
```

Repeat `--blue-reference-cube` or `--red-reference-cube` to use multiple
independent exposures. If the target occupies much of the field, supply a
2-D `--<arm>-reference-object-mask` where object pixels are >0 and sky pixels
are 0. Note that the reference cubes should be from kcwidrp, not later in 
the reduction.

Products per arm:

```text
blue_atmospheric_mask.fits
blue_atmospheric_mask.ecsv
blue_atmospheric_reference.ecsv
blue_atmospheric_mask_qc.png
```

and the corresponding RED files.

Mask semantics:

```text
0 = not excluded by atmospheric mask
1 = exclude from stellar analysis
```

## `validation/03_validate_reduction_package.py`

Checks the final handoff package:

- `i/v/m/e` stacks exist and are internally consistent;
- each atmospheric mask exactly matches its stack's native wavelength grid;
- the mask fraction is within a sanity bound;
- catastrophic negative groups inside the final atmospheric mask are accounted
  for;
- no catastrophic unexplained negative group remains outside the mask;
- upstream stack-integrity failure status is propagated.

Example:

```bash
python validation/03_validate_reduction_package.py \
    --blue-icube /path/blue_kskywizard_icubes.fits \
    --red-icube /path/red_kskywizard_icubes.fits \
    --blue-atmospheric-mask validation_products/02_atmospheric_masks/blue_atmospheric_mask.fits \
    --blue-atmospheric-intervals validation_products/02_atmospheric_masks/blue_atmospheric_mask.ecsv \
    --red-atmospheric-mask validation_products/02_atmospheric_masks/red_atmospheric_mask.fits \
    --red-atmospheric-intervals validation_products/02_atmospheric_masks/red_atmospheric_mask.ecsv \
    --stack-integrity-summary validation_products/01_stack_integrity/stack_integrity_summary.json \
    --atmospheric-mask-summary validation_products/02_atmospheric_masks/atmospheric_mask_summary.json \
    --output-dir validation_products/03_reduction_package
```

A successful run creates:

```text
validation_products/03_reduction_package/CRD_DRP_reduction_manifest.json
```

This manifest is the contractual handoff to CRD_DAP.

## Planned CRD_DAP integration

The next CRD_DAP revision will make Script 01 read
`CRD_DRP_reduction_manifest.json`, verify `validation_status == "PASS"`, load
both arms and both atmospheric masks, and combine those masks with the native
data-quality masks while preserving the components separately for provenance.

Once that is implemented, downstream CRD_DAP scripts will no longer need to
construct their own production atmospheric masks.

## Important design rules

- Do not interpolate over atmospheric wavelengths or overwrite science flux.
- Do not define production masks from pPXF residuals.
- Keep BLUE and RED symmetric at the data-contract level:
  `i/v/m/e + native atmospheric mask + validation`.
- Preserve the interval ECSV as human-readable provenance and the FITS mask as
  the exact machine-readable native-grid product.
