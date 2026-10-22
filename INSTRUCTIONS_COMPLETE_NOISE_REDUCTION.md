# Complete reduction: raw KCWI data through Gaussian stack-noise calibration

The detector-reduction steps below retain the short format of `INSTRUCTIONS.md`. The continuation after them connects the final stack lists to the science and noise replays. Run each reduction in its normal data directory. A repeat of the stack/noise stage needs fresh output directories or deliberate removal of the previous run after preserving anything needed. Commands with `00XXX`, `<night>`, or other placeholders still need your actual frame names.

Check that the installed KSkyWizard and KcwiKit patches described in `patches/README.md` are present before starting. Put the final exposure roots in each stack `.list`: use `_zap` roots for any arm processed by KSkyWizard. The continuation explains how those lists become the science/noise recipes.

---

# CRD_DAP_Reduction
Scripts to reduce raw KCWI data into the needed data products for CRD_DAP. There are a few patches to KSkyWizard and KcwiKit for this reduction. There is also a validation script to check the quality of the reduction.

# Conda environments
We will need three conda environments to run [kcwidrp](https://github.com/Keck-DataReductionPipelines/KCWI_DRP), [kskywizard](https://github.com/zhuyunz/KSkyWizard/tree/main), and [kcwikit](https://github.com/yuguangchen1/KcwiKit). For information on how to initialize these and clone their repos, see the above links.


# ===== RED ARM REDUCTION =====

In kcwidrp conda env, cd to where data are stored. Make sure kcwi.cfg is located in this directory. Run:	
```bash
reduce_kcwi -r -f kr*.fits -g -c kcwi.cfg -st ff
```

Ensure everything ran successfully and geometry .fits files were created. Run:	
```bash
kcwi_gen_log <night>.log -f kr*.fits
```

Make sure the <night>.log is in the intended directory. Run:	
```bash
kcrm_group_frames <night>.log kcrm_<night>.json -c -d -i
```

Open the created .json file and remove the standard star observations. Also ensure that the remaining science exposure have multiple .fits for each exposure length - having only a single exposure will not work. Run:	
```bash
kcrm_create_crmsk kcrm_<night>.json
```

Go back to the kcwir.proc table and remove all of the science frames that had a crmsk made for it. Then rerun the DRP to the end of cube creation to recover masked pixels. Run:	
```bash
reduce_kcwi -r -f kr241226_00XX[X-X].fits -g -c kcwi.cfg -st 2 -k
```

Go back to the kcwir.proc and delete the kr241226_00XX[X-X].fits entries. Also delete the standard star entries. Now, we have to do the standard star WITHOUT the -k to make the invsens file. We will run stage 2 and 3 of the pipeline in a row. Run:	
```bash
reduce_kcwi -r -f kr241226_00XX[X-X].fits -g -c kcwi.cfg -st 2
```

Open kcwir.proc and delete the kr241226_00XX[X-X].fits standard star entries + sky entries. Run:	
```bash
reduce_kcwi -r -f kr241226_00XX[X-X].fits -g -c kcwi.cfg -st 3
```

This creates multiple invsens files, one for each standard star exposure. We only need a single one for KSkyWizard. First, ensure 'compare_standard_invsens.py' is in the directory with all of your data. Then, we will check and compare each one before we choose which to use. Run:	
```bash
python compare_standard_invsens.py redux/kr241226_00XXX_invsens.fits redux/kr241226_00XXX_invsens.fits redux/kr241226_00XXX_invsens.fits
```

Note the one with the lowest QC score, but also check the created plots to confirm. If everything looks decent, choose the one with the lowest QC score. Keep a note of which file that is.

We will SKIP the collapse, flatten, and med filter stage as we are going to let KSkyWizard do the sky subtraction. Then proceed to stage 3 for the science images with the -k keyword. Run:	
```bash
reduce_kcwi -r -f kr241226_00XX[X-X].fits -g -c kcwi.cfg -st 3 -k
```

At this point, we will take the icubed files and run through KSkyWizard to do the sky-subtraction and telluric-correction. To start, ensure the processing.cfg file (located at /Users/maxpiper/miniforge3/envs/kskywizard/lib/python3.11/site-packages/kskywizard) has the correct telgridfile directory (this should be /Users/maxpiper/.pypeit/cache/download/url/5f17ecc1fcc921d6ec01e18d931ec2f8/contents after our edits). We then need the kskywizard env to be activated and then open the kskywizard gui. Run:	
```bash
conda deactivate
conda activate kskywizard
kskywizard
```

From this gui, we need to set the input directory to be where the icubed files and the invsens file live. Set the output directory somewhere you will remember, as this is where we will take the corrected cubes and stack them with KcwiKit. Then, open the DRP invsens file that we selected earlier. We want to exclude any obvious emission or absorption lines that are very sharply peaked and that typically correspond to known elemental wavelengths (like Ha 6564). Hold 'e' and drag cursor over those points. When finished, hit 'f' to re-fit the sensitivity.

Once happy with the cyan line's fit, hit 't' to generate the telluric correction model. This will take a bit and the gui will freeze, so wait until the output updates. Inspect the royal blue line, especially over the regions that have the dark green x's selected for telluric fitting. It's normal for it to follow the sharp spikes that we masked earlier. Once happy, press 'Save updated invsens'.

Now onto sky subtraction. Switch to the 'Science' tab. Enter in the first science frame number, hit enter with the cursor still in that box, and if it's the first time analyzing it, hit 'Load DRP Cube', 'Save Cropped Cube', then 'Load Cropped Cube'. Enter in the target redshift so we can see the spectral lines overlayed. If you are using an off-field sky frame, enter in that frame number. 

If using an in-field sky frame, we need to make it in the terminal. First, open a science exposure in QFitsView to determine the ellipse x and y radii that encompass the galaxy. Then also determine the center pixel of the galaxy. Finally, run this command, making sure to cd into the output directory: 

```bash
cd /Path/to/KSkyWizard/output/directory

cat > kr241226_00172.reg <<'EOF'
# Region file format: DS9 version 4.1
image
ellipse(center_x,center_y,radius_x,radius_y,angle)
EOF
```

Note that, if aligned along PA_kin, angle=0. Then, enter in the ZAP frame number you created, then hit 'Update ZAP mask'. Note that if you did many exposures without dithering, it is valid to use the same ZAP mask for all exposures. Finally, hit 'Process', then type 'q'. This analysis takes a bit and the gui freezes too, so wait until the output updates.

If you are happy with the results, press 'Next-->' and repeat the above sky-subtraction steps for all science frames.


Now onto stacking the sky-subtracted and CR-corrected science frames with KcwiKit. Run:
```bash
cd /Users/maxpiper/KcwiKit
conda activate kcwikit
```

We need to create 2 configuration files that live in the directory of our KSkyWizard datafiles. The first is the .list, which contains a list of the directory to the cubes we want to stack, followed by how many pixels from the bottom and top of the cube we want to pad (can usually just say 3 3). Run:

```bash
cat > /Users/maxpiper/KCWI_Data/redux/kskywizard/red_kskywizard.list <<'EOF'
/Users/maxpiper/KCWI_Data/redux/kskywizard/kr<night>_00XXX_zap 3 3
/Users/maxpiper/KCWI_Data/redux/kskywizard/kr<night>_00XXX_zap 3 3
/Users/maxpiper/KCWI_Data/redux/kskywizard/kr<night>_00XXX_zap 3 3
EOF
```

The second is the .par file which contains the x and y dimensions of the final cube in pixel space (dimension x y) (for large, 67 112; medium, 67 56; small, 134 56) and the pixel scales in arcsec/pix (xpix, ypix) (for large and medium slicer: 0.3; for small slicer: 0.15). Note that we will keep the align_box line empty until later. Run:

```bash
cat > /Users/maxpiper/KCWI_Data/redux/kskywizard/red_kskywizard.par <<'EOF'
dimension x y
align_box x_lower_left y_lower_left x_upper_right y_upper_right
orientation PA
xpix x
ypix y
EOF
```

It it easiest to open a cube in QFitsView to determine how big the cube needs to be as well as what pixel coordinates the alignment box needs. Typically, you can set the orientation to be the set sky PA so it maps cleanly onto MaNGA maps.

We will run this in a Jupyter notebook. Run:	
```bash
jupyter notebook
```

Open the 'kcwikit_alignment_script.ipynb' file and change the beginning directory to where we are keeping all the data, and change 'fn' to be the name of the .list file just created. First run the 'kcwi.kcwi_align(fn, noalign=True)' command, which will create the thum0.fits file. It's important to open this file in QFitsView and enter the pixel coordinates here for the alignment box in the .par file we created earlier, as this is the end-product cube and these coordinates will cross-over correctly. Save that file, then run 'kcwi.kcwi_align(fn)'. Inspect the .pngs and ensure that the cross is in the center of all the images. Then, run the final 'kcwi.kcwi_stack(fn, drizzle=0.9, overwrite=True)' stack line to make your stacked cube!



# ===== Blue ARM REDUCTION =====

In the kcwidrp conda env, cd into where data are stored. Make sure the kcwi.cfg file is in this directory. Run:	
```bash
reduce_kcwi -b -f kb*.fits -g -c kcwi.cfg -st 1
```

We will then make sure the galaxy continuum is not in the sky subtraction. Run:

```bash
cd  /Users/maxpiper/KcwiKit
jupyter notebook
```
Then, open 'create_region_masks.ipynb'. Change the directory to it points to the kcwidrp redux directory, and ensure the frame roots contain all of the blue frames we want. Run the cell, then inspect the output. We want the left plot's ellipse to comfortable encompass the entire galaxy and its faint surrounding region. The right will show the intf detector representation with the actual pixels that will be excluded from the sky calc. Also check the output 'Value IFU area masked' and make sure this is sensible and not ~0% or ~100%. If you are happy with the result, run the next cell, which renames those files to the real _smsk.fits binary mask files. Then, run the final cell to create the kcwi.sky file.

Finally, remove all OBJECT and SKY entries in the kcwib.proc file, and remove the '_sky.fits' and '_intk.fits' files from the /redux/ directory for the science frames. cd back to where the data are stored. Run:	

```bash
reduce_kcwi -b -f kb*.fits -g -c kcwi.cfg -st 2
```

Remove all the science entries from kcwib.proc (including their sky's), then we will run stage three of the pipeline. Run:	

```bash
reduce_kcwi -b -f kb*.fits -g -c kcwi.cfg -st 3
```

We will then combine the standard star frames to improve the flux-calibration accuracy. Create a .list file with the full path of the '_invsens.fits' files to be combined. Run:

```bash
cat > blue.list <<'EOF'
/Path/to/reduced/kbYYMMDD_XXXXX_invsens.fits
/Path/to/reduced/kbYYMMDD_XXXXX_invsens.fits
/Path/to/reduced/kbYYMMDD_XXXXX_invsens.fits
EOF
```

Then, we will combine them. Run:	

```bash
kcwi_combinestd blue.list
```

This will create a new 'blue_invsens.fits' file. Make sure to put this into the /redux/ folder. This will also create a bokeh figure comparing the individual curves of the standards. Inspect this - if there are outliers, remove those entries from the .list file and rerun the above command. Then, open kcwib.proc and remove all of the standard star '_invsens.fits' lines except for one. In this line, replace both kbYYMMDD_XXXXX.fits entries with 'blue.fits'.

Finally, remove all the science frame entries from kcwib.proc and rerun the last stage of the pipeline. Run:	

```bash
reduce_kcwi -b -f kb*.fits -g -c kcwi.cfg -st 3
```

For an observing run where blue frames also need KSkyWizard/ZAP, process each chosen blue cube with the KSkyWizard GUI steps described for red above before making the blue stack list. Use the final blue `_zap_icubes.fits` outputs; the current single-exposure blue run instead uses its final DRP `_icubes.fits`.

Now onto stacking the blue science frames with KcwiKit. In the terminal, activate the kcwikit conda env. 

We need to create 2 configuration files that live in the directory of our KSkyWizard datafiles. The first is the .list, which contains a list of the directory to the cubes we want to stack, followed by how many pixels from the bottom and top of the cube we want to pad (can usually just say 3 3). Run:

```bash
cat > /Users/maxpiper/KCWI_Data/redux/kskywizard/blue_kskywizard.list <<'EOF'
/Users/maxpiper/KCWI_Data/redux/kb<night>_00XXX 3 3
/Users/maxpiper/KCWI_Data/redux/kb<night>_00XXX 3 3
/Users/maxpiper/KCWI_Data/redux/kb<night>_00XXX 3 3
EOF
```

Use one row per final blue exposure. The present blue reduction uses DRP `_icubes.fits` products, so its roots above omit `_zap`. If a later run takes multiple blue exposures through KSkyWizard/ZAP, put those final products in the list instead, with one `_zap` root per exposure. For example, a list root ending in `kb<night>_00XXX_zap` selects `kb<night>_00XXX_zap_icubes.fits`.

The second is the .par file which contains the x and y dimensions of the final cube in pixel space (dimension x y) (for large, 67 112; medium, 67 56; small, 134 56) and the pixel scales in arcsec/pix (xpix, ypix) (for large and medium slicer: 0.3; for small slicer: 0.15). Note that we will keep the align_box line empty until later. Run:

```bash
cat > /Users/maxpiper/KCWI_Data/redux/kskywizard/blue_kskywizard.par <<'EOF'
dimension x y
align_box x_lower_left y_lower_left x_upper_right y_upper_right
orientation PA
xpix x
ypix y
EOF
```

**IT IS IMPORTANT TO USE THE SAME .par ARGUMENTS FOR THE BLUE SIDE AS YOU DID FOR THE RED** The mapping from the blue bins to the red bins relies on both cubes having the same spatial gridding. The easiest way to do this is to open the red.par file, copy all the arguments into the above, then run the terminal command to create the blue.par file.

We will run this in a Jupyter notebook. Run:	

```bash
cd /Users/maxpiper/KcwiKit
jupyter notebook
```

Open the 'kcwikit_alignment_script.ipynb' file and change the beginning directory to where we are keeping all the data, and change 'fn' to be the name of the .list file just created. First run the 'kcwi.kcwi_align(fn, noalign=True)' command, which will create the thum0.fits file. It's important to open this file in QFitsView and enter the pixel coordinates here for the alignment box in the .par file we created earlier, as this is the end-product cube and these coordinates will cross-over correctly. Save that file, then run 'kcwi.kcwi_align(fn)'. Inspect the .pngs and ensure that the cross is in the center of all the images. Then, run the final 'kcwi.kcwi_stack(fn, drizzle=0.9, overwrite=True)' stack line to make your stacked cube!

---

# Fixed science stacks, Gaussian noise, and DAP handoff

The `.list` for each arm is the authority for **which final exposures** are stacked and for their two trim values. The list root must point to a real final `_icubes.fits` exposure; include `_zap` in the root if KSkyWizard/ZAP made the final cube. Add or remove rows as your observing run requires. Keep each arm's `.par` and saved `.shift.list` beside its `.list`, with the same basename. The shift rows must correspond to the listed exposures in the same order.

The interactive alignment steps above are still needed to set the `.par` and measure shifts. If disk space is tight, you can omit their final manual `kcwi.kcwi_stack(...)` call: the `--science-only` command below makes the final science stack used by the DAP. KcwiKit alignment plots and `thum0.fits` remain useful visual checks. These commands do not rerun detector reduction or ZAP.

For the current Mac reduction, both lists live under `/Users/maxpiper/KCWI_Data/redux/kskywizard`. The older list files there still name pre-ZAP or other-volume roots; before using this guide, replace those rows with the actual final inputs. A different output location is fine if the variables below match your reduction.

```bash
export KCWI_DATA=/Users/maxpiper/KCWI_Data
export DRP_REDUX="$KCWI_DATA/redux"
export KSKYWIZARD_OUT="$DRP_REDUX/kskywizard"
export NOISE_ROOT="$KCWI_DATA/noise_calibration"
mkdir -p "$NOISE_ROOT"
```

## 1. Generate the stack recipes from the two lists

`make_recipe.py` reads every non-comment row, resolves its final FITS cube, checks the matching `.par` and `.shift.list`, and writes a recipe with the exact exposure count, order, trims, and fixed drizzle settings. It stops on missing files or mismatched shift rows. The names are symmetric for blue and red; whether a particular exposure was ZAP-processed is encoded in that exposure's list root. Run after both lists and alignments are final:

```bash
/Users/maxpiper/miniforge3/envs/kcwikit/bin/python /Users/maxpiper/CRD_DRP/noise/make_recipe.py \
  --list "$KSKYWIZARD_OUT/red_kskywizard.list" --arm RH3 \
  --output "$NOISE_ROOT/red_recipe.json"
/Users/maxpiper/miniforge3/envs/kcwikit/bin/python /Users/maxpiper/CRD_DRP/noise/make_recipe.py \
  --list "$KSKYWIZARD_OUT/blue_kskywizard.list" --arm BL \
  --output "$NOISE_ROOT/blue_recipe.json"
```

If your manual stack used settings other than drizzle 0.9, sigma clip 10, or `npix_trim=3`, pass the actual `--drizzle`, `--nsigma-clip`, and `--npix-trim` values to **both** recipe commands. The tool will not overwrite an existing recipe; use a new reduction directory for a new run.

Audit the inputs without creating a stack:

```bash
/Users/maxpiper/miniforge3/envs/kcwikit/bin/python /Users/maxpiper/CRD_DRP/noise/stack_noise.py \
  --recipe "$NOISE_ROOT/red_recipe.json" --output-dir "$NOISE_ROOT/red_science" \
  --science-only --audit-only
/Users/maxpiper/miniforge3/envs/kcwikit/bin/python /Users/maxpiper/CRD_DRP/noise/stack_noise.py \
  --recipe "$NOISE_ROOT/blue_recipe.json" --output-dir "$NOISE_ROOT/blue_science" \
  --science-only --audit-only
```

## 2. Make and inspect the science stacks

```bash
/Users/maxpiper/miniforge3/envs/kcwikit/bin/python /Users/maxpiper/CRD_DRP/noise/stack_noise.py \
  --recipe "$NOISE_ROOT/red_recipe.json" --output-dir "$NOISE_ROOT/red_science" --science-only
/Users/maxpiper/miniforge3/envs/kcwikit/bin/python /Users/maxpiper/CRD_DRP/noise/stack_noise.py \
  --recipe "$NOISE_ROOT/blue_recipe.json" --output-dir "$NOISE_ROOT/blue_science" --science-only
```

Each directory contains `science_stack/stack_{i,v,m,e}cubes.fits` and a `recipe_for_noise.json` tied to that **new** science stack. Do not move the science stack after this point: the noise recipe stores its absolute path.

For visual QA, open each `stack_icubes.fits` in QFitsView and collapse a clean continuum range into a white-light image. Confirm the galaxy is centered, no exposure seam dominates, and the spatial footprint agrees with the alignment plots. Open `stack_ecubes.fits` on a central wavelength plane: exposure time should cover the galaxy and decline only where expected near edges. Compare a few wavelength planes in `stack_vcubes.fits` and `stack_mcubes.fits`; variance must be finite and nonnegative wherever exposure is positive, and masks should not unexpectedly remove the galaxy. Browse spectra in central and outer spaxels for abrupt spikes, then read `CRVAL3`, `CD3_3`, `NAXIS3`, `WAVGOOD0`, and `WAVGOOD1` in the FITS headers to confirm the usable wavelength range. The earlier alignment notebook makes alignment plots, and validation 01 below makes minimum-flux and corruption maps. There is currently no separate one-command white-light/coverage QA report; that part is a visual QFitsView check.

Run the stack-integrity check before spending time on noise realizations:

```bash
cd /Users/maxpiper/CRD_DRP
/Users/maxpiper/miniforge3/envs/crd_dap/bin/python validation/01_validate_kcwikit_stacks.py \
  --blue-icube "$NOISE_ROOT/blue_science/science_stack/stack_icubes.fits" \
  --red-icube "$NOISE_ROOT/red_science/science_stack/stack_icubes.fits" \
  --output-dir validation_products/01_stack_integrity
```

Read `stack_integrity_summary.json` and inspect its PNG maps. Investigate any `FAIL`; some narrow sky residuals can trigger the conservative pre-mask scan. Check the candidate wavelengths against the atmospheric mask and review their spatial/spectral behavior before deciding whether the stack is scientifically usable. A `WARN` also needs review of the listed spectral features.

## 3. Generate Gaussian noise-stack realizations

Each draw replaces exposure flux with independent zero-mean Gaussian samples from `UNCERT`, retaining the other FITS extensions. The saved shifts and stack configuration are replayed. The tool first reproduces the science stack, then checks that each draw leaves variance, mask, and exposure products unchanged. Run red and blue as separate commands:

```bash
/Users/maxpiper/miniforge3/envs/kcwikit/bin/python /Users/maxpiper/CRD_DRP/noise/stack_noise.py \
  --recipe "$NOISE_ROOT/red_science/recipe_for_noise.json" \
  --output-dir "$NOISE_ROOT/red_noise" --realizations 32 --seed 20261021 --workers 2
/Users/maxpiper/miniforge3/envs/kcwikit/bin/python /Users/maxpiper/CRD_DRP/noise/stack_noise.py \
  --recipe "$NOISE_ROOT/blue_science/recipe_for_noise.json" \
  --output-dir "$NOISE_ROOT/blue_noise" --realizations 32 --seed 20261022 --workers 2
```

`--workers` runs independent realizations in isolated processes. Start with 2: more workers need more RAM and temporary disk, and may not scale linearly `--workers 1` is the low-memory behavior. Before a 32-draw run, check free space; each retained mock flux cube is roughly the size of its science flux cube. Four realizations (`--realizations 4`) are a quick workflow test and meet Script 2's minimum, but 32 give a less noisy estimate of the per-bin variance. They are a practical target, not a scientifically established magic number. With only four Gaussian draws, the median sample variance used by Script 2 can bias the resulting scale low (about 11% in an idealized constant-variance case); at 32 draws the analogous bias is about 1%. Check scale stability across seeds/counts when interpreting results. A completed `noise_manifest.json` records the exact inputs and output hashes. If a stack fails, inspect the incomplete output and rerun into a new directory.

## 4. Make atmospheric masks and validate the reduction package

Use the science stacks, never the mock flux stacks. The sky-reference cubes andr ed telluric model are chosen separately from the science `.list`: for this 8143-1902 run the references below are the recorded selections. Change them for a new target or observing run; a blue stack can include multiple listed exposures. The mask output grid follows the new science stacks.

```bash
cd /Users/maxpiper/CRD_DRP
/Users/maxpiper/miniforge3/envs/crd_dap/bin/python validation/02_build_atmospheric_masks.py \
  --blue-icube "$NOISE_ROOT/blue_science/science_stack/stack_icubes.fits" \
  --red-icube "$NOISE_ROOT/red_science/science_stack/stack_icubes.fits" \
  --blue-reference-cube "$DRP_REDUX/kb241226_00118_icubes.fits" \
  --blue-reference-cube "$DRP_REDUX/kb241226_00119_icubes.fits" \
  --red-reference-cube "$DRP_REDUX/kr241226_00172_icubed.fits" \
  --red-reference-cube "$DRP_REDUX/kr241226_00173_icubed.fits" \
  --red-reference-cube "$DRP_REDUX/kr241226_00174_icubed.fits" \
  --red-reference-cube "$DRP_REDUX/kr241226_00175_icubed.fits" \
  --red-telluric-spectrum "$KSKYWIZARD_OUT/kr241226_00110_standard_spec1d_tellmodel.fits" \
  --red-telluric-medium air \
  --known-line-list data/known_sky_lines.ecsv \
  --output-dir validation_products/02_atmospheric_masks

/Users/maxpiper/miniforge3/envs/crd_dap/bin/python validation/03_validate_reduction_package.py \
  --blue-icube "$NOISE_ROOT/blue_science/science_stack/stack_icubes.fits" \
  --red-icube "$NOISE_ROOT/red_science/science_stack/stack_icubes.fits" \
  --blue-atmospheric-mask validation_products/02_atmospheric_masks/blue_atmospheric_mask.fits \
  --blue-atmospheric-intervals validation_products/02_atmospheric_masks/blue_atmospheric_mask.ecsv \
  --red-atmospheric-mask validation_products/02_atmospheric_masks/red_atmospheric_mask.fits \
  --red-atmospheric-intervals validation_products/02_atmospheric_masks/red_atmospheric_mask.ecsv \
  --atmospheric-mask-summary validation_products/02_atmospheric_masks/atmospheric_mask_summary.json \
  --output-dir validation_products/03_reduction_package
```

Inspect `reduction_package_validation.json`. This post-mask check tests whether extreme negative features remain outside the atmospheric mask; it does not replace your review of the separate pre-mask integrity summary. Only proceed to DAP with a passing reduction package and a resolved explanation for any pre-mask failure. The handoff file is
`/Users/maxpiper/CRD_DRP/validation_products/03_reduction_package/CRD_DRP_reduction_manifest.json`.

## 5. Point the DAP at the science and noise products

In `/Users/maxpiper/CRD_DAP/config/<target>.py`, set the four `BL_*CUBE` paths to `blue_science/science_stack/stack_{i,v,m,e}cubes.fits`, the four `RH3_*CUBE` paths to the red equivalents, and the two `*_NOISE_STACK_MANIFEST` paths to `blue_noise/noise_manifest.json` and `red_noise/noise_manifest.json`. Update the target's master arcs and ancillary inputs when those change. Use the new reduction-package manifest in Script 1; then run Scripts 2 and 3 with new run names. Script 2 uses the complete noise ensembles to estimate a scalar for every actual bin. Do not resume DAP checkpoints from a different stack.

For the current target, after editing `config/8143-1902.py` to point to the
new products, the handoff commands are:

```bash
conda deactivate
conda activate crd_dap
cd /Users/maxpiper/CRD_DAP
/Users/maxpiper/miniforge3/envs/crd_dap/bin/python scripts/01_prepare_and_register_cubes.py \
  --config config/8143-1902.py \
  --reduction-manifest /Users/maxpiper/CRD_DRP/validation_products/03_reduction_package_red_zap/CRD_DRP_reduction_manifest.json \
  --run-name 8143-1902_scalar_redzap_S01

/Users/maxpiper/miniforge3/envs/crd_dap/bin/python scripts/02_make_master_BL_bins.py \
  --config config/8143-1902.py \
  --script1-run /Users/maxpiper/CRD_DAP/runs/8143-1902_scalar_redzap_S01 \
  --run-name 8143-1902_scalar_redzap_S02
  
/Users/maxpiper/miniforge3/envs/crd_dap/bin/python scripts/03_build_RH3_likelihood_cubes.py \
  --config config/8143-1902.py \
  --script1-run /Users/maxpiper/CRD_DAP/runs/8143-1902_scalar_redzap_S01 \
  --script2-run /Users/maxpiper/CRD_DAP/runs/8143-1902_scalar_redzap_S02 \
  --run-name 8143-1902_scalar_redzap_S03 --workers 3
```

Use new run names if these directories already exist, and substitute the new target's config and run names for another observing run.

These replays measure reduced-exposure `UNCERT` propagated through a fixed KcwiKit stack. They do not simulate uncertainty in detector reduction, the fitted ZAP sky model, telluric correction, or alignment fit. The current red RL development cubes do not cover the full planned Ca II triplet range.
