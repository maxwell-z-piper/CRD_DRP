# CRD_DAP_Reduction
Scripts to reduce raw KCWI data into the needed data products for CRD_DAP. This also contains validation scripts to check the quality of the reduction.

# Conda environments
We will need three conda environments to run [kcwidrp](https://github.com/Keck-DataReductionPipelines/KCWI_DRP), [kskywizard](https://github.com/zhuyunz/KSkyWizard/tree/main), and [kcwikit](https://github.com/yuguangchen1/KcwiKit). For information on how to initialize these and clone their repos, see the above links.

=============================
===== RED ARM REDUCTION =====
=============================

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
cat > red_kskywizard.list <<'EOF'
/Users/maxpiper/Desktop/CRD_Decomposition/KCWI_TestData/kr<night>_00XXX 3 3
/Users/maxpiper/Desktop/CRD_Decomposition/KCWI_TestData/kr<night>_00XXX 3 3
/Users/maxpiper/Desktop/CRD_Decomposition/KCWI_TestData/kr<night>_00XXX 3 3
EOF
```

The second is the .par file which contains the x and y dimensions of the final cube in pixel space (dimension x y) (for large, 67 112; medium, 67 56; small, 134 56) and the pixel scales in arcsec/pix (xpix, ypix) (for large and medium slicer: 0.3; for small slicer: 0.15). Note that we will keep the align_box line empty until later. Run:

```bash
cat > red_kskywizard.par <<'EOF'
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



==============================
===== Blue ARM REDUCTION =====
==============================

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
/Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kbYYMMDD_XXXXX_invsens.fits
/Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kbYYMMDD_XXXXX_invsens.fits
/Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kbYYMMDD_XXXXX_invsens.fits
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


Now onto stacking the blue science frames with KcwiKit. In the terminal, activate the kcwikit conda env. 

We need to create 2 configuration files that live in the directory of our KSkyWizard datafiles. The first is the .list, which contains a list of the directory to the cubes we want to stack, followed by how many pixels from the bottom and top of the cube we want to pad (can usually just say 3 3). Run:

```bash
cat > blue_kskywizard.list <<'EOF'
/Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kb<night>_00XXX 3 3
/Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kb<night>_00XXX 3 3
/Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kb<night>_00XXX 3 3
EOF
```

The second is the .par file which contains the x and y dimensions of the final cube in pixel space (dimension x y) (for large, 67 112; medium, 67 56; small, 134 56) and the pixel scales in arcsec/pix (xpix, ypix) (for large and medium slicer: 0.3; for small slicer: 0.15). Note that we will keep the align_box line empty until later. Run:

```bash
cat > blue_kskywizard.par <<'EOF'
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




