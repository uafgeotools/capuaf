#!/bin/bash
#==================================================================================
# Illinois example
#
# $eid = 20080418093700
# $model = cus (central US)
# $CAPHOME - Directory where cap.c is located
# $CAPRUN/inv/cus - Directory from where you will be running cap
# $CAPHOME/EXAMPLES/20080418093700 - data directory (this will be copied (rsync) to $CAPRUN/inv/cus
# $CAPHOME/EXAMPLES/20080418093700_check - Results (plots and output files) provided for comparison
# $ExampleDir - Directory where you want to generate new plots and output files
#
# Runs and save in the $ExampleDir
# MIGHT want to rename the previous $ExampleDir and save it (just in case)
# Search $ExampleDir and make sure it does not overwrite anything
# You will be comparing the results between $ExampleDir and $CAPHOME/EXAMPLES/20080418093700_check
#==================================================================================
eid=20080418093700
model=cus 

#!!!!!!!! MAKE SURE!!!!!!!!!!!!!
# > cd $CAPHOME
# > git checkout cap.c
# > open cap.c and set: 
#    FTC_data=1
#    FTC_green=0 
# then recompile 
# > make all
# FTC means Filter-then-Cut (the data and the green)

#==================================================================================
# --[EXAMPLE 0]---[PREPARE GENERATE GREEN'S FUNCTIONS and DATA]---------------
# This is same as the example in ../README_EXAMPLE
# modify cap.pl to indicate that you want to compute your own Green's functions
# Uncomment the following line in cap.pl
# $green = "$caprun/models";                # user testing

# 1. run the frequency-wavenumber code fk
# $CAPRUN is the directory where you will be running CAP
cd $CAPRUN/models/cus
cp $SUTIL/grp-utils/fk3.0/cus .
# Generate green's functions for cus model with double-couple source at 15 km depth and recievers at 140 145 205 230 260 275 295 410km. 
# If you are performing depth test you will need to 
# For more info
# > fk.pl
echo "fk.pl -Mcus/15 -N512/0.4/2 140 145 205 230 260 275 295 410"
fk.pl -Mcus/15 -N512/0.4/2 140 145 205 230 260 275 295 410
# If you are performing depth test you will need to uncomment the following lines
# fk.pl -Mcus/5 -N512/0.4/2 140 145 205 230 260 275 295 410
# fk.pl -Mcus/10 -N512/0.4/2 140 145 205 230 260 275 295 410
# fk.pl -Mcus/20 -N512/0.4/2 140 145 205 230 260 275 295 410
# fk.pl -Mcus/25 -N512/0.4/2 140 145 205 230 260 275 295 410
# fk.pl -Mcus/30 -N512/0.4/2 140 145 205 230 260 275 295 410

# note: this uses the "cus" 1D model for "central U.S."
# note: You do not need to run fk if you have a pre-computed library of Green's functions.
# Just Uncomment the following line in cap.pl
# $green = "/store/wf/FK_synthetics";        # standard models at UAF

# 2. run fk again to generate isotropic components
fk.pl -Mcus/15 -N512/0.4/2 -S0 140 145 205 230 260 275 295 410
# If you are performing depth test you will need to uncomment the following lines
# fk.pl -Mcus/5 -N512/0.4/2 -S0 140 145 205 230 260 275 295 410
# fk.pl -Mcus/10 -N512/0.4/2 -S0 140 145 205 230 260 275 295 410
# fk.pl -Mcus/20 -N512/0.4/2 -S0 140 145 205 230 260 275 295 410
# fk.pl -Mcus/25 -N512/0.4/2 -S0 140 145 205 230 260 275 295 410
# fk.pl -Mcus/30 -N512/0.4/2 -S0 140 145 205 230 260 275 295 410

# This will create 3 additional green functions for each distance of the form dist.grn.(a-c),
# e.g., 140.grn.a, 140.grn.b, 140.grn.c.

# 3. prepare data directory and directory for storing the results ($ExampleDir). Copy the the data to directory where you will be running cap
# All commands will be run form here (Alternatively, one could run from $CAPHOME/EXAMPLES/ without copying the data)
# (CAUTION: ExampleDir will be overwritten)
ExampleDir=$CAPRUN/inv/cus/20080418093700_check2
mkdir -p $ExampleDir

cd $CAPRUN/inv/cus
rsync -av $CAPHOME/EXAMPLES/20080418093700 .

#Check that you have new weight.dat file with 12 columns (3 extra for Pnl_window, surf_pick, surf_window).
#The old weight file is under the name weight_orig.dat
more 20080418093700/weight.dat

# 4. Run cap 
echo "cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15/5.0 -E0 -K0 -Y1 20080418093700"
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15/5.0 -E0 -K0 -Y2 20080418093700
# save results
cp 20080418093700/cus_015.ps $ExampleDir/cus_015.ps
cp 20080418093700/cus_015.out $ExampleDir/cus_015.out

# 5. Make sure both the results are same
# gv $CAPHOME/EXAMPLES/20080418093700_check/cus_015.ps &
# gv $ExampleDir/cus_015.ps
echo "diff $ExampleDir/cus_015.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015.out"
diff $ExampleDir/cus_015.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015.out

#====================================================================================
# From here on out in README_EXAMPLE, we assume that the Green's functions are pre-computed.
# To run these examples green's functions needs to be pre-computed (see README_EXAMPLE for how to generate green's function)
# Make this change in cap.pl: (uncomment #standard models at UAF" and comment the "#user testing")
# $green = "/store/wf/FK_synthetics";               # UAF linux network
# # $green = "$caprun/models";                      # user testing (generate your own green's function)
# If you running CAP on cluster (pacman) uncomment the following:
# $green = "/import/c/d/ERTHQUAK/FK_synthetics ";   # UAF cluster

###################################################################################################

# -------------[EXAMPLES BEGIN]----------------------------------
## EXAMPLE 1: Illinois event, double couple (DC)

# perform the depth test
# 1. First make sure there are not other *.out files
rm $eid/cus*.ps
rm $eid/cus*.out

# 2. Run cap for different depths by setting the range and increment for which you want to run depth test using -A flag (-Adep_min/dep_max/dep_inc).
# One output (.out) file should be generated for each depth.
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15/5.0 -E0 -K0 -Y2 -A5/30/5 20080418093700

# 3. perform depth test
depth_test 20080418093700 cus 

# 4. save the output files for comparison
mv $eid/"$model"* $ExampleDir/
mv dep_20080418093700.ps $ExampleDir/depth_dc.ps
mv junk1.out $ExampleDir/junk1_dc.out
mv junk2.out $ExampleDir/junk2_dc.out

# 5. Check results
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $ExampleDir/depth_dc.ps $CAPHOME/EXAMPLES/20080418093700_check/depth_dc.ps
echo "$ExampleDir/junk1_dc.out $CAPHOME/EXAMPLES/20080418093700_check/junk1_dc.out"
diff $ExampleDir/junk1_dc.out $CAPHOME/EXAMPLES/20080418093700_check/junk1_dc.out 

#==================================================================================
## EXAMPLE 2: Illinois event, DC with P wave polarities in weight file
# 0. Follow same step 1-3 EXAMPLE 0

# NOTE: skip_zero_weights=1 (default) will not plot the polarity information but it will be used
# In cap.c skip_zero_weights=0 and compile
#-----------------------------------------------------------------------------------
# NOTE: skip_zero_weights=1 is the default; CAP will not read an input seismogram with all-zero weights
# NOTE: skip_zero_weights=0 allows the use of polarities even when all of a station's weights =0 (in the input file);
#       CAP still needs at least 1 waveform for the inversion
# KEY NOTE: skip_zero_weights=0 every entry in the weight file needs to have either a polarity or a wiggle

# FUTURE: Whether skip_zero_weights=0 or 1, what we want is that if there is a station with all-zero weights and no polarity,
#        then we want CAP to simply skip over that station, as if it wasn't in the file.
#-----------------------------------------------------------------------------------
# > make cap

# Run cap with :
# a.  available polarities
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -E1 -K1 -Y2 -Mcus_15/5.0 20080418093700 -Zweight_pfm_all.dat
cp 20080418093700/cus_015.ps $ExampleDir/cus_015_DC_pol.ps
cp 20080418093700/cus_015.out $ExampleDir/cus_015_DC_pol.out

#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $CAPHOME/EXAMPLES/20080418093700_check/cus_015_DC_pol.ps &
# gv $ExampleDir/cus_015_DC_pol.ps &
diff $ExampleDir/cus_015_DC_pol.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_DC_pol.out
echo "diff $ExampleDir/cus_015_DC_pol.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_DC_pol.out"
#-----------------------------------------------------------------------------------

# b. run cap with only one polarity measurement
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -E1 -K1 -Y2 -Mcus_15/5.0 20080418093700 -Zweight_pfm_one.dat 
# gv ./20080418093700/cus_015.ps &

# c. run cap with all polarities but one flipped: PVMO
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -E1 -K1 -Y2 -Mcus_15/5.0 20080418093700 -Zweight_pfm_badPVMO.dat 
# gv ./20080418093700/cus_015.ps &
#==> solution is slightly worse, as PVMO was nearly nodal and is now within a colored quadrant

# d. run cap with all polarities but one flipped: WVT
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -E1 -K1 -Y2 -Mcus_15/5.0 20080418093700 -Zweight_pfm_badWVT.dat 
# gv ./20080418093700/cus_015.ps &
# ==> solution is significantly worse, since WVT was not nodal and is now within a colored quadrant

# e. run cap with all polarities but one flipped: BLO
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -E1 -K1 -Y2 -Mcus_15/5.0 20080418093700 -Zweight_pfm_badBLO.dat 
# gv ./20080418093700/cus_015.ps &
# ==> no solution is possible

# note: the command "-I10/0.1 -J0/0/0/0" could be added but will not affect the results

#==================================================================================
## EXAMPLE 3: Illinois event, DC with P wave polarities in weight file
# Nothing to check
# Given bad polarity observations (eg PVMO flipped from (correct) negative to (incorrect) positive),
# cap may still allow a solution, within error threshold, even if it conflicts with the polarity observation. 
# KEY: threshold should be NEGATIVE if polarities are allowed to conflict expected polarity.

# This example runs 6 cases for different first-motion threshold (fm_thr) with values:
# fm_thr=-1.0, -0.2, -0.1, 0, 0.4, 1.0 (see cap flag -F to set thresholds).

# 0. Follow same step 3 EXAMPLE 1

# 1. Run cap with with 

# a. no MT grid point is rejected, so the solution is same as EXAMPLE 1 (PVMO obs does not match PVMO pred)
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F-1.0 -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -E1 -K1 -Y2 -Mcus_15/5.0 20080418093700 -Zweight_pfm_badPVMO.dat 

# b. a relatively small number of MT grid points are rejected that have large-amplitude mismatch of polarities;
#    this solution, with PVMO not matching PVMO pred, is not rejected beacause the PVMO is low amplitude
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F-0.2 -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -E1 -K1 -Y2 -Mcus_15/5.0 20080418093700 -Zweight_pfm_badPVMO.dat 

# c. PVMO is "just" inside the colored quadrant
#    note: PVMO could be "just" inside the white quadrant as well
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F-0.1 -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -E1 -K1 -Y2 -Mcus_15/5.0 20080418093700 -Zweight_pfm_badPVMO.dat 

# d. PVMO is "just" inside the colored quadrant (PVMO obs does match PVMO pred)
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F0.0  -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -E1 -K1 -Y2 -Mcus_15/5.0 20080418093700 -Zweight_pfm_badPVMO.dat 

# e. the solution seeks to maximize the minimum amplitude associated with any (correct) polarity (PVMO obs does match PVMO pred)
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F0.4  -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -E1 -K1 -Y2 -Mcus_15/5.0 20080418093700 -Zweight_pfm_badPVMO.dat 

# f. no solution, since every MT grid point is rejected with some station having too small an amplitude
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F1.0  -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -E1 -K1 -Y2 -Mcus_15/5.0 20080418093700 -Zweight_pfm_badPVMO.dat 
#==================================================================================
## EXAMPLE 4: Illinois event, FMT with Zhu line search
#  Follow [STEP 0]

# 0. run cap with the -J flag
# ISO and CLVD will be searched for in their general range: (-1 < ISO < +1) and (-0.5 < CLVD < 0.25).
# However you can change the increment by using -J flag in the input command line.
# depth search range for the Illinois event (FMT) - Zhu and BenZion parameterization
# note: -J0/0.1/0/0.1 means line search from (0,0) in increments of 0.1 in each direction

# perform the depth test
# 1. First make sure there are not other *.out files
rm $eid/cus*.ps
rm $eid/cus*.out

# 2. Run cap for different depths by setting the range and increment for which you want to run depth test using -A flag (-Adep_min/dep_max/dep_inc).
# One output (.out) file should be generated for each depth.
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15/5.0 -E0 -K0 -J0/0.1/0/0.1 -A5/30/5 20080418093700

# 3. perform depth test
depth_test 20080418093700 cus 

# 4. Copy only the best solution file
cp 20080418093700/cus_015_fmt.ps $ExampleDir/cus_015_fmt_line.ps
cp 20080418093700/cus_015.out $ExampleDir/cus_015_fmt_line.out
cp 20080418093700/cus_015_beach_fmt.ps $ExampleDir/cus_015_beach_fmt_line.ps
# Copy the depth result files
cp dep_20080418093700.ps $ExampleDir/depth_fmt_line.ps
cp junk1.out $ExampleDir/junk1_fmt_line.out
cp junk2.out $ExampleDir/junk2_fmt_line.out

# 5. Check results
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_line.ps &
# gv $ExampleDir/cus_015_fmt_line.ps &
diff $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_line.out $ExampleDir/cus_015_fmt_line.out
echo "diff $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_line.out $ExampleDir/cus_015_fmt_line.out"
# gv $CAPHOME/EXAMPLES/20080418093700_check/depth_fmt_line.ps &
# gv $ExampleDir/depth_fmt_line.ps &
# diff $CAPHOME/EXAMPLES/20080418093700_check/junk1_fmt_line.out $ExampleDir/junk1_fmt_line.out
diff $CAPHOME/EXAMPLES/20080418093700_check/junk1_fmt_line.out $ExampleDir/junk1_fmt_line.out
echo "diff $CAPHOME/EXAMPLES/20080418093700_check/junk1_fmt_line.out $ExampleDir/junk2_fmt_line.out"

#==================================================================================
## EXAMPLE 5: Illinois event, FMT with fixed focal mech (strike/dip/rake/iso/clvd/mag)
#  Follow [STEP 0]

cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15/5.21 -E0 -K0 -Y2 -R296/296/84/84/6/6 -J0.16/0/0.04/0 -I10/0 20080418093700
cp 20080418093700/cus_015_fmt.ps $ExampleDir/cus_015_fmt_fixed_zhu.ps
cp 20080418093700/cus_015.out $ExampleDir/cus_015_fmt_fixed_zhu.out
cp 20080418093700/cus_015_beach_fmt.ps $ExampleDir/cus_015_beach_fmt_fixed_zhu.ps

#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_fixed_zhu.ps &
# gv $ExampleDir/cus_015_fmt_fixed_zhu.ps &
diff $ExampleDir/cus_015_fmt_fixed_zhu.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_fixed_zhu.out
echo "diff $ExampleDir/cus_015_fmt_fixed_zhu.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_fixed_zhu.out"

#==================================================================================
## EXAMPLE 6: Illinois event, DC, full grid search (no line search)
# This is like EXAMPLE 1, but here we use a single, uniform grid search with no line searches or refined grid searches.
# This is aimed at performing a secondary analysis of uncertainties outside of CAP.

# 1. Run cap 
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15/5.0 -I10/0.1 -E1 -K1 -Y2 -J0/0/0/0 20080418093700
cp 20080418093700/cus_015_fmt.ps $ExampleDir/cus_015_dc_grid.ps 
cp 20080418093700/cus_015.out $ExampleDir/cus_015_dc_grid.out
cp 20080418093700/cus_015_beach_fmt.ps $ExampleDir/cus_015_beach_dc_grid.ps 
# note: -J0/0/0/0 means fixed ISO and CLVD components -- the command can be omitted and the result will be the same
# note: -I here will use a 0.1 interval of Mw

# 2. Check results
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $CAPHOME/EXAMPLES/20080418093700_check/cus_015_dc_grid.ps &
# gv $ExampleDir/cus_015_dc_grid.ps &
diff $ExampleDir/cus_015_dc_grid.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_dc_grid.out
echo "diff $ExampleDir/cus_015_dc_grid.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_dc_grid.out"

#==================================================================================
## EXAMPLE 7: same as EXAMPLE 6 but with fixed magnitude
#  Follow [STEP 0]

# 1. Run cap 
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15/5.2 -I10/0 -E1 -K1 -Y2 -J0/0/0/0 20080418093700
cp 20080418093700/cus_015_fmt.ps $ExampleDir/cus_015_dc_fixed.ps
cp 20080418093700/cus_015.out $ExampleDir/cus_015_dc_fixed.out
cp 20080418093700/cus_015_beach_fmt.ps $ExampleDir/cus_015_beach_dc_fixed.ps

# 2. Check results
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $CAPHOME/EXAMPLES/20080418093700_check/cus_015_dc_fixed.ps &
# gv $ExampleDir/cus_015_dc_fixed.ps &
diff $ExampleDir/cus_015_dc_fixed.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_dc_fixed.out
echo "diff $ExampleDir/cus_015_dc_grid.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_dc_grid.out"

#==================================================================================
## EXAMPLE 8: Illinois event, FMT with lune parameterization, full grid search (no line search)
# 0. Follow same step 3 EXAMPLE 6
#  Follow [STEP 0]

# 1.  run cap (this will take a couple minutes)
# ISO and CLVD will be searched for in their general range (-90 < ISO < +90 ; -30 < CLVD < +30) (lune parameters).
# However you can change the increment by using -J flag in the input command line.
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15/5.0 -I10/0.2/10 -E1 -K1 -Y2 -J-90/90/-30/30 20080418093700
cp 20080418093700/cus_015_fmt.ps $ExampleDir/cus_015_fmt_grid.ps
cp 20080418093700/cus_015.out $ExampleDir/cus_015_fmt_grid.out
cp 20080418093700/cus_015_beach_fmt.ps $ExampleDir/cus_015_beach_fmt_grid.ps
# note: -J0/5/0/5 means 5 deg increment in ISO (delta) and 5 deg increment in CLVD (longitude); the 0's are ignored in this mode
# note: -I here will use a 0.2 interval of Mw and 10 degree for ISO and CLVD lune parameters

# 2. Check results
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $CAPHOME/EXAMPLES/20080418093700_check/cus_015_dc_fixed.ps &
# gv $ExampleDir/cus_015_fmt_grid.ps &
diff $ExampleDir/cus_015_fmt_grid.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_grid.out
echo "diff $ExampleDir/cus_015_fmt_grid.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_grid.out"

#==================================================================================
## EXAMPLE 9: Illinois event, Direct search FMT (LUNE parameterization) with fixed focal mech (strike/dip/rake/iso/clvd/mag)
# 0. Follow same step 3 EXAMPLE 6
#  Follow [STEP 0]

# 1. Run cap 
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15/5.21 -E1 -K1 -Y2 -R296/296/84/84/6/6 -J-30/30/10/10 -I10/0/0 20080418093700
cp 20080418093700/cus_015_fmt.ps $ExampleDir/cus_015_fmt_fixed_lune.ps
cp 20080418093700/cus_015.out $ExampleDir/cus_015_fmt_fixed_lune.out
cp 20080418093700/cus_015_beach_fmt.ps $ExampleDir/cus_015_beach_fmt_fixed_lune.ps
# note: the solution is slightly better than the result in EXAMPLE 4 (why? precision?)

# 2. Check results
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_fixed_lune.ps &
# gv $ExampleDir/cus_015_fmt_fixed_lune.ps &
diff $ExampleDir/cus_015_fmt_fixed_lune.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_fixed_lune.out
echo "diff $ExampleDir/cus_015_fmt_fixed_lune.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_fixed_lune.out"

#==================================================================================
## EXAMPLE 10: Illinois event, DC with random sample search (optional: grid-search over Mw)
#  Follow [STEP 0]

# 1. Run cap 
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15/5.2 -I2/0 -E1 -K2 -Y2 -J0/0/0/0 20080418093700
cp 20080418093700/cus_015_fmt.ps $ExampleDir/cus_015_dc_rand.ps
cp 20080418093700/cus_015.out $ExampleDir/cus_015_dc_rand.out
cp 20080418093700/cus_015_beach_fmt.ps $ExampleDir/cus_015_beach_dc_rand.ps
# NOTES
# - the number of samples is based on the orientation increment (-I) which is 2 deg here; decrease this for more samples
# - If running only at one particular magnitude, set -I magnitude increment to 0 (e.g., -I2/0)
# - To search over magnitude too, use nonzero magnitude increment
# - For any other ISO or CLVD value, set it using -J flag. Just remember to keep the increment=0 (e.g. -J10/0/10/0)
# - the -J0/0/0/0 can be omitted, and the output will be identical

# 2. Check results
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $CAPHOME/EXAMPLES/20080418093700_check/cus_015_dc_rand.ps &
# gv $ExampleDir/cus_015_dc_rand.ps &
diff $ExampleDir/cus_015_dc_rand.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_dc_rand.out
echo "diff $ExampleDir/cus_015_fmt_fixed_lune.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_fixed_lune.out"

#==================================================================================
## EXAMPLE 11: Illinois event, FMT-fixed-magnitude with random sample search (optional: grid-search over Mw)
#  Follow [STEP 0]

# 1. Run cap 
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15/5.2 -I2/0/10 -E1 -K2 -Y2 -J0/1/0/1 20080418093700
cp 20080418093700/cus_015_fmt.ps $ExampleDir/cus_015_fmt_rand.ps
cp 20080418093700/cus_015.out $ExampleDir/cus_015_fmt_rand.out
cp 20080418093700/cus_015_beach_fmt.ps $ExampleDir/cus_015_beach_fmt_rand.ps
# NOTES
# - the number of samples is based on the orientation increment (-I) which is 2 deg here; decrease this for more samples
# - If running only at one particular magnitude, set -I magnitude increment to 0 (e.g., -I2/0)
# - To search over magnitude too, use nonzero magnitude increment
# - the ISO and CLVD increments need to be set to any nonzero number (e.g., -J10/1/10/2 would perform the same task)

# 2. Check results
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_rand.ps &
# gv $ExampleDir/cus_015_fmt_rand.ps &
diff $ExampleDir/cus_015_fmt_rand.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_rand.out
echo "diff $ExampleDir/cus_015_fmt_rand.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_rand.out"

#==================================================================================
## EXAMPLE 12: Illinois event, Fixed FMT (direct search at any ISO, CLVD point)
#  Follow [STEP 0]

# 1. Run cap 
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15/5.2 -I2/0/10 -E1 -K1 -Y1 -J10/10/5/5 20080418093700
cp 20080418093700/cus_015_fmt.ps $ExampleDir/cus_015_fmt_fixed.ps
cp 20080418093700/cus_015.out $ExampleDir/cus_015_fmt_fixed.out
cp 20080418093700/cus_015_beach_fmt.ps $ExampleDir/cus_015_beach_fmt_fixed.ps

# 2. Check results
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_fixed.ps &
# gv $ExampleDir/cus_015_fmt_fixed.ps &
diff $ExampleDir/cus_015_fmt_fixed.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_fixed.out
echo "diff $ExampleDir/cus_015_fmt_fixed.out $CAPHOME/EXAMPLES/20080418093700_check/cus_015_fmt_fixed.out"