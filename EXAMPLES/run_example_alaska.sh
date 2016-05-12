#!/bin/bash
#==================================================================================
# Alaska examples 
#
# The figures generated here were used in following papers (See these papers for more details):
# - SilwalTape2016 paper (JGR)
# - Tape2013 (EPSL)
#
# $eid = 20090407201255351 (Anchorage event) - SilwalTape2016 paper (JGR)
# $eid = 20120411092157444 (Nenana triggered event) - Tape2013 (EPSL)
# $CAPHOME - Directory where cap.c is located
# $CAPRUN/inv/cus - Directory from where you will be running cap
# $ExampleDir - Directory where you want to generate new plots and output files
#
# Runs and save in the $ExampleDir
# MIGHT want to rename the previous $ExampleDir and save it (just in case)
# Search $ExampleDir and make sure it does not overwrite anything
# You will be comparing the results between $ExampleDir and $CAPHOME/EXAMPLES/20080418093700_check
#====================================================================================

# EXAMPLE 1: Running Alaska Example (20090407201255351)
# These figures were used in SilwalTape2016 and the supplementary data stored in ScholarWorks@UA collections

# 1. Prepare cap.c 
# cap.c is in $CAPRUN (see README_EXAMPLE) OR
# > which cap.c

#!!!!!!!! MAKE CHANGES IN cap.c !!!!!!!!!!!!!
# > cd $CAPHOME
# > git checkout cap.c
# > open cap.c and set: 
#    FTC_data=1
#    FTC_green=0 
# then recompile 
# > make all
# FTC means Filter-then-Cut (the data and the green)

# 2. Prepare the data directory and directory for storing the results ($ExampleDir)
eid=20090407201255351
model=scak
mwa=4.0
mwb=5.0
dep=39

# Create a new results directory (CAUTION: ExampleDir will be overwritten)
ExampleDir=$CAPRUN/inv/"$model"/${eid}_check5
mkdir -p $ExampleDir

# Prepare data directory. Copy the the data to directory where you will be running cap
# All commands will be run form here (Alternatively, one could run from $CAPHOME/EXAMPLES/ without copying the data)
cd $CAPRUN/inv/"$model"
rsync -av $CAPHOME/EXAMPLES/$eid .

#========main paper figure (subset of stations)======================
# 3. Run cap
cap.pl -H0.02 -P1/15/60 -p1 -S2/10/0 -T15/120 -D1/1/0.5 -C0.25/0.6667/0.025/0.0625 -W1 -M"$model"_$dep -m$mwa/$mwb/0.1 -I1/1/36/10/19 -R0/0 -Zweight111_subset.dat -Y1 $eid

# 4. Save output files for comparison
cp ${eid}/scak_039.ps $ExampleDir/${eid}_main.ps
cp ${eid}/scak_039.out $ExampleDir/${eid}_main.out
cp ${eid}/scak_039_beach.ps $ExampleDir/${eid}_beach_main.ps

# 5. Compare results
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $ExampleDir/${eid}_main.ps $CAPHOME/EXAMPLES/20090407201255351_check/20090407201255351_main.ps &
echo "diff $ExampleDir/${eid}_main.out $CAPHOME/EXAMPLES/20090407201255351_check/20090407201255351_main.out"
diff $ExampleDir/${eid}_main.out $CAPHOME/EXAMPLES/20090407201255351_check/20090407201255351_main.out

#===== Depth Search - Using different norms and weight files ====================
# 6. Run cap for different depths, different norms (L1 and L2) and weights
# Eight different weight files options: weightijk.dat (Last three numbers are for (i)-body wave, (j)-surface waves, (k)-station choice)
#   i = 1 or 0 for whether body waves are used or not
#   j = 1 or 0 for whether surface waves are used or not
#   k = 0 for all stations, 1 for selected stations, 2 for AEC stations
#
# -> weight011.dat - only surface waves at selected good stations (M011)
# -> weight012.dat - only surface waves at AEC stations (M012)
# -> weight101.dat - only body waves at selected good stations (M101)
# -> weight110.dat - both body and surface wave at all available stations (M110)
# -> weight111.dat - both body and surface wave at selected good stations (M111)
# -> weight112.dat - both body and surface wave at AEIC stations (M112)
# -> az_sort_weight110.dat - same as weight110.dat but azimuthally sorted stations
# -> az_sort_weight111.dat - same as weight111.dat but azimuthally sorted stations

wts=(101 011 112 012 110 111) # weight files identifier
#wts=(111)
norms=(1 2)
#norms=(1)
dep_min=15
dep_max=75
dep_inc=2

# Run for all weights and norms
for ii in ${wts[@]}
do
    for norm in ${norms[@]}
    do 
	rm "$eid"/"$model"*.out
	rm "$eid"/"$model"*.ps
	mkdir -p L$norm/M$ii
	cap.pl -H0.02 -P1/15/60 -p1 -S2/10/0 -T15/120 -D1/1/0.5 -C0.25/0.6667/0.025/0.0625 -W1 -M"$model"_41 -m$mwa/$mwb/0.1 -I1/1/36/10/19 -R0/0 -Zweight$ii.dat -Y$norm -A$dep_min/$dep_max/$dep_inc $eid
	depth_test "$eid" "$model"
	mv "$eid"/"$model"*.out L$norm/M$ii/
	mv "$eid"/"$model"*.ps L$norm/M$ii/
	mv dep_"$eid".ps L$norm/M$ii/
	mv junk1.out L$norm/M$ii/
	mv junk2.out L$norm/M$ii/
	# Find the best solution (minimum misfit depth and fault plane solution) and move the output files for that depth to the ExampleDir (also rename: see mindep.pl)
	$CAPHOME/EXAMPLES/mindep.pl $eid $norm $ii $ExampleDir
    done
done

# 7. Check results
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $ExampleDir/20090407201255351_L1_M111.ps $CAPHOME/EXAMPLES/20090407201255351_check/20090407201255351_L1_M111.ps &
echo "diff $ExampleDir/20090407201255351_L1_M111.out $CAPHOME/EXAMPLES/20090407201255351_check/20090407201255351_L1_M111.out"
diff $ExampleDir/20090407201255351_L1_M111.out $CAPHOME/EXAMPLES/20090407201255351_check/20090407201255351_L1_M111.out
# In the same way you can check all other files (different combinations of norm (L1 or L2) and weights (see $wts above)) 
#====================================================================================

# 8. Plot waveforms sorted by azimuth
# Azimuthal soritng example
cap.pl -H0.02 -P1/20/60 -p1 -S3/10/0 -T15/120 -D1/1/0.5 -C0.2/0.5/0.025/0.06 -W1 -Mscak_41 -m4.5 -I1/1/36/10/19 -R0/0 -Zaz_sort_weight111.dat -Y1 $eid
cp "$eid"/scak_041.ps $ExampleDir/20090407201255351_L1_M111_az.ps 
cp "$eid"/scak_041.out $ExampleDir/20090407201255351_L1_M111_az.out
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $ExampleDir/20090407201255351_L1_M111_az.ps $CAPHOME/EXAMPLES/20090407201255351_check/20090407201255351_L1_M111_az.ps
echo "diff $ExampleDir/20090407201255351_L1_M111_az.out $CAPHOME/EXAMPLES/20090407201255351_check/20090407201255351_L1_M111_az.out"
diff $ExampleDir/20090407201255351_L1_M111_az.out $CAPHOME/EXAMPLES/20090407201255351_check/20090407201255351_L1_M111_az.out
#====================================================================================

# 9. Random search
cap.pl -H0.02 -P1/20/60 -p1 -S3/10/0 -T15/120 -D1/1/0.5 -C0.2/0.5/0.025/0.06 -W1 -Mscak_41 -m4.5 -I100000 -R0/0 -Zweight111.dat -Y1 $eid
cp "$eid"/scak_041.ps $ExampleDir/20090407201255351_L1_M111_rand.ps 
cp "$eid"/scak_041.out $ExampleDir/20090407201255351_L1_M111_rand.out
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $ExampleDir/20090407201255351_L1_M111_rand.ps $CAPHOME/EXAMPLES/20090407201255351_check/20090407201255351_L1_M111_rand.ps
echo "diff $ExampleDir/20090407201255351_L1_M111_rand.out $CAPHOME/EXAMPLES/20090407201255351_check/20090407201255351_L1_M111_rand.out"
diff $ExampleDir/20090407201255351_L1_M111_rand.out $CAPHOME/EXAMPLES/20090407201255351_check/20090407201255351_L1_M111_rand.out

# 10. Input source function
# run source_example.m
cap.pl -H0.02 -P1/20/60 -p1 -S3/10/0 -T15/120 -D1/1/0.5 -C0.2/0.5/0.025/0.06 -W1 -Mscak_41 -m4.5 -I100000 -Zweight111_subset.dat -Y1 -R0/0 -L./sin_source.sac 20090407201255351

# EXAMPLE 2: Running Nenana triggering example (20120411092157444)
# This was used in Tape2013 (Nenana eq triggering paper) and Tape2015 (MFSZ paper)
# Plots are slightly different because older version (not documented which one) was used
# for making the paper figures.
# Figures here are improvement over the paper figures
# Note: Plots match with figures in the supplement of Tape2015 (MFSZ paper)

# 1. Prepare cap.c 
#!!!!!!!! MAKE CHANGES IN CAP.C!!!!!!!!!!!!!
# > cd $CAPHOME
# > git checkout cap.c
# > open cap.c and set: 
#    FTC_data=1
#    FTC_green=0 
# then recompile 
# > make all
# FTC means Filter-then-Cut (the data and the green)

# 2. Prepare the data directory and directory for storing the results ($ExampleDir)
eid=20120411092157444
model=tactmod

ExampleDir=$CAPRUN/inv/"$model"/${eid}_check5
mkdir -p $ExampleDir

cd $CAPRUN/inv/"$model"
rsync -av $CAPHOME/EXAMPLES/$eid .

rm -f $eid/"$model"*ps
rm -f $eid/"$model"*out
# 3. Run cap
cap.pl -H0.02 -P1/1/15 -p1 -S0.5/2/0 -T2/40 -D1/1/0.5 -C1/10/0.3/0.6 -W1 -M"$model"_16 -m3/4/0.1 -I1/1/36/10/19 -R0/0 -Zweight111.dat -Y1 -A10/30/1 $eid 
depth_test $eid "$model"

# 4. Save output files for comparison
cp $eid/"$model"*ps $ExampleDir
cp $eid/"$model"*out $ExampleDir
cp dep_$eid.ps $ExampleDir

# 5. Compare results
# Minimum depth occurs at 16 km
# Also check SAC header for AEIC catalog depth results
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $ExampleDir/tactmod_016.ps $CAPHOME/EXAMPLES/20120411092157444_check/tactmod_016.ps &
echo "diff $ExampleDir/tactmod_016.out $CAPHOME/EXAMPLES/20120411092157444_check/tactmod_016.out &"
diff $ExampleDir/tactmod_016.out $CAPHOME/EXAMPLES/20120411092157444_check/tactmod_016.out &