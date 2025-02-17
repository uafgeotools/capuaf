#!/bin/bash
# full moment tensor inversions of the Uturuncu main event

# flags to set if running this file as a shell script.
# 0 = FALSE, 1 = TRUE
run_nopol=0  # run or not inversions without first motion polarities (3.5 hrs to complete, per event)
run_depth=0  # run or not depth test (30 min to complete) 
replace=0    # replace or not figures in 20100516063454464_check 

# ==============================================================================
# Part 1. check that cap.c is compiled with the right flags, verify weight files

# 1. go to cap directory and make sure cap is set to its default condition
cd $CAPHOME

# 2. Make sure cap.c (git checkout cap.c) is compiled with the following flags, then compile.
# NOTE The following examples will run with these set of flags, except for first
# motion polarity example and runs for event 20110622023324299
# FTC_data=1
# FTC_green=0
# skip_zero_weights=0
make cap

# 3. go to run directory
cd $CAPHOME/EXAMPLES

# 4. check values in the file w_utuhalf_P01_V10_R01_S10.dat
# NOTE 
# File name reflects weights for each waveform component, eg: 
#   (PV - PR - SurfV - SurfR - SurfT)  = 10 - 1 - 10 - 10 - 10
# All P arrivals times were obtained from analyst picks
# Surface wave time shifts are set to 0.0 for utuhalf
more 20100516063454464/w_utuhalf_P01_V10_R01_S10.dat

# ==============================================================================
# Part 2. run cap inversion with polarities + body waves + surface waves
# NOTE 
# ISO will be gridded by considering equal spacing in sin(iso) space
# DIP will be gridded by considering equal spacing in cos(dip) space
# All weights should be set as 10 1 10 10 10, which gives higher weights to PV component and surface waves
# All P arrivals times were obtained from analyst picks

# 1. go to the data directory and run cap. 
# NOTE fits for body waves will be scaled by an absolute amplitude, then multiplied by a factor of 3. (this is done with flag "-P1/...")
cd $CAPHOME/EXAMPLES
cap.pl -H0.01 -P0.58/1/30 -p0.67 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -I5/0.1/5 -J-90/90/-30/30 -E1 -K1 -Y1 -Mutuhalf_4/2.4 -Zw_utuhalf_P01_V10_R01_S10.dat 20100516063454464

# 2. compare results
if [ $replace -eq 1 ] ; then
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S10_fmt_amp_abs.ps
fi
gv 20100516063454464/utuhalf_004_fmt.ps &
gv 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S10_fmt_amp_abs.ps &

# 4. run cap
cap.pl -H0.01 -P0/1/30 -p0 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -I5/0.1/5 -J-90/90/-30/30 -E1 -K1 -Y1 -Mutuhalf_4/2.4 -Zw_utuhalf_P01_V10_R01_S10.dat 20100516063454464 

# 5. compare results and verify .out file
if [ $replace -eq 1 ] ; then
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S10_fmt.ps
cp 20100516063454464/utuhalf_004.out 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S10.out
fi
gv 20100516063454464/utuhalf_004_fmt.ps &
gv 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S10_fmt.ps
diff 20100516063454464/utuhalf_004.out 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S10.out

# NOTE 
# 1. Flag -P0/... normalizes seismograms in each window
# 2. Lune grid will NOT do a refined search, so pick -I accordingly (here 5 deg increments on the lune, 0.1 increments in Mw)
# 3. All P arrivals are hand-picked (to within 0.01 s), we set the allowable time shift for P to 0.1 s
# 4. Flag -L gives the half-duration of the synthetics. This is a key command to get high-frequency synthetics in P (but not too high)
# 5. This result is an improvement of the BSSA2014 poster (same event)
# 6. surface wave time shifts are allowed to be between 0.0 s and 4.0 s, which reflects our belief that either:
# (a) the model is systematically too slow
# (b) the depth is wrong
# (c) the origin time is wrong
# 7. We are able to specify the 0.0 to 4.0 s (*) shifts with two changes:
# (1) the relative surface time shifts are +/- 2.0, which is set as -S0.2/2.0/0
# (2) the systematic S_shifts in weight file (last column) are set to 2.0 
# (*) actual time shifts are then 0.0 s (-2.0 + 2.0) and 4.0 s (+2.0 + 2.0)
 
# -----------------------------------------------------------
# 6. Run cap but using different weights for  waveform components
# NOTE Weights are set as 1 1 60 60 60, which gives most weight to surface waves
cap.pl -H0.01 -P0.58/1/30 -p0.67 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -I5/0.1/5 -J-90/90/-30/30 -E1 -K1 -Y1 -Mutuhalf_4/2.4 -Zw_utuhalf_P01_V01_R01_S60.dat 20100516063454464 
if [ $replace -eq 1 ] ; then 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V01_R01_S60_fmt.ps
cp 20100516063454464/utuhalf_004.out 20100516063454464_check/20100516063454464_utuhalf_004_P01_V01_R01_S60.out
fi

# 7. compare results
gv 20100516063454464/utuhalf_004_fmt.ps &
gv 20100516063454464_check/20100516063454464_utuhalf_004_P01_V01_R01_S60_fmt.ps &

# -----------------------------------------------------------
# Part 3. run depth test using all waveform components and first motion polarities

# 1. Run cap inversions at each depth
# NOTE this task takes about 30min to complete
if [ $run_depth -eq 1 ] ; then 
for depth in 1 2 3 4 5 6 7 8 9 10; do cap.pl -W1 -Y1 -D1/1/0.5 -H0.01 -C2/10/0.25/0.5 -F0 -L0.1 -S0.05/2.0/0 -P0.58/1/30 -p0.67 -T1.5/60 -E1 -K1 -I5/0.1/5 -J-90/90/-30/30 -Mutuhalf_$depth/2.8 -Zw_utuhalf_P01_V10_R01_S10.dat 20100516063454464 ; done

# 2. Plot the result
# NOTE to reproduce figure in the Uturuncu FMT paper open script depth.pl and 
# set flag "$onlydc = 0" to plot beachballs as full moment tensors. Then
# uncomment lines #288, and comment line #337 and uncomment #338 (or follow instruction there)
depth_test 20100516063454464 utuhalf

# 3. Compare results
gv dep_20100516063454464.ps &
gv 20100516063454464_check/depth_test_20100516063454464_utuhalf_L_dt.eps &
fi

# -----------------------------------------------------------
# Part 4. run cap inversions with and without polarities, with body waves. Then with body waves + surface waves.
# NOTE runs without polarities take about 3.5 hours to complete.

if [ $run_nopol -eq 1 ] ; then
# 1. Run cap with body waves + surface waves. No polarities. Component weights 10-1-10
# NOTE runs without polarities take about 3.5 hours to complete.
cap.pl -H0.01 -P0.58/1/30 -p0.67 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -E1 -K1 -Y1 -I5/0.1/5 -J-90/90/-30/30 -Mutuhalf_4/2.9 -Zw_utuhalf_P00_V10_R01_S10.dat 20100516063454464 
if [ $replace -eq 1 ] ; then 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P00_V10_R01_S10_fmt.ps
cp 20100516063454464/utuhalf_004.out 20100516063454464_check/20100516063454464_utuhalf_004_P00_V10_R01_S10.out
fi

# 2. Compare results
gv 20100516063454464/utuhalf_004_fmt.ps &
gv 20100516063454464_check/20100516063454464_utuhalf_004_P00_V10_R01_S10_fmt.ps
fi

# 3. Run cap with polarities + body waves. No surface waves. Component weights 10-1-0
cap.pl -H0.01 -P0.58/1/30 -p0.67 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -E1 -K1 -Y1 -I5/0.1/5 -J-90/90/-30/30 -Mutuhalf_4/2.9 -Zw_utuhalf_P01_V10_R01_S00.dat 20100516063454464 
if [ $replace -eq 1 ] ; then 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00_fmt.ps
cp 20100516063454464/utuhalf_004.out 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00.out
fi

# 4. Compare results
gv 20100516063454464/utuhalf_004_fmt.ps &
gv 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00_fmt.ps &

# 5. Run with body waves only. No surface waves. No polarities. Component weights 10-1-0
# NOTE runs without polarities take about 3.5 hours to complete.
if [ $run_nopol -eq 1 ] ; then
cap.pl -H0.01 -P0.58/1/30 -p0.67 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -E1 -K1 -Y1 -I5/0.1/5 -J-90/90/-30/30 -Mutuhalf_4/2.9 -Zw_utuhalf_P00_V10_R01_S00.dat 20100516063454464 
if [ $replace -eq 1 ] ; then 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P00_V10_R01_S00_fmt.ps
cp 20100516063454464/utuhalf_004.out 20100516063454464_check/20100516063454464_utuhalf_004_P00_V10_R01_S00.out
fi

# 6. Compare results
gv 20100516063454464/utuhalf_004_fmt.ps &
gv 20100516063454464_check/20100516063454464_utuhalf_004_P00_V10_R01_S00_fmt.ps
fi

# -----------------------------------------------------------
# Part 5. run cap inversions with individual components. All inversions with polarities.
# 1. Run cap with component PV only
cap.pl -H0.01 -P0.58/1/30 -p0.67 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -E1 -K1 -Y1 -I5/0.1/5 -J-90/90/-30/30 -Mutuhalf_4/2.9 -Zw_utuhalf_P01_V10_R00_S00.dat 20100516063454464 
if [ $replace -eq 1 ] ; then 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R00_S00_fmt.ps
cp 20100516063454464/utuhalf_004.out 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R00_S00.out
fi

# 2. Compare results
gv 20100516063454464/utuhalf_004_fmt.ps &
gv 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R00_S00_fmt.ps &

# 3. Run cap with component PR only
cap.pl -H0.01 -P0.58/1/30 -p0.67 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -E1 -K1 -Y1 -I5/0.1/5 -J-90/90/-30/30 -Mutuhalf_4/2.4 -Zw_utuhalf_P01_V00_R01_S00.dat 20100516063454464 
if [ $replace -eq 1 ] ; then 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V00_R01_S00_fmt.ps
cp 20100516063454464/utuhalf_004.out 20100516063454464_check/20100516063454464_utuhalf_004_P01_V00_R01_S00.out
fi

# 4. Compare results
gv 20100516063454464/utuhalf_004_fmt.ps &
gv 20100516063454464_check/20100516063454464_utuhalf_004_P01_V00_R01_S00_fmt.ps &

# 5. Run cap with component Surf only
cap.pl -H0.01 -P0.58/1/30 -p0.67 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -E1 -K1 -Y1 -I5/0.1/5 -J-90/90/-30/30 -Mutuhalf_4/2.9 -Zw_utuhalf_P01_V00_R00_S10.dat 20100516063454464 
if [ $replace -eq 1 ] ; then 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V00_R00_S10_fmt.ps
cp 20100516063454464/utuhalf_004.out 20100516063454464_check/20100516063454464_utuhalf_004_P01_V00_R00_S10.out
fi

# 6. Compare results
gv 20100516063454464/utuhalf_004_fmt.ps &
gv 20100516063454464_check/20100516063454464_utuhalf_004_P01_V00_R00_S10_fmt.ps &

# -----------------------------------------------------------
# Part 6. run cap inversions with combinations of components. All inversions with polarities.
# 1. polarities + body waves (PV + PR)
# NOTE this is the same inversion as in Part 4, step 3: Run cap with body waves. No surface waves. No polarities. Component weights 10-1-0
# The example is repeated here to have a complete set of tests.
cap.pl -H0.01 -P0.58/1/30 -p0.67 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -E1 -K1 -Y1 -I5/0.1/5 -J-90/90/-30/30 -Mutuhalf_4/2.9 -Zw_utuhalf_P01_V10_R01_S00.dat 20100516063454464

# 2. compare results
gv 20100516063454464/utuhalf_004_fmt.ps &
gv 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00_fmt.ps &

# 3. polarities + Surface waves
cap.pl -H0.01 -P0.58/1/30 -p0.67 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -I5/0.1/5 -J-90/90/-30/30 -E1 -K1 -Y1 -Mutuhalf_4/2.9 -Zw_utuhalf_P01_V00_R00_S10.dat 20100516063454464
if [ $replace -eq 1 ] ; then 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V00_R00_S10_fmt.ps
cp 20100516063454464/utuhalf_004.out 20100516063454464_check/20100516063454464_utuhalf_004_P01_V00_R00_S10.out
fi

# 4. compare results
gv 20100516063454464/utuhalf_004_fmt.ps &
gv 20100516063454464_check/20100516063454464_utuhalf_004_P01_V00_R00_S10_fmt.ps &

# 5. PV + PR + Surf. No polarities
# NOTE runs without polarities take about 3.5 hours to complete.
# NOTE this is the same inversion as in Part 4, step 1: Run cap body waves + surface waves. No polarities. Component weights 10-1-10
# The example is repeated here to have a complete set of tests.
if [ $run_nopol -eq 1 ] ; then
cap.pl -H0.01 -P0.58/1/30 -p0.67 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -I5/0.1/5 -J-90/90/-30/30 -E1 -K1 -Y1 -Mutuhalf_4/2.4 -Zw_utuhalf_P00_V10_R01_S10.dat 20100516063454464
fi

# 6. compare results
gv 20100516063454464/utuhalf_004_fmt.ps &
gv 20100516063454464_check/utuhalf_004_P00_V10_R01_S10_fmt.ps &

# -----------------------------------------------------------
# Part 7. run cap inversions for specific physical models:  double couple, opening crack, closing crack
# Inversions are with polarities + body waves.

# 1. polarities + PV + PR. No surface waves. Constrained to DC.
cap.pl -H0.01 -P0.58/1/30 -p0.67 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -E1 -K1 -Y1 -I5/0.1 -Mutuhalf_4/2.9 -Zw_utuhalf_P01_V10_R01_S00.dat 20100516063454464 
if [ $replace -eq 1 ] ; then 
cp 20100516063454464/utuhalf_004.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00_dc.ps
cp 20100516063454464/utuhalf_004.out 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00_dc.out
fi

# 2. Compare results
gv 20100516063454464/utuhalf_004.ps &
gv 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00_dc.ps &

# 3. Run cap for opening crack (oc)
# NOTE start magnitude is 2.4
cap.pl -H0.01 -P0.58/1/30 -p0.67 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -E1 -K1 -Y1 -I5/0.1/5 -J-90/90/-30/-30 -Mutuhalf_4/2.4 -Zw_utuhalf_P01_V10_R01_S00.dat 20100516063454464 
if [ $replace -eq 1 ] ; then 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00_co.ps
cp 20100516063454464/utuhalf_004.out 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00_co.out
fi

# 4. Compare results
gv 20100516063454464/utuhalf_004_fmt.ps &
gv 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00_co.ps &

# 5. Run cap for closing crack (cc)
cap.pl -H0.01 -P0.58/1/30 -p0.67 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -E1 -K1 -Y1 -I5/0.1/5 -J-90/90/30/30 -Mutuhalf_4/2.9 -Zw_utuhalf_P01_V10_R01_S00.dat 20100516063454464 
if [ $replace -eq 1 ] ; then 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00_cc.ps
cp 20100516063454464/utuhalf_004.out 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00_cc.out
fi

# 6. Compare results
gv 20100516063454464/utuhalf_004_fmt.ps &
gv 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00_cc.ps &

# ==============================================================================
