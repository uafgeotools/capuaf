#!/bin/bash
# ==============================================================================
# shows filter improvement for small signals and short time windoes when using flag filter-then-cut,
# as opposed to the default CAP which first cuts the data and then filters it (cut-then-filter)
#
# Uturuncu event 20110622023324299
#
# NOTE this should NOT be run as a script since each example requires to compile CAP with different flags
#
# NOTE FTC_data=0 (original CAP) first cuts the waveform windows, then filters each windows (cut-then-filter).
# NOTE FTC_data=1 first filters the waveform, then cuts the waveform into windows (filter-then-cut).
#==============================================================================

# 0. go to cap directory and make sure cap is set to its default condition
cd $CAPHOME
git checkout cap.c

# 1. Use option cut-then-filter. In cap.c set the flags below, then compile and run
# FTC_data=0
# FTC_green=0
# skip_zero_weights=0
make cap
cd EXAMPLES
cap.pl -H0.01 -P1/1/30 -p1 -S0.1/2.0/0 -T1.5/60.0 -D1/1/0.5 -C2/10/0.1/0.4 -W1 -Y1 -F0 -L0.05 -E1 -K1 -I5/0/5 -J33.748989/33.748989/30/30 -R35/35/49.864567/49.864567/40/40 -Mutuhalf_9/0.5 -Zweight.dat 20110622023324299
if [ $replace -eq 1] ; then
cp 20110622023324299/utuhalf_009_fmt.ps 20110622023324299_check/utuhalf_009_fmt_cut_then_filter.ps
fi

# 2. compare results
gv 20110622023324299/utuhalf_009_fmt.ps &
gv 20110622023324299_check/utuhalf_009_fmt_cut_then_filter.ps &

# 3. Use option filter-then-cut. In cap.c set the flags below, then compile and run
# FTC_data=1
# FTC_green=0
# skip_zero_weights=0
make cap
cd $CAPHOME/EXAMPLES
cap.pl -H0.01 -P1/1/30 -p1 -S0.1/2.0/0 -T1.5/60.0 -D1/1/0.5 -C2/10/0.1/0.4 -W1 -Y1 -F0 -L0.05 -E1 -K1 -I5/0/5 -J33.748989/33.748989/30/30 -R35/35/49.864567/49.864567/40/40 -Mutuhalf_9/0.5 -Zweight.dat 20110622023324299
if [ $replace -eq 1] ; then
cp 20110622023324299/utuhalf_009_fmt.ps 20110622023324299_check/utuhalf_009_fmt_filter_then_cut_obs.ps
fi

# 4. compare results
gv 20110622023324299/utuhalf_009_fmt.ps &
gv 20110622023324299_check/utuhalf_009_fmt_filter_then_cut_obs.ps &

# 5. Use option filter-then-cut for both data and greens. In cap.c set the flags below, then compile and run 
# NOTE plotting flag -P and -p are set to 0 which normalize seismograms in each window
# FTC_data=1
# FTC_green=1
# skip_zero_weights=0
make cap
cd $CAPHOME/EXAMPLES
cap.pl -H0.01 -P0/1/30 -p0 -S0.1/2.0/0 -T1.5/60.0 -D1/1/0.5 -C2/10/0.1/0.4 -W1 -Y1 -F0 -L0.05 -E1 -K1 -I5/0/5 -J33.748989/33.748989/30/30 -Mutuhalf_9/0.5 -Zweight.dat 20110622023324299
if [ $replace -eq 1] ; then
cp 20110622023324299/utuhalf_009_fmt.ps 20110622023324299_check/utuhalf_009_fmt_filter_then_cut_obs_greens.ps
fi

# 6. compare results
gv 20110622023324299/utuhalf_009_fmt.ps &
gv 20110622023324299_check/utuhalf_009_fmt_filter_then_cut_obs_greens.ps &

# ==============================================================================

