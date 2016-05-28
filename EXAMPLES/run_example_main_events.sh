#!/bin/sh
#
# run_example_main_events.sh
# This script runs CAP for the main events for 
# Illinois(Default CAP example),
# Uturuncu volcano (Alvizuri & Tape, 2016),
# Alaska (Silwal & Tape, 2016)
#
# This script should be kept simple. The goal is to easily and quickly check
# whether CAP is producing sensible results.
#
#
# 20160519 celso alvizuri - cralvizuri@alaska.edu 
#----------------------------------------------------------

# compile cap
cd ~/REPOSITORIES/cap/
make cap

# main event examples
cd /REPOSITORIES/cap/EXAMPLES/

# Uturuncu
# This example will perform a full grid search at constant magnitude.
# NOTE This example searches through 20M+ solutions. 
#      Timings on eagle are <17 min without openMP, <5 min with openMP

sampFMT="-I13/35/73/18/37"  # GRID sampling. dense, similar to Uturuncu
sampFMT="-I7/17/35/9/35"    # GRID sampling. not dense, produces similar results to dense
sampFMT="-I22121190"        # RAND sampling. (T = 17m3.815s) dense
sampFMT="-I1311975"         # RAND sampling. (T = 1m2.073s) not dense, produces similar results to dense
cap.pl -H0.01 -P0.58/1/30 -p0.67 -S0.05/2.0/0 -T1.5/60 -D1/1/0.5 -C2/10/0.25/0.5 -F0.0 -L0.1 -W1 -Y1 -Mutuhalf_4 -m2.8 -Zw_utuhalf_P01_V10_R01_S10.dat 20100516063454464 $sampFMT

# check result
diff 20100516063454464_check/20100516063454464_utuhalf_004_P00_V10_R01_S10.out OUTPUT_DIR/20100516063454464_utuhalf_004.out

# Alaska DC
cap.pl -H0.02 -P1/15/60 -p1 -S2/10/0 -T15/120 -D1/1/0.5 -C0.25/0.6667/0.025/0.0625 -W1 -Mscak_39 -m4.0/5.0/0.1 -I1/1/36/10/19 -R0/0 -Zweight111_subset.dat -Y1 20090407201255351
# check result
diff 20090407201255351_check/20090407201255351_main.out OUTPUT_DIR/20090407201255351_scak_039.out

# Alaska FMT
sampFMT="-I13/35/73/18/37"  # GRID sampling. dense, similar to Uturuncu
sampFMT="-I7/17/35/9/35"    # GRID sampling. not dense, produces similar results to dense
sampFMT="-I22121190"        # RAND sampling. (T = 17m3.815s) dense
sampFMT="-I1311975"         # RAND sampling. (T = 1m2.073s) not dense
cap.pl -H0.02 -P1/15/60 -p1 -S2/10/0 -T15/120 -D1/1/0.5 -C0.25/0.6667/0.025/0.0625 -W1 -Mscak_39 -m4.0/5.0/0.1 $sampFMT -Zweight111_subset.dat -Y1 20090407201255351

# Illinois
cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15 -m5.2 -Y1 20080418093700 -I1/1/73/18/37

# check result
diff 20080418093700_check/cus_015.out OUTPUT_DIR/20080418093700_cus_015.out

