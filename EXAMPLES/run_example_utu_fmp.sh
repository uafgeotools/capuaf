#!/bin/bash
# ==============================================================================
# first motion polarities for the Uturuncu main event.
#
# NOTE
#   The CAP run in this example creates file out.misfit.fmp_
#
#   This file (out.misfit.fmp_) is 700+ MB for a grid spacing of 5 degrees (lune
#   and orientation) and can be deleted after computing summary
# ==============================================================================

# 1. Set cap.c to its original condition
cd $CAPHOME
git checkout cap.c 

# 2. Make sure cap is compiled with the following flags, then compile
#
# only_first_motion=1
# misfit_on_lune=0
# FTC_data=1
# FTC_green=0
# skip_zero_weights=0

# 3. Compile cap
make cap

# 4. Run cap
# NOTE this run produces file out.misfit.fmp_ 
cd EXAMPLES
cap.pl -E1 -K1 -Y1 -J-90/90/-30/30 -I5/0/5 -Mutuhalf_4/2.8 -Zw_utuhalf_P01_V10_R01_S10.dat 20100516063454464

# 5. Compute summary of first motion polarities using the output file computed above.
# NOTE when finished python will create the summary file:
# out_misfit_fmp_summary, this file 
python $CAPHOME/UTILS/celso/get_fmt_misfit_fmp.py out.misfit.fmp_

# 4. Plot polarity misfit on the lune with file out.misfit_summary
# PENDING

# ==============================================================================
