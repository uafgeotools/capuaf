#!/bin/bash
# 
# This README_EXAMPLE_utu_fmp runs CAP in FMP mode and computes a summary of
# first motion polarities for the Uturuncu main event.
#
# NOTE
#   The CAP run in this example creates file out.misfit.fmp_
#
#   This file (out.misfit.fmp_) is 700+ MB for a grid spacing of 5 degrees (lune
#   and orientation) and can be deleted after computing summary
#
#   Other reference scripts
#       README_EXAMPLE_utu_fk   -- compute library of Greens functions for Uturuncu main event
#       README_EXAMPLE_utu_FTC  -- shows improvement of waveform filtering when applying filter filter-then-cut
# 
# For more details see
#
# @article{AlvizuriTape2016,
#      AUTHOR = {C. Alvizuri and C. Tape},
#      TITLE = {{Full moment tensors for small events ($\mw < 3$) at Uturuncu volcano, Bolivia}},
#      JOURNAL = {Geophys.~J.~Int. \rm(in prep.)},
#      PAGES = {},
#      VOLUME = {},
#      NUMBER = {},
#      EID = {},
#      DOI = {},
#      YEAR = {2015}
# }
# 
# 20151110 celso alvizuri - cralvizuri@alaska.edu 
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

