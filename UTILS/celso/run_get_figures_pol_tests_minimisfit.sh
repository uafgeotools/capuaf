#!/bin/bash
#
# Commands to generate figure 10.
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
#
# 20151013 celso alvizuri - cralvizuri@alaska.edu 
#-----------------------------------------------------------

# NOTE this set of tests require to change and run fmt_misfit_mini.gmt twice
#
# without polarities
# datadir_main="/home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_${pol}_V10_R01_S10_L_dt"
#

dir0="/home/alvizuri/shared/plutons/inversions/inv20151028_utupaper_FMP_effects"
dirf="/home/alvizuri/manuscripts/2015/fmt_uturuncu/figures/fmp_influence"

#-----------------------------------------------------------
# figure 10
sh fmt_misfit_mini.gmt 20110528061823523 utuhalf 9 0 10 01 10
sh fmt_misfit_mini.gmt 20110528061823523 utuhalf 9 1 10 01 10 nocpt

sh fmt_misfit_mini.gmt 20110708145821435 utuhalf 9 0 10 01 10
sh fmt_misfit_mini.gmt 20110708145821435 utuhalf 9 1 10 01 10 nocpt

sh fmt_misfit_mini.gmt 20110823072642794 utuhalf 9 0 10 01 10
sh fmt_misfit_mini.gmt 20110823072642794 utuhalf 9 1 10 01 10 nocpt

sh fmt_misfit_mini.gmt 20110823072700783 utuhalf 9 0 10 01 10
sh fmt_misfit_mini.gmt 20110823072700783 utuhalf 9 1 10 01 10 nocpt

# move results to manuscript dir
evid=20110528061823523; cp ${dir0}/job_${evid}_P00_V10_R01_S10/${evid}/fmt/${evid}_utuhalf_009_P00_V10_R01_S10_misfit_mini.eps ${dirf}/${evid}_utuhalf_009_P00_V10_R01_S10_misfit_mini.eps
evid=20110528061823523; cp ${dir0}/job_${evid}_P01_V10_R01_S10/${evid}/fmt/${evid}_utuhalf_009_P01_V10_R01_S10_misfit_mini.eps ${dirf}/${evid}_utuhalf_009_P01_V10_R01_S10_misfit_mini.eps
evid=20110708145821435; cp ${dir0}/job_${evid}_P00_V10_R01_S10/${evid}/fmt/${evid}_utuhalf_009_P00_V10_R01_S10_misfit_mini.eps ${dirf}/${evid}_utuhalf_009_P00_V10_R01_S10_misfit_mini.eps
evid=20110708145821435; cp ${dir0}/job_${evid}_P01_V10_R01_S10/${evid}/fmt/${evid}_utuhalf_009_P01_V10_R01_S10_misfit_mini.eps ${dirf}/${evid}_utuhalf_009_P01_V10_R01_S10_misfit_mini.eps
evid=20110823072642794; cp ${dir0}/job_${evid}_P00_V10_R01_S10/${evid}/fmt/${evid}_utuhalf_009_P00_V10_R01_S10_misfit_mini.eps ${dirf}/${evid}_utuhalf_009_P00_V10_R01_S10_misfit_mini.eps
evid=20110823072642794; cp ${dir0}/job_${evid}_P01_V10_R01_S10/${evid}/fmt/${evid}_utuhalf_009_P01_V10_R01_S10_misfit_mini.eps ${dirf}/${evid}_utuhalf_009_P01_V10_R01_S10_misfit_mini.eps
evid=20110823072700783; cp ${dir0}/job_${evid}_P00_V10_R01_S10/${evid}/fmt/${evid}_utuhalf_009_P00_V10_R01_S10_misfit_mini.eps ${dirf}/${evid}_utuhalf_009_P00_V10_R01_S10_misfit_mini.eps
evid=20110823072700783; cp ${dir0}/job_${evid}_P01_V10_R01_S10/${evid}/fmt/${evid}_utuhalf_009_P01_V10_R01_S10_misfit_mini.eps ${dirf}/${evid}_utuhalf_009_P01_V10_R01_S10_misfit_mini.eps


