#!/bin/bash
#
# Commands to generate figures 6, and Supplement figures S14 -- S19.
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

# pick figures to generate

# set data directories
#from_dir="/home/alvizuri/projects/source_inversion_new/analysis/dw_20150604_results_fmtpaper_test_components_utuhalf_004"
#to_dir="/home/alvizuri/manuscripts/2015/fmt_uturuncu/data/20100516063454464/fmt_misfit"
evid_main="20100516063454464"

# NOTE set path in fmt_misfit.gmt to 
# datadir_main="/home/alvizuri/shared/plutons/inversions/inv20151024_mainevent_tests" 

sh fmt_misfit.gmt ${evid_main} utuhalf 4 1 10 1 10     # F6, main event
sh fmt_misfit.gmt ${evid_main} utuhalf 4 0 10 1 10     # S14
sh fmt_misfit.gmt ${evid_main} utuhalf 4 0 10 1 0      # S15
sh fmt_misfit.gmt ${evid_main} utuhalf 4 1 10 1 0      # S16
sh fmt_misfit.gmt ${evid_main} utuhalf 4 1 10 0 0      # S17
sh fmt_misfit.gmt ${evid_main} utuhalf 4 1  0 1 0      # S18
sh fmt_misfit.gmt ${evid_main} utuhalf 4 1  0 0 10     # S19

# move results to manuscript dir
dir0="/home/alvizuri/shared/plutons/inversions/inv20151024_mainevent_tests"
dirf="/home/alvizuri/manuscripts/2015/fmt_uturuncu/figures/mainevent"
cp ${dir0}/job_${evid_main}_P01_V10_R01_S10/${evid_main}/fmt/${evid_main}_utuhalf_004_misfit.eps ${dirf}/${evid_main}_utuhalf_004_P01_V10_R01_S10_L_dt_misfit.eps
cp ${dir0}/job_${evid_main}_P00_V10_R01_S10/${evid_main}/fmt/${evid_main}_utuhalf_004_misfit.eps ${dirf}/${evid_main}_utuhalf_004_P00_V10_R01_S10_L_dt_misfit.eps
cp ${dir0}/job_${evid_main}_P00_V10_R01_S00/${evid_main}/fmt/${evid_main}_utuhalf_004_misfit.eps ${dirf}/${evid_main}_utuhalf_004_P00_V10_R01_S00_L_dt_misfit.eps
cp ${dir0}/job_${evid_main}_P01_V10_R01_S00/${evid_main}/fmt/${evid_main}_utuhalf_004_misfit.eps ${dirf}/${evid_main}_utuhalf_004_P01_V10_R01_S00_L_dt_misfit.eps
cp ${dir0}/job_${evid_main}_P01_V10_R00_S00/${evid_main}/fmt/${evid_main}_utuhalf_004_misfit.eps ${dirf}/${evid_main}_utuhalf_004_P01_V10_R00_S00_L_dt_misfit.eps
cp ${dir0}/job_${evid_main}_P01_V00_R01_S00/${evid_main}/fmt/${evid_main}_utuhalf_004_misfit.eps ${dirf}/${evid_main}_utuhalf_004_P01_V00_R01_S00_L_dt_misfit.eps
cp ${dir0}/job_${evid_main}_P01_V00_R00_S10/${evid_main}/fmt/${evid_main}_utuhalf_004_misfit.eps ${dirf}/${evid_main}_utuhalf_004_P01_V00_R00_S10_L_dt_misfit.eps

