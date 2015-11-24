#!/bin/bash
#
# Commands to generate figures 7, 8, and Supplement figure S23.
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

# NOTE set path in fmt_misfit_mini.gmt to 
# datadir_main="/home/alvizuri/shared/plutons/inversions/inv20151024_mainevent_tests" 
exit    # check path above then comment this line

evid_main="20100516063454464"

dir0="/home/alvizuri/shared/plutons/inversions/inv20151024_mainevent_tests"
dirf="/home/alvizuri/manuscripts/2015/fmt_uturuncu/figures/mainevent"

#-----------------------------------------------------------
# figure 7
# 1. get cpt files
#    run main event NO POLARITIES to get cpt files
sh fmt_misfit_mini.gmt 20100516063454464 utuhalf 4 0 10 1 10 
sh fmt_misfit_mini.gmt 20100516063454464 utuhalf 4 1 10 1 10 nocpt
sh fmt_misfit_mini.gmt 20100516063454464 utuhalf 4 0 10 1 0
sh fmt_misfit_mini.gmt 20100516063454464 utuhalf 4 1 10 1 0 nocpt

# move results to manuscript dir
cp ${dir0}/job_${evid_main}_P00_V10_R01_S10/${evid_main}/fmt/${evid_main}_utuhalf_004_P00_V10_R01_S10_misfit_mini.eps ${dirf}/${evid_main}_utuhalf_004_P00_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_${evid_main}_P01_V10_R01_S10/${evid_main}/fmt/${evid_main}_utuhalf_004_P01_V10_R01_S10_misfit_mini.eps ${dirf}/${evid_main}_utuhalf_004_P01_V10_R01_S10_L_dt_misfit_mini_const_cpt.eps
cp ${dir0}/job_${evid_main}_P00_V10_R01_S00/${evid_main}/fmt/${evid_main}_utuhalf_004_P00_V10_R01_S00_misfit_mini.eps ${dirf}/${evid_main}_utuhalf_004_P00_V10_R01_S00_L_dt_misfit_mini.eps
cp ${dir0}/job_${evid_main}_P01_V10_R01_S00/${evid_main}/fmt/${evid_main}_utuhalf_004_P01_V10_R01_S00_misfit_mini.eps ${dirf}/${evid_main}_utuhalf_004_P01_V10_R01_S00_L_dt_misfit_mini.eps


#-----------------------------------------------------------
# figure 8
# each generates its own cpt
sh fmt_misfit_mini.gmt 20100516063454464 utuhalf 4 1 10 0  0
sh fmt_misfit_mini.gmt 20100516063454464 utuhalf 4 1  0 1  0
sh fmt_misfit_mini.gmt 20100516063454464 utuhalf 4 1  0 0 10
sh fmt_misfit_mini.gmt 20100516063454464 utuhalf 4 1 10 1 10

# move results to manuscript dir
cp ${dir0}/job_${evid_main}_P01_V00_R00_S10/${evid_main}/fmt/${evid_main}_utuhalf_004_P01_V00_R00_S10_misfit_mini.eps ${dirf}/${evid_main}_utuhalf_004_P01_V00_R00_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_${evid_main}_P01_V00_R01_S00/${evid_main}/fmt/${evid_main}_utuhalf_004_P01_V00_R01_S00_misfit_mini.eps ${dirf}/${evid_main}_utuhalf_004_P01_V00_R01_S00_L_dt_misfit_mini.eps
cp ${dir0}/job_${evid_main}_P01_V10_R00_S00/${evid_main}/fmt/${evid_main}_utuhalf_004_P01_V10_R00_S00_misfit_mini.eps ${dirf}/${evid_main}_utuhalf_004_P01_V10_R00_S00_L_dt_misfit_mini.eps
cp ${dir0}/job_${evid_main}_P01_V10_R01_S10/${evid_main}/fmt/${evid_main}_utuhalf_004_P01_V10_R01_S10_misfit_mini.eps ${dirf}/${evid_main}_utuhalf_004_P01_V10_R01_S10_L_dt_misfit_mini.eps

#-----------------------------------------------------------
# figure S23
# get cpt from results 1 10-1-0

# get the other results
sh fmt_misfit_mini.gmt ${evid_main} utuhalf 4 1 10 1 0
sh fmt_misfit_mini.gmt ${evid_main} utuhalf 4 1 10 1 0 dc
sh fmt_misfit_mini.gmt ${evid_main} utuhalf 4 1 10 1 0 co
sh fmt_misfit_mini.gmt ${evid_main} utuhalf 4 1 10 1 0 cc

# move results to manuscript dir
cp ${dir0}/job_${evid_main}_P01_V10_R01_S00/${evid_main}/fmt/${evid_main}_utuhalf_004_P01_V10_R01_S00_misfit_mini.eps ${dirf}/${evid_main}_utuhalf_004_P01_V10_R01_S00_L_dt_misfit_mini.eps
cp ${dir0}/job_${evid_main}_P01_V10_R01_S00_cc/${evid_main}/fmt/${evid_main}_utuhalf_004_P01_V10_R01_S00_misfit_mini.eps ${dirf}/${evid_main}_utuhalf_004_P01_V10_R01_S00_L_dt_cc_misfit_mini.eps
cp ${dir0}/job_${evid_main}_P01_V10_R01_S00_co/${evid_main}/fmt/${evid_main}_utuhalf_004_P01_V10_R01_S00_misfit_mini.eps ${dirf}/${evid_main}_utuhalf_004_P01_V10_R01_S00_L_dt_co_misfit_mini.eps
cp ${dir0}/job_${evid_main}_P01_V10_R01_S00_dc/${evid_main}/fmt/${evid_main}_utuhalf_004_P01_V10_R01_S00_misfit_mini.eps ${dirf}/${evid_main}_utuhalf_004_P01_V10_R01_S00_L_dt_dc_misfit_mini.eps


