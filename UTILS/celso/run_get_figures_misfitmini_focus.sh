#!/bin/bash
#
# Commands to generate figures 12, 13, 14.
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
# datadir_main="/home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V10_R01_S10_L_dt"

# figures 12, 13, 14
#dir0="dw_20150618_utuhalf_P01_V10_R01_S10_L_dt"
dir0="/home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V10_R01_S10_L_dt"
dirf="/home/alvizuri/manuscripts/2015/fmt_uturuncu/figures"

sh fmt_misfit_mini.gmt 20100930095232848 utuhalf  9
sh fmt_misfit_mini.gmt 20100930095413030 utuhalf  9
sh fmt_misfit_mini.gmt 20100930095612039 utuhalf  9
sh fmt_misfit_mini.gmt 20100930102501621 utuhalf 10
sh fmt_misfit_mini.gmt 20100930114531523 utuhalf  9
sh fmt_misfit_mini.gmt 20100930155702647 utuhalf  9
sh fmt_misfit_mini.gmt 20100930163913172 utuhalf  9
sh fmt_misfit_mini.gmt 20100930172333129 utuhalf  9
sh fmt_misfit_mini.gmt 20100930183517772 utuhalf  9
sh fmt_misfit_mini.gmt 20100930215912435 utuhalf 10
sh fmt_misfit_mini.gmt 20101001111055009 utuhalf  9
sh fmt_misfit_mini.gmt 20101001111242832 utuhalf 10
sh fmt_misfit_mini.gmt 20101005032124419 utuhalf  9
sh fmt_misfit_mini.gmt 20111028072530002 utuhalf 10
sh fmt_misfit_mini.gmt 20110425035005648 utuhalf 08
sh fmt_misfit_mini.gmt 20110429191551352 utuhalf  9
sh fmt_misfit_mini.gmt 20111222231621860 utuhalf  9
sh fmt_misfit_mini.gmt 20111224064311431 utuhalf  9
sh fmt_misfit_mini.gmt 20110708145821435 utuhalf  9
sh fmt_misfit_mini.gmt 20110430155647801 utuhalf  9
sh fmt_misfit_mini.gmt 20110809032109852 utuhalf 10
sh fmt_misfit_mini.gmt 20110823072642794 utuhalf  9
sh fmt_misfit_mini.gmt 20111224035332713 utuhalf 10
sh fmt_misfit_mini.gmt 20111224064906473 utuhalf  9
sh fmt_misfit_mini.gmt 20111224172941331 utuhalf 10

cp ${dir0}/job_20110430155647801/20110430155647801/fmt/20110430155647801_utuhalf_009_misfit_mini.eps ${dirf}/iso/20110430155647801_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20110809032109852/20110809032109852/fmt/20110809032109852_utuhalf_010_misfit_mini.eps ${dirf}/iso/20110809032109852_utuhalf_010_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20110823072642794/20110823072642794/fmt/20110823072642794_utuhalf_009_misfit_mini.eps ${dirf}/iso/20110823072642794_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20111224035332713/20111224035332713/fmt/20111224035332713_utuhalf_010_misfit_mini.eps ${dirf}/iso/20111224035332713_utuhalf_010_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20111224064906473/20111224064906473/fmt/20111224064906473_utuhalf_009_misfit_mini.eps ${dirf}/iso/20111224064906473_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20111224172941331/20111224172941331/fmt/20111224172941331_utuhalf_010_misfit_mini.eps ${dirf}/iso/20111224172941331_utuhalf_010_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20110425035005648/20110425035005648/fmt/20110425035005648_utuhalf_008_misfit_mini.eps ${dirf}/cracks/20110425035005648_utuhalf_008_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20110429191551352/20110429191551352/fmt/20110429191551352_utuhalf_009_misfit_mini.eps ${dirf}/cracks/20110429191551352_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20111222231621860/20111222231621860/fmt/20111222231621860_utuhalf_009_misfit_mini.eps ${dirf}/cracks/20111222231621860_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20111224064311431/20111224064311431/fmt/20111224064311431_utuhalf_009_misfit_mini.eps ${dirf}/cracks/20111224064311431_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20110708145821435/20110708145821435/fmt/20110708145821435_utuhalf_009_misfit_mini.eps ${dirf}/cracks/20110708145821435_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20100930095232848/20100930095232848/fmt/20100930095232848_utuhalf_009_misfit_mini.eps ${dirf}/cluster/20100930095232848_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20100930095413030/20100930095413030/fmt/20100930095413030_utuhalf_009_misfit_mini.eps ${dirf}/cluster/20100930095413030_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20100930095612039/20100930095612039/fmt/20100930095612039_utuhalf_009_misfit_mini.eps ${dirf}/cluster/20100930095612039_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20100930102501621/20100930102501621/fmt/20100930102501621_utuhalf_010_misfit_mini.eps ${dirf}/cluster/20100930102501621_utuhalf_010_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20100930114531523/20100930114531523/fmt/20100930114531523_utuhalf_009_misfit_mini.eps ${dirf}/cluster/20100930114531523_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20100930155702647/20100930155702647/fmt/20100930155702647_utuhalf_009_misfit_mini.eps ${dirf}/cluster/20100930155702647_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20100930163913172/20100930163913172/fmt/20100930163913172_utuhalf_009_misfit_mini.eps ${dirf}/cluster/20100930163913172_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20100930172333129/20100930172333129/fmt/20100930172333129_utuhalf_009_misfit_mini.eps ${dirf}/cluster/20100930172333129_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20100930183517772/20100930183517772/fmt/20100930183517772_utuhalf_009_misfit_mini.eps ${dirf}/cluster/20100930183517772_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20100930215912435/20100930215912435/fmt/20100930215912435_utuhalf_010_misfit_mini.eps ${dirf}/cluster/20100930215912435_utuhalf_010_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20101001111055009/20101001111055009/fmt/20101001111055009_utuhalf_009_misfit_mini.eps ${dirf}/cluster/20101001111055009_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20101001111242832/20101001111242832/fmt/20101001111242832_utuhalf_010_misfit_mini.eps ${dirf}/cluster/20101001111242832_utuhalf_010_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20101005032124419/20101005032124419/fmt/20101005032124419_utuhalf_009_misfit_mini.eps ${dirf}/cluster/20101005032124419_utuhalf_009_P01_V10_R01_S10_L_dt_misfit_mini.eps
cp ${dir0}/job_20111028072530002/20111028072530002/fmt/20111028072530002_utuhalf_010_misfit_mini.eps ${dirf}/cluster/20111028072530002_utuhalf_010_P01_V10_R01_S10_L_dt_misfit_mini.eps

