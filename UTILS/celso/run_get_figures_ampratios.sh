#!/bin/bash
#
# Commands to generate histograms and scatter plots. Figures 3 and 9, and
# Supplement figures S2, S3, S4. 
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
# 20151115 celso alvizuri - cralvizuri@alaska.edu 
#-----------------------------------------------------------

# histograms of amplitude ratios
# figure 3a, ln(Vobs/Vsyn), ln(Robs, Rsyn)
# NOTE this will generate histogram and scatter plot
# uturuncu fmt, weights 1-1-1
python get_stats_ampratios_composite_v2.py /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V01_R01_S01_L_dt/utuhalf_P01_V01_R01_S01_L_dt_amps_body_v2v_all /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V01_R01_S01_L_dt/utuhalf_P01_V01_R01_S01_L_dt_amps_body_r2r_all a2a

mv amp_ratios_plot_histogram.eps /home/alvizuri/manuscripts/2015/fmt_uturuncu/figures/amp_ratios/amp_ratios_histo_utuhalf_P01_V01_R01_S01_L_dt.eps

# figure 3b, ln(Vobs/Vsyn), ln(Robs, Rsyn)
# NOTE this will generate histogram and scatter plot
# uturuncu fmt, weights 10-1-10
python get_stats_ampratios_composite_v2.py /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V10_R01_S10_L_dt/utuhalf_P01_V10_R01_S10_L_dt_amps_body_v2v_all /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V10_R01_S10_L_dt/utuhalf_P01_V10_R01_S10_L_dt_amps_body_r2r_all a2a

mv amp_ratios_plot_histogram.eps /home/alvizuri/manuscripts/2015/fmt_uturuncu/figures/amp_ratios/amp_ratios_histo_utuhalf_P01_V10_R01_S10_L_dt.eps

# supplement figure S3
python get_stats_ampratios_composite_v2.py /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V01_R01_S01_L_dt/utuhalf_P01_V01_R01_S01_L_dt_amps_body_v2v_all /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V01_R01_S01_L_dt/utuhalf_P01_V01_R01_S01_L_dt_amps_body_r2r_all v2r

mv amp_ratios_plot_histogram.eps  /home/alvizuri/manuscripts/2015/fmt_uturuncu/figures/amp_ratios/amp_ratios_histo_v2r_utuhalf_P01_V01_R01_S01_L_dt.eps
mv amp_ratios_plot_scatter.eps  /home/alvizuri/manuscripts/2015/fmt_uturuncu/figures/amp_ratios/amp_ratios_scatter_v2r_utuhalf_P01_V01_R01_S01_L_dt.eps

#-----------------------------------------------------------
# histogram+scatter comparing influence of first motion polarities on moment tensor inversions

# figure 9
# NOTE adjust plot limits in the python script: Lines 51, 90
python utu63_scatter_histo.py /home/alvizuri/REPOSITORIES/manuscripts/alvizuri/papers/2014fmt/data/theta_dvr_eid_utuhalf
mv output.eps /home/alvizuri/manuscripts/2015/fmt_uturuncu/figures/scatter_histo_theta_dvr_utuhalf_V10_R01_S10.eps
# figure S2a
# NOTE adjust plot limits in the python script: Lines 51, 90
python utu63_scatter_histo.py /home/alvizuri/REPOSITORIES/manuscripts/alvizuri/papers/2014fmt/data/theta_dvr_eid_utu1d
mv output.eps /home/alvizuri/manuscripts/2015/fmt_uturuncu/figures/scatter_histo_theta_dvr_utu1d_V10_R01_S10.eps

exit
# figure S2b 
# NOTE adjust plot limits in the python script: Lines 52, 95
python utu63_scatter_histo.py /home/alvizuri/REPOSITORIES/manuscripts/alvizuri/papers/2014fmt/data/theta_dvr_eid_utuhalf_utu1d_pol
mv output.eps /home/alvizuri/manuscripts/2015/fmt_uturuncu/figures/scatter_histo_theta_dvr_utuhalf_utu1d_P01_V10_R01_S10.eps


#-----------------------------------------------------------

# figure S4
sh plot_map_amplitude_ratios.gmt /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V01_R01_S01_L_dt/utuhalf_P01_V01_R01_S01_L_dt_amps_body_v2v_stats
sh plot_map_amplitude_ratios.gmt /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V01_R01_S01_L_dt/utuhalf_P01_V01_R01_S01_L_dt_amps_body_r2r_stats
sh plot_map_amplitude_ratios.gmt /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V10_R01_S10_L_dt/utuhalf_P01_V10_R01_S10_L_dt_amps_body_v2v_stats
sh plot_map_amplitude_ratios.gmt /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V10_R01_S10_L_dt/utuhalf_P01_V10_R01_S10_L_dt_amps_body_r2r_stats



