#!/bin/bash
# 
# generate histograms of amplitude ratios
# 
# 20151115 celso alvizuri - cralvizuri@alaska.edu 
#-----------------------------------------------------------



# histograms of amplitude ratios
# figure 3a, ln(Vobs/Vsyn), ln(Robs, Rsyn)
# NOTE this will generate histogram and scatter plot
# uturuncu fmt, weights 1-1-1
python get_stats_ampratios_composite_v2.py /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V01_R01_S01_L_dt/utuhalf_P01_V01_R01_S01_L_dt_amps_body_v2v_all /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V01_R01_S01_L_dt/utuhalf_P01_V01_R01_S01_L_dt_amps_body_r2r_all a2a

#mv amp_ratios_plot_histogram.eps
#mv amp_ratios_plot_scatter.eps

# figure 3b, ln(Vobs/Vsyn), ln(Robs, Rsyn)
# NOTE this will generate histogram and scatter plot
# uturuncu fmt, weights 10-1-10
python get_stats_ampratios_composite_v2.py /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V10_R01_S10_L_dt/utuhalf_P01_V10_R01_S10_L_dt_amps_body_v2v_all /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V10_R01_S10_L_dt/utuhalf_P01_V10_R01_S10_L_dt_amps_body_r2r_all a2a


# supplement figure S3
python get_stats_ampratios_composite_v2.py /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V01_R01_S01_L_dt/utuhalf_P01_V01_R01_S01_L_dt_amps_body_v2v_all /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V01_R01_S01_L_dt/utuhalf_P01_V01_R01_S01_L_dt_amps_body_r2r_all v2r

#-----------------------------------------------------------
# histogram+scatter comparing influence of first motion polarities on moment tensor inversions

# figure 9
# be sure to 
python utu63_scatter_histo.py /home/alvizuri/REPOSITORIES/GEOTOOLS/shared/alvizuri/papers/2014fmt/data/theta_dvr_eid_utuhalf

# figure S2a
python utu63_scatter_histo.py /home/alvizuri/REPOSITORIES/GEOTOOLS/shared/alvizuri/papers/2014fmt/data/theta_dvr_eid_utu1d

# figure S2b (adjust python limits)
python utu63_scatter_histo.py /home/alvizuri/REPOSITORIES/GEOTOOLS/shared/alvizuri/papers/2014fmt/data/theta_dvr_eid_utuhalf_utu1d_pol

#-----------------------------------------------------------

# figure S4
sh plot_map_amplitude_ratios.gmt /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V01_R01_S01_L_dt/utuhalf_P01_V01_R01_S01_L_dt_amps_body_v2v_stats
sh plot_map_amplitude_ratios.gmt /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V01_R01_S01_L_dt/utuhalf_P01_V01_R01_S01_L_dt_amps_body_r2r_stats
sh plot_map_amplitude_ratios.gmt /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V10_R01_S10_L_dt/utuhalf_P01_V10_R01_S10_L_dt_amps_body_v2v_stats
sh plot_map_amplitude_ratios.gmt /home/alvizuri/shared/plutons/inversions/inv20150618_utuhalf_P01_V10_R01_S10_L_dt/utuhalf_P01_V10_R01_S10_L_dt_amps_body_r2r_stats



