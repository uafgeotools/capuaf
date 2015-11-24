#!/bin/bash
#
# Commands to generate waveform fit plots for the main event.
# Figures 4 and Supplement figures S7--S13, S20, S22.
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
# 20151118 celso alvizuri - cralvizuri@alaska.edu 
#-----------------------------------------------------------


# main event figure 4 -- normalized amplitudes
# Event 20100516063454464 Model utuhalf_004 FM 150 76.443817 -75 Mw 2.80 rms 1.790e-07 2013 ERR 0 0 0 ISO -9.594069 0.00 CLVD 10.00 0.00 VR 22.1 data2 2.027e-07
# NOTE uncomment line 254 cap_plt.pl : #$stams = $stamb;
cap_wf.pl -P0.5e+0.5/1/30 -W1 -E1 -K1 -Y1 -I5/0/5 -J-9.594069/-9.594069/10/10 -R150/150/76.443817/76.443817/-75/-75 -D1/1/0.5 -H0.01 -C2/10/0.25/0.5 -S0.05/2.0/0 -T1.5/60 -L0.10 -F0 -Mutuhalf_4/2.8 -Zw_utuhalf_P01_V10_R01_S10.dat 20100516063454464 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S10_fmt.ps

# main event figure S7 -- abs amplitudes
cap_wf.pl -P7e-5/1/30 -p0.2 -W1 -E1 -K1 -Y1 -I5/0/5 -J-9.594069/-9.594069/10/10 -R150/150/76.443817/76.443817/-75/-75 -D1/1/0.5 -H0.01 -C2/10/0.25/0.5 -S0.05/2.0/0 -T1.5/60 -L0.10 -F0 -Mutuhalf_4/2.8 -Zw_utuhalf_P01_V10_R01_S10.dat 20100516063454464 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S10_fmt_amp_abs.ps

# main event figure S8 -- P00 V10 R01 S10
# Event 20100516063454464 Model utuhalf_004 FM  145 72.962448  -70 Mw 2.80 rms 1.786e-07  2013 ERR   0   0   0 ISO -12.839588 0.00 CLVD 5.00 0.00 VR 22.4 data2 2.027e-07
cap_wf.pl -P7e-5/1/30 -p0.2 -W1 -E1 -K1 -Y1 -I5/0/5 -J-12.839588/-12.839588/5/5 -R145/145/72.962448/72.962448/-70/-70 -D1/1/0.5 -H0.01 -C2/10/0.25/0.5 -S0.05/2.0/0 -T1.5/60 -L0.10 -F0 -Mutuhalf_4/2.8 -Zw_utuhalf_P00_V10_R01_S10.dat 20100516063454464 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P00_V10_R01_S10_fmt.ps

# main event figure S9 -- P00 V10 R01 S00
# Event 20100516063454464 Model utuhalf_004 FM  130 58.170238  -80 Mw 3.00 rms 1.541e-07    33 ERR   0   0   0 ISO -26.387798 0.00 CLVD 20.00 0.00 VR 32.8 data2 1.879e-07
cap_wf.pl -P7e-5/1/30 -p0.2 -W1 -E1 -K1 -Y1 -I5/0/5 -J-26.387798/-26.387798/20/20 -R130/130/58.170238/58.170238/-80/-80 -D1/1/0.5 -H0.01 -C2/10/0.25/0.5 -S0.05/2.0/0 -T1.5/60 -L0.10 -F0 -Mutuhalf_4/3.0 -Zw_utuhalf_P00_V10_R01_S00.dat 20100516063454464 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P00_V10_R01_S00_fmt.ps

# main event figure S10 -- P01 V10 R01 S00
# Event 20100516063454464 Model utuhalf_004 FM  145 65.782738  -80 Mw 2.90 rms 1.562e-07    33 ERR   0   0   0 ISO -22.885382 0.00 CLVD 20.00 0.00 VR 30.9 data2 1.879e-07
cap_wf.pl -P7e-5/1/30 -p0.2 -W1 -E1 -K1 -Y1 -I5/0/5 -J-22.885382/-22.885382/20/20 -R145/145/65.782738/65.782738/-80/-80 -D1/1/0.5 -H0.01 -C2/10/0.25/0.5 -S0.05/2.0/0 -T1.5/60 -L0.10 -F0 -Mutuhalf_4/2.9 -Zw_utuhalf_P01_V10_R01_S00.dat 20100516063454464 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00_fmt.ps

# main event figure S11 -- P01 V10 R00 S00
# Event 20100516063454464 Model utuhalf_004 FM  145 65.782738  -80 Mw 3.00 rms 1.404e-07    16 ERR   0   0   0 ISO -22.885382 0.00 CLVD 20.00 0.00 VR 47.6 data2 1.940e-07
cap_wf.pl -P7e-5/1/30 -p0.2 -W1 -E1 -K1 -Y1 -I5/0/5 -J-22.885382/-22.885382/20/20 -R145/145/65.782738/65.782738/-80/-80 -D1/1/0.5 -H0.01 -C2/10/0.25/0.5 -S0.05/2.0/0 -T1.5/60 -L0.10 -F0 -Mutuhalf_4/3.0 -Zw_utuhalf_P01_V10_R00_S00.dat 20100516063454464 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R00_S00_fmt.ps

# main event figure S12 -- P01 V00 R01 S00
# Event 20100516063454464 Model utuhalf_004 FM   10  5.000000  -75 Mw 2.20 rms 1.170e-07    16 ERR   0   0   0 ISO   0.000000 0.00 CLVD 0.00 0.00 VR 15.2 data2 1.271e-07
cap_wf.pl -P7e-5/1/30 -p0.2 -W1 -E1 -K1 -Y1 -I5/0/5 -J0/0/0/0 -R10/10/5/5/-75/-75 -D1/1/0.5 -H0.01 -C2/10/0.25/0.5 -S0.05/2.0/0 -T1.5/60 -L0.10 -F0 -Mutuhalf_4/2.2 -Zw_utuhalf_P01_V00_R01_S00.dat 20100516063454464 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V00_R01_S00_fmt.ps

# main event figure S13 -- P01 V00 R00 S10
# Event 20100516063454464 Model utuhalf_004 FM  150 72.962448  -55 Mw 2.80 rms 1.842e-07  1980 ERR   0   0   0 ISO  -3.184737 0.00 CLVD 15.00 0.00 VR 21.7 data2 2.082e-07
cap_wf.pl -P7e-5/1/30 -p0.2 -W1 -E1 -K1 -Y1 -I5/0/5 -J-3.184737/-3.184737/15/15 -R150/150/72.962448/72.962448/-55/-55 -D1/1/0.5 -H0.01 -C2/10/0.25/0.5 -S0.05/2.0/0 -T1.5/60 -L0.10 -F0 -Mutuhalf_4/2.8 -Zw_utuhalf_P01_V00_R01_S00.dat 20100516063454464 
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V00_R00_S10_fmt.ps

# main event figure S20 -- DC, P01 V10 R01 S00
# Event 20100516063454464 Model utuhalf_004 FM  150 72.962448  -55 Mw 2.80 rms 1.842e-07  1980 ERR   0   0   0 ISO  -3.184737 0.00 CLVD 15.00 0.00 VR 21.7 data2 2.082e-07
cap_wf.pl -P7e-5/1/30 -p0.2 -W1 -E1 -K1 -Y1 -I5/0/0 -D1/1/0.5 -H0.01 -C2/10/0.25/0.5 -S0.05/2.0/0 -T1.5/60 -L0.10 -F0 -Mutuhalf_4/2.7 -Zw_utuhalf_P01_V10_R01_S00.dat 20100516063454464
cp 20100516063454464/utuhalf_004.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00_dc_fmt.ps

# main event figure S22 -- CO, P01 V10 R01 S00
# Event 20100516063454464 Model utuhalf_004 FM  345 34.875584  -90 Mw 3.10 rms 1.636e-07    33 ERR   0   0   0 ISO -30.000000 0.00 CLVD 30.00 0.00 VR 24.2 data2 1.879e-07
cap_wf.pl -P7e-5/1/30 -p0.2 -W1 -E1 -K1 -Y1 -I5/0/5 -J-30/-30/30/30 -R345/345/34.875584/34.875584/-90/-90 -D1/1/0.5 -H0.01 -C2/10/0.25/0.5 -S0.05/2.0/0 -T1.5/60 -L0.10 -F0 -Mutuhalf_4/3.1 -Zw_utuhalf_P01_V10_R01_S00.dat 20100516063454464
cp 20100516063454464/utuhalf_004_fmt.ps 20100516063454464_check/20100516063454464_utuhalf_004_P01_V10_R01_S00_cc_fmt.ps



