#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
This code takes in a file name with P arrivals and makes a different weight file with the P arrivals.  The new weight file can be a body weight file, surface wieight file or using everything.
Kyle Smith
Made January 19, 2018 
"""

import os, sys
import obspy
from read_cap_weight import *
import math
import string

# name of file that contains the P arrivals
file_w_Parrivals = "/home/ksmith/REPOSITORIES/capuaf/MTs/Dec10MT/V1/20171210122855089/weight_body_150.dat"

# name of the output file
outfile = "/home/ksmith/REPOSITORIES/capuaf/MTs/Dec10MT/V1/20171210122855089/weight_body_150_sc.dat"

# read in weight file
stnm,edist,PV_wt,PR_wt,SV_wt,SR_wt,ST_wt,P_arrival,P_len,S_arrival,S_len,waveform_shft = read_cap_weight(file_w_Parrivals)

# make new file
f = open(outfile, "w")
stfmt = '%34s %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.4f %4.0f %4.0f %4.0f %4.0f \n';
scalewts = 100 
filenumber = 2

# make copy
if filenumber == 1:
    for linenum in range(0,len(stnm)):
        xx = stfmt % (stnm[linenum],edist[linenum],PV_wt[linenum],PR_wt[linenum],SV_wt[linenum],SR_wt[linenum],ST_wt[linenum],P_arrival[linenum],P_len[linenum], S_arrival[linenum],S_len[linenum],waveform_shft[linenum])
        f.write(xx)

# scale file
if filenumber == 2:
    for linenum in range(0,len(stnm)):
        xx = stfmt % (stnm[linenum],edist[linenum],PV_wt[linenum]*scalewts,PR_wt[linenum]*scalewts,SV_wt[linenum]*scalewts,SR_wt[linenum]*scalewts,ST_wt[linenum]*scalewts,P_arrival[linenum],P_len[linenum], S_arrival[linenum],S_len[linenum],waveform_shft[linenum])
        f.write(xx)

# make standard weight file (all 1s)
if filenumber == 3:
    for linenum in range(0,len(stnm)):
        xx = stfmt % (stnm[linenum],edist[linenum],scalewts,scalewts,scalewts,scalewts,scalewts,P_arrival[linenum],P_len[linenum], S_arrival[linenum],S_len[linenum],waveform_shft[linenum])
        f.write(xx)

# make body waves weight file
if filenumber == 4:
    for linenum in range(0,len(stnm)):
        xx = stfmt % (stnm[linenum],edist[linenum],scalewts,scalewts,0,0,0,P_arrival[linenum],P_len[linenum], S_arrival[linenum],S_len[linenum],waveform_shft[linenum])
        f.write(xx)

# make surface waves weight file
if filenumber == 5:
    for linenum in range(0,len(stnm)):
        xx = stfmt % (stnm[linenum],edist[linenum],0,0,scalewts,scalewts,scalewts,P_arrival[linenum],P_len[linenum], S_arrival[linenum],S_len[linenum],waveform_shft[linenum])
        f.write(xx)

f.close()
