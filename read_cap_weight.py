#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to read information from weight files in cap
"""

def read_cap_weight(filename):
    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    count = 0
    stnm = [];edist = [];PV_wt = [];PR_wt = [];SV_wt = [];SR_wt = [];ST_wt = []; P_arrival = []; P_len =[]; S_arrival = []; S_len = []; waveform_shft = [];

    for line in lines:
        line_elements = line.split()
        stnm.append(line_elements[0])
        edist.append(float(line_elements[1]))
        PV_wt.append(float(line_elements[2]))
        PR_wt.append(float(line_elements[3]))
        SV_wt.append(float(line_elements[4]))
        SR_wt.append(float(line_elements[5]))
        ST_wt.append(float(line_elements[6]))
        P_arrival.append(float(line_elements[7]))
        P_len.append(float(line_elements[8]))
        S_arrival.append(float(line_elements[9]))
        S_len.append(float(line_elements[10]))
        waveform_shft.append(float(line_elements[11]))
        
    return stnm,edist,PV_wt,PR_wt,SV_wt,SR_wt,ST_wt,P_arrival,P_len,S_arrival,S_len,waveform_shft
    
