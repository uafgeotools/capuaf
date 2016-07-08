#!/usr/bin/env python
'''
See B. Savage script for reference and for reading RESP values.

This script converts response formats from RESP to Antelope for use with the LLNL database.
This script checks for some constants and outputs a response file for each sensor epoch.

20160504 celso alvizuri - cralvizuri@alaska.edu 
'''

import sys
import os
from math import *
from numpy import pi

if len(sys.argv) < 2 :
    print "Usage: " + os.path.basename(sys.argv[0]) + " RESP-file"
    print "      Converts a RESP file into an Antelope \"paz\" file"
    exit(-1)

file = sys.argv[1]

fp = open( file )

hash = { 
    "B050F03":    { "name": "station",  'id': 2 },
    "B050F16":    { "name": "network",  'id': 2 },
    "B052F04":    { "name": "channel",  'id': 2 },
    "B052F03":    { "name": "location", 'id': 2 },
    "B052F22":    { "name": "start",    'id': 3 },  # test
    "B052F23":    { "name": "end",      'id': 3 },  # test
    "B053F03":    { "name": "type",     "id": 4 },
    "B053F04":    { "name": "stage",    "id": 4 },
    "B053F05":    { "name": "units",    "id": 5 },
    "B053F07":    { "name": "A0",       "id": 4 },
    "B053F08":    { "name": "fn",       "id": 3 },
    "B053F14":    { "name": "npoles",   "id": 4 },
    "B053F09":    { "name": "nzeros",   "id": 4 },
    "B053F15-18": { "name": "poles",    "id": range(2,4) },
    "B053F10-13": { "name": "zeros",    "id": range(2,4) },
    "B058F03":    { "name": "stage",    "id": 4 },
    "B058F04":    { "name": "sd",       "id": 2 },
    "B058F05":    { "name": "fs",       "id": 4 },
    
    }


stages = list()
data = dict()

block = ""

# check A0 values
# function parameters are made explicit for now but they can be cleaned up
def getA0_theo(fn, fs, poles, zeros, npoles, nzeros, A0):
    # Add in the additional zeros, converting the output to displacement in meters
    for i in range(gamma):
        zeros.append(0.0)
        zeros.append(0.0)
    nzeros = nzeros + gamma

    # Save Poles and Zeros in a list of Complex Numbers
    P = list()
    for i in range(0, npoles*2, 2) :
        P.append( complex(poles[i], poles[i+1]) )
    Z = list()
    for i in range(0, nzeros*2, 2) :
        Z.append( complex(zeros[i], zeros[i+1]) )

    # If there are no poles or zeros, set the calculated A0 to defined A0
    if npoles == 0 or nzeros == 0:
        calc_A0 = A0
    else : # Compute A0 from the poles and Zeros
        f0 = complex(0.0, 2.0 * pi * fs)
        denom = f0 - Z[0]
        for i in range(1,len(Z)):    denom = denom * ( f0 - Z[i])
        numer = f0 - P[0]
        for i in range(1,len(P)):    numer = numer * ( f0 - P[i])
        calc_A0 = abs( numer / denom ) * 2 * pi * fs
        #calc_A0 = abs( numer / denom )

    # Handle cases where frequencies do not match
    if abs(fs - fn) > 1e-4 :
        print "Warning: [",station,network,channel,"]",
        print "Sensitivity and Normalization Frequencies not equal"
        print "      Will Calcuate Normalization Constant (rdseed/evalresp default)"
        print "      Sensitivity Frequency   ", fs
        print "      Normalization Frequency ", fn
        A0 = calc_A0

    # Handle cases where the A0 values do not match
    if abs(A0 - calc_A0)/calc_A0 > 0.005 :
        print "Warning: [",station,network,channel,"]",
        print "Calculated and Defined A0 do not match"
        print "      Defined   ", A0, "  from collected RESP file values"
        print "      Calculated", calc_A0, "  from poles and zeros"
        print "      Using Calculated Normalization Constant (rdseed/evalresp default)"
        A0 = calc_A0

    # Find the final normalization constant
    # (doesn't produce the right value. I'm disabling this part)
    #A0 = A0 / A0d
    return A0, calc_A0

# output to file
def pz2file(poles, zeros, npoles, nzeros, A0, epoch):
    # Create a file name out save the polezero file
    #out = "SAC_PZs_" + network + "_" + station + "_" + channel + "_" + location + "_from_" + file
    #out = "antpz_" + network + "_" + station + "_" + channel + "_" + location + "_" + epoch
    out = network + "_" + station + "_" + channel + "_" + epoch
    fp = open(out, 'w')

#    print >> fp, "ZEROS",nzeros
#    for i in range(0,nzeros*2,2) :
#        if zeros[i] != 0.0 or zeros[i+1] != 0.0:
#            print >> fp, "%.4f  %.4f" % ( zeros[i], zeros[i+1] )
#  
#    print >> fp, "POLES",npoles
#    for i in range(0,npoles*2,2) :
#        print >> fp, "%.4f  %.4f" % ( poles[i], poles[i+1] )
#    print >>fp, "CONSTANT %.6e" % ( A0 * sd )

    print >> fp, "# This response file was derived from: %s" % file
    print >> fp, "# Epoch %s" % epoch
    print >> fp, "# Format: Antelope \"paz\""
    print >> fp, "# See man pages for description of the format"
    print >> fp, "# Antelope version 5.4"
    print >> fp, "# "
    print >> fp, "theoretical 1 complete paz LLNL/UAF"
    print >> fp, "%s" % (A0)
    print >> fp, "%s" % (npoles)
    for i in range(0,npoles*2,2) :
        print >> fp, ("%13.6e  %13.6e" % (poles[i], poles[i+1]))
    print >> fp, ("%s" % (nzeros))
    for i in range(0,nzeros*2,2) :
            print >> fp, ("%13.6e  %13.6e" % (zeros[i], zeros[i+1]))

    fp.close()

# process file
for line in fp :
    line = line.replace("\n", "")
    v = line.split()
    key = v[0]
    c = key[0] 
    if c == "#":
        continue
    cur = key[0:4]
    if cur != block:
        if 'block' in data :
            if (data['block'] != "B058") or (data['block'] == "B058" and data['stage'] == 0) :
                stages.append(data.copy())
                data.clear()

    block = cur
    if key in hash :
        data['block'] = block
        name = hash[key]['name']
        id   = hash[key]['id']

        if name == "npoles" or name == "nzeros" :
            data[name] = int(v[id])
        elif name == "poles" or name == "zeros" :
            if not name in data:
                data[name] = list()
            for i in id :
                data[name].append( float(v[ i ]))
        elif name == "A0" :
            data[name] = float(v[id])
        elif name == "sd" or name == "fs" :
            if data['stage'] == 0 :
                data[name] = float(v[id])
        elif name == "stage" :
            data[name] = int(v[id])
        elif name == "fn" :
            if data['stage'] == 1 :
                data[name] = float(v[id])
        elif name == "units" :
            if data['stage'] == 1 :
                data[name] = v[id]
        else :
            data[ name ] = v[ id ]

fp.close()

if (data['block'] != "B058") or (data['block'] == "B058" and data['stage'] == 0) :
    stages.append(data)

# get values
A0     = 1
npoles = 0
nzeros = 0
gamma  = 0
nepoch = 0

start = 0
end = 0

#poles  = list()
#zeros  = list()

#for blocks in stages:
#    print(blocks)
#    print("yes")

pz = {}

# NOTE there are more stages than epochs!
for s in stages:
    if 'start' in s:
        y0 = s['start'].split(",")[0]
        d0 = s['start'].split(",")[1]
        nepoch = nepoch + 1
    if 'end' in s:
        end = s['end']
        y1 = s['end'].split(",")[0]
        d1 = s['end'].split(",")[1]
        epoch = ("%s%s_%s%s" % (y0, d0, y1, d1 ))
        pz[nepoch, 'epoch'] = epoch
#        print("%s%s_%s%s" % (y0, d0, y1, d1 ))
#    print("nepoch %s %s  to  %s" % (nepoch, start, end))

    # Convert units to Displacement by adding additional Zeros
    if 'units' in s:
        if s['units'] == "M/S":      gamma = 1
        if s['units'] == "M/S**2":   gamma = 2
     # Compute the Product of all A0 values (converting to correct units)
    if 'A0' in s :
        if nepoch > 0:
            A0 = 1
        A0 = A0 * s['A0']
        pz[(nepoch, 'A0')] = A0

        if s['type'] == "B" :
            A0 = A0 * (2 * pi)**( s['npoles'] - s['nzeros'] )

    # Compute additional normalization factor         
    if 'fn' in s :
        fn = s['fn']
        pz[nepoch, 'fn'] = fn
        A0d = (2 * pi * fn)**gamma
    else:
        print("Warning. fn not found. setting to null")
        fn = ''
    if 'sd' in s :            sd = s['sd'] * (2 * pi * s['fs'])**gamma
        
    # Save Needed Values
    if 'fs' in s :            
        fs = s['fs']
        pz[nepoch, 'fs'] = fs
    else:
        print("Warning. fs not found. setting to null")
        fs = ''

    if 'station' in s:        station = s['station']
    if 'network' in s:        network = s['network']
    if 'channel' in s:        channel = s['channel']
    if 'location' in s:
        location = s['location']
        if location == "??":
            location = ""

    # Determine how many total poles and zeros exist
    if 'nzeros' in s :        
        nzeros = nzeros + s['nzeros']
    if 'npoles' in s :        
        npoles = npoles + s['npoles']
 
    # Save all poles and zeros (converting to correct units)
    if 'poles' in s:
        poles  = list()
        for i in range(s['npoles'] * 2) :
            if s['type'] == "B" :
                s['poles'][i] = s['poles'][i] * 2 * pi
            poles.append(s['poles'][i])
        pz[nepoch,'p'] = poles, npoles
        npoles = 0
#   else:
#       print("warning. no poles, set dict value to null")
#       pz[nepoch,'p'] = '', 0
#        print("poles: %s" % poles)

    if 'zeros' in s:
        zeros  = list()
        for i in range(s['nzeros'] * 2) :
            if s['type'] == "B" :
                s['zeros'][i] = s['zeros'][i] * 2 * pi
            zeros.append(s['zeros'][i])
        pz[nepoch,'z'] = zeros, nzeros
        nzeros = 0
#    else:
#        print("warning. no zeros, set dict value to null")
#        pz[nepoch,'z'] = '', 0
#        print("zeros: %s" % zeros)

for iepoch in range(1, nepoch+1):
    # check A0 values
    getA0_theo(pz[iepoch, 'fn'], pz[iepoch, 'fs'], pz[(iepoch, 'p')][0], pz[(iepoch, 'z')][0],
            pz[(iepoch, 'p')][1], pz[(iepoch, 'z')][1], pz[iepoch, 'A0'])

    # output to file
    pz2file(pz[(iepoch, 'p')][0], pz[(iepoch, 'z')][0], pz[(iepoch, 'p')][1],
            pz[(iepoch, 'z')][1], pz[iepoch, 'A0'], pz[iepoch, 'epoch'])

