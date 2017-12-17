#!/usr/bin/env python
"""
convert tf_NN.pz files to antelope pole-zero format


20160505 cralvizuri <cralvizuri@alaska.edu>
"""

import sys
import os

if len(sys.argv) < 2 :
    print "Usage: " + os.path.basename(sys.argv[0]) + " RESP-file"
    print "      Converts a RESP file into a SAC PoleZero file"
    print "      SAC Polezero file is built from station/net/channel"
    exit(-1)

filename = sys.argv[1]
fid = open(filename)

P = list()
Z = list()
l4 = list()
np = 0
nz = 0
for line in fid:
    l2 = line.replace("\n", "")
    l3 = l2.split()
    l4.append(l3)

    if "POLES" in l3:
        np = int(l3[1])
    elif "ZEROS" in l3:
        nz = int(l3[1])
    elif "CONSTANT" in l3:
        constant = float(l3[1])
fid.close()

for i in range(len(l4)):
    if "POLES" in l4[i]:
        for ip in range(i+1, i+np+1):
            P.append(l4[ip])
    if "ZEROS" in l4[i]:
        for iz in range(i+1, i+nz+1):
            Z.append(l4[iz])


def pz2file(P, Z, np, nz, filename):
    out = "antpz_" + filename
    global fp
    fp = open(out, 'w')
    print >> fp, ("# This response file was derived from: %s" % filename)
    print >> fp, ("# Epoch: not available for this file type")
    print >> fp, ("# Format: Antelope \"paz\"")
    print >> fp, ("# See man pages for description of the format")
    print >> fp, ("# Antelope version 5.4")
    print >> fp, ("#")
    print >> fp, ("theoretical 1 complete paz LLNL/UAF")
    print >> fp, ("%12.6e" % (constant))

    printpz(P, np, fp)
    printpz(Z, nz, fp)

    fp.close()

def printpz(PZ, npz, fid):
    print >> fp, ("%s" % npz)
    for i in range(npz):
        print >> fp, ("%13.6e  %13.6e" % (float(PZ[i][0]), float(PZ[i][1])))

pz2file(P, Z, np, nz, filename)

