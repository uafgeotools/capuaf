#!/bin/bash
#
# COMMANDS TO RECREATE INVERSION RESULTS
# Uturuncu event 20100516063454464 
# Velocity model utuhalf
# FMT with lune parameterization
# Full grid search (no line search)
# ==============================================================================

# Part 0. generate Green's functions 
# NOTE "utuhalf" is a homogeneous half-space
# NOTE You do not need to run fk if you have a pre-computed library of Green's functions

# 1. run the frequency-wavenumber code fk
mkdir -p $CAPRUN/models/utuhalf
cd $CAPRUN/models/utuhalf/
cp $CAPHOME/EXAMPLES/20100516063454464/utuhalf .
fk.pl -Mutuhalf/4 -N16384/0.01 11 14 18 19 21 23 31 33 37 50

# 2. run fk again to generate isotropic components
# This will create 3 additional green functions for each distance of the form dist.grn.(a-c),
# eg. 11.grn.a, 11.grn.b, 11.grn.c
fk.pl -Mutuhalf/4 -N16384/0.01 -S0 11 14 18 19 21 23 31 33 37 50
