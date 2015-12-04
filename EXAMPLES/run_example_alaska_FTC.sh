#!/bin/bash
#====================================================================================
# tests results of FTC (Filter-then-Cut)
#
# $eid = 20080918194353069 (southern Alaska event for FTC testing)
# $CAPHOME - Directory where cap.c is located
# $CAPRUN/inv/cus - Directory from where you will be running cap
# $ExampleDir - Directory where you want to generate new plots and output files
#
# Runs and save in the $ExampleDir
# MIGHT want to rename the previous $ExampleDir and save it (just in case)
# Search $ExampleDir and make sure it does not overwrite anything
# You will be comparing the results between $ExampleDir and $CAPHOME/EXAMPLES/20080418093700_check
#====================================================================================
# EXAMPLE 1: Running Alaska Example (20080918194353069) - Checking FTC data and green's outputs

# 1. Prepare the data directory and directory for storing the results ($ExampleDir)
eid=20080918194353069
model=scak

ExampleDir=$CAPRUN/inv/"$model"/${eid}_check2
mkdir -p $ExampleDir

# Prepare data directory. Copy the the data to directory where you will be running cap
# All commands will be run form here (Alternatively, one could run from $CAPHOME/EXAMPLES/ without copying the data)
cd $CAPRUN/inv/"$model"
rsync -av $CAPHOME/EXAMPLES/$eid .

# 2. Prepare cap.c 
# NOTE: To run the example you will have to copy paste commands from here and make changes to cap.c after
#!!!!!!!! MAKE CHANGES IN CAP.C !!!!!!!!!!!!!
# > cd $CAPHOME
# > git checkout cap.c
# > open cap.c and set: 
#    FTC_data=0
#    FTC_green=0 
# then recompile 
# > make all

# 3. Run cap
cap.pl -H0.02 -P1/20/60 -p2 -S3/10/0 -T15/120 -D1/1/0.5 -C0.2/0.5/0.025/0.06 -W1 -M"$model"_72/4.8 -I10/0.1 -Zweight101.dat -E1 -K1 -Y1 $eid
cp $eid/scak_072.ps $ExampleDir/FTC_off.ps
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $ExampleDir/FTC_off.ps $CAPHOME/EXAMPLES/20080918194353069_check/FTC_off.ps &


# 4. Prepare cap.c AGAIN
#!!!!!!!! MAKE CHANGES IN CAP.C !!!!!!!!!!!!!! 
# > cd $CAPHOME
# > git checkout cap.c
# > open cap.c and set: 
#    FTC_data=1
#    FTC_green=1 
# then recompile 
# > make all

# 5. Run cap AGAIN
cap.pl -H0.02 -P1/20/60 -p2 -S3/10/0 -T15/120 -D1/1/0.5 -C0.2/0.5/0.025/0.06 -W1 -M"$model"_72/4.8 -I10/0.1 -Zweight101.dat -E1 -K1 -Y1 $eid
cp $eid/scak_072.ps $ExampleDir/FTC_on.ps
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $ExampleDir/FTC_on.ps $CAPHOME/EXAMPLES/20080918194353069_check/FTC_on.ps &

# 6. Prepare cap.c AGAIN
#!!!!!!!! MAKE CHANGES IN CAP.C !!!!!!!!!!!!!
# > cd $CAPHOME
# > git checkout cap.c
# > open cap.c and set: 
#    FTC_data=1
#    FTC_green=0 
# then recompile 
# > make all

# 7. Run cap AGAIN
cap.pl -H0.02 -P1/20/60 -p2 -S3/10/0 -T15/120 -D1/1/0.5 -C0.2/0.5/0.025/0.06 -W1 -M"$model"_72/4.8 -I10/0.1 -Zweight101.dat -E1 -K1 -Y1 $eid
cp $eid/scak_072.ps $ExampleDir/FTC_on_off.ps
#-------- To check the result make sure both files are same (if running from the terminal, you will have to replace $ExampleDir with the path above (search $ExampleDir)-----------
# gv $ExampleDir/FTC_on_off.ps $CAPHOME/EXAMPLES/20080918194353069_check/FTC_on_off.ps &

#====================================================================================
