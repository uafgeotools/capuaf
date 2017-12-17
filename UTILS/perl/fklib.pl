#!/usr/bin/perl -w

#---------------------------------------------
# script to generate commands for computing a full
# Green's function library for a given layered model (for "fk")
#
# our database: /store/wf/FK_synthetics
# 1. set parameters here, then run fklib.pl
# 2. see README_fklib
#
#---------------------------------------------

$pwd = $ENV{PWD};
$geotools = $ENV{GEOTOOLS};
$home = $ENV{HOME};

# base directory (where to calculate the green's functions)
$basedir = "/store/wf/FK_synthetics/";

# specify what sets of Green's functions you want
$idc = 1;    # Green's functions for double couple
$iiso = 1;   # additional Green's functions for full moment tensors

# range of depths -- INCLUSIVE
# note: the depth is computed every 1 km; depinc determines how many different version of fk you want to run at any time
$depmax = 200; $depinc = 10;  # depinc = number of csh files = number of run directories (example: directories inside tactmod_all)
# I

# range of distances -- INCLUSIVE
# WARNING: max number of distances is set in FK in model.h as ndis
$distmin = 1; $distmax = 500; $distinc = 1; $nsamp=8192; $samplerate = .02;

# FOR TESTING ONLY (Comment this line for full run)
$idc=1; $iiso=1; $depmax = 100; $depinc = 100; $distmax = 101; $distinc = 100; $nsamp=512; $samplerate = .02;

# input model
$iAK = 1; # Alaska specific model files are saved as ak_$model
$model = "scak";
$model = "tactmod";
#$model = "aleut";
#$model = "northak";
#$model = "northak";
#$model = "prem_cont_no_water"; $iAK = 0; # PREM continental without water layer; and set iAK=0 since its not Alaska specific

# Base directory for precreating the green's function directories (example: tactmod_all and scak_all) and csh file
$basedir = "$home/Downloads/tempx";

print "pwd = $pwd\n";
print "GEOTOOLS = $geotools\n";
print "basedir = $basedir\n";
print "rundir = $basedir/${model}_all\n";

# write the C-shell script to file
$cshfilemain = "fklib.csh";
open(CSHMAIN,">$cshfilemain");

# NEEDS CLARIFICATION
print CSHMAIN "mkdir $basedir\n";
print CSHMAIN "cd $basedir\n";
print CSHMAIN "mkdir ${model}_all\n";

# loop over the minimum starting depth
for ($i = 0; $i < $depinc-1; $i = $i+1) {

  $depmin = $i;

  # write the C-shell script to file
  $cshfile = sprintf("fklib_depmin%3.3i_depinc%3.3i.csh",$depmin,$depinc);
  open(CSH,">$cshfile");

  $dep = $depmin;
  while ($dep <= $depmax) {
    $sdep = sprintf("%.0f",$dep);

    # KEY COMMAND
    if($idc==1) {
       $cmd = "-M${model}/$sdep -N${nsamp}/${samplerate}";
       print CSH "echo $cmd\n";
       print CSH "fk.pl $cmd";
       $dist = $distmin;
       while ($dist <= $distmax) {
         print CSH " $dist";
         $dist = $dist + $distinc;
       }
       print CSH " >& ${model}_${sdep}_ofile\n";
    }
    if($iiso==1) {
       $cmd = "-M${model}/$sdep -N${nsamp}/${samplerate} -S0";
       print CSH "echo $cmd\n";
       print CSH "fk.pl $cmd";
       $dist = $distmin;
       while ($dist <= $distmax) {
         print CSH " $dist";
         $dist = $dist + $distinc;
       }
       print CSH " >& ${model}_${sdep}_iso_ofile\n";
    }
    # next depth
    $dep = $dep + $depinc;
}
  close CSH;

  # run these lines as admin from /store/wf/FK_synthetics/${model}_all/
  $slab = sprintf("%2.2i",$i);
  print CSHMAIN "mkdir ${model}_all/${model}_${slab}\n";
  print CSHMAIN "cp $pwd/$cshfile ${model}_all/${model}_${slab}\n";
  # pre-compiled fk code
  print CSHMAIN "cp /usr/local/seismo/bin/fk ${model}_all/${model}_${slab}\n";
  # 1D model -- THIS WILL NEED TO BE CHANGED
  if ($iAK==1) {
      print CSHMAIN "cp $geotools/tomo_util/fkmodels/ak_${model} ${model}_all/${model}_${slab}/${model}\n";}
  else {
      print CSHMAIN "cp $geotools/tomo_util/fkmodels/${model} ${model}_all/${model}_${slab}/${model}\n";}
}

print "done writing $depinc csh files, one for each set of depths\n";

close CSHMAIN;

print "see $cshfilemain\n";

#================================================
