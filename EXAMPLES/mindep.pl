#!/bin/env perl

# script for doing specific jobs at best depth. example : saving the output files into another directory, creating log files at best depth and magnitude

($eid,$norm,$wts,$resultdir) = @ARGV;

$depfile = "./L$norm/M$wts/junk1.out";
system("more $depfile");
open(IN,$depfile); @deplines=<IN>;
$nlines=@deplines;
#$resultdir="RESULTS";

$min0 = -1.0e+19;	# impossibly large misfit value

for ($i=0;$i<$nlines;$i++){
  (undef,undef,undef,$smodeldep,undef,$strike,$dip,$rake,undef,$Mw,undef,$err,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,$VR,undef,undef) = split(" ",$deplines[$i]);
  ($model,$depth) = split("_",$smodeldep);
  if ($VR > $min0){
    $min0=$VR;
    $mindep=$depth;
    $bestMw=$Mw;
    $beststrike=$strike;
    $bestdip=$dip;
    $bestrake=$rake;
    $bestVR=$VR;
  }
}

# 1. Move figures and files for best depth solution into another directory
if (1){
    #system("rm -fr $resultdir");
    system("mkdir -p $resultdir");
    system("cp L${norm}/M${wts}/scak_${mindep}.ps $resultdir/${eid}_L${norm}_M${wts}.ps");
    system("cp L${norm}/M${wts}/scak_${mindep}.out $resultdir/${eid}_L${norm}_M${wts}.out");
    system("cp L${norm}/M${wts}/dep_${eid}.ps $resultdir/dep_${eid}_L${norm}_M${wts}.ps");
}
print "$eid/L$norm/M$wts/scak_$mindep.out\n";

# 2. Create summary of solutions
if (0){
    $summary_file = sprintf("$resultdir/$eid/summary.dat",$eid);
    print "$summary_file\n";
    open(my $sum,'>>',$summary_file);
    say $sum "$norm $wts $beststrike $bestdip $bestrake $bestMw $mindep $bestVR";
    close $sum;
}
