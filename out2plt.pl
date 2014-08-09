#!/usr/bin/env perl
#
# script for plotting the CAP results directly from output file

# these are the only things one need to change based on the site installation
$home = $ENV{HOME};                        # my home directory
$caphome = $ENV{CAPHOME};                  # CAP home directory
$caprun = $ENV{CAPRUN};                    # run directory

$debug =0;

require "$caphome/cap_plt.pl";             # include plot script

if (@ARGV < 2) {die("Usage: cap_sol eid\n")}
($eid, $outfile) = @ARGV;

chdir($eid);

open(FFF,$outfile);
@rslt = <FFF>;
$nrow = @rslt-1;
(undef,$md_dep,undef,$m1,undef,$m2,undef,$amplify,undef,$ampfact,undef,$ncom,undef,$spib,undef,$spis,$filterB,$filterS,$fmt_flag)=split(" ",@rslt[$nrow]);

if (debug) {
print "$nrow";
print "$md_dep, $m1, $m2, $amplify, $ampfact, $ncom, $spib, $spis, $filterB, $filterS, $fmt_flag";

$filterBand = $filterB." ".$filterS;
print "$filterBand";
}
&plot($md_dep, $m1, $m2, $amplify, $ampfact, $ncom, $spib, $spis, $filterBand, $fmt_flag);
chdir("../");
