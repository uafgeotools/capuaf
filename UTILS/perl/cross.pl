#!/usr/bin/perl -w

if (@ARGV < 1) {die("Usage: cap_sol eid\n")}
($eid) = @ARGV;

#==========Read files=================
$datadir = "/home/vipul/gmt/data/cap/MOOS/${eid}/L1/M111/";
# read in event info
$capevent = "$datadir/${eid}_event.dat";
if (not -f $capevent) {die("Check if capevent $capevent exist or not\n");}
(undef,$elat,$elon,$edep)  = split(" ",`grep $eid $capevent`);

# read output file
$capout = "$datadir/${eid}.out";
if (not -f $capout) {die("Missing File: Output file $capout missing");}
open(IN,$capout); @linescap = <IN>; 
(undef,undef,undef,$smodeldep,undef,$strike,$dip,$rake,undef,$Mw,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef) = split(" ",$linescap[0]);
($model,$depth) = split("_",$smodeldep);

#==========I/O Files=================
$cptfile = "mt.cpt";   # color palette table file
$cshfile = "cap_cross.csh";   # color shape file
$psfile = "${eid}_sec.ps";     # output postscript file

#==========Remove older Files=================
system("rm -f $cshfile");
system("rm -f $psfile");

open(CSH,">$cshfile");
# ========GMT commands===========
print CSH "gmtset PAPER_MEDIA Custom_4ix5i MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 0.2c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT 1 ANOT_FONT 0.5 LABEL_FONT 1 LABEL_OFFSET 0.1c HEADER_FONT_SIZE 15 FRAME_PEN 1p TICK_PEN 1p PAGE_ORIENTATION landscape CHAR_ENCODING symbol\n";

#===================== plot XYZ =====================================
$scale=100;
$rak_min=-90/$scale;
$rak_max=90/$scale;
$dip_min=0/$scale;
$dip_max=90/$scale;
$stk_min=0/$scale;
$stk_max=360/$scale;
$rak_len=$rak_max-$rak_min;
$dip_len=$dip_max-$dip_min;
$stk_len=$stk_max-$stk_min;
$J1 = "-JX$rak_len/$dip_len";
$J2 = "-JX$rak_len/$stk_len";
$J3 = "-JX$dip_len/$stk_len";
$XY1 = "-Xa12c -Ya2.5c";
$XY2 = "-Xa12c -Ya5.5c";
$XY3 = "-Xa17.5c -Ya5.5c";
$XY4 = "-Xa16.5c -Ya0.5c";

# Different color palette files for misfit and omega
# Manually find min and max by: example:
#    >> minmax /home/vipul/gmt/data/cap/MOOS/20071010180326301/L1/M111//20071010180326301_omega_misfit_outline.dat 
#    >> N = 122	<0/180>	<29.522876/39.63632>
# And then subtract and add 1 to min and max
print CSH "makecpt -Cpanoply -T29.5/40.6/.01 -D -Z> $cptfile \n";
$omegacpt = "color_omega.cpt";
$Tomega = "-T0/180/15 -D -I";
print CSH "makecpt -Csplit $Tomega > $omegacpt\n";

$stkfile = "$datadir/${eid}_stkfile.dat";
$dipfile = "$datadir/${eid}_dipfile.dat";
$rakfile = "$datadir/${eid}_rakfile.dat";
$stkcont = "$datadir/${eid}_stkcont.dat";
$dipcont = "$datadir/${eid}_dipcont.dat";
$rakcont = "$datadir/${eid}_rakcont.dat";

$E = "-E-70/90";
# Put (strike dip rake) value here
#$strike=220; $dip=90; $rake=0;
$strike=40; $dip=86; $rake=5;
$strikem=309.6; $dipm=85; $rakem=-4;

$dip=sprintf("%d",$dip);

#print CSH "makecpt -Cwysiwyg -T0/20/1 -D -Z > $cptfile \n";
#print CSH "makecpt -Cpanoply -T0/16/.01 -D -Z> $cptfile \n";
print CSH "psbasemap $J1 -R-90/90/0/90 -Ba45WeSn $XY1 -K -V>>$psfile\n";
print CSH "psxy $stkfile $J1 -R-90/90/0/90 -Ba45:rake:/a45:dip:WeSn -C$cptfile -Ss17p $XY1 -K -O -V>>$psfile\n";
print CSH "pscontour $stkcont -C$omegacpt $J1 -R-90/90/0/90 -Gn1 -W1p $XY1 -K -O -V>>$psfile\n";
print CSH "psxy $J1 -R-90/90/0/90 -SR -W0.5p,0/0/0,- $XY1 -K -O -V >>$psfile<<EOF\n$rake 0\n$rake 90\nEOF\n";
print CSH "psxy $J1 -R-90/90/0/90 -W0.5p,0/0/0,- $XY1 -K -O -V >>$psfile<<EOF\n-90 $dip\n90 $dip\nEOF\n";
print CSH "psxy $J1 -R-90/90/0/90 -Ba45:rake:/a45:dip:WeSn -Sa12p -G255/255/255 -W0/0/0 $XY1 -K -O -V>>$psfile<<EOF\n$rake $dip\nEOF\n";
print CSH "psmeca $J1 -R-90/90/0/90 -Sa0.2 -Gred -L1p/0/0/0 -K -O -N -V $XY1 >>$psfile<<EOF\n$rakem $dipm 10 $strikem $dipm $rakem $Mw 0 0\nEOF\n";
print CSH "psmeca $J1 -R-90/90/0/90 -Sa0.2 -G100 -L1p/0/0/0 -K -O -N -V $XY1 >>$psfile<<EOF\n$rake $dip 10 $strike $dip $rake $Mw 0 0\nEOF\n";
print CSH "pstext $J1 -R-90/90/0/90 -W255/255/255,o/5p $XY1 -K -N -O -V>> $psfile<<EOF\n-90 90 10 0 1 LT strike=$strike\nEOF\n";

#print CSH "pscontour $xyzfile -JX5c/5c -R-90/90/0/90 -Ba45:rake:/a45:dip:WeSn -C$cptfile -Lthinnest -I -Xa2c -Ya2.5c -W0.5p -K -O -V>>$psfile\n";
# print CSH "xyz2grd $xyzfile -G$stkgrd -I20+/10+ -R-90/90/10/90 -V\n";
# print CSH "grdimage $stkgrd -C$cptfile -R-90/90/0/90 -JX5c/5c -Ba45WeSn -Xa2c -Ya2.5c -K -O -V>>$psfile\n";
print CSH "psbasemap $J2 -R-90/90/0/360 -Ba45Wesn $XY2 -K -O -V>>$psfile\n";
print CSH "psxy $dipfile $J2 -R-90/90/0/360 -Ba45:rake:/a45:strike:Wesn -C$cptfile -Ss17p $XY2 -K -O -V>>$psfile\n";
print CSH "pscontour $dipcont -C$omegacpt $J2 -R-90/90/0/360 -Gn1 -W1p $XY2 -K -O -V>>$psfile\n";
print CSH "psxy $J2 -R-90/90/0/360 -W0.5p,0/0/0,- $XY2 -K -O -V >>$psfile<<EOF\n$rake 0\n$rake 360\nEOF\n";
print CSH "psxy $J2 -R-90/90/0/360 -W0.5p,0/0/0,- $XY2 -K -O -V >>$psfile<<EOF\n-90 $strike\n90 $strike\nEOF\n";
print CSH "psxy $J2 -R-90/90/0/360 -Ba45:rake:/a45:strike:Wesn -Sa12p -G255/255/255 -W0 $XY2 -K -O -V>>$psfile<<EOF\n$rake $strike\nEOF\n";
print CSH "psmeca $J2 -R-90/90/0/360 -Sa0.2 -Gred -L1p/0/0/0 -K -O -N -V $XY2>>$psfile<<EOF\n$rakem $strikem 10 $strikem $dipm $rakem $Mw 0 0\nEOF\n";
print CSH "psmeca $J2 -R-90/90/0/360 -Sa0.2 -G100 -L1p/0/0/0 -K -O -N -V $XY2>>$psfile<<EOF\n$rake $strike 10 $strike $dip $rake $Mw 0 0\nEOF\n";
print CSH "pstext $J2 -R-90/90/0/360 -W255/255/255,o/5p $XY2 -K -N -O -V>> $psfile<<EOF\n-90 360 10 0 1 LT dip=$dip\nEOF\n";
#print CSH "pstext -JM7i -R0/1/0/1 -W255/255/255,o/5p -Xa12.1c -Ya14.3c -K -N -O -V>> $psfile<<EOF\n0 0 10 0 1 LB Dip=$dip\nEOF\n";

print CSH "psbasemap $J3 -R0/90/0/360 -Ba45wESn $XY3 -Ya5.5c  -K -O -V>>$psfile\n";
print CSH "psxy $rakfile $J3 -R0/90/0/360 -Ba45:dip:/a45::wESn -C$cptfile -Ss17p -: $XY3 -K -O -V>>$psfile\n";  # -: toggles between (x,y) and (y,x)
print CSH "pscontour $rakcont -C$omegacpt $J3 -R0/90/0/360 -Gn1 -W1p $XY3 -: -K -O -V>>$psfile\n";
print CSH "psxy $J3 -R0/90/0/360 -W0.5p,0/0/0,- $XY3 -K -O -V >>$psfile<<EOF\n$dip 0\n$dip 360\nEOF\n";
print CSH "psxy $J3 -R0/90/0/360 -W0.5p,0/0/0,- $XY3 -K -O -V >>$psfile<<EOF\n0 $strike\n90 $strike\nEOF\n";
print CSH "psxy $J3 -R0/90/0/360 -Ba45:dip:/a45::wESn -Sa12p -G255/255/255 -W0 -: $XY3 -K -O -V>>$psfile<<EOF\n$strike $dip\nEOF\n";
print CSH "psmeca $J3 -R0/90/0/360 -Sa0.2 -Gred -L1p/0/0/0 -K -O -N -V $XY3 >>$psfile<<EOF\n$dipm $strikem 10 $strikem $dipm $rakem $Mw 1 1\nEOF\n";
print CSH "psmeca $J3 -R0/90/0/360 -Sa0.2 -G100 -L1p/0/0/0 -K -O -N -V $XY3 >>$psfile<<EOF\n$dip $strike 10 $strike $dip $rake $Mw 1 1\nEOF\n";
print CSH "pstext $J3 -R0/90/0/360 -W255/255/255,o/5p $XY3 -K -N -O -V>> $psfile<<EOF\n0 360 10 0 1 LT rake=$rake\nEOF\n";
#print CSH "pstext -JM7i -R0/1/0/1 -W255/255/255,o/5p -Xa17.6c -Ya14.3c -K -N -O -V>> $psfile<<EOF\n0 0 10 0 1 LB Rake=$rake\nEOF\n";

#print CSH "pstext -JM7i -R0/1/0/1 -Xa12c -Ya15.2c -K -N -O -V>> $psfile<<EOF\n0 0 12 0 1 LB (f)\nEOF\n";
print CSH "psscale -C$cptfile -D19.75c/3.25c/2.5c/0.2c -K -O -Ac -B4::/:: >> $psfile\n"; # relative to the basemap

$J = "-JX4.5c/4.5c";
$R = "-R-2/2/-2/2";
$X = "-X1c";
$Y = "-Y19.5c";

#print CSH "psbasemap -JX3c/3c -R-1/1/-1/1 -Ba30/a1wesn -X9c -Y22c -K -O -V>>$psfile\n";
print CSH "psmeca $J $R -Sa2.5c -G150/150/150 -L1p/0/0/0 -N $XY4 -K -O -V >>$psfile<<EOF\n0 0 $edep $strike $dip $rake 5 0 0\nEOF\n";

#========================================================================================
# Always use
print CSH "pstext -JM7i -R0/1/0/1 $Xr $Yr -N -O -V>> $psfile<<EOF\n0 0 14 0 1 LB \nEOF\n";

close (CSH);
system("csh -f $cshfile");
#system("gv $psfile &");
