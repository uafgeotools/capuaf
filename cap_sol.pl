#!/usr/bin/perl -w

if (@ARGV < 1) {die("Usage: cap_sol eid\n")}
($eid) = @ARGV;

#==========Read files=================
$datadir = "/home/vipul/gmt/data/${eid}";
# read in event info
$capevent = "$datadir/${eid}_event.dat";
if (not -f $capevent) {die("Check if capevent $capevent exist or not\n");}
(undef,$elat,$elon,$edep)  = split(" ",`grep $eid $capevent`);

# read station info
$capstation = "$datadir/${eid}_station.dat";  
if (not -f $capstation) {die("Missing File: Couldn't find $capstation in the running directory");}
open(IN,$capstation); @capstationlines = <IN>; $nsta = @capstationlines;

# read output file
$capout = "$datadir/${eid}.out";
if (not -f $capout) {die("Missing File: Output file $capout missing");}
open(IN,$capout); @linescap = <IN>; 
(undef,undef,undef,$smodeldep,undef,$strike,$dip,$rake,undef,$Mw,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef) = split(" ",$linescap[0]);
(undef,undef,undef,$elat,undef,$elon,undef,$edep)  =split(" ",$linescap[1]);
($model,$depth) = split("_",$smodeldep);
$Mw=sprintf("%.1f",$Mw);
$depth=sprintf("%d",$depth);
@capcols = (14,22,13,17,21,6,5,9);  # column of data to plot
@wcols = (11,19,11,15,19,3,3,7);    # columns with weight factors

$histogramfile = "$datadir/${eid}_histo.dat";
if (not -f $histogramfile) {die("Missing File: Histogram file $histogramfile missing");}

$samplesfile = "$datadir/${eid}_samples.dat";
#$samplesfile = "/home/vipul/gmt/data/samples.dat";
if (not -f $samplesfile) {die("Missing File: Moment tensor file $samplesfile missing");}
open(IN,$samplesfile); @sampleslines = <IN>; $nsamples = @sampleslines;
#($omega,$misfit,$stk_samples,$dip_samples,$rak_samples) = split(" ",$sampleslines[0]);

$histogramfile2 = "$datadir/${eid}_histo2.dat";
open(IN,$histogramfile2); @histo2 = <IN>;
($strike2,$dip2,$rake2)= split(" ",$histo2[0]);

#==========I/O Files=================
$cptfile = "map_alaska.cpt";   # color palette table file
$cshfile = "cap_solution.csh";   # color shape file
$psfile = "${eid}_solution.ps";     # output postscript file
$grdfile = "/home/admin/share/datalib/topography/GLOBAL/ETOPO1/ETOPO1_Bed_g.grd";
$gradfile = "/home/admin/share/datalib/topography/GLOBAL/ETOPO1/ETOPO1_Bed_g.grad";
$faultfile = "/home/admin/share/datalib/faults/alaska/DGGS2012/ak_fault_database_faults.gmt";
#$histogramfile = "/home/vipul/gmt/data/${eid}/${eid}_histo.dat";
$histogramfile2 = "$datadir/${eid}_histo2.dat";

# ==========Plotting specifications=====
$xmin = -155;
$xmax = -145;
$ymin = 57;
$ymax = 65;
$scale = 7500000;
$xcen = ($xmin + $xmax)/2;
$ycen = ($ymin + $ymax)/2;
$J = "-Jb$xcen/$ycen/$ymin/$ymax/1:$scale";
$R = "-R$xmin/$xmax/$ymin/$ymax";
#$B = "-Ba4f2:.Southern_Alaska:WesN";
$B = "-Ba4f2WesN";
$X = "-X14c";
$Y = "-Y16c";
$Xr = "-X-14c";
$Yr = "-Y-16c";
$title = "eid: $eid";
$stainfoa= "-Si5p -W3p,0/0/0";
$stainfob = "-Si5p -W1p,255/255/255";

# OPtions
$iinset = 0;
$iomega = 1;
$iminus = 0;

# $alaska = "-175/-135/50/70";
# $Ra = "-R$alaska";
# $Ja = "-Jb$region/0.5c";
# $Ba = "-Ba10f5";

open(CSH,">$cshfile");

# ========GMT commands===========
print CSH "gmtset PAPER_MEDIA letter MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 0.2c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT 1 ANOT_FONT 0.5 LABEL_FONT 1 LABEL_OFFSET 0.1c HEADER_FONT_SIZE 15 FRAME_PEN 1p TICK_PEN 1p PAGE_ORIENTATION portrait CHAR_ENCODING symbol\n";
# cut section of grid
#print CSH "grdcut $grdfile -G$gridfile $R -V\n";
# make color pallete for plotting map
print CSH "makecpt -Cglobe -T-3000/3000/200 -D -Z > $cptfile \n";
#print CSH "grd2cpt $grdfile -Cglobe -Z >> $cptfile \n";
# create basemap
#print CSH "psbasemap $R $J −Ba2f1d:" ":/a2f1d:" ":WeSn -K -V -X3c  -Y3c > $psfile \n";
print CSH "psbasemap $R $J $B $X $Y -K  > $psfile \n";
# plot topography
print CSH "grdimage $grdfile $R $J $B -E100 -C$cptfile -I$gradfile -V -K -O >> $psfile \n";
# plot faultline
print CSH "psxy $faultfile $J $R -m -W/255/0/0 -V -K -O >> $psfile\n";
# plot coastline
print CSH "pscoast $R $J -A0 -Df -W -N1 -L-148/58/58.5/200+p1p,0/0/0,solid+f255/255/255 -O -K -V >> $psfile\n";
# print CSH "psxy $J $R -W1p,255/255/255 -Si10p -V -K -O >> $psfile<<EOF\n-150 60\nEOF\n";
# print CSH "pstext $J $R -G0/0/0 -W255/255/255 -O -V -K>> $psfile<<EOF\n$elon $elat  14 0 1 CM $eid\nEOF\n";
print CSH "pstext -JM7i -R0/1/0/1 -Xa0c -Ya11c -K -N -O -V>> $psfile<<EOF\n0 0 12 0 1 LB (c)\nEOF\n";

# scale bar (????)
# print CSH "pscoast $J $R -L-152/57.5/57.5/100+p0.5p,0/0/0,solid+f255/255/255 -N1 -K -O -V >> $psfile";

#=================== plot stations=======================
for ($xx=1; $xx <= 8; $xx++){
for ($i = 1; $i <= $nsta; $i++) {
   ($sta,$net,$lat,$lon,$dist,$az) = split(" ",$capstationlines[$i-1]);
   # plot an open triangle at each station in the output file
   print CSH "psxy $J $R $stainfoa -K -O -V >>$psfile<<EOF\n $lon $lat\nEOF\n";
   print CSH "psxy $J $R $stainfob -K -O -V >>$psfile<<EOF\n $lon $lat\nEOF\n";
   @capout = split(" ",`grep $sta $capout`);
   $iplot = $capcols[$xx-1];
   $fplot = $capout[$iplot-1]; 
   if($capout[0]) {
     #print "-- @capout -- $Tdt $Rdt $Zdt\n";
     $iplot2 = $wcols[$xx-1];
     $wplot = $capout[$iplot2-1];
     #print "\n-- xx = $xx -- iplot2 = $iplot2 -- wplot = $wplot\n";
     # only plot the colored ray path if the station was used in the inversion (nonzero weight)
     if ($wplot > 0) {
       print CSH "psxy $J $R -W0.5p,0/0/0 -K -O -V >>$psfile<<EOF\n$elon $elat\n$lon $lat\nEOF\n";
       #print CSH "psxy $J $R -W -m -K -O -V >>$psfile<<EOF\n> -Z$fplot\n$elon $elat\n$lon $lat\nEOF\n";
     }
   }
 }
}

# plot moment tensor beachball
print CSH "psmeca $J $R -Sa0.2 -G150/150/150 -L1p/0/0/0 -K -O -V >>$psfile<<EOF\n$elon $elat $edep $strike $dip $rake $Mw 0 0\nEOF\n";

# Plot cities
print CSH "pstext $R $J -G0/0/0 -S2p,255/255/255 -K -N -O -V>> $psfile<<EOF\n-150.02 61.17 12 0 1 LB A\nEOF\n"; # Anchorage
print CSH "pstext $R $J -G0/0/0 -S2p,255/255/255 -K -N -O -V>> $psfile<<EOF\n-147.87 64.82 12 0 1 LB F\nEOF\n"; # Fairbanks

#print CSH "psscale -C$cptfile -D3.5c/-1c/6c/0.25ch -O -Ac -B1500:Elevation:/:m: -E >> $psfile\n"; # relative to the basemap

# ========Specifications for inset========

if ($iinset==1) {
  $xmina = -175;
  $xmaxa = -135;
  $ymina = 52;
  $ymaxa = 72;
  $scalea = 65000000;
  $xcena = ($xmina + $xmaxa)/2;
  $ycena = ($ymina + $ymaxa)/2;
  $Ja = "-Jb$xcena/$ycena/$ymina/$ymaxa/1:$scalea";
  $Ra = "-R$xmina/$xmaxa/$ymina/$ymaxa";
  $Ba = " -B100wesn";
  $Xa = "-Xa3c";
  $Ya = "-Ya-0.7c";
  #$gridfilea = "alaska_grid.grd";

  #print CSH "grdcut $grdfile -G$gridfilea $Ra -V\n";
  #print CSH "makecpt -Ccopper -T0/1/.0001 -D -Z > $cptfilea \n";
  print CSH "psbasemap $Ra $Ja $Ba $Xa $Ya -O -K -V >> $psfile \n";
  #print CSH "grdimage $grdfile -C$cptfilea $Ra $Ja $Ba $Xa $Ya  -O >> $psfile \n";
  # plot coastline
  print CSH "pscoast $Ra $Ja -A500/0/1 -G150 -Dl -S220/255/255 -C220/255/255 -W -N1/1p $Xa $Ya -K -O -V >> $psfile\n";
  # section of inset
  print CSH "psxy $Ja $Ra -W1p,0/0/0 $Xa $Ya -O -V -K >> $psfile<<EOF\n$xmin $ymin\n$xmin $ymax\n$xmax $ymax\n$xmax $ymin\n$xmin $ymin\nEOF\n";
  # plot faultline
  print CSH "psxy $faultfile $Ja $Ra -m -W/255/0/0 $Xa $Ya -V -O -K >> $psfile\n";
  #print CSH "psscale -C$cptfile -D3.5c/-1c/6c/0.25ch -O -Ac -B1500:Elevation:/:m: -E >> $psfile\n"; # relative to the basemap
  #print CSH "pstext $Ja $Ra -G0/0/0 -W255/255/255 $Xa $Ya -O -V>> $psfile<<EOF\n-150 60  14 0 1 CM \nEOF\n";
}

# Origin Marker
print CSH "pstext -JM7i -R0/1/0/1 -N $Xr $Yr -O -K -V>> $psfile<<EOF\n0 0 14 0 1 LB \nEOF\n";

# title
print CSH "pstext -JM7i -R0/1/0/1 -N -Xa1c -Ya26.75c -O -K -V>> $psfile<<EOF\n0 0 14 0 1 LB $title\nEOF\n";
print CSH "pstext -JM7i -R0/1/0/1 -N -Xa1c -Ya26.25c -O -K -V>> $psfile<<EOF\n0 0 10 0 1 LB Strike $strike Dip $dip Rake $rake\nEOF\n";
print CSH "pstext -JM7i -R0/1/0/1 -N -Xa1c -Ya25.75c -O -K -V>> $psfile<<EOF\n0 0 10 0 1 LB Strike $strike2 Dip $dip2 Rake $rake2\nEOF\n";
print CSH "pstext -JM7i -R0/1/0/1 -N -Xa1c -Ya25.25c -O -K -V>> $psfile<<EOF\n0 0 10 0 1 LB Mw $Mw Depth $depth km\nEOF\n";

# ====================plotting histograms====================
print CSH "gmtset LABEL_OFFSET 0.03c\n";
print CSH "pshistogram $histogramfile2 -JX3c/3c -R0/360/0/50 -G0/0/0 -L0.25p,0/0/0 -Z1 -W10  -Ba90:strike:/a10::WeSn -Xa6.5c -Ya2c -T0 -K -O -V >>$psfile\n";
print CSH "pshistogram $histogramfile -JX3c/3c -R0/360/0/50 -G0/255/0 -L0.25p,0/0/0 -Z1 -W10  -Ba90f45:strike:/a10::WeSn -Xa6.5c -Ya2c -T0 -K -O -V >>$psfile\n";
print CSH "pshistogram $histogramfile2 -JX3c/3c -R0/90/0/50 -G0/0/0 -L0.25p,0/0/0 -Z1 -W5  -Ba30:dip:/a10::WeSn -Xa6.5c -Ya6.9c -T1 -K -O -V >>$psfile\n";
print CSH "pshistogram $histogramfile -JX3c/3c -R0/90/0/50 -G0/255/0 -L0.25p,0/0/0 -Z1 -W5  -Ba30:dip:/a10::WeSn -Xa6.5c -Ya6.9c -T1 -K -O -V >>$psfile\n";
print CSH "pshistogram $histogramfile2 -JX3c/3c -R-180/180/0/50 -G0/0/0 -L0.25p,0/0/0 -Z1 -W10  -Ba90f45:rake:/a10::WeSn -Xa6.5c -Ya11.8c -T2 -K -O -V >>$psfile\n";
print CSH "pshistogram $histogramfile -JX3c/3c -R-180/180/0/50 -G0/255/0 -L0.25p,0/0/0 -Z1 -W10  -Ba90f45:rake:/a10::WeSn -Xa6.5c -Ya11.8c -T2 -K -O -V >>$psfile\n";
print CSH "pstext -JM7i -R0/1/0/1 -Xa6.5c -Ya15.2c -K -N -O -V>> $psfile<<EOF\n0 0 12 0 1 LB (e)\nEOF\n";
#============================ plot beachballs===================
$Xm = "-X1c";
$Ym = "-Y1c";
$Xr = "-X-1c";
$Yr = "-Y-1c";
$Bm = "-BWeSn";
$Jm = "-JX4c/17c";
$Rm = "-R0/4/0/16.5";

print CSH "makecpt -Cpanoply -T0/16/.01 -D -Z> $cptfile \n";
#print CSH "makecpt -Cwysiwyg -T0/20/1 -D -Z > $cptfile \n";
#print CSH "gmtset BASEMAP_FRAME_RGB 255/0/0 \n";
print CSH "psbasemap $Jm $Rm $Bm $Xm $Ym -K -O -V>>$psfile\n";
print CSH "gmtset BASEMAP_FRAME_RGB 0/0/0 \n";
#print CSH "psmeca $Jm $Rm -Sa0.5 -G0/0/255 -L1p/0/0/0 -K -O -V >>$psfile<<EOF\n1 3 $edep $strike $dip $rake $Mw 0 0 0\nEOF\n";

# HORIZONTAL
#for ($i = 1; $i <= $nsamples; $i++) {
#  $xloc = 2*$i-1;
#  ($omega,$misfit,$stk_sample,$dip_sample,$rak_sample,$xi) = split(" ",$sampleslines[$i-1]);
#  if ($i <= 10) {
#    print CSH "psmeca $Jm $Rm -Sa0.4 -M -G0/0/255 -L1p/0/0/0 -Z$cptfile -K -O -V >>$psfile<<EOF\n$xloc 3 $misfit $stk_sample $dip_sample $rak_sample $Mw 0 0 $omega\nEOF\n";
#  } else {
#    $xloc = 2*($i-10)-1;
#    print CSH "psmeca $Jm $Rm -Sa0.4 -M -G0/0/255 -L1p/0/0/0 -Z$cptfile -K -O -V >>$psfile<<EOF\n$xloc 1 $misfit $stk_sample $dip_sample $rak_sample $Mw 0 0 $omega\nEOF\n";
#  }
#}

#VERTICAL
for ($i = $nsamples; $i >= 1; $i--) {
  ($omega,$misfit,$stk_sample,$dip_sample,$rak_sample,$xi) = split(" ",$sampleslines[$nsamples-$i]);
  $XI=sprintf("%.0f",$omega);
  if ($i%2 != 1) {
    $yloc = ($i-1)/1.25;  
    print CSH "psmeca $Jm $Rm -Sa0.4 -M -G0/0/255 -L1p/0/0/0 -Z$cptfile -K -O -V >>$psfile<<EOF\n1 $yloc $misfit $stk_sample $dip_sample $rak_sample $Mw 0 0 $XI\nEOF\n";
  } else {
    $yloc = $i/1.25;
    print CSH "psmeca $Jm $Rm -Sa0.4 -M -G0/0/255 -L1p/0/0/0 -Z$cptfile -K -O -V >>$psfile<<EOF\n3 $yloc $misfit $stk_sample $dip_sample $rak_sample $Mw 0 0 $XI\nEOF\n";
  }
}
#print CSH "psbasemap -JX2c/2c -R0/1/0/1 -Bwesn -Xa0c -Ya15c -K -O -V>>$psfile\n";  # for making a box for first sample (in case it is the best solution)
print CSH "psscale -C$cptfile -D18.5c/1.5c/3c/0.2c -K -O -Ac -B4::/:: -Ef >> $psfile\n"; # relative to the basemap
print CSH "pstext -JM7i -R0/1/0/1 -Xa20c -Ya1.2c -K -N -O -V>> $psfile<<EOF\n0 0 10 90 1 LB misfit\nEOF\n";
print CSH "pstext -JM7i -R0/1/0/1 -Xa0c -Ya17.25c -K -N -O -V>> $psfile<<EOF\n0 0 12 0 1 LB (d)\nEOF\n";
print CSH "pstext -JM7i -R0/1/0/1 -N $Xr $Yr -K -O -V>> $psfile<<EOF\n0 0 14 0 1 LB \nEOF\n";
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

$omegacpt = "color_omega.cpt";
$Tomega = "-T0/180/15 -D -I";
print CSH "makecpt -Csplit $Tomega > $omegacpt\n";
$stkfile = "$datadir/${eid}_stkfile.dat";
$dipfile = "$datadir/${eid}_dipfile.dat";
$rakfile = "$datadir/${eid}_rakfile.dat";
$stkcont = "$datadir/${eid}_stkcont.dat";
$dipcont = "$datadir/${eid}_dipcont.dat";
$rakcont = "$datadir/${eid}_rakcont.dat";

# plot corss-sections of opposite solution instead ; Also plot opposite moment tensor
if ($iminus==1){
$oppfile="$datadir/${eid}_opp_solution.dat";
open(IN,$oppfile); @opp = <IN>; 
($strike,$dip,$rake,$err_opp,$omega_opp,$xi_opp)=split(" ",$opp[0]);
$stkfile = "$datadir/${eid}_stkfile_opp.dat";
$dipfile = "$datadir/${eid}_dipfile_opp.dat";
$rakfile = "$datadir/${eid}_rakfile_opp.dat";
$stkcont = "$datadir/${eid}_stkcont_opp.dat";
$dipcont = "$datadir/${eid}_dipcont_opp.dat";
$rakcont = "$datadir/${eid}_rakcont_opp.dat";
}

#print CSH "makecpt -Cwysiwyg -T0/20/1 -D -Z > $cptfile \n";
print CSH "makecpt -Cpanoply -T0/16/.01 -D -Z> $cptfile \n";
print CSH "psbasemap $J1 -R-90/90/0/90 -Ba45WeSn $XY1 -K -O -V>>$psfile\n";
print CSH "psxy $stkfile $J1 -R-90/90/0/90 -Ba45:rake:/a45:dip:WeSn -C$cptfile -Ss17p $XY1 -K -O -V>>$psfile\n";
print CSH "pscontour $stkcont -C$omegacpt $J1 -R-90/90/0/90 -Gn1 -W1p $XY1 -K -O -V>>$psfile\n";
print CSH "psxy $J1 -R-90/90/0/90 -W0.5p,0/0/0,- $XY1 -K -O -V >>$psfile<<EOF\n$rake 0\n$rake 90\nEOF\n";
print CSH "psxy $J1 -R-90/90/0/90 -W0.5p,0/0/0,- $XY1 -K -O -V >>$psfile<<EOF\n-90 $dip\n90 $dip\nEOF\n";
print CSH "psxy $J1 -R-90/90/0/90 -Ba45:rake:/a45:dip:WeSn -Sa12p -G255/255/255 -W0/0/0 $XY1 -K -O -V>>$psfile<<EOF\n$rake $dip\nEOF\n";
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
print CSH "pstext $J2 -R-90/90/0/360 -W255/255/255,o/5p $XY2 -K -N -O -V>> $psfile<<EOF\n-90 360 10 0 1 LT dip=$dip\nEOF\n";
#print CSH "pstext -JM7i -R0/1/0/1 -W255/255/255,o/5p -Xa12.1c -Ya14.3c -K -N -O -V>> $psfile<<EOF\n0 0 10 0 1 LB Dip=$dip\nEOF\n";

print CSH "psbasemap $J3 -R0/90/0/360 -Ba45wESn $XY3-Xa17.5c -Ya5.5c  -K -O -V>>$psfile\n";
print CSH "psxy $rakfile $J3 -R0/90/0/360 -Ba45:dip:/a45::wESn -C$cptfile -Ss17p -: $XY3 -K -O -V>>$psfile\n";  # -: toggles between (x,y) and (y,x)
print CSH "pscontour $rakcont -C$omegacpt $J3 -R0/90/0/360 -Gn1 -W1p $XY3 -: -K -O -V>>$psfile\n";
print CSH "psxy $J3 -R0/90/0/360 -W0.5p,0/0/0,- $XY3 -K -O -V >>$psfile<<EOF\n$dip 0\n$dip 360\nEOF\n";
print CSH "psxy $J3 -R0/90/0/360 -W0.5p,0/0/0,- $XY3 -K -O -V >>$psfile<<EOF\n0 $strike\n90 $strike\nEOF\n";
print CSH "psxy $J3 -R0/90/0/360 -Ba45:dip:/a45::wESn -Sa12p -G255/255/255 -W0 -: $XY3 -K -O -V>>$psfile<<EOF\n$strike $dip\nEOF\n";
print CSH "pstext $J3 -R0/90/0/360 -W255/255/255,o/5p $XY3 -K -N -O -V>> $psfile<<EOF\n0 360 10 0 1 LT rake=$rake\nEOF\n";
#print CSH "pstext -JM7i -R0/1/0/1 -W255/255/255,o/5p -Xa17.6c -Ya14.3c -K -N -O -V>> $psfile<<EOF\n0 0 10 0 1 LB Rake=$rake\nEOF\n";

print CSH "pstext -JM7i -R0/1/0/1 -Xa12c -Ya15.2c -K -N -O -V>> $psfile<<EOF\n0 0 12 0 1 LB (f)\nEOF\n";
# ========================omega distribution==========================
$Xo = "-X7.5c";
$Yo = "-Y16.5c";
$Xr = "-X-7.5c";
$Yr = "-Y-16.5c";
$Jo = "-JX5c/4c";
$J1 = "-JX5c/3c";

if ($iomega==1) {
$omegaerr = "$datadir/${eid}_omega_err.dat";
$omegaerrpost = "$datadir/${eid}_omega_err_post.dat";
$omegaerrdis = "$datadir/${eid}_omega_err_dis.dat";

$Bo = "-Ba45:omega,\@~w\@~:/a4WeSn";
$Ro = "-R0/180/0/20";
$max = 3.1416;

# Plot all random samples, posterior and random samples from posterior
print CSH "psbasemap $Jo $Ro $Bo $Xo $Yo -K -O -V>>$psfile\n";
#print CSH "psxy $omegaerr $Jo $Ro $Bo -Sc1p -C$cptfile $Xo $Yo -K -O -V>>$psfile\n";   # for coloring samples based on their misfit value
#print CSH "psxy $omegaerr $Jo $Ro $Bo -Sc1p -G0/0/255 -K -O -V>>$psfile\n"; # homogeneously genertated random samples (blue)
print CSH "psxy $omegaerrpost $Jo $Ro $Bo -Sc2p -G0/255/0 -K -O -V>>$psfile\n"; # posterior samples (green)
print CSH "psxy $Jo $Ro $Bo -Sc5p -G255/255/0 -W0,1p -K -O -V>>$psfile<<EOF\n$omega_opp $err_opp\nEOF\n"; # opposite solution (-M i.e. omega--180)
# print CSH "psxy $samplesfile $Jo $Ro $Bo -Sc2p -G255/0/0 -K -O -V>>$psfile\n"; (randomly picked 20 samples from posterior - red)
#print CSH "psscale -C$cptfile -D6c/1.6c/3.2c/0.25c -K -O -Ac -B0::/:: -Ef >> $psfile\n";
print CSH "pstext -JM7i -R0/1/0/1 -Xa-1c -Ya1.5c -K -N -O -V>> $psfile<<EOF\n0 0 10 90 1 LB misfit\nEOF\n";

# Plot step curve for posterior and homogeneous distribution of Kagan angles
$omegapdf =  "$datadir/${eid}_omega_homo_pdf.dat";
$omegapostpdf =  "$datadir/${eid}_omega_post_pdf.dat";
print CSH "psbasemap $J1 -R0/180/0/3 -Ba45/a1wEsn -Xa0c -Ya4c -K -O -V>>$psfile\n";
print CSH "psxy $omegapdf $J1 -R0/$max/0/3 -Ba45/a1wEsn -Xa0c -Ya4c -W1p,0/0/255 -K -O -V>>$psfile\n";
print CSH "psxy $omegapostpdf $J1 -R0/$max/0/3 -Ba45/a1wEsn -Xa0c -Ya4c -W1p,0/255/0 -K -O -V>>$psfile\n";
print CSH "pstext -JM7i -R0/1/0/1 -Xa-1c -Ya4.5c -G0/255/0 -K -N -O -V>> $psfile<<EOF\n0 0 10 90 1 LB \@~s\@~(\@~w\@~)\nEOF\n";
print CSH "pstext -JM7i -R0/1/0/1 -Xa-1c -Ya5.5c -G0/0/255 -K -N -O -V>> $psfile<<EOF\n0 0 10 90 1 LB \@~m\@~(\@~w\@~)\nEOF\n";

# Plot the likelihood function
print CSH "psbasemap $J1 -R0/180/0/4 -Ba45/a1WesN -Xa0c -Ya7c -K -O -V>>$psfile\n";
print CSH "psxy $omegaerrdis $J1 -R0/$max/0/4 -Ba45/a1Wesn -Xa0c -Ya7c -W1p,255/0/0 -K -O -V>>$psfile\n";
#print CSH "psxy -JX6c/4c -R0/3.1416/0/4 -Ba45/a1Wesn -Xa0c -Ya7c -W1p,255/0/0 -K -O -V>>$psfile<<EOF\n0 5\n180 5\nEOF\n";
print CSH "pstext -JM7i -R0/1/0/1 -Xa-1c -Ya8.5c -K -N -O -V>> $psfile<<EOF\n0 0 10 90 1 LB L(\@~w\@~)\nEOF\n";

# Write the value of uncertainity
$OMEGAerr =  "$datadir/${eid}_OMEGA.dat";
open(IN,$OMEGAerr); @OMEGAval = <IN>; 
($OMEGA) = split(" ",$OMEGAval[0]);
#$OMEGA =$OMEGAval[0];
print CSH "pstext -JM7i -R0/1/0/1 -Xa3c -Ya9c -W255/255/255,o/5p -K -N -O -V>> $psfile<<EOF\n0 0 14 0 1 LB \@~W\@~ = $OMEGA\@+o\@+\nEOF\n";
#print CSH "pstext $OMEGAerr -JM7i -R0/1/0/1 -Xa0c -Ya12c -K -N -O -V>> $psfile \n";
print CSH "psbasemap -JX5c/10c -R0/180/0/4 -K -O -V>>$psfile\n";
print CSH "psxy -JX5c/10c -R0/180/0/4 -W1p,0/0/0,- -K -O -V >>$psfile<<EOF\n$OMEGA 0\n$OMEGA 4\nEOF\n";
print CSH "pstext -JM7i -R0/1/0/1 -Xa-0.75c -Ya10.5c -K -N -O -V>> $psfile<<EOF\n0 0 12 0 1 LB (b)\nEOF\n";
}
############ ELSE XI 
else {
$xierr = "$datadir/${eid}_xi_err.dat";
$xierrpost = "$datadir/${eid}_xi_err_post.dat";
$xierrdis = "$datadir/${eid}_xi_err_dis.dat";

$Bo = "-Ba30:\@~x\@-0\@-\@~:/a4WeSn";
$Ro = "-R0/120/0/20";
$max = 2.0944;
print CSH "psbasemap $Jo $Ro $Bo $Xo $Yo -K -O -V>>$psfile\n";
#print CSH "psxy $omegaerr $Jo $Ro $Bo -Sc1p -C$cptfile $Xo $Yo -K -O -V>>$psfile\n";
print CSH "psxy $xierr $Jo $Ro $Bo -Sc1p -G0/0/255 -K -O -V>>$psfile\n";
print CSH "psxy $xierrpost $Jo $Ro $Bo -Sc2p -G0/255/0 -K -O -V>>$psfile\n";
print CSH "psxy $Jo $Ro $Bo -Sc5p -G255/255/0 -K -O -V>>$psfile<<EOF\n$xi_opp $err_opp\nEOF\n";
for ($i = 1; $i <= $nsamples; $i++) {
  ($omega,$misfit,$stk_sample,$dip_sample,$rak_sample,$xi) = split(" ",$sampleslines[$i-1]);
  #print CSH "psxy $Jo $Ro $Bo -Sc2p -G255/0/0 -K -O -V>>$psfile<<EOF\n$xi $misfit\nEOF\n";   # plots 20 random samples from posterior space
}
#print CSH "psscale -C$cptfile -D6c/1.6c/3.2c/0.25c -K -O -Ac -B0::/:: -Ef >> $psfile\n";
print CSH "pstext -JM7i -R0/1/0/1 -Xa-1c -Ya1.5c -K -N -O -V>> $psfile<<EOF\n0 0 10 90 1 LB misfit\nEOF\n";

$xipdf =  "$datadir/${eid}_xi_homo_pdf.dat";
$xipostpdf =  "$datadir/${eid}_xi_post_pdf.dat";
print CSH "psbasemap $J1 -R0/120/0/4 -Ba30/a1wEsn -Xa0c -Ya4c -K -O -V>>$psfile\n";
print CSH "psxy $xipdf $J1 -R0/$max/0/4 -Ba30/a1wEsn -Xa0c -Ya4c -W1p,0/0/255 -K -O -V>>$psfile\n";
print CSH "psxy $xipostpdf $J1 -R0/$max/0/4 -Ba30/a1wEsn -Xa0c -Ya4c -W1p,0/255/0 -K -O -V>>$psfile\n";  # 120 degree = 2.0944 radians
print CSH "pstext -JM7i -R0/1/0/1 -Xa-1c -Ya4.5c -G0/255/0 -K -N -O -V>> $psfile<<EOF\n0 0 10 90 1 LB \@~s\@~(m)\nEOF\n";
print CSH "pstext -JM7i -R0/1/0/1 -Xa-1c -Ya5.5c -G0/0/255 -K -N -O -V>> $psfile<<EOF\n0 0 10 90 1 LB \@~m\@~(m)\nEOF\n";

print CSH "psbasemap $J1 -R0/120/0/6 -Ba30/a1WesN -Xa0c -Ya7c -K -O -V>>$psfile\n";
print CSH "psxy $xierrdis $J1 -R0/$max/0/6 -Ba30/a1Wesn -Xa0c -Ya7c -W1p,255/0/0 -K -O -V>>$psfile\n";
#print CSH "psxy -JX6c/4c -R0/3.1416/0/4 -Ba45/a1Wesn -Xa0c -Ya7c -W1p,255/0/0 -K -O -V>>$psfile<<EOF\n0 5\n180 5\nEOF\n";
print CSH "pstext -JM7i -R0/1/0/1 -Xa-1c -Ya8.5c -K -N -O -V>> $psfile<<EOF\n0 0 10 90 1 LB L(m)\nEOF\n";

$XIerr =  "$datadir/${eid}_XI.dat";
open(IN,$XIerr); @XIval = <IN>; 
($XI) = split(" ",$XIval[0]);
print CSH "pstext -JM7i -R0/1/0/1 -Xa4c -Ya9c -W255/255/255,o/5p -K -N -O -V>> $psfile<<EOF\n0 0 12 0 1 RB \@~P\@~ = $XI\@+o\@+\nEOF\n";
#print CSH "pstext $OMEGAerr -JM7i -R0/1/0/1 -Xa0c -Ya12c -K -N -O -V>> $psfile \n";

print CSH "psbasemap -JX5c/10c -R0/120/0/4 -K -O -V>>$psfile\n";
print CSH "psxy -JX5c/10c -R0/120/0/4 -W1p,0/0/0,- -K -O -V >>$psfile<<EOF\n$XI 0\n$XI 4\nEOF\n";
print CSH "pstext -JM7i -R0/1/0/1 -Xa-0.75c -Ya10.5c -K -N -O -V>> $psfile<<EOF\n0 0 12 0 1 LB (b)\nEOF\n";
}

print CSH "pstext -JM7i -R0/1/0/1 $Xr $Yr -K -N -O -V>> $psfile<<EOF\n0 0 14 0 1 LB \nEOF\n";
#========================================================================================

#Plot the P-T axis on the beachball
$pxsfile="$datadir/${eid}_pxs.dat";
open(IN,$pxsfile); @psf = <IN>; 
($pxs,$pys)=split(" ",$psf[0]);
$txsfile="$datadir/${eid}_txs.dat"; 
open(IN,$txsfile); @tsf = <IN>; 
($txs,$tys)=split(" ",$tsf[0]);

$J = "-JX4.5c/4.5c";
$R = "-R-1/1/-1/1";
$X = "-X1c";
$Y = "-Y19.5c";
$Xr = "-X-1c";
$Yr = "-Y-19.5c";

#print CSH "psbasemap -JX3c/3c -R-1/1/-1/1 -Ba30/a1wesn -X9c -Y22c -K -O -V>>$psfile\n";
print CSH "psmeca $J $R -Sa4.5c -G150/150/150 -L1p/0/0/0 -N $X $Y -K -O -V >>$psfile<<EOF\n0 0 $edep $strike $dip $rake 5 0 0\nEOF\n";
print CSH "psxy $pxsfile $J $R -Sc1p -G0/0/255 -K -O -V>>$psfile\n";
print CSH "psxy $txsfile $J $R -Sc1p -G255/0/0 -K -O -V>>$psfile\n";
print CSH "psxy $J $R -Sc5p -W1p/0/0/0 -G150/150/255 -K -N -O -V>>$psfile<<EOF\n$pxs $pys\nEOF\n";
print CSH "psxy $J $R -Sc5p -W1p/0/0/0 -G255/150/150 -K -N -O -V>>$psfile<<EOF\n$txs $tys\nEOF\n";
print CSH "pstext -JM7i -R0/1/0/1 -Ya4.75c -K -N -O -V>> $psfile<<EOF\n0 0 12 0 1 LB (a)\nEOF\n";

print CSH "pstext -JM7i -R0/1/0/1 $Xr $Yr -K -N -O -V>> $psfile<<EOF\n0 0 14 0 1 LB \nEOF\n";

#========================================================================================
# Always use
print CSH "pstext -JM7i -R0/1/0/1 $Xr $Yr -N -O -V>> $psfile<<EOF\n0 0 14 0 1 LB \nEOF\n";

close (CSH);
system("csh -f $cshfile");
#system("gv $psfile &");
