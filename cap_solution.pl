#!/usr/bin/perl -w

if (@ARGV < 1) {die("Usage: map_vipul eid\n")}
($eid) = @ARGV; 

#==========Read files=================
# read in event info
$capevent = "/home/vipul/gmt/data/${eid}_event.dat";
if (not -f $capevent) {die("Check if capevent $capevent exist or not\n");}
(undef,$elat,$elon,$edep)  = split(" ",`grep $eid $capevent`);

# read station info
$capstation = "/home/vipul/gmt/data/${eid}_station.dat";  
if (not -f $capstation) {die("Missing File: Couldn't find $capstation in the running directory");}
open(IN,$capstation); @capstationlines = <IN>; $nsta = @capstationlines;

# read output file
$outputfile = "/home/vipul/gmt/data/${eid}.out";
if (not -f $outputfile) {die("Missing File: Output file $outputfile missing");}
open(IN,$outputfile); @outputlines = <IN>; 
(undef,undef,undef,$smodeldep,undef,$strike,$dip,$rake,undef,$Mw,undef,$rms,$val1,undef,$sigstrike,$sigdip,$sigrake,undef,$iso1,$iso2,undef,$clvd1,$clvd2) = split(" ",$outputlines[0]);
($model,$depth) = split("_",$smodeldep);

$histogramfile = "/home/vipul/gmt/data/${eid}_histo.dat";
if (not -f $histogramfile) {die("Missing File: Histogram file $histogramfile missing");}

$samplesfile = "/home/vipul/gmt/data/${eid}_samples.dat";
#$samplesfile = "/home/vipul/gmt/data/samples.dat";
if (not -f $samplesfile) {die("Missing File: Moment tensor file $samplesfile missing");}
open(IN,$samplesfile); @sampleslines = <IN>; $nsamples = @sampleslines;
#($omega,$misfit,$stk_samples,$dip_samples,$rak_samples) = split(" ",$sampleslines[0]);

#==========I/O Files=================
$cptfile = "map_alaska.cpt";   # color palette table file
$cshfile = "cap_solution.csh";   # color shape file
$psfile = "${eid}_solution.ps";     # output postscript file
$grdfile = "/home/admin/share/datalib/topography/GLOBAL/ETOPO1/ETOPO1_Bed_g.grd";
$gradfile = "/home/admin/share/datalib/topography/GLOBAL/ETOPO1/ETOPO1_Bed_g.grad";
$faultfile = "/home/admin/share/datalib/faults/alaska/DGGS2012/ak_fault_database_faults.gmt";
$histogramfile = "/home/vipul/gmt/data/${eid}_histo.dat";

# ==========Plotting specifications=====
$xmin = -155;
$xmax = -145;
$ymin = 56;
$ymax = 64;
$scale = 9000000;
$xcen = ($xmin + $xmax)/2;
$ycen = ($ymin + $ymax)/2;
$J = "-Jb$xcen/$ycen/$ymin/$ymax/1:$scale";
$R = "-R$xmin/$xmax/$ymin/$ymax";
#$B = "-Ba4f2:.Southern_Alaska:WesN";
$B = "-Ba4f2WesN";
$X = "-X14c";
$Y = "-Y15c";
$Xr = "-X-14c";
$Yr = "-Y-15c";
$title = "Southern Alaska (eid: $eid)";

$stainfoa= "-Si5p -W3p,0/0/0";
$stainfob = "-Si5p -W1p,255/255/255";

# $alaska = "-175/-135/50/70";
# $Ra = "-R$alaska";
# $Ja = "-Jb$region/0.5c";
# $Ba = "-Ba10f5";

open(CSH,">$cshfile");

# ========GMT commands===========

print CSH "gmtset PAPER_MEDIA letter MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 0.2c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT 1 ANOT_FONT 0.5 LABEL_FONT 1 HEADER_FONT_SIZE 15 FRAME_PEN 1p TICK_PEN 1p PAGE_ORIENTATION portrait\n";
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
# plot coastline
print CSH "pscoast $R $J -A0 -Ia/0/0/255 -Df -W -N1 -O -K -V >> $psfile\n";
# plot faultline
print CSH "psxy $faultfile $J $R -m -W/255/0/0 -V -K -O >> $psfile\n";
# print CSH "psxy $J $R -W1p,255/255/255 -Si10p -V -K -O >> $psfile<<EOF\n-150 60\nEOF\n";
#print CSH "pstext $J $R -G0/0/0 -W255/255/255 -O -V -K>> $psfile<<EOF\n$elon $elat  14 0 1 CM $eid\nEOF\n";
print CSH "pstext -JM7i -R0/1/0/1 -N -Xa7c -Ya12c -O -K -V>> $psfile<<EOF\n0 0 14 0 1 RM $title\nEOF\n";
print CSH "pstext -JM7i -R0/1/0/1 -N -Xa7c -Ya11.5c -O -K -V>> $psfile<<EOF\n0 0 10 0 1 RM Strike: $strike Dip: $dip Rake: $rake Mw: $Mw Depth: $depth km\nEOF\n";
# scale bar (????)
# print CSH "pscoast $J $R -L-152/57.5/57.5/100+p0.5p,0/0/0,solid+f255/255/255 -N1 -K -O -V >> $psfile";

#=================== plot stations=======================
for ($i = 1; $i <= $nsta; $i++) {
   ($sta,$net,$lat,$lon,$dist,$az) = split(" ",$capstationlines[$i-1]);
   # plot an open triangle at each station in the output file
   print CSH "psxy $J $R $stainfoa -K -O -V >>$psfile<<EOF\n $lon $lat\nEOF\n";
   print CSH "psxy $J $R $stainfob -K -O -V >>$psfile<<EOF\n $lon $lat\nEOF\n";
 }

# plot moment tensor beachball
print CSH "psmeca $J $R -Sa0.2 -G0/0/255 -L1p/0/0/0 -K -O -V >>$psfile<<EOF\n$elon $elat $edep $strike $dip $rake $Mw 0 0\nEOF\n";

#print CSH "psscale -C$cptfile -D3.5c/-1c/6c/0.25ch -O -Ac -B1500:Elevation:/:m: -E >> $psfile\n"; # relative to the basemap

# ========Specifications for inset========
$xmina = -175;
$xmaxa = -135;
$ymina = 52;
$ymaxa = 72;
$scalea = 50000000;
$xcena = ($xmina + $xmaxa)/2;
$ycena = ($ymina + $ymaxa)/2;
$Ja = "-Jb$xcena/$ycena/$ymina/$ymaxa/1:$scalea";
$Ra = "-R$xmina/$xmaxa/$ymina/$ymaxa";
$Ba = " -B100wesn";
$Xa = "-Xa2c";
$Ya = "-Ya-1c";
$gridfilea = "alaska_grid.grd";

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


# Origin Marker
print CSH "pstext -JM7i -R0/1/0/1 -N $Xr $Yr -O -K -V>> $psfile<<EOF\n0 0 14 0 1 LB vipul\nEOF\n";

# ====================plotting histograms====================
print CSH "pshistogram $histogramfile -JX3c/3c -R0/360/0/20 -G0/255/0 -L0.25p,0/0/0 -Z1 -W10  -Ba90:strike:/a5::WeSn -Xa1.5c -Ya24.3c -T0 -K -O -V >>$psfile\n";
print CSH "pshistogram $histogramfile -JX3c/3c -R0/90/0/50 -G0/255/0 -L0.25p,0/0/0 -Z1 -W10  -Ba20:dip:/a10::WeSn -Xa1.5c -Ya19.4c -T1 -K -O -V >>$psfile\n";
print CSH "pshistogram $histogramfile -JX3c/3c -R-90/90/0/50 -G0/255/0 -L0.25p,0/0/0 -Z1 -W10  -Ba30:rake:/a10::WeSn -Xa1.5c -Ya14.5c -T2 -K -O -V >>$psfile\n";

#============================ plot beachballs===================
$Xm = "-X0.5c";
$Ym = "-Y8.5c";
$Xr = "-X-0.5c";
$Yr = "-Y-8.5c";
$Bm = "-BWeSn";
$Jm = "-JX20c/4c";
$Rm = "-R0/20/0/4";

print CSH "makecpt -Cpanoply -T0/16/.01 -D -Z> $cptfile \n";
#print CSH "makecpt -Cwysiwyg -T0/20/1 -D -Z > $cptfile \n";
#print CSH "gmtset BASEMAP_FRAME_RGB 255/0/0 \n";
print CSH "psbasemap $Jm $Rm $Bm $Xm $Ym -K -O -V>>$psfile\n";
print CSH "gmtset BASEMAP_FRAME_RGB 0/0/0 \n";
#print CSH "psmeca $Jm $Rm -Sa0.5 -G0/0/255 -L1p/0/0/0 -K -O -V >>$psfile<<EOF\n1 3 $edep $strike $dip $rake $Mw 0 0 0\nEOF\n";

for ($i = 1; $i <= $nsamples; $i++) {
  $xloc = 2*$i-1;
  ($omega,$misfit,undef,$stk_sample,$dip_sample,$rak_sample) = split(" ",$sampleslines[$i-1]);
  if ($i <= 10) {
    print CSH "psmeca $Jm $Rm -Sa0.5 -G0/0/255 -L1p/0/0/0 -Z$cptfile -K -O -V >>$psfile<<EOF\n$xloc 3 $misfit $stk_sample $dip_sample $rak_sample $Mw 0 0 $omega\nEOF\n";
  } else {
    $xloc = 2*($i-10)-1;
    print CSH "psmeca $Jm $Rm -Sa0.5 -G0/0/255 -L1p/0/0/0 -Z$cptfile -K -O -V >>$psfile<<EOF\n$xloc 1 $misfit $stk_sample $dip_sample $rak_sample $Mw 0 0 $omega\nEOF\n";
  }
}

print CSH "psscale -C$cptfile -D19.75c/2c/4c/0.25c -K -O -Ac -B2::/:err: -Ef >> $psfile\n"; # relative to the basemap
print CSH "pstext -JM7i -R0/1/0/1 -N $Xr $Yr -K -O -V>> $psfile<<EOF\n0 0 14 0 1 LB \nEOF\n";

#===================== plot XYZ =====================================
$xyzfile = "/home/vipul/gmt/data/319605_stkfile.dat";
#$stkgrd = "/home/vipul/gmt/data/319605_stkgrd.grd";
#print CSH "makecpt -Cwysiwyg -T0/20/1 -D -Z > $cptfile \n";
print CSH "makecpt -Cpanoply -T0/16/.01 -D -Z> $cptfile \n";
print CSH "psbasemap -JX5c/5c -R-90/90/0/90 -Ba45WeSn -Xa2c -Ya2.5c  -K -O -V>>$psfile\n";
print CSH "psxy $xyzfile -JX5c/5c -R-90/90/0/90 -Ba45:rake:/a45:dip:WeSn -C$cptfile -Ss8p  -Xa2c -Ya2.5c  -K -O -V>>$psfile\n";
print CSH "psxy -JX5c/5c -R-90/90/0/90 -Ba45:rake:/a45:dip:WeSn -Sx8p -W255/255/255 -Xa2c -Ya2.5c -K -O -V>>$psfile<<EOF\n$rake $dip\nEOF\n";

#print CSH "pscontour $xyzfile -JX5c/5c -R-90/90/0/90 -Ba45:rake:/a45:dip:WeSn -C$cptfile -Lthinnest -I -Xa2c -Ya2.5c -W0.5p -K -O -V>>$psfile\n";
# print CSH "xyz2grd $xyzfile -G$stkgrd -I20+/10+ -R-90/90/10/90 -V\n";
# print CSH "grdimage $stkgrd -C$cptfile -R-90/90/0/90 -JX5c/5c -Ba45WeSn -Xa2c -Ya2.5c -K -O -V>>$psfile\n";

$xyzfile = "/home/vipul/gmt/data/319605_dipfile.dat";
#$dipgrd = "/home/vipul/gmt/data/319605_dipgrd.grd";
print CSH "psbasemap -JX5c/5c -R-90/90/0/360 -Ba45WeSn -Xa9c -Ya2.5c  -K -O -V>>$psfile\n";
print CSH "psxy $xyzfile -JX5c/5c -R-90/90/0/360 -Ba45:rake:/a45:strike:WeSn -C$cptfile -Ss5p  -Xa9c -Ya2.5c -K -O -V>>$psfile\n";
print CSH "psxy -JX5c/5c -R-90/90/0/360 -Ba45:rake:/a45:strike:WeSn -Sx8p -W255/255/255 -Xa9c -Ya2.5c -K -O -V>>$psfile<<EOF\n$rake $strike\nEOF\n";

$xyzfile = "/home/vipul/gmt/data/319605_rakfile.dat";
print CSH "psbasemap -JX5c/5c -R0/360/0/90 -Ba45WeSn -Xa16c -Ya2.5c  -K -O -V>>$psfile\n";
print CSH "psxy $xyzfile -JX5c/5c -R0/360/0/90 -Ba45:strike:/a45:dip:WeSn -C$cptfile -Ss8p  -Xa16c -Ya2.5c -K -O -V>>$psfile\n";
print CSH "psxy -JX5c/5c -R0/360/0/90 -Ba45:strike:/a45:dip:WeSn -Sx8p -W255/255/255 -Xa16c -Ya2.5c -K -O -V>>$psfile<<EOF\n$strike $dip\nEOF\n";

# ========================omega distribution==========================
$omegaerr = "/home/vipul/gmt/data/${eid}_omega_err.dat";
$omegaerrpost = "/home/vipul/gmt/data/${eid}_omega_err_post.dat";
$omegaerrdis = "/home/vipul/gmt/data/${eid}_omega_err_dis.dat";
$Xo = "-X6.5c";
$Yo = "-Y14.5c";
$Xr = "-X-6.5c";
$Yr = "-Y-14.5c";
$Bo = "-Ba45:omega:/a2:error:WeSn";
$Jo = "-JX6c/4c";
$Ro = "-R0/180/0/16";
print CSH "psbasemap $Jo $Ro $Bo $Xo $Yo -K -O -V>>$psfile\n";
#print CSH "psxy $omegaerr $Jo $Ro $Bo -Sc1p -C$cptfile $Xo $Yo -K -O -V>>$psfile\n";
print CSH "psxy $omegaerr $Jo $Ro $Bo -Sc1p -G0/0/255 -K -O -V>>$psfile\n";
print CSH "psxy $omegaerrpost $Jo $Ro $Bo -Sc2p -G0/255/0 -K -O -V>>$psfile\n";
print CSH "psxy $samplesfile $Jo $Ro $Bo -Sc2p -G255/0/0 -K -O -V>>$psfile\n";
print CSH "psscale -C$cptfile -D6c/2c/4c/0.25c -K -O -Ac -B0::/:: -Ef >> $psfile\n";

#print CSH "psbasemap -JX6c/3c -R0/180/0/30 -Bf45/a10wEsN -Xa0c -Ya4c -K -O -V>>$psfile\n";
#print CSH "pshistogram $omegaerr -JX6c/3c -R0/180/0/30 -G0/0/255 -L0.25p,0/0/0 -Z1 -W5 -B00wesn -Xa0c -Ya4c -T0 -K -O -V >>$psfile\n";
#print CSH "pshistogram $omegaerrpost -JX6c/3c -R0/180/0/30 -G0/255/0 -L0.25p,0/0/0 -Z1 -W5 -B00wesn -Xa0c -Ya4c -T0 -K -O -V >>$psfile\n";
$omegapdf =  "/home/vipul/gmt/data/${eid}_omega_homo_pdf.dat";
$omegapostpdf =  "/home/vipul/gmt/data/${eid}_omega_post_pdf.dat";
print CSH "psbasemap -JX6c/3c -R0/3.1416/0/3 -Ba1/a1wEsn -Xa0c -Ya4c -K -O -V>>$psfile\n";
print CSH "psxy $omegapdf -JX6c/3c -R0/3.1416/0/3 -Ba1/a1wEsn -Xa0c -Ya4c -W1p,0/0/255 -K -O -V>>$psfile\n";
print CSH "psxy $omegapostpdf -JX6c/3c -R0/3.1416/0/3 -Ba1/a1wEsn -Xa0c -Ya4c -W1p,0/255/0 -K -O -V>>$psfile\n";

print CSH "psbasemap -JX6c/3c -R0/3.1416/0/3 -Ba1/a1WeN -Xa0c -Ya7c -K -O -V>>$psfile\n";
print CSH "psxy $omegaerrdis -JX6c/3c -R0/3.1416/0/3 -Ba1/a1WeN -Xa0c -Ya7c -W1p,0/255/255 -K -O -V>>$psfile\n";
print CSH "psxy -JX6c/3c -R0/3.1416/0/3 -Ba1/a1WeN -Xa0c -Ya7c -W1p,255/0/0 -K -O -V>>$psfile<<EOF\n0 5\n180 5\nEOF\n";

print CSH "pstext -JM7i -R0/1/0/1 $Xr $Yr -K -N -O -V>> $psfile<<EOF\n0 0 14 0 1 LB \nEOF\n";
#========================================================================================
# Always use
print CSH "pstext -JM7i -R0/1/0/1 $Xr $Yr -N -O -V>> $psfile<<EOF\n0 0 14 0 1 LB \nEOF\n";

close (CSH);
system("csh -f $cshfile");
system("gv $psfile &");
