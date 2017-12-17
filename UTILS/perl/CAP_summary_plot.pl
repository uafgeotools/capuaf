#!/usr/bin/perl -w
# 
# This script makes the summary plots for CAP uncertainty analysis 
# Uncertainty analysis is done in MATLAB (see CAP_unc.m, run_CAP_unc.m and MTbrick_Csection.m)
#
# INPUT: Set the $datadir (where files from run_CAP_unc.m are generated) first
#        Set the $output dir too (default is ./)
#
# Usage: CAP_summary_plot.pl eid OR 
#        CAP_summary_plot.pl eid norm weight
# Example: CAP_summary_plot.pl 20070911234634153  OR
#          CAP_summary_plot.pl 20070911234634153 L1 M111
#          CAP_summary_plot.pl 20090407201255351 L1 M111  # mtvipul paper example
# 
# OUTPUT: ./eid_norm_weight.ps
# 
# Parent script: cap_sol.pl
# Other sister scipts: unc_theory.pl (for generating explanantion figure of P_av estimation);
#                      cross.pl      (cross-section in  strike-dip-rake space)
# Vipul Silwal
# (vipulsilwal@gmail.com)
# 
# 5 Aug, 2015
# 
# 

($eid,$norm,$wts) = @ARGV;

#if (@ARGV < 1) {die("Usage: cap_sol eid\n")}

# Default parameters (See SilwalTape2016 for details of moment tensor naming convention)
# We may not need this convention in future 
# Use better output name
if (@ARGV eq 1) {
    print "-------------------------------------\n";
    print "norm and weight not specified\n Using default: L1 norm and M111\n"; 
    print "-------------------------------------\n";
    $norm = "L1"; $wts = "M111"; # these are only used for using proper input and output directory; and output figure name
}

#==========set dat directory info========
$label="MOOS";
#$label="MISC";
#$label="NW";
#$label='test';

#==========Read files=================
# Change these to where your data file are located
#$datadir = "/home/vipul/CAP/inv/scak/20090407201255351_gmt_data";
#$datadir = "/home/vipul/gmt/data/cap/$label/${eid}/$norm/$wts";
#$datadir = "/home/vipul/REPOSITORIES/cap/EXAMPLES/RESULTS_20090407201255351/gmt/data"; # For the EXAMPLE run
#$datadir = "/home/vipul/gmt/data/cap/test"; # for testing
#$datadir = "/home/vipul/CAP/inv/scak/20090407201255351_gmt_data"; # for testing mtvipul example
#$datadir = "/home/vipul/REPOSITORIES/cap/20090407201253480_gmt_data";
#$Nst=2;
#$datadir = "/home/vipul/gmt/data/cap/test_stn${Nst}"; # for testing effect of station
#$datadir = "/home/vipul/gmt/data/cap/test_stn${Nst}_SH"; # for testing effect of adding SH waves (one station only)
#$datadir = "/home/vipul/gmt/data/cap/test_stn${Nst}_PT"; # for testing effect of fixing PT axis (two stations only)
#---celso thesis
#$datadir = '/home/vipul/celso_thesis/20090407201255351/gmt_data/';
#$datadir = '/home/vipul/celso_thesis/Little_Skull_Main/gmt_data/'; $iweird = 1;
#$datadir = '/home/vipul/CAP/inv/scak/MOOS_updated/waveforms_alaska/20070911234634153/gmt_data/'; # TEST CASE
#$datadir = '/home/vipul/CAP/inv/scak/20090407201253480_save/20090407201253480_gmt_data_51/';
#$datadir = '/home/vipul/CAP/inv/scak/NEHRP/beluga/20080205035142446/gmt_data/';
#$datadir = '/home/vipul/CAP/inv/scak/NEHRP/beluga/20080126042942584/gmt_data/';
$datadir = './gmt_data/';

# New way of pointing to the datadir - seems much more reasonable
#$datadir = '/home/vipul/CAP/inv/scak/MOOS_updated/waveforms_alaska';
#$datadir = $datadir.'/'.$eid.'/gmt_data/';

$dlib = "/home/admin/share/datalib";

if (@ARGV eq 2) {
    $datadir = ${norm}; # Assume that the second input parameter is the $datadir
    $norm = "L1"; $wts = "M111"; #  default 
}

print "$datadir";

# Options (Input parameter)
$iinset = 0; # Plot inset
$idenali = 0; # for events around Denali
$half_page = 0; # half-page version
$islab_extent = 1;  # Plot slab extent (informative in case of southern Alaska events)
# if you want to plot only some of the subfigures

$figA = 1;
$figB = 1;
$figC = 1;
$figG = 1;
$figH = 1;
$figI = 1;
if ($half_page) {
    $figB = 0;
    $figG = 0;
    $figI = 0;}
$ORIGIN="";
$dmisfit = 4; # for misfit color scale and (misfit vs omega plot yTicks)

#==========Parameters=================

$ballinfo = "-G200 -L1p,0/0/0";
$paxiscolor = "0";    # 255/150/150
$taxiscolor = "0";    # 150/150/255
$axisinfo = "-a0.3c/cc -e$paxiscolor -p1p,255 -g$taxiscolor -t1p,255";
$tickinc = 30;

# read station info
$capstation = "$datadir/${eid}_station_list_ALL.dat";
if (not -e $capstation) {die("Missing File: Couldn't find $capstation in the running directory");}

open(IN,$capstation); @capstationlines = <IN>; $nsta = @capstationlines;

# read output file
$capout = "$datadir/${eid}.out";
if (not -f $capout) {die("Missing File: Output file $capout missing");}
open(IN,$capout); @linescap = <IN>; 
(undef,undef,undef,$smodeldep,undef,$strike,$dip,$rake,undef,$Mw,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef) = split(" ",$linescap[0]);
(undef,undef,undef,$elat,undef,$elon,undef,$edep)  =split(" ",$linescap[1]);
($evid,$smodel,$depth) = split("_",$smodeldep);
if ($iweird == 1) { # only for LSM earthquake
$smodel = 'wes';
$depth = 10;}
print "$dip $depth";
$Mw=sprintf("%.1f",$Mw);
$depth=sprintf("%d",$depth);
@capcols = (14,22,13,17,21,6,5,9);  # column of data to plot
@wcols = (11,19,11,15,19,3,3,7);    # columns with weight factors
# integer values for labels
$fstrike=sprintf("%d",$strike);
$fdip=sprintf("%d",$dip);
$frake=sprintf("%d",$rake);

#==========I/O Files=================== (OR PERHAPS DEFINE FILES WHERE THEY ARE USED)
$outputdir = "./";
#$outputdir = "OUTPUT/${eid}";
#$outputdir = "./NW/";
system("mkdir -p $outputdir");
$cptfile = "map_alaska.cpt";   # color palette table file
$cshfile = "cap_solution.csh";   # color shape file
$psfile = "./$outputdir/${eid}_${norm}_${wts}.ps";     # output postscript file
$grdfile = "/home/admin/share/datalib/topography/GLOBAL/ETOPO1/ETOPO1_Bed_g.grd";
$gradfile = "/home/admin/share/datalib/topography/GLOBAL/ETOPO1/ETOPO1_Bed_g.grad";
$faultfile = "/home/admin/share/datalib/faults/alaska/DGGS2012/ak_fault_database_faults.gmt";

#=================================================
open(CSH,">$cshfile");

# ========GMT commands============================
print CSH "gmtset PAPER_MEDIA letter MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 0.2c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT 1 ANOT_FONT 0.5 LABEL_FONT 1 LABEL_OFFSET 0.1c HEADER_FONT_SIZE 15 FRAME_PEN 1p TICK_PEN 0.5p PAGE_ORIENTATION portrait CHAR_ENCODING ISOLatin1+\n";

# Always use (shift orgin if needed - moves all the figures)
print CSH "pstext -JM7i -R0/1/0/1 -X-0.5c -Y-0.5c -K> $psfile<<EOF\n0 0 14 0 1 LB $ORIGIN\nEOF\n";

#===========FIGURE A=============================
# The map
if ($figA){
    # search ifkmod below to specify additional parameters for each region
    if($smodel eq "scak") {
	$ifkmod = 1;
	print "\nifkmod = $ifkmod: southern Alaska\n";
	$xtick1 = 2; $ytick1 = 1; $xtick2 = 1; $ytick2 = 0.5;
    }
    if($smodel eq "tactmod") {
	$ifkmod = 2;
	print "\nifkmod = $ifkmod: central Alaska\n";
	$xtick1 = 1; $ytick1 = 1; $xtick2 = 0.5; $ytick2 = $xtick2;
    }
    if($label eq "NW") {
	$ifkmod = 3;
	print "\nifkmod = $ifkmod: central Alaska\n";
	$xtick1 = 4; $ytick1 = 1; $xtick2 = 2; $ytick2 = 0.5;
    }
    if($smodel eq "wes") {
	$ifkmod = 4;
	print "\nifkmod = $ifkmod: central Alaska\n";
	$xtick1 = 4; $ytick1 = 1; $xtick2 = 2; $ytick2 = 0.5;
    }
    
    #--------------becase different sized lat-lon box (interior or southern Alaska) requires different paramters
    if ($ifkmod==1) { # SOUTHERN ALASKA
	$xmin = -155.0; $xmax = -145.0; $ymin = 58; $ymax = 65.0;
	$sbarinfo = "-L-147/58.5/58.5/100+p1.5p,0/0/0,solid+f255/255/255";
	$Jbscale = 8000000;
	$Y = "-Y16c";
	$Yr = "-Y-16c";
    }
    if ($ifkmod==2) {   # INTERIOR ALASKA (NENANA)
	$xmin = -152.5; $xmax = -143.5; $ymin = 60.9; $ymax = 67.1;
	$sbarinfo = "-L-146/66.5/63.3/100+p1.5p,0/0/0,solid+f255/255/255";
	$xtick1 = 2; $ytick1 = 1; $xtick2 = 1; $ytick2 = 0.5;
	$Jbscale = 7000000;
	$Y = "-Y16c";
	$Yr = "-Y-16c";
    }
    if ($idenali==1) {  # INTERIOR ALASKA (DENALI) - USE ONLY WHEN NEEDED
	$xmin = -154.5; $xmax = -144; $ymin = 58; $ymax = 66; $Jbscale = 8500000;
	$sbarinfo = "-L-151.8/58.7/58.7/200+p1.5p,0/0/0,solid+f255/255/255";
	$Y = "-Y16c";
	$Yr = "-Y-16c";
	$title0 = "Denali, April 16, 2014";
    }
    if ($ifkmod==3) {  # NORTH_WESTERN ALASKA - USE ONLY WHEN NEEDED
	$xmin = -164; $xmax = -147; $ymin = 59.5; $ymax = 69; $Jbscale = 13000000;
	$sbarinfo = "-L-160/61/61/200+p1.5p,0/0/0,solid+f255/255/255";
	$xtick1 = 4; $ytick1 = 1; $xtick2 = 1; $ytick2 = 0.5;
	$Y = "-Y18c";
	$Yr = "-Y-18c";
	$inoatak = 1;
    }
    if ($ifkmod==4) {  # NORTH_WESTERN ALASKA - USE ONLY WHEN NEEDED
	$xmin = -125; $xmax = -115; $ymin = 30; $ymax = 40; $Jbscale = 13000000;
	$sbarinfo = "-L-160/61/61/200+p1.5p,0/0/0,solid+f255/255/255";
	$xtick1 = 2; $ytick1 = 2; $xtick2 = 1; $ytick2 = 1;
	$Y = "-Y18c";
	$Yr = "-Y-18c";
	$inoatak = 1;
    }

    #------------------------------------------------------
    $xcen = ($xmin + $xmax)/2;
    $ycen = ($ymin + $ymax)/2;
    $J = "-Jb$xcen/$ycen/$ymin/$ymax/1:$Jbscale";
    $R = "-R$xmin/$xmax/$ymin/$ymax";
    #$B = "-Ba4f2:.Southern_Alaska:WesN";
    #$B = "-Ba4f2WesN";
    $B = sprintf("-Ba%3.3ff%3.3fd:\" \":/a%3.3ff%3.3fd:\" \"::.\" \":WesN",$xtick1,$xtick2,$ytick1,$ytick2);
    $xshift = 1.5;
    $X = "-X${xshift}c";
    $Xr = "-X-${xshift}c";
    $stainfoa= "-Si5p -W3p,0/0/0";
    $stainfob = "-Si5p -W1p,255/255/255";
    
    # Position for the figure
    print CSH "pstext -JM7i -R0/1/0/1 -N $X $Y -O -K -V>> $psfile<<EOF\n0 0 14 0 1 LB \nEOF\n";
    # cut section of grid
    #print CSH "grdcut $grdfile -G$gridfile $R -V\n";
    # make color pallete for plotting map
    print CSH "makecpt -Cglobe -T-3000/3000/200 -D -Z > $cptfile \n";
    #print CSH "grd2cpt $grdfile -Cglobe -Z >> $cptfile \n";
    # create basemap
    #print CSH "psbasemap $R $J −Ba2f1d:" ":/a2f1d:" ":WeSn -K -V -X3c  -Y3c > $psfile \n";
    print CSH "psbasemap $R $J $B -K -O -V >> $psfile \n";
    # plot topography
    print CSH "grdimage $grdfile $R $J $B -E100 -C$cptfile -I$gradfile -V -K -O >> $psfile \n";
    # plot faultline
    print CSH "psxy $faultfile $J $R -m -W1p,255/0/0 -V -K -O >> $psfile\n";
    # plot coastline
    # print CSH "pscoast $R $J -A0 -Df -W -N1 -L-150/59.5/64/200+p1p,0/0/0,solid+f255/255/255 -O -K -V >> $psfile\n";
    print CSH "pscoast $R $J -A10/0/1 -Df -W -N1 $sbarinfo -O -K -V >> $psfile\n";
    # print CSH "psxy $J $R -W1p,255/255/255 -Si10p -V -K -O >> $psfile<<EOF\n-150 60\nEOF\n";
    # print CSH "pstext $J $R -G0/0/0 -W255/255/255 -O -V -K>> $psfile<<EOF\n$elon $elat  14 0 1 CM $eid\nEOF\n";
    $xtag = $xshift - 2.25;
    if ($half_page==0){
    if ($inoatak){
	print CSH "pstext -JM7i -R0/1/0/1 -Xa${xtag}c -Ya9c -K -N -O -V>> $psfile<<EOF\n0 0 10 0 1 LB (a)\nEOF\n";}
    else{
	print CSH "pstext -JM7i -R0/1/0/1 -Xa${xtag}c -Ya11c -K -N -O -V>> $psfile<<EOF\n0 0 10 0 1 LB (a)\nEOF\n";}
}
    # scale bar (????)
    # print CSH "pscoast $J $R -L-152/57.5/57.5/100+p0.5p,0/0/0,solid+f255/255/255 -N1 -K -O -V >> $psfile";
    
    if ($islab_extent==1) {
    $slabextent = "$dlib/seismicity/regional/alaska/slabextent/slab_extent_all.lonlat";  # with Wrangell
    if (not -f $slabextent) {
	die("Check if slabextent $slabextent exist or not\n");
    }
    print CSH "psxy $slabextent $J $R -m -W2p,0/0/255 -K -V -O >> $psfile\n";
}
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
    # note: magnitude is fixed for plotting
    print CSH "psmeca $J $R -Sa0.2 -G150/150/150 -L0.5p,0/0/0 -K -O -V >>$psfile<<EOF\n$elon $elat $edep $strike $dip $rake 5 0 0\nEOF\n";
    
    # Plot cities
    print CSH "pstext $R $J -G0/0/0 -S2p,255/255/255 -K -N -O -V>> $psfile<<EOF\n-150.02 61.17 12 0 1 LB A\nEOF\n"; # Anchorage
    print CSH "pstext $R $J -G0/0/0 -S2p,255/255/255 -K -N -O -V>> $psfile<<EOF\n-147.87 64.82 12 0 1 LB F\nEOF\n"; # Fairbanks
    # ========Specifications for inset========
    
if ($iinset==1) {
    $xmina = -175;
    $xmaxa = -135;
    $ymina = 52;
    $ymaxa = 72;
    $scalea = 70000000;
    $xcena = ($xmina + $xmaxa)/2;
    $ycena = ($ymina + $ymaxa)/2;
    $Ja = "-Jb$xcena/$ycena/$ymina/$ymaxa/1:$scalea";
    $Ra = "-R$xmina/$xmaxa/$ymina/$ymaxa";
    $Ba = " -B100wesn";
    $Xa = "-Xa3.9c";
    $Ya = "-Ya-2.75c";
    if ($idenali==1){$Ya = "-Ya-1.25c"; }# for denali event
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
}
    
    # Origin Marker (Move origin back)
    print CSH "pstext -JM7i -R0/1/0/1 -N $Xr $Yr -O -K -V>> $psfile<<EOF\n0 0 14 0 1 LB \nEOF\n";
}

#===========FIGURE C=============================
# Uncertainty analysis
if ($figC){
    $xshift = 10.5;
    $yshift = 24;
    $X = "-X${xshift}c";
    $Y = "-Y${yshift}c";
    $Xr = "-X-${xshift}c";
    $Yr = "-Y-${yshift}c";

    # heights of the 3 connected subplots
    $height1 = 2.75;
    $height2 = $height1;
    $height3 = 1.80;
    # absolute y-position of base of each subplot
    $y1 = 0;
    $y2 = $y1 - $height2;
    $y3 = $y2 - $height3 - 0.3;
    $y4 = $y3 - 3.5;

    # Get outline of all posterior points and also find the max and min
    $omegaoutline = "$datadir/${eid}_omega_misfit_outline.dat";
    open(IN,$omegaoutline); @omega_misfit = <IN>;
    $min_misfit = 10000;
    $max_misfit = 0;
    for ($ii=0;$ii<@omega_misfit;$ii++){
	($omega1,$misfit1) = split(" ",$omega_misfit[$ii]);
	if ($misfit1 < $min_misfit){
	    $min_misfit = $misfit1;}
	if ($misfit1 > $max_misfit){
	    $max_misfit = $misfit1;}
    }
    $max_err = sprintf("%1.0f",$max_misfit+1);
    $min_err = sprintf("%1.0f",$min_misfit-1);
    
    # Read file with multiple Mrefs
    $iOtherMrefs = 0;
    if ($iOtherMrefs){
	$Mref_pts = "$datadir/${eid}_Mref_pts.dat";
	open(IN,$Mref_pts); @Mpts = <IN>;
    }
    #------------------- Plot all random samples (just their outline - to reduce file size), and random subset of posterior samples
    $J = "-JX5c/${height1}c";
    $wtick = "a60f${tickinc}";
    #$B = "-B$wtick:'\@~w\@~':/a5::WesN";
    $B = "-B$wtick/a${dmisfit}::WesN";
    $R = "-R0/180/$min_err/$max_err";
    print CSH "psbasemap $J $R $B $X $Y -K -O -V>>$psfile\n";
    # outline of homogeneous samples
    print CSH "psxy $omegaoutline $J $R -B -G0/0/255 -L -W -K -O -V>>$psfile\n"; # homogeneously genertated random samples (blue)
    # ------------------posterior samples
    # Get subset of posterior samples (Green points)
    $omegaerrpost = "$datadir/${eid}_omega_err_post.dat";
    open(IN,$omegaerrpost); @omegaerrpostlines=<IN>;
    $post_samples = @omegaerrpostlines; # subset of posterior samples (400 samples)
    for ($ii=0;$ii<$post_samples;$ii++) {
	($omega[$ii],$misfit[$ii],$pxs[$ii],$pys[$ii],$txs[$ii],$tys[$ii]) = split(" ",$omegaerrpostlines[$ii]);
	print CSH "psxy $J $R -B -Sc2p -G124/252/0 -K -O -V>>$psfile<<EOF\n$omega[$ii] $misfit[$ii]\nEOF\n"; # posterior samples (green)
    }
    if ($iOtherMrefs){ # those red points (Mref for another runs)
	for ($ii=0;$ii<@Mpts;$ii++){
	($Mref_omega, $Mref_misfit, $Mref_label) = split(" ",$Mpts[$ii]);
	print CSH "psxy $J $R -Sc3p -G255/0/255 -Lf -K -O -V>>$psfile<<EOF\n$Mref_omega $Mref_misfit\nEOF\n";
	print CSH "pstext $J $R -Xa0c -Ya0.2c -G255/0/255 -K -N -O -V>> $psfile<<EOF\n$Mref_omega $Mref_misfit 8 0 1 LB $Mref_label\nEOF\n";}
    }
    print CSH "pstext -JM7i -R0/1/0/1 -Xa-1c   -Ya1.c  -K -N -O -V>> $psfile<<EOF\n0 0 10 90 1 LB misfit\nEOF\n";
    if ($half_page==0){print CSH "pstext -JM7i -R0/1/0/1 -Xa-1.3c -Ya3c -K -N -O -V>> $psfile<<EOF\n0 0 10 0 1 LB (c)\nEOF\n";}
    print CSH "pstext -JM7i -R0/1/0/1 -Xa2.35c -Ya3.3c -K -N -O -V>> $psfile<<EOF\n0 0 12 0 1 LB \@~w\@~\nEOF\n";

    #----------------- Plot the smoothed pdf for posterior samples and the mesa
    # QUESTION: How do we know that the axes are equal (assuming radians for the x-axis)?
    $mesapdf = "$datadir/${eid}_mesa_pdf.dat";
    $omegapostpdf = "$datadir/${eid}_post_pdf.dat";
    $J = "-JX5c/${height2}c";
    $B = "-B$wtick/a1::Wes";
    $R = "-R0/180/0/2.3";
    $Y = "-Ya${y2}c";
    $X = "-Xa0c";
    print CSH "psbasemap $J $R $B $X $Y -K -O -V>>$psfile\n";
    #print CSH "psbasemap $J $R -B/a10e $X $Y -K -O -V>>$psfile\n";
    print CSH "psxy $mesapdf $J $R $X $Y -W1p,0/0/255 -K -O -V>>$psfile\n"; # blue
    print CSH "psxy $omegapostpdf $J $R $X $Y -W1p,124/252/0 -K -O -V>>$psfile\n"; # green

    $y1x = 3.3; $y2x = 4.1; $y3x = 5.1;
    $yshift = -5.5;
    $y1m = $y1x + $yshift; $y2m = $y2x + $yshift; $y3m = $y3x + $yshift;
    print CSH "pstext -JM7i -R0/1/0/1 -Xa-1c -Ya${y1m}c -G50/205/50 -K -N -O -V>> $psfile<<EOF\n0 0 8 90 1 LB P'(\@~w\@~)\nEOF\n";
    print CSH "pstext -JM7i -R0/1/0/1 -Xa-1.15c -Ya${y1m}c -G50/205/50 -K -N -O -V>> $psfile<<EOF\n0 0 8 90 1 LB \136 \nEOF\n";
    print CSH "pstext -JM7i -R0/1/0/1 -Xa-1c -Ya${y2m}c -G0/0/255 -K -N -O -V>> $psfile<<EOF\n0 0 8 90 1 LB V'(\@~w\@~)\nEOF\n";
    print CSH "pstext -JM7i -R0/1/0/1 -Xa-1.15c -Ya${y2m}c -G0/0/255 -K -N -O -V>> $psfile<<EOF\n0 0 8 90 1 LB \136 \nEOF\n";
    if ($half_page==0){print CSH "pstext -JM7i -R0/1/0/1 -Xa-1.3c -Ya${y3m}c -K -N -O -V>> $psfile<<EOF\n0 0 10 0 1 LB (d)\nEOF\n";}

    #----------------- Plot the smoothed cdf (cumulative) for posterior samples and the mesa
    # QUESTION: How do we know that the axes are equal (assuming radians for the x-axis)?
    $mesacdf = "$datadir/${eid}_mesa_cdf.dat";
    $omegapostcdf = "$datadir/${eid}_post_cdf.dat";
    $result_file = "$datadir/${eid}_result.dat";
    open(IN,$result_file); @result=<IN>;
    ($P,$vP,$OMEGA) = split(" ",$result[0]);
    $J = "-JX5c/${height3}c";
    $B = "-B$wtick:'\@~w\@~':/a1::Wes";
    $R = "-R0/180/0/1";
    $Y = "-Ya${y3}c";
    $X = "-Xa0c";
    print CSH "psbasemap $J $R $B $X $Y -K -O -V>>$psfile\n";
    print CSH "psbasemap $J $R -Ba200/a1ne $X $Y -K -O -V>>$psfile\n";
    #print CSH "psxy $J $R $X $Y -W0.5p,0/0/0, -K -O -V>>$psfile<<EOF\n0 $P\n$OMEGA $P\nEOF\n";
    #print CSH "psxy $J $R $X $Y -W0.5p,0/0/0, -K -O -V>>$psfile<<EOF\n0 $vP\n$OMEGA $vP\nEOF\n";
    #print CSH "psxy $J $R $X $Y -W0.5p,0/0/0, -K -O -V>>$psfile<<EOF\n$OMEGA $vP\n$OMEGA $P\nEOF\n";
    #print CSH "psxy $J $R $X $Y -W0.5p,0/0/0,- -K -O -V>>$psfile<<EOF\n$OMEGA 0\n$OMEGA $P\nEOF\n";
    print CSH "psxy $mesacdf $J $R $X $Y -W1p,0/0/255 -K -O -V>>$psfile\n"; # blue
    print CSH "psxy $omegapostcdf $J $R $X $Y -W1p,124/252/0 -K -O -V>>$psfile\n"; # green

    $y1x = 5.7; $y2x = 6.4; $y3x = 7.25;
    $yshift = -10.5;
    $y1m = $y1x + $yshift; $y2m = $y2x + $yshift; $y3m = $y3x + $yshift;
    print CSH "pstext -JM7i -R0/1/0/1 -Xa-1c -Ya${y1m}c -G50/205/50 -K -N -O -V>> $psfile<<EOF\n0 0 8 90 1 LB P(\@~w\@~)\nEOF\n";
    print CSH "pstext -JM7i -R0/1/0/1 -Xa-1.15c -Ya${y1m}c -G50/205/50 -K -N -O -V>> $psfile<<EOF\n0 0 8 90 1 LB \136 \nEOF\n";
    print CSH "pstext -JM7i -R0/1/0/1 -Xa-1c -Ya${y2m}c -G0/0/255 -K -N -O -V>> $psfile<<EOF\n0 0 8 90 1 LB V(\@~w\@~)\nEOF\n";
    print CSH "pstext -JM7i -R0/1/0/1 -Xa-1.15c -Ya${y2m}c -G0/0/255 -K -N -O -V>> $psfile<<EOF\n0 0 8 90 1 LB \136 \nEOF\n";
    if ($half_page==0){print CSH "pstext -JM7i -R0/1/0/1 -Xa-1.3c -Ya${y3m}c -K -N -O -V>> $psfile<<EOF\n0 0 10 0 1 LB (e)\nEOF\n";}

    # write OMEGA
    $omega_position = $OMEGA + 5;
    #print CSH "pstext $J $R $X $Y -W255/255/255 -K -N -O -V>> $psfile<<EOF\n$omega_position 0.07 10 0 1 LB \@~W\@~ = $OMEGA\nEOF\n";
    
    #-------------------- Plot the P(V) curve
    $vPfile = "$datadir/${eid}_vP.dat";
    $vPfile_avg = "$datadir/${eid}_vP_avg.dat";
    open(IN,$vPfile); @VP=<IN>;
    $vP_avg = 0;
    for ($ii=0;$ii<@VP;$ii++){
	($p[$ii],$vp[$ii]) = split(" ",$VP[$ii]);
	$vP_avg = $vP_avg + $vp[$ii];}
    $vP_avg = sprintf("%.2f",$vP_avg/(@VP));
    $J = "-JX3c/3c";
    $B = "-Ba1f0.1/a1f0.1::WS";
    $R = "-R0/1.0/-0./1.0";
    $Y = "-Ya${y4}c";
    $xshift = 2;
    $X = "-Xa${xshift}c";
    $ixy = 1; # to plot X-Y dashed line
    
    print CSH "psbasemap $J $R $B $X $Y -K -O -V>>$psfile\n";
    print CSH "psbasemap $J $R -Ba5/a5en $X $Y -K -O -V>>$psfile\n";
    #print CSH "psxy $J $R $X $Y -W0.5p,0/0/0, -K -O -V>>$psfile<<EOF\n$P 0\n$P $vP\nEOF\n";
    #print CSH "psxy $J $R $X $Y -W0.5p,0/0/0, -K -O -V>>$psfile<<EOF\n0 $vP\n$P $vP\nEOF\n";
    print CSH "psxy $vPfile_avg $J $R $X $Y -G211/211/211 -K -O -V>>$psfile\n";
    print CSH "psxy $vPfile $J $R $X $Y -N -W1.1p,255/0/0, -K -O -V>>$psfile\n";
    if ($ixy){
	print CSH "psxy $J $R $X $Y -W0.5p,0/0/0,- -K -O -V>>$psfile<<EOF\n0 0\n1 1\nEOF\n";
    }
    print CSH "pstext -JM7i -R0/1/0/1 -K -O -V>> $psfile<<EOF\n0 0 14 0 1 LB $ORIGIN\nEOF\n";

    $y1x = 9.35; $y2x = 7.7; $y3x = 8.5; $y4x = 11.4; $y5x = 10.75;
    $y11x = 9.75; $y33x = 8.6;
    $yshift = -16.7;
    $y1m = $y1x + $yshift; $y2m = $y2x + $yshift; $y3m = $y3x + $yshift; $y4m = $y4x + $yshift; $y5m = $y5x + $yshift;
    $y11m = $y11x + $yshift; $y33m = $y33x + $yshift;
    $x1m = $xshift -0.35; $x2m = $xshift + 1.4; $x3m = $xshift + 1.15; $x5m = $xshift -1.5;
    $x33m = $xshift + 1.5;

    #print CSH "gmtset CHAR_ENCODING Standard\n";
    print CSH "pstext -JM7i -R0/1/0/1 -Xa${x1m}c -Ya${y1m}c -G0/0/0 -K -N -O -V>> $psfile<<EOF\n0 0 12 90 12 LB \303\nEOF\n"; # P(V) Y axis label
    print CSH "pstext -JM7i -R0/1/0/1 -Xa${x1m}c -Ya${y11m}c -G0/0/0 -K -N -O -V>> $psfile<<EOF\n0 0 12 90 0 LB   (V)\nEOF\n"; # P(V) Y axis label
    print CSH "pstext -JM7i -R0/1/0/1 -Xa${x2m}c -Ya${y2m}c -G0/0/0 -K -N -O -V>> $psfile<<EOF\n0 0 10 0 0 LB V\nEOF\n"; # V X axis label
    print CSH "pstext -JM7i -R0/1/0/1 -Xa${x3m}c -Ya${y33m}c -K -N -O -V>> $psfile<<EOF\n0 0 14 0 12 LB \303 \nEOF\n";  #Pav (confidence measure)
    print CSH "pstext -JM7i -R0/1/0/1 -Xa${x33m}c -Ya${y3m}c -K -N -O -V>> $psfile<<EOF\n0 0 10 0 0 LB \@-AV\@- = $vP_avg\nEOF\n";  #Pav (confidence measure)
    #print CSH "pstext -JM7i -R0/1/0/1 -Xa3c -Ya${y4m}c -K -N -O -V>> $psfile<<EOF\n0 0 10 0 1 LB v($P) = $vP\nEOF\n";
    #print CSH "gmtset CHAR_ENCODING ISOLatin1+\n";

    if ($half_page==0){print CSH "pstext -JM7i -R0/1/0/1 -Xa${x5m}c -Ya${y5m}c -K -N -O -V>> $psfile<<EOF\n0 0 10 0 1 LB (f)\nEOF\n";}
    # move origin back
    print CSH "pstext -JM7i -R0/1/0/1 $Xr $Yr -K -O -V>> $psfile<<EOF\n0 0 14 0 1 LB $ORIGIN\nEOF\n";
}

#===========FIGURE H==============================
# Beachaball and the title
if ($figH){
    # Get the information for auxilary plane (first element of auxilary fault plane samples files - also used making histograms)
    $Mreffile = "$datadir/${eid}_Mref.dat";
    open(IN,$Mreffile); @Mref = <IN>;
    ($strike,$dip,$rake,$strike2,$dip2,$rake2)= split(" ",$Mref[0]);
    # integer values for labels
    $fstrike=sprintf("%d",$strike);
    $fdip=sprintf("%d",$dip);
    $frake=sprintf("%d",$rake);
    # Get subset of posterior samples (used for plotting P-T axis points)
    $omegaerrpost = "$datadir/${eid}_omega_err_post.dat";
    open(IN,$omegaerrpost); @omegaerrpostlines=<IN>;
    $post_samples = @omegaerrpostlines; # subset of posterior samples (400 samples)
    
    $J = "-JX4.5c/4.5c";
    $R = "-R-1/1/-1/1";
    # move origin
    $xshift = 16.5;
    $yshift = 19.5;
    $X = "-X${xshift}c";
    $Y = "-Y${yshift}c";
    # move back
    $Xr = "-X-${xshift}c";
    $Yr = "-Y-${yshift}c";
    
    #Plot beachball and P-T axes
    #NOTE: It does not appear to be possible to (easily) plot the P-T axes alone, on top of the sets of P-T axes.
    #      They have to be plotted along with the beachball, it seems.
    print CSH "psmeca $J $R -Sa4.5c $ballinfo $axisinfo $X $Y -N -K -O -V >>$psfile<<EOF\n0 0 $edep $strike $dip $rake 5 0 0\nEOF\n";

    #Plot the sets of P-T axes on the beachball
    #QUESTION: Wouldn't it be safer to pass the strike-dip-rake values of the posterior MTs into psmeca?
    #ANSWER:   Could not figure out how to just plot P-T axis and not the quadrant boundraies with psmeca
    if (1){
	for ($ii=0;$ii<$post_samples;$ii++){
	    ($omega[$ii],$misfit[$ii],$pxs[$ii],$pys[$ii],$txs[$ii],$tys[$ii]) = split(" ",$omegaerrpostlines[$ii]);
	    $pxsi = $pxs[$ii]; $pysi = $pys[$ii];$txsi = $txs[$ii];$tysi = $tys[$ii];
	    print CSH "psxy $J $R -Sc1.5p -G0/0/255 -K -O -V>>$psfile<<EOF\n$pxsi $pysi\nEOF\n";
	    print CSH "psxy $J $R -Sc1.5p -G255/0/0 -K -O -V>>$psfile<<EOF\n$txsi $tysi\nEOF\n";
	}
    }

    if ($half_page==0){print CSH "pstext -JM7i -R0/1/0/1 -Xa0c -Ya4.75c -N -K -O -V>> $psfile<<EOF\n0 0 10 0 1 LB (h)\nEOF\n";}
    
    # change to nearest integer
    $strike=sprintf("%d",$strike);
    $dip=sprintf("%d",$dip);
    $rake=sprintf("%d",$rake);
    $strike2=sprintf("%d",$strike2);
    $dip2=sprintf("%d",$dip2);
    $rake2=sprintf("%d",$rake2);
    # Move origin back
    print CSH "pstext -JM7i -R0/1/0/1 $Xr $Yr -K -O -V>> $psfile<<EOF\n0 0 14 0 1 LB $ORIGIN\nEOF\n";

    # Write the title
    $y1m = $yshift + 6.0; $dy = 0.5; $y2m = $y1m + $dy; $y3m = $y2m + $dy; $y4m = $y3m + $dy; $y5m = $y4m + $dy;
    $title = "${eid}_${norm}_${wts} ";
    $OptTitle = "";
    #$OptTitle = "Mmin=M\@-0\@-; Mref=M\@-1\@-; \@~w\@~(Mmin,Mref) = 85"; % for mtvipul paper example
    #$OptTitle = "Mmin=M1; Mref=M2; \@~w\@~(Mmin,Mref) = 101"; # optional title
    # print CSH "pstext -JM7i -R0/1/0/1 -N -Xa1c -Ya27.25c -O -K -V>> $psfile<<EOF\n0 0 12 0 1 LB $title0 \nEOF\n";
    print CSH "pstext -JM7i -R0/1/0/1 -N -Xa${xshift}c -Ya${y5m}c -O -K -V>> $psfile<<EOF\n0 0 10 0 1 LB $OptTitle\nEOF\n";
    print CSH "pstext -JM7i -R0/1/0/1 -N -Xa${xshift}c -Ya${y4m}c -O -K -V>> $psfile<<EOF\n0 0 10 0 1 LB $title\nEOF\n";
    print CSH "pstext -JM7i -R0/1/0/1 -N -Xa${xshift}c -Ya${y3m}c -O -K -V>> $psfile<<EOF\n0 0 10 0 1 LB Strike $strike Dip $dip Rake $rake\nEOF\n";
    print CSH "pstext -JM7i -R0/1/0/1 -N -Xa${xshift}c -Ya${y2m}c -O -K -V>> $psfile<<EOF\n0 0 10 0 1 LB Strike $strike2 Dip $dip2 Rake $rake2\nEOF\n";
    print CSH "pstext -JM7i -R0/1/0/1 -N -Xa${xshift}c -Ya${y1m}c -O -K -V>> $psfile<<EOF\n0 0 10 0 1 LB Mw $Mw Depth $depth km\nEOF\n";
}
#===========FIGURE B=============================
# XYZ cross-sections
if ($figB){
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
    $xshift1 = 2.25; $xshift2 = 7.75;
    $yshift1 = 2.25; $yshift2 = 5.25;
    $XY1 = "-Xa${xshift1}c -Ya${yshift1}c";
    $XY2 = "-Xa${xshift1}c -Ya${yshift2}c";
    $XY3 = "-Xa${xshift2}c -Ya${yshift2}c";
    $R1 = "-R-90/90/0/90";
    $R2 = "-R-90/90/0/360";
    $R3 = "-R0/90/0/360";

    print CSH "makecpt -Cpanoply -T$min_err/$max_err/.01 -D -Z> $cptfile \n";
    $omegacpt = "color_omega.cpt";
    $Tomega = "-T0/180/15 -D -I";
    print CSH "makecpt -Csplit $Tomega > $omegacpt\n";
    # cross-sections in misfit space (colored base)
    $stkfile = "$datadir/${eid}_stkfile.dat";
    $dipfile = "$datadir/${eid}_dipfile.dat";
    $rakfile = "$datadir/${eid}_rakfile.dat";
    # cross-sections in omega (for contours)
    $stkcont = "$datadir/${eid}_stkcont.dat";
    $dipcont = "$datadir/${eid}_dipcont.dat";
    $rakcont = "$datadir/${eid}_rakcont.dat";

    $raketick = "a90f${tickinc}";
    $striketick = "a180f${tickinc}";
    $diptick = "a90f${tickinc}";

    #print CSH "makecpt -Cwysiwyg -T0/20/1 -D -Z > $cptfile \n";
    #print CSH "makecpt -Cpanoply -T0/16/.01 -D -Z> $cptfile \n";
    # fixed-strike cross-section (bottom left)
    print CSH "psbasemap $J1 $R1 -B${raketick}:rake:/${diptick}:dip:WeSn $XY1 -K -O -V>>$psfile\n";
    print CSH "psxy $stkfile $J1 $R1 -B -C$cptfile -Ss17p $XY1 -K -O -V>>$psfile\n";
    print CSH "pscontour $stkcont -C$omegacpt $J1 $R1 -Gn1 -W1p $XY1 -K -O -V>>$psfile\n";
    print CSH "psxy $J1 $R1 -W0.5p,0/0/0,- $XY1 -K -O -V >>$psfile<<EOF\n$rake 0\n$rake 90\nEOF\n";
    print CSH "psxy $J1 $R1 -W0.5p,0/0/0,- $XY1 -K -O -V >>$psfile<<EOF\n-90 $dip\n90 $dip\nEOF\n";
    print CSH "psxy $J1 $R1 -N -B -Sa12p -G255/255/255 -W0/0/0 $XY1 -K -O -V>>$psfile<<EOF\n$rake $dip\nEOF\n";
    print CSH "pstext $J1 $R1 -W255/255/255,o/5p $XY1 -K -N -O -V>> $psfile<<EOF\n-90 90 10 0 1 LT strike=$fstrike\nEOF\n";
    
    # fixed-dip cross-section (top left)
    print CSH "psbasemap $J2 $R2 -B${raketick}:rake:/${striketick}:strike:Wesn $XY2 -K -O -V>>$psfile\n";
    print CSH "psxy $dipfile $J2 $R2 -B -C$cptfile -Ss17p $XY2 -K -O -V>>$psfile\n";
    print CSH "pscontour $dipcont -C$omegacpt $J2 $R2 -Gn1 -W1p $XY2 -K -O -V>>$psfile\n";
    print CSH "psxy $J2 $R2 -W0.5p,0/0/0,- $XY2 -K -O -V >>$psfile<<EOF\n$rake 0\n$rake 360\nEOF\n";
    print CSH "psxy $J2 $R2 -W0.5p,0/0/0,- $XY2 -K -O -V >>$psfile<<EOF\n-90 $strike\n90 $strike\nEOF\n";
    print CSH "psxy $J2 $R2 -N -B -Sa12p -G255/255/255 -W0 $XY2 -K -O -V>>$psfile<<EOF\n$rake $strike\nEOF\n";
    print CSH "pstext $J2 $R2 -W255/255/255,o/5p $XY2 -K -N -O -V>> $psfile<<EOF\n-90 360 10 0 1 LT dip=$fdip\nEOF\n";
    
    # fixed-rake cross section (top right)
    print CSH "psbasemap $J3 $R3 -B${diptick}:dip:/${striketick}::wESn $XY3 -K -O -V>>$psfile\n";
    print CSH "psxy $rakfile $J3 $R3 -B -C$cptfile -Ss17p -: $XY3 -K -O -V>>$psfile\n";  # -: toggles between (x,y) and (y,x)
    print CSH "pscontour $rakcont -C$omegacpt $J3 $R3 -Gn1 -W1p $XY3 -: -K -O -V>>$psfile\n";
    print CSH "psxy $J3 $R3 -W0.5p,0/0/0,- $XY3 -K -O -V >>$psfile<<EOF\n$dip 0\n$dip 360\nEOF\n";
    print CSH "psxy $J3 $R3 -W0.5p,0/0/0,- $XY3 -K -O -V >>$psfile<<EOF\n0 $strike\n90 $strike\nEOF\n";
    print CSH "psxy $J3 $R3 -N -B -Sa12p -G255/255/255 -W0 -: $XY3 -K -O -V>>$psfile<<EOF\n$strike $dip\nEOF\n";
    print CSH "pstext $J3 $R3 -W255/255/255,o/5p $XY3 -K -N -O -V>> $psfile<<EOF\n0 360 10 0 1 LT rake=$frake\nEOF\n";
    
    $x1m = $xshift1 + 7.5; $x2m = $xshift1 + 7.25;
    print CSH "psscale -C$cptfile -D${x1m}c/3c/2.5c/0.2c -K -O -Ac -B${dmisfit}::/:: >> $psfile\n"; # relative to the basemap
    print CSH "pstext -JM7i -R0/1/0/1 -Xa${x2m}c -Ya2.7c -K -N -O -V>> $psfile<<EOF\n0 0 11 90 1 LB misfit\nEOF\n";
    
    $xtag = $xshift1 - 1.5;
    if ($half_page==0){
    print CSH "pstext -JM7i -R0/1/0/1 -Xa${xtag}c -Ya14.75c -K -N -O -V>> $psfile<<EOF\n0 0 10 0 1 LB (b)\nEOF\n";}
}
#===========FIGURE G=============================
# Plot Histograms
if ($figG){
    $histogramfile = "$datadir/${eid}_histo.dat"; # fault plane
    $histogramfile2 = "$datadir/${eid}_histo2.dat"; # auxiliary fault plane
    
    print CSH "gmtset LABEL_OFFSET 0.03c\n";
    $histo_ymax = 30;
    $J = "-JX3c/3c";
    $R1 = "-R0/360/0/$histo_ymax"; # strike 
    $R2 = "-R0/90/0/$histo_ymax";  # dip
    $R3 = "-R-180/180/0/$histo_ymax"; # rake
    
    $xshift = 12.5;
    $X = "-Xa${xshift}c";
    $Yinit = 2.25;
    $Yshft = 4.5;
    $Y = "-Ya${Yinit}c";
    print CSH "pshistogram $histogramfile2 $J $R1 -G0/0/0 -L0.25p,0/0/0 -Z1 -W10  -Ba180f${tickinc}:strike:/a10::WeSn $X $Y -T0 -K -O -V >>$psfile\n";
    print CSH "pshistogram $histogramfile $J $R1 -G0/255/0 -L0.25p,0/0/0 -Z1 -W10  -Ba180f${tickinc}:strike:/a10::WeSn $X $Y -T0 -K -O -V >>$psfile\n";
    $Yinit = $Yinit+$Yshft;
    $Y = "-Ya${Yinit}c";
    print CSH "pshistogram $histogramfile2 $J $R2 -G0/0/0 -L0.25p,0/0/0 -Z1 -W5  -Ba${tickinc}:dip:/a10::WeSn $X $Y -T1 -K -O -V >>$psfile\n";
    print CSH "pshistogram $histogramfile $J $R2 -G0/255/0 -L0.25p,0/0/0 -Z1 -W5  -Ba${tickinc}:dip:/a10::WeSn $X $Y -T1 -K -O -V >>$psfile\n";
    $Yinit = $Yinit+$Yshft;
    $Y = "-Ya${Yinit}c";
    print CSH "pshistogram $histogramfile2 $J $R3 -G0/0/0 -L0.25p,0/0/0 -Z1 -W10  -Ba180f${tickinc}:rake:/a10::WeSn $X $Y -T2 -K -O -V >>$psfile\n";
    print CSH "pshistogram $histogramfile $J $R3 -G0/255/0 -L0.25p,0/0/0 -Z1 -W10  -Ba180f${tickinc}:rake:/a10::WeSn $X $Y -T2 -K -O -V >>$psfile\n";
    
    $xtag = $xshift - 1.5;
    print CSH "pstext -JM7i -R0/1/0/1 -Xa${xtag}c -Ya14.75c -K -N -O -V>> $psfile<<EOF\n0 0 10 0 1 LB (g)\nEOF\n";
}

#===========FIGURE I=============================
# Plots sample beachballs of the posterior samples
if ($figI){
    $xshift=17;
    $yshift=1;
    $X = "-X${xshift}c";
    $Y = "-Y${yshift}c";
    $Xr = "-X-${xshift}c";
    $Yr = "-Y-${yshift}c";
    $B = "-BWeSn";
    $J = "-JX4c/17c";
    $R = "-R0/4/0/16.5";
    
    print CSH "makecpt -Cpanoply -T$min_err/$max_err/.01 -D -Z> $cptfile \n";
    print CSH "psbasemap $J $R $B $X $Y -K -O -V>>$psfile\n";
    print CSH "gmtset BASEMAP_FRAME_RGB 0/0/0 \n";
    
    $samplesfile = "$datadir/${eid}_samples.dat";
    if (not -f $samplesfile) {die("Missing File: Moment tensor file $samplesfile missing");}
    open(IN,$samplesfile); @sampleslines = <IN>; $nsamples = @sampleslines;
    #($omega,$misfit,$stk_samples,$dip_samples,$rak_samples) = split(" ",$sampleslines[0]);

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
    
    # VERTICAL
    for ($i = $nsamples; $i >= 1; $i--) {
	($omega,$misfit,$stk_sample,$dip_sample,$rak_sample) = split(" ",$sampleslines[$nsamples-$i]);
	$omega=sprintf("%.0f",$omega);
	if ($i%2 != 1) {
	    $yloc = ($i-1)/1.25;
	    print CSH "psmeca $J $R -Sa0.4 -M -G0/0/255 -L1p/0/0/0 -Z$cptfile -K -O -V >>$psfile<<EOF\n1 $yloc $misfit $stk_sample $dip_sample $rak_sample $Mw 0 0 $omega\nEOF\n";
	} else {
	    $yloc = $i/1.25;
	    print CSH "psmeca $J $R -Sa0.4 -M -G0/0/255 -L1p/0/0/0 -Z$cptfile -K -O -V >>$psfile<<EOF\n3 $yloc $misfit $stk_sample $dip_sample $rak_sample $Mw 0 0 $omega\nEOF\n";
	}
    }
    print CSH "pstext -JM7i -R0/1/0/1 -Xa-0.5c -Ya17.25c -K -N -O -V>> $psfile<<EOF\n0 0 10 0 1 LB (i)\nEOF\n";
    # Move origin back
    print CSH "pstext -JM7i -R0/1/0/1 -N $Xr $Yr -K -O -V>> $psfile<<EOF\n0 0 14 0 1 LB $ORIGIN\nEOF\n";
}


#========================================================================================
# Always use
print CSH "pstext -JM7i -R0/1/0/1 -N -O -V>> $psfile<<EOF\n0 0 14 0 1 LB \nEOF\n";

close (CSH);
system("csh -f $cshfile");
#system("mv 20090407201255351_L1_M111.ps /home/vipul/CAP/inv/scak/MOOS/20090407201255351/stn${Nst}");
#system("mkdir -p ./SCAK/${norm}/${wts}/");
#system("mv $psfile ./SCAK/${norm}/${wts}/");
#return;
#system("gv $psfile &");
