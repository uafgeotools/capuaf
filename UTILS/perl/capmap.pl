#!/usr/bin/perl -w

#==========================================================
#  capmap.pl
#  plot maps of misfit for CAP moment tensor inversion results
#  location of this script: $CAPHOME/UTILS/perl/capmap.pl
#
# EXAMPLES:
# > cd $CAPHOME/UTILS/perl/examples/20090407201255351
# > $CAPHOME/UTILS/perl/capmap.pl 20090407201255351
# > cd $CAPHOME/UTILS/perl/examples/HOYA
# > $CAPHOME/UTILS/perl/capmap.pl HOYA
#
# FOR A NEW EXAMPLE:
# 1. You need the following files in your current working directory:
#    20090407201255351_station_list_ALL.dat -- from pysep
#    20090407201255351.out -- from cap output (needs to be renamed)
# 2. Search capmap.pl for "KEY COMMAND" to see where to make changes ($USECUSTOM = 1)
#
# OBSOLETE EXAMPLES:
# > cd $CAPHOME/UTILS/perl/examples/old_examples
# > capmap.pl 319605             # crustal MOOS event
# > capmap.pl 289317             # slab MOOS event
# > capmap.pl 20120411091837     # Nenana interior
# > capmap.pl 266672             # crustal interior
# > capmap.pl 294837             # small MOOS event
# > capmap.pl 206                # Uturuncu event
# > capmap.pl 20100516063454     # Uturuncu event
#
# FUTURE WORK
#   This script is a monster and probably ought to be rewritten from scratch.
#   Have this script automatically run when cap is run.
#   Eliminate the need to rename the cap out file.
#   Need to have capmap.pl automatically identify if .out examples are for FMT or DC.
#
#==========================================================

if (@ARGV < 1) {die("Usage: capmap.pl eid\n")}
($eid) = @ARGV;

# KEY COMMANDS
$USECUSTOM = 0;
# srad_km:  max distance to stations
# Jbscale:  size of basemap (increase number to decrease size)
# color scale for time shifts for surface waves (dtminS, dtmaxS) and P (dtminP, dtmaxP)
$srad_km = 450; $Jbscale = 7500000; $dtminS = -6; $dtmaxS = 14; $dtminP = -2; $dtmaxP = 2;  # 20190326
#$srad_km = 300; $Jbscale = 5000000; $dtminS = -1; $dtmaxS = 4; $dtminP = -1; $dtmaxP = 3;  # 20181027

# read event information. This is only for old output files.
$capevent = "${eid}_event_info.dat";
open(IN,$capevent); @capeventlines = <IN>;
(undef,$elat,$elon,$edep)  = split(" ",$capeventlines[0]);

# read in a list of all available waveforms (station, network, lat, lon)
$capstation = "${eid}_station_list_ALL.dat";
if (not -f $capstation) {die("Check if capstation $capstation exist or not\n");}
open(IN,$capstation); @capstationlines = <IN>; $nsta = @capstationlines;

# read in cap output file
$capout = "${eid}.out";
if (not -f $capout) {die("Check if capout $capout exist or not\n");}
open(IN,$capout); @linescap = <IN>;

# from cap.c, conversion between M0 and Mw
# amp = pow(10.,1.5*m0+16.1-20);    # Mw to M0 in CAP
(undef,undef,undef,$smodeldep,undef,$strike,$dip,$rake,undef,$Mw,undef,$rms,undef,undef,$gamma,undef,$delta,undef,$vr,undef,undef) = split(" ",$linescap[0]);
if (not -f $capsevent){(undef,undef,undef,$elat,undef,$elon,undef,$edep) = split(" ",$linescap[1]);}
(undef,undef,undef,$M0dyne,$ma,$mb,$mc,$md,$me,$mf) = split(" ",$linescap[2]);
(undef,$smodel,$edepcap) = split("_",$smodeldep);
@Mcmt = ($mf,$ma,$md,$mc,-$me,-$mb);  # convert from CAP(AkiRichards) to psmeca(GCMT)
$M0 = $M0dyne*1e-7;  # N-m

#print "\n -- @linescap[0] -- \n";
#print "\n -- $smodel -- $strike $dip $rake -- $Mw -- $rms -- $sigstrike $sigdip $sigrake --\n";

$topocorr = 0;

# default values
$bmin = -7; $bmax = 0;        # basement
$dmin = 0; $dmax = 160; $depinc = 40; $deptick = 10;
$itopocolor = 1;                # 1-7 (default 2 or 7; 5-6=Jeff)

#$origin_inset = "-Xa3.25 -Ya0.1";
$Jbscalei = 50000000;

$iportrait = 1;
$otitle1 = "-Xa0.0 -Ya7.10";
$otitle2 = "-Xa0.0 -Ya6.75";
$obar = "-Xa6.4 -Ya0";  # scale bar
$obar2 = "-Xa6.4 -Ya2.5";  # scale bar
$obar3 = "-Xa2.2 -Ya5.5";  # scale bar
#$obar3 = "-Xa2.2 -Ya4.3"; # publication (Nenana)

$Eres = "-E100";
$coast_res = " -Di -A0/0/1 -W0.75p,0/0/0";  # -Df
$coastinfo = "-N1/2p/0/0/0 -N2/1p/0/0/0 $coast_res -S150/255/255 -C150/255/255";

$iB = 8;
$stag = "Alaska";
# These are default locations (perhaps will change for different regions)
# Set region specific values 
$origin = "-X0.6 -Y1";
$origin_pp2 = "-X4.2 -Y0";
$origin_inset = "-Xa4.5 -Ya6.2";

$plot_unused_stations = 1;

# THE 1D MODEL USED IN THE CAP INVERSION WILL DETERMINE THE PLOTTING REGION

# search ifkmod below to specify additional parameters for each region
if($smodel eq "scak") {
   $ifkmod = 1;
   $pmax = 1;
   print "\nifkmod = $ifkmod: southern Alaska\n";
   $xtick1 = 2; $ytick1 = 1; $xtick2 = 1; $ytick2 = 0.5;
   $emax = 3000; $emin = -$emax; $itopocolor = 4;
   $ifmt = 0; $cmtsize = 0.4;

   #$coastinfo = "-N1/1p/0/0/0 $coast_res -S150/255/255 -C150/255/255 -Ia/1.5p,150/255/255";   # rivers and covered bathymetry
   $coastinfo = "-N1/1p/0/0/0 $coast_res -S150/255/255";
   $coastinfo = "-N1/1p/0/0/0 $coast_res";

} elsif($smodel eq "tactmod") {
   $ifkmod = 2;
   $pmax = 1;
   print "\nifkmod = $ifkmod: central Alaska\n";
   $xtick1 = 1; $ytick1 = 1; $xtick2 = 0.5; $ytick2 = $xtick2;
   $emax = 3000; $emin = -$emax; $itopocolor = 4;
   $ifmt = 0; $cmtsize = 0.4;

} elsif($smodel eq "utu06") {
   $ifkmod = 3;
   $pmax = 1;
   print "\nifkmod = $ifkmod: PLUTONS - Bolivia\n";
   $xtick1 = 0.2; $ytick1 = 0.2; $xtick2 = 0.1; $ytick2 = $xtick2;
   $emax = 6000; $emin = 0; $itopocolor = 1;
   $ifmt = 1; $cmtsize = 0.7;
   $topocorr = 4;

} elsif($smodel eq "wes") {
   $ifkmod = 4;
   $pmax = 1;
   print "\nifkmod = $ifkmod: Western U.S.\n";
   $xtick1 = 5; $ytick1 = 5; $xtick2 = 1; $ytick2 = $xtick2;
   $emax = 4161; $emin = -5097; $itopocolor = 4;
   $ifmt = 1; $cmtsize = 0.7;
   $topocorr = 0;
   $plot_unused_stations = 0;

}  else {
  print "\n smodel: $smodel\n";
  die("MODEL IS NOT ALLOWED");
}

if($USECUSTOM==1) {$ifkmod=5;}

# correction for high elevation stations
$edepcap = $edepcap - $topocorr;

# for title (note: ISOLatin1+ character encoding for the +/- symbol)
# $captag = "model $smodeldep (strike $strike \\261 $sigstrike, dip $dip \\261 $sigdip, rake $rake \\261 $sigrake) Mw $Mw VR $vr";
$gamma = sprintf("%.1f",$gamma);
$delta = sprintf("%.1f",$delta);
$dip = sprintf("%d",$dip);
$vr = sprintf("%.1f",$vr);

#----------------

# psmeca format for complete AEC MT catalog (with AECmt eid labels)
#$cmtfile = "/home/carltape/gmt/CMT/alaska/MT_AEC_rsub00_tsub01_12_psmeca_eid";

# get source information from AEC MT catalog (psmeca format)
#print "\n-- grep $eid $cmtfile--\n";
#@cmtline = split(" ",`grep $eid $cmtfile`);
##chomp(@cmtline);
#pop(@cmtline);
#print "\n -- @cmtline --\n";
#$elon = $cmtline[0];
#$elat = $cmtline[1];
#$edep = $cmtline[2];  # depth from AEC MT inversion
#$etag = sprintf("AEC MT %s, depth %.1f km",$eid,$edep);

#==================================================================================

$dir0 = "/home/carltape/gmt";
$cshfile = "capmap.csh";
$name = "capmap";

# DATA FILES
$dlib = "/home/admin/share/datalib";
$platefile  = "$dir0/plates/plate_boundaries/bird_boundaries";
#$plate_labs = "$dir0/alaska_2005/alaska_plate_labs.xyz";
$file_trench = "$dlib/plate_boundaries/Tape_sub_zones/cht_grav_low_outer_RUL";
$iRUL = 1;
if (not -f $platefile) {die("Check if platefile $platefile exist or not\n");}
if (not -f $file_trench) {die("Check if file_trench $file_trench exist or not\n");}

# topo files
$grdfile0 = "$dlib/topography/GLOBAL/ETOPO1/ETOPO1_Ice_g.grd";
$gradfile0 = "$dlib/topography/GLOBAL/ETOPO1/ETOPO1_Ice_g.grad";
if (not -f $grdfile0) {die("Check if grdfile $grdfile0 exist or not\n");}
if (not -f $gradfile0) {die("Check if gradfile $gradfile0 exist or not\n");}

# CMT sources
$dir_tomo    = "/home/carltape/ADJOINT_TOMO/ADJOINT_TOMO_OUTPUT";
$dir_sources = "/home/carltape/results/SOURCES/alaska_16";

# cities
$cities = "$dir0/cities/alaska_cities";

#==================================================================================

if($iportrait == 1){$orient = "-P"; $rotangle = 0}
else {$orient = " "; $rotangle = 90}

# plotting specifications
$fsize0 = "24";
$fsize1 = "18";
$fsize2 = "12";
$fsize3 = "10";
$fontno = "1";    # 1 or 4
$tick   = "0.3c";
#$tick   = "0.0c";  # temp
$fpen   = "2p";
$tpen   = "2p";

# plotting specificiations
# NOTE: pen attribute (-W) does not work for psvelo
#$coast_info     = "-A5000 -Dl -W1.0p -S220/255/255 -C220/255/255 -G200";  # A : smallest feature plotted, in km^2; D : resolution

#$coast_info     = "${coast_res} -W1.5p -Na/1.0p -S220/255/255 -C220/255/255 -G200";  # -I1/1.0p for rivers
#$coast_info2    = "${coast_res} -W1.0p -Na/1.0p";
#$coast_info2    = "${coast_res} -W1.0p -Na/1.5p";
#$coast_info3    = "${coast_res} -W1.5p -Na/1.0p -S220/255/255 -C220/255/255 -G200";
#$coast_info4    = "${coast_res} -Na/0.5p,255/255/255 -S220/255/255 -C220/255/255 -G200";
$coast_info_bath = "${coast_res} -W1p,0 -Na/1.0p,255/255/255 -G200";
$coast_info_inset  = "-W0.5p -Dl -A500/0/1 -N1/1p -S220/255/255 -C220/255/255 -G150";

$plate_info     = "-m -W2p,255/0/0";
$plate_infoK    = "-m -W2p,0/0/0";
$plate_infoY    = "-m -W2p,255/255/0";
$plate_infoR    = "-m -W2p,255/0/0";
#$plate_infoR    = "-m -W4p,255/0/0";

$fault_info     = "-m -W1.5p,255/0/0";
$fault_infoK    = "-m -W1.5p,0/0/0";

$stainfo  = "-Si12p -W1p/0/0/0";
#$stainfo  = "-Si18p -W1p/0/0/0";  # MOOSrs
$stainfo0 = "$stainfo -G255/255/255";
$station_info   = "-Sc5p -W0.5p/255/255/255 -G255/0/0";
$station_info2  = "-Si10p -W0.5p/0/0/0 -G255/0/0";

$textinfo       = "-G255 -S1.5p";
$textinfo2      = "-G0/0/0 -S2p,255/255/255";
$textinfo3      = "-G0/0/0 -S2p,255/255/0";
$inset_box_info = "-W3.0p,0/0/255 -N";

# BOOLEAN: WHICH FIGURES TO PLOT
$ixv    = 1;
$ipdf   = 0;

# topography/bathymetry
$Dlen1 = 2;
$Dx = 0;
$Dy = $Dlen1/2;
$Dwid = 0.20;
$Dscale_topo = "-D$Dx/$Dy/$Dlen1/$Dwid";
$Bscale_topo = "-B1000f500:\"Elevation (m)\":";

$Dscale_grav = "-D$Dx/$Dy/$Dlen1/$Dwid";
$Bscale_grav = "-B20f10:\"Isostatic Gravity Anomaly (mgal)\": -E10p";

$Dscale_base = "-D$Dx/$Dy/$Dlen1/$Dwid";
$Bscale_base = "-B1f0.5:\"Basement surface (km)\": -Eb10p";

# CMT depth scale
$Dlen2 = 2;
$Dx = 0;
$Dy = $Dlen2/2;
#$depinc = 40;
$Dscale_dep = "-D$Dx/$Dy/$Dlen2/$Dwid";
$Bscale_dep = "-Ba${depinc}f${deptick}:\"Earthquake depth-km\": -Eb10p";

# plot title
$J_title = "-JM7i";  # -JM7i
$R_title = "-R0/1/0/1";
$fsize_title1 = 16;
$fsize_title2 = 14;

# open the shell script
open(CSH,">$cshfile");
print CSH "gmtset PAPER_MEDIA letter CHAR_ENCODING ISOLatin1+ MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2 HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 FRAME_PEN $fpen TICK_PEN $tpen\n";

#==================================================================================
# PLOTTING FEATURES

$iscaletopo = 0;
$iscalebar = 1;                 # default 1
$ibasements = 0; $iscalebase = 0;
$islab = 0; $iscaleslab = 0;
$ivolc = 0;
$ifaults = 1;                   # default 1
$ifaults_extra = 0;             # default 0
$iplates = 1;                   # default 1
$itrench = 0;
$ibox = 0;
$iinset = 0;

$istation_seis = 0;
$istationlabels = 0;
@isall = (1,0,0,0,1,0,0,0,0,0,0);  # flags for numerous stations files

$islab_extent = 1;      # default 1
$iroads = 0;
$iplatelabels = 0;
$ifaultlabels = 0;
$itopolabels = 0;
$icities = 0; $icitylabels = 0;

#$title1 = sprintf("CAP results for eid %s, depth %.0f km (catalog %.1f km)",$eid,$edepcap,$edep);
#$title2 = "model $smodeldep (gamma $gamma, delta $delta, strike $strike, dip $dip, rake $rake) Mw $Mw VR $vr";
$title1 = "Event $eid";
$title2 = sprintf("model %s, depth %.0f km (catalog %.1f km)",$smodeldep,$edepcap,$edep);

# USER CHANGE THESE
$islabc = 0;
$itopo = 1;

#==================================================================================
# LOOP OVER DIFFERENT SCALAR QUANTITIES FROM CAP

# 23 columns of CAP: weights in columns 3,7,11,15,19
# NOTE: TIME SHIFTS FOR Z AND R ARE FIXED -- (BOTH PNL AND SURF)

@caplabs = ("ZRdt","Tdt","Zcc","Rcc","Tcc","ZRPnldt","ZPnlcc","RPnlcc","PVratio","PRratio","SVratio","SRratio","STratio");
@capvals = ("Rayleigh(Z,R) time shift, s","Love time shift, s","Rayleigh(Z) CC","Rayleigh(R) CC","Love CC","P(Z,R) time shift, s","P(Z) CC","P(R) CC","ln(Aobs/Asyn), P(Z)","ln(Aobs/Asyn), P(R)","ln(Aobs/Asyn), Rayleigh(Z)","ln(Aobs/Asyn), Rayleigh(R)","ln(Aobs/Asyn), Love");

# example ruler from a CAP output file
#  1                                2       3     4   5    6      7      8        9   10  11    12   13    14      15       16  17  18   19    20   21       22       23  24    25 26  27
# 20090407201255351.YV.KASH..BH  49.0/-0.00 0   0.00 63  2.86  3.08 8.27e-04 3.79e-05 0   0.00  5  2.86  3.41 9.47e-04 3.12e-05 1   0.58 91 -0.58  0.25 2.88e-04 2.25e-04 1   1.40 46 -0.58  
# 28     29       30     31   32  33  34    35       36         37     38 39
# 0.40 3.26e-04 2.18e-04 1   0.41 85 -0.30  0.10 1.79e-04   1.61e-04    0 0
# the CAP output files do not seem to match these values. Pending.
@capcols = (20,34,19,26,33,6,5,12,7,14,21,28,35);  # column of data to plot
@wcols = (17,31,17,24,31,3,3,10,3,10,17,24,31); # columns with weight factors
@icc = (0,0,1,1,1,0,1,1,2,2,2,2,2);           # whether CC (1), tshift (0), data-syn ratio(2)

# color limits
$ccmin = 30;
$ccmax = 100;
$ccinc = 10;

# Surface wave min-max shift
#$dtmaxsurf = 4.5;
$dtmaxsurf = 8.5;
$dtminsurf = -$dtmaxsurf;
#$dtminsurf = 0; $dtmaxsurf = 5;
# P wave min-max shift
#dtmaxpnl = 2.5;
$dtmaxpnl = 2.2;
$dtminpnl = -$dtmaxpnl;
if($USECUSTOM) {
  $dtminsurf = $dtminS;
  $dtmaxsurf = $dtmaxS;
  $dtminpnl = $dtminP;
  $dtmaxpnl = $dtmaxP;
}
$dtinc = 1;
$P_amp_ratio_min = -3;
$P_amp_ratio_max = -$P_amp_ratio_min;
$S_amp_ratio_min = -2;
$S_amp_ratio_max = -$S_amp_ratio_min;
$Pinc = 1;
$P_amp_ratio_inc = 1;
$S_amp_ratio_inc = 1;
@cmins = ($dtminsurf,$dtminsurf,$ccmin,$ccmin,$ccmin,$dtminpnl,$ccmin,$ccmin,$P_amp_ratio_min,$P_amp_ratio_min,$S_amp_ratio_min,$S_amp_ratio_min,$S_amp_ratio_min);
@cmaxs = ($dtmaxsurf,$dtmaxsurf,$ccmax,$ccmax,$ccmax,$dtmaxpnl,$ccmax,$ccmax,$P_amp_ratio_max,$P_amp_ratio_max,$S_amp_ratio_max,$S_amp_ratio_max,$S_amp_ratio_max);
@cticks = ($dtinc,$dtinc,$ccinc,$ccinc,$ccinc,$dtinc,$ccinc,$ccinc,$P_amp_ratio_inc,$P_amp_ratio_inc,$S_amp_ratio_inc,$S_amp_ratio_inc,$S_amp_ratio_inc);

#==================================================================================
# LOOP OVER DIFFERENT SCALAR QUANTITIES TO PLOT

# KEY COMMAND: min and max indices for plotting the maps listed above
#$xpmin = 1; $xpmax = @caplabs;  # for full set of figures
$xpmin = 1; $xpmax = $xpmin;    # for testing

print "\n PLOTTING MAPS FROM $xpmin TO $xpmax\n";

for ($xx = $xpmin; $xx <= $xpmax; $xx++) {
#==================================================================================

$clab = $caplabs[$xx-1];
$fname = "${eid}_${xx}_$clab";
$psfile = "$fname.ps";

#==================================================================================
# LOOP OVER TWO DIFFERENT MAP REGIONS (FOR EACH PAGE)
#$pmax = 2;
for ($pp = 1; $pp <= $pmax; $pp++) {
#==================================================================================

    if ($ifkmod==1) {   # SOUTHERN ALASKA + INTERIOR ALASKA
	$icities = 0;
	$iroads = 1;
	if ($elat >= 62) {   # INTERIOR ALASKA (DENALI) - only one map ($pmax=1)
	    $xmin = -157.0; $xmax = -140.0; $ymin = 58; $ymax = 68.0;
	    $sbarinfo = "-L-144/58.7/58.7/200+p1.5p,0/0/0,solid+f255/255/255"; $Jbscale = 6000000;
	    $ibasementc = 0;
	    $ititle = 1;	# publication
	    $iscalecap = 1;
	    $pmax = 1;
	    $otitle1 = "-Xa0.0 -Ya9.1";
	    $otitle2 = "-Xa0.0 -Ya8.7";
	    $obar3 = "-Xa2.2 -Ya7.5"; # scale bar
	    $orient = "-P"; $rotangle = 0;
	    $iinset = 1;
	} elsif ($pp==1) {   # SOUTHERN ALASKA - map_1
	    $xmin = -155.0; $xmax = -142.0; $ymin = 57; $ymax = 66.0;
	    $sbarinfo = "-L-146/57.5/57.5/100+p1.5p,0/0/0,solid+f255/255/255"; $Jbscale = 8000000;
	    # For offshore Kodiak earthquake
	    #$xmin = -165.0; $xmax = -132.5; $ymin = 53; $ymax = 66.0;
	    #$sbarinfo = "-L-140/55/55/300+p1.5p,0/0/0,solid+f255/255/255"; $Jbscale = 11000000;
	    # smaller region for NEHRP
	    #$xmin = -154.0; $xmax = -146.0; $ymin = 59; $ymax = 64.0;
	    #$sbarinfo = "-L-149/59.5/59.5/100+p1.5p,0/0/0,solid+f255/255/255"; $Jbscale = 4000000;
	    $ibasementc = 0;
	    $ititle = 1;	# publication
	    $iscalecap = 1;
	    $origin = "-X0.6 -Y1";
	    $origin_inset = "-Xa4.5 -Ya6.2";
	    $orient = '-P'
	} else {   # SOUTHERN ALASKA ZOOMED into MOOS - map_2
	    $xmin = -151.5; $xmax = -147.75; $ymin = 59.5; $ymax = 62.0;
	    $sbarinfo = "-L-149/59.7/59.7/100+p1.5p,0/0/0,solid+f255/255/255";
	    $Jbscale = 2500000;
	    $ibasementc = 1;
	    $ititle = 0;
	    $iscalecap = 0;
	    $origin_pp2 = "-X4.2 -Y0";
	}
    } elsif ($ifkmod==2){   # INTERIOR ALASKA (NENANA)
	$icities = 1;
	$iroads = 1;
	if ($pp==1) {
	    $xmin = -153.0; $xmax = -144.0; $ymin = 63; $ymax = 66.5;
	    $sbarinfo = "-L-147/63.3/63.3/100+p1.5p,0/0/0,solid+f255/255/255";
	    $Jbscale = 3000000;
	    $ibasementc = 0;
	    $ititle = 1;	# publication
	    $iscalecap = 1;
	}
	if ($pp==2) {
	    $xmin = -150.2; $xmax = -148.5; $ymin = 64.5; $ymax = 65;
	    $sbarinfo = "-L-148.5/64.6/64.6/10+p1.5p,0/0/0,solid+f255/255/255";
	    $Jbscale = 800000;
	    $ibasementc = 0;
	    $ititle = 0;	# publication
	    $iscalecap = 1;
	}
    } elsif ($ifkmod==3){   # UTURUNCU BOLIVIA
	$xmin = -67.8; $xmax = -66.7; $ymin = -22.8; $ymax = -21.8;
	$sbarinfo = "-L-67/-22.7/-22.7/20+p1.5p,0/0/0,solid+f255/255/255";
	$Jbscale = 800000;
	$ibasementc = 0;
	$ititle = 1;		# publication
	$iscalecap = 1;

    } elsif ($ifkmod==4) {
	#$xmin = -125.5; $xmax = -104; $ymin = 30; $ymax = 47;
	$xmin = -124; $xmax = -110; $ymin = 32; $ymax = 42;
	#$sbarinfo = "-L-107.0/46.0/40/100+p1.5p,0/0/0,solid+f255/255/255";
	$sbarinfo = "-L-122.5/33.8/33.8/100+p1.5p,0/0/0,solid+f255/255/255";
	$Jbscale = 8000000;
	$ibasementc = 0;
	$ititle = 1;		# publication
	$iscalecap = 1;
	$orient = "-P";
	$rotangle = 0;
	$iplates = 0;
    }

  # centered on epicenter
  if ($USECUSTOM) {
    print "\n CUSTOM REGION CENTERED ON EPICENTER\n";
    $dy = $srad_km/100;   # degrees
    $fac = cos($elat*3.14159/180.0);
    $xmin = $elon - $dy/$fac;
    $xmax = $elon + $dy/$fac;
    $ymin = $elat - $dy;
    $ymax = $elat + $dy;

    $xran = $xmax - $xmin;
    $yran = $ymax - $ymin;
    $xL = $xmin + 0.9*$xran;
    $yL = $ymin + 0.1*$yran;
    $sbarinfo = "-L$xL/$yL/$yL/100+p1.5p,0/0/0,solid+f255/255/255";
    $icities = 0;
    $ititle = 1;
    $iscalecap = 1;
    $ibasementc = 0;
    $islab_extent = 0;

    $xtick1 = 2; $ytick1 = 1; $xtick2 = 0.5; $ytick2 = $xtick2;
    $emax = 3000; $emin = -$emax; $itopocolor = 4;
    $ifmt = 0; $cmtsize = 0.4;
  }

  $xcen = ($xmin + $xmax)/2;
  $ycen = ($ymin + $ymax)/2;
  $J = "-Jb$xcen/$ycen/$ymin/$ymax/1:$Jbscale";
  $R = "-R$xmin/$xmax/$ymin/$ymax";

#==================================================================================
# BASEMAP OPTIONS

# scale for plots
@Bopts = ("WESN","Wesn","wesN","wEsn","weSn","WesN","wEsN","wESn","WeSn","WEsn","weSN","wESN","WESn","WeSN","WEsN","wesn");
$B0 = sprintf("-Ba%3.3ff%3.3fd:\" \":/a%3.3ff%3.3fd:\" \"::.\" \":",$xtick1,$xtick2,$ytick1,$ytick2);
#$B0 = "-Ba${xtick}f${xtick}g${xtick}:.\" \":";   # with gridlines

$B = "$B0".$Bopts[$iB];
#$B = "$B0"."wesn";   # temp

if ($pp==1) {
    print CSH "psbasemap $J $R $B -K -V $orient $origin > $psfile\n";
} else {
  print CSH "psbasemap $J $R $B -K -O -V $origin_pp2 >> $psfile\n";
}

if ($itopo==1) {
  $colordir = "/home/carltape/gmt/color_maps/";
  @cfile = (0,0,0,0,1,1,1);  # =1 if using a color file directly
  @topocolor = ("globe","topo","relief","gray","world_bright.cpt","continent.cpt","ghayes.cpt");
  $icfile = $cfile[$itopocolor-1];
  
  if ($icfile==0) {
    # make colorpoint file (-Z for continuous topography)
    $cptfile = "color_topo.cpt";
    $dc = ($emax-$emin)/20;
    $Ttopo = "-T$emin/$emax/$dc -Z -D";
    $colormap = $topocolor[$itopocolor-1];
    print CSH "makecpt -C$colormap $Ttopo > $cptfile\n";
  } else {
    $cptfile = "$colordir/$topocolor[$itopocolor-1]";
  }

  #print CSH "grdimage $grdfile0 -C$cptfile $J $R -I$gradfile0 $Eres -K -O -V >> $psfile\n";
  print CSH "grdimage $grdfile0 -C$cptfile $J $R -I$gradfile0 $Eres -K -O -V >> $psfile\n";

  if ($iscaletopo==1) {
    print CSH "psscale -C$cptfile $Dscale_topo $Bscale_topo $obar -V -K -O >> $psfile\n";
  }
}

# slab contours
if ($islab==1 || $islabc==1) {
  # make colorpoint file
  $cptfile = "color_slab.cpt";
  $smin = -160; $smax = 0;
  #$ds = ($smax-$smin)/20;
  $ds = $depinc;
  $Tslab = "-T$smin/$smax/$ds -D -I";  # -Z for continuous
  print CSH "makecpt -Cseis $Tslab > $cptfile\n";

  $slabfile = "$dlib/subduction_zones/slab/alu_slab1.0_clip.grd";
  if (not -f $slabfile) {die("Check if slabfile $slabfile exist or not\n");}
  $clipfile = "$dlib/subduction_zones/slab/alu_slab1.0.clip";
  if (not -f $clipfile) {die("Check if clipfile $clipfile exist or not\n");}
}

# plot coastlines
print CSH "pscoast $J $R $coastinfo -K -O -V >> $psfile\n";

# slab contours
if ($islabc==1) {
  # make colorpoint file
  $cptfile = "color_con.cpt";
  $ds = 10;
  $Tslab = "-T$smin/$smax/$ds -D -I";
  print CSH "makecpt -Cseis $Tslab > $cptfile\n";
  print CSH "grdcontour $slabfile -C$cptfile -A- -W+1.5p $J $R -K -O -V >> $psfile\n";

  if ($iscaleslab == 1) {
    $Bscale_slab = sprintf("-B40f20:\"Depth to subduction interface, km\": -Eb10p",$pwr);
    print CSH "psscale -C$cptfile $Dscale_topo $Bscale_slab $obar2 -V -K -O >> $psfile\n";
  }
}

if ($ibasementc==1 || $ibasements==1) {
  # make colorpoint file (-Z for continuous topography)
  $cptfile = "color.cpt";
  $db = ($bmax-$bmin)/20;
  $Tbase = "-T$bmin/$bmax/$db -Z -D";
  print CSH "makecpt -Cseis $Tbase > $cptfile\n";

  $basefile = "$dlib/basins/cook_inlet/cook_basement_shellenbaum2010_xyz.dat";
  if (not -f $basefile) {die("Check if basefile $basefile exist or not\n");}
  if ($ibasements==1) {
     print CSH "awk '{print \$1,\$2,\$3/-1000}' $basefile | psxy -C$cptfile -Sc4p $J $R -K -O -V >> $psfile\n";
   }

  if ($iscalebase==1) {
    print CSH "psscale -C$cptfile $Dscale_base $Bscale_base $obar -V -K -O >> $psfile\n";
  }

  if ($ibasementc == 1) {
    $cptfile = "color_bcon.cpt";
    $Tbase = "-T$bmin/$bmax/1 -D";
    print CSH "makecpt -Cseis $Tbase > $cptfile\n";
    print CSH "awk '{print \$1,\$2,\$3/-1000}' $basefile | pscontour -C$cptfile -W0.75p,0/0/0,-- -A- $J $R -K -O -V >> $psfile\n";
  }
}

# plot roads
if ($iroads==1) {
  $roads = "$dlib/cities_roads_etc/roads/alaska_roads.xy";
  if (not -f $roads) {
    die("Check if roads $roads exist or not\n");
  }
  $roadinfo = "-W1p,0/0/0";
  print CSH "psxy $roads $J $R $roadinfo -m -K -O -V >> $psfile\n";
}

if ($iplates==1) {
  print CSH "psxy ${platefile} $J $R $plate_infoK -K -V -O >> $psfile\n";
}
if ($itrench==1) {
  # plot trenches (with teeth)
  $plateinfo = "-m -W0.75p,0/0/0";
  if ($iRUL==0) {
    $trenchinfo = "$plateinfo -G0/0/0 -Sf0.2/0.08rt";
  } else {
    $trenchinfo = "$plateinfo -G0/0/0 -Sf0.2/0.08lt";
  }
  print CSH "psxy $file_trench $J $R $trenchinfo -O -K -V >> $psfile\n";
}

if ($ifaults==1 || $iinset==1) {
  #$faultfile = "/home/carltape/gmt/faults/alaska/alaska_faults_new2.gmtlin";
  #$faultfile = "/home/carltape/gmt/faults/alaska/alaska_faults_new2.gmtlin_mod";
  $faultfile = "/home/admin/share/datalib/faults/alaska/DGGS2012/ak_fault_database_faults.gmt";
  if (not -f $faultfile) {
    die("Check if faultfile $faultfile exist or not\n");
  }
}

if ($ifaults==1) {
  $faultcol = "255/0/0";
  $faultinfo = "-m -W1.5p,$faultcol";
  $faultinfo2 = "-m -W0.75p,$faultcol";
  print CSH "psxy $faultfile $J $R $faultinfo -K -V -O >> $psfile\n";
}

# extent of the slab
if ($islab_extent==1) {
  $slabextent = "/home/admin/share/datalib/seismicity/regional/alaska/slabextent/slab_extent_all.lonlat";
  if (not -f $slabextent) {
    die("Check if slabextent $slabextent exist or not\n");
  }
  print CSH "psxy $slabextent $J $R -W3p,0/0/255 -m -K -V -O >> $psfile\n";
}

if ($ivolc==1) {
  # global volcanoes
  $volcfile = "$dlib/volcanoes/global_volcs_lonlat.dat";
  if (not -f $volcfile) {
    die("Check if volcfile $volcfile exist or not\n");
  }
  $volcinfo = "-W0.5p,0/0/0 -St10p -G255/0/0";
  print CSH "psxy $volcfile $J $R $volcinfo -K -V -O >> $psfile\n";

  # Alaska volcanoes
  $volcfile = "$dlib/volcanoes/alaska/AKvolc_lonlat.dat";
  if (not -f $volcfile) {
    die("Check if volcfile $volcfile exist or not\n");
  }
  #$volcinfo = "-W0.5p,0/0/0 -St10p -G255/0/0";
  print CSH "psxy $volcfile $J $R $volcinfo -K -V -O >> $psfile\n";
}

  # symbols for various station files
  $stainfok = "-Si12p -W3p,0/0/0";
  #$stainfo0a = "-Si12p -W1p,255/255/255 -G0/0/0"; $stainfo0b = $stainfo0a; 

  $textinfoL = "-G0/0/0 -W255/255/255";
  $flab = 6;

print CSH "psbasemap $J $R $B -K -V -O >> $psfile\n";

#==================================================================================
# PLOTTING CAP RESULTS

 # make colorpoint file
 $cmin = $cmins[$xx-1];
 $cmax = $cmaxs[$xx-1];
 $ctick = $cticks[$xx-1];
 $ctick2 = $ctick/2;

  $cpt_cap = "color_cap.cpt";
  $colormap = "seis";
  if($icc[$xx-1] == 1) {          # Cross-correlation
    $dc = ($cmax-$cmin)/12;
    $Tcap = "-T$cmin/$cmax/$dc";
    print CSH "makecpt -C$colormap $Tcap > $cpt_cap\n";

  } elsif($icc[$xx-1] == 0) {     # time-shift
    $dc = ($cmax-$cmin)/9;
    $Tcap = "-T$cmin/$cmax/$dc -D";   # -Z for continuous, -D
    print CSH "makecpt -C$colormap $Tcap > $cpt_cap\n";

  } elsif($icc[$xx-1] == 2) {     # log(data/syn)
    $dc = ($cmax-$cmin)/9;
    $Tcap = "-T$cmin/$cmax/$dc -D";   # -Z for continuous, -D
    print CSH "makecpt -C$colormap $Tcap > $cpt_cap\n";
}
# PLOT COLORED RAY PATHS
 for ($i = 1; $i <= $nsta; $i++) {
    ($sta,$net,$lat,$lon,$dist,$az) = split(" ",$capstationlines[$i-1]);
    #print "-- STA = $sta; NET = $net\n"; # for debugging
    # WARNING: THIS COULD RETURN MULTIPLE LINES
    #@capout1 = split(" ",`grep $sta $capout`);   # get ALL stations
    #$number_of_same_name_stations = @capout1;
    #print("============@capout1  $number_of_same_name_stations\n");
    # THIS WILL RETURN MULTIPLE LINES AS ARRAYS
    @capout1 = `grep $sta $capout`;
    $number_of_same_name_stations = @capout1;
    #print("=======> @capout1  $number_of_same_name_stations\n"); # for debugging

    # LOOP OVER MULTIPLE LINES (MULTIPLE STATIONS WITH SAME NAME)
    for ($j = 0; $j < $number_of_same_name_stations; $j++) {
	@capout_stn = split(" ", @capout1[$j]);

	$number_of_elements = @capout_stn;

	#($t1,$t2,$t3,$t4,$t5,$t6,$t7,$t8,$t9,$t10,$t11,$t12,$Zcc,$Zdt,$t15,$t16,$Rcc,$Rdt,$t19,$t20,$Tcc,$Tdt) = split(" ",`grep $sta $capout`);
	$iplot = $capcols[$xx-1];
	$fplot = $capout_stn[$iplot-1];
	#print "-- iplot = $iplot, fplot = $fplot --\n @capout_stn\n"; # for debugging
	if($capout_stn[0]) {
	    $iplot2 = $wcols[$xx-1];
	    $wplot = $capout_stn[$iplot2-1];
	    #print "-- xx = $xx -- iplot2 = $iplot2 -- wplot = $wplot\n"; # for debugging
	    # only plot the colored ray path if the station was used in the inversion (nonzero weight)
	    if($wplot > 0) {
		print CSH "psxy $J $R -W3p,0/0/0 -K -O -V >>$psfile<<EOF\n$elon $elat\n$lon $lat\nEOF\n";
		print CSH "psxy $J $R -W-2p -m -C$cpt_cap -K -O -V >>$psfile<<EOF\n> -Z$fplot\n$elon $elat\n$lon $lat\nEOF\n";
	    }
	}
    }
}

# PLOT COLORED STATIONS (WITH LABELS)
# If station has waveforms but is not used, then plot an open triangle.
# If station is used, then plot a colored triangle (CC or time shift).
for ($i = 1; $i <= $nsta; $i++) {

   ($sta,$net,$lat,$lon,$dist,$az) = split(" ",$capstationlines[$i-1]);
   # plot an open triangle at each station in the output file
   $stainfo4a = $stainfok; $stainfo4b = "-Si12p -W1p,255/255/255";
   if($plot_unused_stations) {
       print CSH "psxy $J $R $stainfo4a -K -O -V >>$psfile<<EOF\n $lon $lat\nEOF\n";
       print CSH "psxy $J $R $stainfo4b -K -O -V >>$psfile<<EOF\n $lon $lat\nEOF\n";
   }

   # for waveforms used in the inversion, plot OVER with a colored symbol
   @capout = split(" ",`grep $sta $capout`);
   @capout1 = `grep $sta $capout`;
   $number_of_same_name_stations = @capout1;
   for ($j = 0; $j < $number_of_same_name_stations; $j++) {
       @capout_stn = split(" ", @capout1[$j]);
       if(@capout_stn) {
	   $iplot = $capcols[$xx-1];
	   $fplot = $capout_stn[$iplot-1];
	   $iplot2 = $wcols[$xx-1];
	   $wplot = $capout_stn[$iplot2-1];
	   # only plot the colored label if the station was used in the inversion (nonzero weight)
	   if($wplot > 0) {
	       print CSH "psxy $J $R -Si24p -W0.5p/0/0/0 -C$cpt_cap -K -O -V >>$psfile<<EOF\n $lon $lat $fplot\nEOF\n";
	       print CSH "pstext $J $R $textinfoL -K -O -V >>$psfile<<EOF\n $lon $lat $flab 0 1 RB $sta\nEOF\n";
	   }
       }
       if($plot_unused_stations) {
	   print CSH "pstext $J $R $textinfoL -K -O -V >>$psfile<<EOF\n $lon $lat $flab 0 1 RB $sta\nEOF\n";
       }
   }
}

  # make colorpoint file
  $cptcmt = "color_cmt1.cpt";
  $cptcmt2 = "color_cmt2.cpt";
  $colorbar = "seis";
  $dinc = $depinc;
  print CSH "makecpt -C$colorbar -T$dmin/$dmax/$dinc -D > $cptcmt\n";
  #print CSH "makecpt -C$colorbar -T-$dmax/-$dmin/$dinc -I -D > $cptcmt2\n";

# plot moment tensor
if($ifmt==0) {
  # plot moment tensor using the strike/dip/rake value (-Sa)
  $cmtinfo = "-Sa${cmtsize} -L1p/0/0/0";
  print CSH "psmeca $J $R $cmtinfo -Z$cptcmt -K -O -V >>$psfile<<EOF\n$elon $elat $edep $strike $dip $rake $Mw 0 0\nEOF\n";
  #print CSH "psmeca $J $R $cmtinfo -Z$cptcmt -K -O -V >>$psfile<<EOF\n$elon $elat $edep $strike $dip $rake $Mw 0 0 $eid\nEOF\n";

} else {
  # plot moment tensor using the Mij values
  # note: the iexp for psmeca here is arbitrary and trades off with SmX.X
  # note: if you want real magnitude, then recalculate M and iexp to consider the moment (M0)
  $cmtinfo = "-Sm0.4 -L1p/0/0/0";
  print CSH "psmeca $J $R $cmtinfo -Z$cptcmt -K -O -V >> $psfile<<EOF\n$elon $elat $edep @Mcmt 24@\nEOF\n";

  #$cmtinfo = "-Sm0.4 -L1p/0/0/0";
  #print CSH "psmeca $J $R $cmtinfo -Z$cptcmt -K -O -V >> $psfile<<EOF\n@cmtline\nEOF\n";
  #print CSH "psscale -C$cptcmt2 $Dscale_dep -Ef10p $Bscale_dep $obar -V -K -O >> $psfile\n";
}

  if ($iscalecap==1) {
    $Dscale_dt = "-D$Dx/$Dy/3.5/${Dwid}h";
    if ($icc[$xx-1] == 1) {
      $E = "-Eb10p";
    } else {
      $E = "-E10p";
    }
    $Bscale_dt = "-B${ctick}f${ctick2}:\"$capvals[$xx-1]\": $E";
    print CSH "psscale -C$cpt_cap $Dscale_dt $Bscale_dt $obar3 -V -K -O >> $psfile\n";
  }

#die("testing");

#==================================================================================
# TEXT LABELS

# plot box
if ($ibox==1) {
   $boxinfo = "-W2p,0/0/0";
   #$xmin = 200; $xmax = 226; $ymin = 56; $ymax = 68;
   #$boxinfo = "-W2p,0/0/255 -Ap";  # WHY DOES THE -A COMMAND NO LONGER WORK
   #print CSH "psxy $J $R $boxinfo -K -O -V >>$psfile<<EOF\n $xmin $ymin\n$xmax $ymin\n$xmax $ymax\n$xmin $ymax\n$xmin $ymin\nEOF\n";
   #$boxfile = "alaska_inset_points.dat";
   #print CSH "psxy $boxfile $J $R $boxinfo -K -O -V >>$psfile\n";

   #$boxfile = "interior_points.dat";
   #$boxfile = "ak_south_points.dat";
   $boxfile = "nenana_basin_points.dat";
   #$boxfile = "cook_inlet_basin_points.dat";
   print CSH "psxy $boxfile $J $R $boxinfo -K -O -V >>$psfile\n";
}

# plot cities
if ($icitylabels==1 || $icities==1) {
  if (not -f $cities) {
    print "try: cp /home/carltape/gmt/cities/alaska_cities_sub $cities\n";
    die("cities $cities does not exist\n");
  }

  $fsizecity = 12;
  $textinfo3 = "-G0/0/0 -S2p,255/255/255"; # -N or not
  $city_info = "-Sc7p -W1.0p/0/0/0 -G255/255/0";
  if (0==1) {
    print CSH "psxy $J $R $city_info -N -K -O -V >>$psfile<<EOF\n-147.87 64.82\nEOF\n";	# Fairbanks dot
    print CSH "pstext $textinfo3 -N -D0.05/0.05 $J $R -K -V -O >> $psfile<<EOF\n-147.87 64.82 $fsizecity 0 1 LB Fairbanks\nEOF\n"; # Fairbanks label
    print CSH "psxy $J $R $city_info -N -K -O -V >>$psfile<<EOF\n-118.40 33.93\nEOF\n";
    print CSH "pstext $textinfo3 -N -D0.05/0.05 $J $R -K -V -O >> $psfile<<EOF\n-118.40 33.93 $fsizecity 0 1 LB Los Angeles\nEOF\n";
    print CSH "psxy $J $R $city_info -N -K -O -V >>$psfile<<EOF\n-122.3 47.45\nEOF\n";
    print CSH "pstext $textinfo3 -N -D0.05/0.05 $J $R -K -V -O >> $psfile<<EOF\n-122.3 47.45 $fsizecity 0 1 LB Seattle\nEOF\n";
  }
  if ($icities==1) {
    print CSH "awk '{print \$2,\$1}' $cities | psxy ${city_info} $J $R -K -V -O >> $psfile\n";
  }
  if ($icitylabels==1) {
    print CSH "awk '{print \$2,\$1,$fsizecity,0,1,\"LB\",\$3}' $cities | pstext $textinfo3 -D0.05/0.05 $J $R -K -V -O >> $psfile\n";
  }
}

if ($itopolabels==1) {
  if (not -f ${topo_labs}) { die("Check if topo_labs ${topo_labs} exist or not\n") }
  print CSH "pstext ${topo_labs} $textinfo $J $R -K -V -O >> $psfile\n";
}
if ($ifaultlabels==1) {
  $fsizefault = 14;
  if (not -f ${fault_labs}) {
      die("fault_labs ${fault_labs} does not exist\n")
  }
  print CSH "pstext ${fault_labs} $textinfo $J $R -K -V -O >> $psfile\n";
  print CSH "awk '{print \$1,\$2,$fsizefault,\$3,1,\"LB\",\$5}' ${fault_labs} | pstext $textinfo -D0.05/0.05 $J $R -K -V -O >> $psfile\n";

}
if ($iplatelabels==1) {
  $text_info = "-G255 -S1p";
  if (not -f ${plate_labs}) {
      print "try: cp $dlib/plate_models/Peter_Bird_02/plate_labs2a_bird.xyz ${plate_labs}\n";
      die("plate_labs ${plate_labs} does not exist\n")
  }
  print CSH "pstext $plate_labs $J $R $text_info -K -O -V >> $psfile\n";
}

# inset map
if ($iinset==1) {
  # inset map
  # WHY IS THE INSET BOX NOT DRAWN WITH GREAT CIRCLE ARCS?
  $xmini = 188; $xmaxi = 228; $ymini = 52; $ymaxi = 72;
  $Rinset = "-R$xmini/$xmaxi/$ymini/$ymaxi";
  $Binset = "-B1000wesn"; 
  $xceni = ($xmini + $xmaxi)/2;
  $yceni = ($ymini + $ymaxi)/2;
  $Jinset = "-Jb$xceni/$yceni/$ymini/$ymaxi/1:$Jbscalei";
  $inset_box_info = "-W2p,0/0/255";

  print CSH "psbasemap $Jinset $Rinset $Binset $origin_inset -K -O -V >> $psfile\n";
  print CSH "pscoast $Jinset $Rinset $coast_info_inset $origin_inset -K -O -V >> $psfile\n";
  print CSH "psxy $platefile $Jinset $Rinset -W1p,255/0/0 -m -K -V -O $origin_inset >> $psfile\n";
  print CSH "psxy $faultfile $Jinset $Rinset -W1p,255/0/0 -m -K -V -O $origin_inset >> $psfile\n";
  $cities = "$dir0/cities/alaska_cities_sub"; $city_info = "-Sc7p -W1.0p/0/0/0 -G255/255/255";
  print CSH "awk '{print \$2,\$1}' $cities | psxy ${city_info} $Jinset $Rinset -K -V -O $origin_inset >> $psfile\n";
  print CSH "psxy $Jinset $Rinset $inset_box_info -K -O -V $origin_inset >>$psfile<<EOF\n $xmin $ymin\n$xmax $ymin\n$xmax $ymax\n$xmin $ymax\n$xmin $ymin\nEOF\n";
}

# title and subtitle
# NOTE: previous errors with pstext may cause title NOT to be plotted outside map region
if ($ititle==1) {
  #print CSH "pstext -N $R_title $J_title -K -O -V >>$psfile<<EOF\n $x_title $y_title $fsize_title1 0 $fontno CM $title1\nEOF\n";
  print CSH "pstext -N $R_title $J_title $otitle1 -K -O -V >>$psfile<<EOF\n 0 0 $fsize_title1 0 $fontno LM $title1\nEOF\n";
  print CSH "pstext -N $R_title $J_title $otitle2 -K -O -V >>$psfile<<EOF\n 0 0 $fsize_title2 0 $fontno LM $title2\nEOF\n";
}

# scale bar in km
# Must specify at least one of -C, -G, -S, -I, -N, and -W
if ($iscalebar==1) {
  print CSH "pscoast $J $R $sbarinfo -N1 -K -O -V >> $psfile\n";  # add -N1 to plot borders
}

#==================================================
}   # pp loop
#==================================================

# print nothing to finish
print CSH "pstext -N $R_title $J_title $otitle1 -O -V >>$psfile<<EOF\n 0 0 $fsize_title1 0 $fontno CM\nEOF\n"; # FINISH
print CSH "echo done with $psfile\n";

#==================================================
}   # xx loop
#==================================================

  close (CSH);
  system("csh -f $cshfile");

  #print "convert $psfile $jpgfile\n";
  #system("convert $psfile -rotate $rotangle $jpgfile");
  #if ($ipdf==1) {system("ps2pdf $psfile")}
  #if ($ixv==1) {system("gv $psfile &")}

#==================================================

# NOTE: make sure that the executable is there
#@pdcat = "/home/carltape/bin/pdcat -r";
@pdcat = "/usr/local/bin/pdcat -r";
$k = 1;

for ($xx = $xpmin; $xx <= $xpmax; $xx++) {
  $clab = $caplabs[$xx-1];
  $fname = "${eid}_${xx}_$clab";
  $psfile = "$fname.ps";
  system("ps2pdf $psfile");

  $pdffile = "$fname.pdf";
  $pdcat[$k] = $pdffile; $k = $k+1;
}

$ofile = "allcap_${eid}.pdf";
print "output file is $ofile\n";
$pdcat[$k] = "./$ofile";
print "@pdcat\n";
#die("TESTING");
`@pdcat`;			# execute

#==================================================
