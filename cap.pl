#!/usr/bin/env perl

# if on pacman uncomment the following
#use v5.10;
#
# A user-friendly PERL interface to the CAP source inversion code cap
#
# written by Lupei Zhu, 3/6/1998, Caltech
# 
# revision history
#	6/18/2001	add usage and documentation.
#	11/05/2012	add isotropic and CLVD search (-J).
#

use File::Copy;

# these are the only things one need to change based on the site installation
$home = $ENV{HOME};                        # my home directory
$caphome = $ENV{CAPHOME};                  # CAP home directory
$caprun = $ENV{CAPRUN};                    # run directory

require "$caphome/cap_plt.pl";             # include plot script
require "$caphome/sub_read_parameter_file.pl";  #  read parameter file

#================defaults======================================
$cmd = "cap";

$inp_cmd = "inp_cmd";

# green's function location
# within each of these directories are subdirectories of models (cus, scak, utuhalf, etc)
#$green = "$home/data/models/Glib";               # original
$green = "/store/wf/FK_synthetics";               # UAF linux network
#$green = "/import/c/d/ERTHQUAK/FK_synthetics ";  # UAF cluster
#$green = "$caprun/models";                       # user testing
#$green = " /home/alvizuri/PROJECTS/CAP/models";  # western US model (temporary location: move to /store)

$repeat = 0;
$fm_thr = 0.01;
$dirct='';
$disp=1;
$mltp=0;
$pol_wt=999;
$weight="weight.dat";
$fmt_flag="false";

# plotting
$plot = 0;
$amplify = 1;
$spib = 40; # sec per inch, body waves
$spis = 45; # spi, surface waves
$keep = 0;
$rise = 0.5;
$ampfact = 1;

# filters and window lengths
($f1_pnl, $f2_pnl, $f1_sw, $f2_sw, $m1, $m2) = (0.02,0.2,0.02,0.1,35,70);

# max. shifts
$max_shft1=1;		# max. shift for Pnl
$max_shft2=5;		# max. shift for surface wave
$tie = 0.5;		# tie between SV and SH

# weights between different portions
$weight_of_pnl=2;		# weight for pnl portions
$power_of_body=1;		# distance scaling power for pnl waves
$power_of_surf=0.5;

# apparent velocities
#($vp, $love, $rayleigh) = (7.8, 3.5, 3.1);
($vp, $love, $rayleigh) = (-1, -1, -1);

# search types
$grid_type = -1.0;

# minimization (norm)
$norm = 1;

# for sorting the output file by distance or azimuth
$isort = 0;

# Use parameter file instead
$parameter_file = '';

#----------------------------------------------------------- 
# DEFAULT VALUES FOR SEARCH PARAMETERS
use Math::Trig ':pi';
$deg2rad = pi / 180.0;

# RANGE LIMITS PER PARAMETER
# magnitude 
($deg, $dm, $dlune) = (10, 0.1, 0.);
$nmw = 1;

# orientation
$k1 = 0;    $k2 = 360;
$h1 = 0;    $h2 = 1;
$s1 = -90;  $s2 = 90;

# lune --> (v,w) 
# v = [-1/3, 1/3] w = [-3pi/8, 3pi/8]
$w2 = (3.0 * pi / 8.0); $w1 = -$w2;
$v2 = (1.0 / 3.0);      $v1 = -$v2;

# NUMBER OF POINTS PER PARAMETER

# search type = GRID -- specify number of points for each parameter
# case 1. full moment tensor. use nv/nw/nstr/ndip/nrak
# case 2. fixed lambda. use nstr/ndip/nrak
# default number of points for each parameter

# for nX, X>1 (at least one parameter value has to exist!)
$nv = 1; $dv = 10;    # gamma
$nw = 1; $dw = 10;    # delta
$nk = 1; $dk = 10;    # strike
$nh = 1; $dh = 10;    # dip
$ns = 1; $ds = 10;    # rake

# magnitude default is to run a single point
$nmw = 1;
$dmw = 0;

# search type = RAND -- specify number of solutions to generate
$nsol    =  100000;     # full moment tensor
$nsol_fixlam =  100000;     # fixed lambda (includes DC)

# number of entries for flag R
$nR = 0;    # if not specified then it's the full range

# flag for old grid
$oldgrid = 0;   # default is oldgrid=0 (ie do not use old grid)

#----------------------------------------------------------- 

# number of freedom per sample for estimating uncertainty
$nof = 0.01;
$dt = 0.1;

# rms thresholds for discarding bad traces
@thrshd = (10., 10., 10., 10., 10.);

# paramters for running cap for multiple depths
$dep_inc = 0;

# command line input, [] means optional, see () for default value
$usage = 
" ===== CAP source inversion using seismic waveforms ====
	Ref: Zhu and Helmberger, 1996, BSSA 86, 1645-1641

  Data preparation:
     Put all three-component waveforms station.[r,t,z] of the event in
  a single directory named with the event ID. The data should be velocity
  in cm/s or displacement in cm in the SAC format, with the reference
  time set at the origin time and epicentral distance and azimuth
  set in the SAC header. There should be another file called $weight
  in the same directory, in the following format:
	station_name dist w1 w2 w3 w4 w5 tp tp_w ts ts_w ti
  where dist specifies the names of Green functions (dist.grn.?) to be used.
  w1 to w5 are the weights for 5 segments of waveforms: PnlZ, PnlR, Z, R, T.
  tp is first P arrival time if it's set to a positive value. tp_w is the body wave window.
  ts is the arrival time for surface waves. ts_w is the surface wave window. ti is the initial
  time shift for the surface waves, positive means that the data is delayed w.r.t. the model.
  If w2 is set to -1, it indicates that the station is at teleseimic distances and only
  the P (PnlZ) and SH (T) are used. In this case, ts is the S arrival time when it is positive.

  The Greens function library:
     The Greens functions are computed using FK, named as xxx.grn.[0-8] where
  xxx is the distance. All Greens functions from one source depth are placed
  in a single directory named as model_depth. They are in SAC format with
  two time marks set: t1 for the first P arrival and t2 for the first S arrival.
  If first-motion data are to be used in the inversion, the Greens functions
  need to have user1 and user2 set as the P and S take-off angles (in degrees from down).

  Time window determination:
     The inversion breaks the whole record into two windows, the Pnl window
  and the surface wave window. These windows are determined in following way:
    1) If the SAC head has time mark t1, t2 set, the code will use them for
       the Pnl window. The same is true for the surface wave window (using t3 and t4).
    Otherwise,
    2) If positive apparent velocities are given to the code (see -V below), it will use
       them to calculate the time windows:
	  t1 = dist/vp - 0.3*m1, t2 = ts + 0.2*m1
	  t3 = dist/vLove - 0.3*m2, t4 = dist/vRayleigh + 0.7*m2
    Otherwise,
    3) Using the tp, ts in the Green function header
 	  t1 = tp - 0.2*m1,  t2 = t1+m1
	  t3 = ts - 0.3*m2,  t4 = t3+m2
    Here m1, m2 are the maximum lengths for the Pnl and surface waves windows
    (see the -T options below).

=====================================================================================================
  Usage: cap.pl -Mmodel_depth/mag [-A<dep_min/dep_max/dep_inc>] [-C<f1_pnl/f2_pnl/f1_sw/f2_sw>]
                  [-D<w1/p1/p2>] [-F<thr>] [-Ggreen] [-Hdt] 
                  [-I<nsol> OR -I<Nv/Nw/Nstrike/Ndip/Nrake>]
                  [-k1 (old grid setup)] [-L<tau>]
                  [-M$model_$dep] [-m$mw OR -m<mw1>/<mw2>/<dmw> ] [-N<n>]
                  [-O] [-P[<Yscale[/Xscale_b[/Xscale_s[/k]]]]>] [-Qnof]
                  [-R<v0/w0/strike0/dip0/rake0> OR -R<v1/v2/w1/w1/strike1/strike2/dip1/dip2/rake1/rake2>] 
                  [-S<s1/s2[/tie]>] [-T<m1/m2>]
                  [-Udirct] [-V<vp/vl/vr>] [-Wi] [-Y<norm>] [-Zstring] event_dirs

    -A  run cap for different depths. (dep_min/dep_max/dep_inc).
    -B  FLAG NOT IN USE.
    -C  filters for Pnl and surface waves, specified by the corner
        frequencies of the band-pass filter. ($f1_pnl/$f2_pnl/$f1_sw/$f2_sw).
    -D	weight for Pnl (w1) and distance scaling powers for Pnl (p1) and surface
        waves (p2). If p1 or p2 is negative, all traces will be normalized. ($weight_of_pnl/$power_of_body/$power_of_surf).
    -E  FLAG NOT IN USE.
    -F	include first-motion data in the search. thr is the threshold ($fm_thr).
    	The first motion data are specified in $weight. The polarities
        can be specified using +-1 for P, +-2 for SV, and +-3 for SH after
        the station name, e.g. LHSA/+1/-3 means that P is up and SH is CCW.
        The Green functions need to have take-off angles stored in the SAC
        header.
        threshold should be NEGATIVE if polarities are allowed to conflict expected polarity (as mentioned in weight file).
    -G  Green's function library location ($green).
    -H  dt ($dt).
    -I  specify number of solutions (random mode) OR number of points per parameter.
        RAND: -I<nsol>  e.g. -I10000  --- will generate 10,000 random solutions.
        GRID: -I<Nv>/<Nw>/<Nstrike>/<Ndip>/<Nrake> where Nx = number of poits for parameter x
    -J  FLAG NOT IN USE.
    -K  sort output file by distance (=1) or azimuth (=2). default (=0) does nothing.
    -k  Specify k1 to build a lune grid, not a UV grid.
        Use only when LUNE_GRID_INSTEAD_OF_UV = 0 (and recompile cap)
    -L  source duration (estimate from mw, can put a sac file name here).
    -M	specify the model and source depth.
    -m	specify point magnitude: -m<mw0> OR magnitude range: -n<mw1>/<mw2>/<dmw>
    -N  repeat the inversion n times and discard bad traces ($repeat).
    -O  output CAP input (off).
    -P	generate waveform-fit plot with plotting scale. ([-P<Yscale>/<XscaleBody>/<XscaleSurf>](/k))
    	Yscale: amplitude in inch for the first trace of the page ($amplify).
        # scale down: -P1e-5, -P1e-6,... seismograms scaled by this amplitude (smaller the number larger the waveform)
        # scale up:   -P1, -P2,...       seismograms also scaled by amplitude but then enlarged by factor 1, 2, (anything greater than 1)... (larger the number larger the waveform)
        # scale each window:  -P0.5e+0.5 seismograms will be scaled for each component window (all same size) - Normalized scaling
        Xscale: seconds per inch. (body: $spib, surface:$spis).
        append k if one wants to keep those waveforms.
        -p  (small p) For amplitude scaling of surface waves; example: If set to 2 surface waves amplitude will be multipled by twice 
    -Q  number of freedom per sample ($nof)
    -R	For double couple use -R0/0.
        For point solution use -Rv/w/k/h/s
        For search range use -Rv1/v2/w1/w2/k1/k2/h1/h2/s1/s2
        Note: This should come after -I flag in the command or it crashes sometimes! 
    -S	max. time shifts in sec for Pnl and surface waves ($max_shft1/$max_shift2) and
        tie between SH shift and SV shift:
        tie=0 		shift SV and SH independently,
        tie=0.5 	force the same shift for SH and SV ($tie).
    -T	max. time window lengths for Pnl and surface waves ($m1/$m2).
    -U  directivity, specify rupture direction on the fault plane (off).
    -V	apparent velocities for Pnl, Love, and Rayleigh waves (off).
    -W  Integration.
        W$disp => (default) Do not integrate. Uses data in its original form
        W1 => Integrate
    -X  weight for normalized polarity misfit [0,1). This is combined with the waveform misfit. Don't use -X1 since there could be multiple solutions that could fit the observed polarity, use -X.99 instead in order to include atleast some waveform measure. Suggested value = 0.5 (Bug: Crashes when using -X flag and no polarity is available in the weight file)
    -Y  specify norm (1 - L1 norm; 2 - L2 norm)
    -Z  specify a different weight file name ($weight).

-----------------------------------------------------------------------------------------------------
Weight file description:
Column No:        1                2          3           4            5            6         7             8                   9               10                  11               12
Description:   Station_Name    Distance   PV_weight   PR_weight    SV_weight   SR_weight   SH_weight   P_arrival_time    P_window_length   S_arrival_time    S_window_length   shift_synthetics

Convention: Postive time-shift means synthetics is arriving earlier (faster velocity model) and it needs to be shifted in the positive direction in order to match it with the data.

=====================================================================================================
Examples:
RANDOM SEARCH
Double Couple:
> cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15 -m4.5/5.5/0.1 -R0/0 -Y1 -I100000 20080418093700
FMT search:
> cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15 -m4.5/5.5/0.1 -R0/0 -Y1 -I100000 20080418093700

GRID SEARCH
Double Couple:
> cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15 -m4.5/5.5/0.1 -R0/0 -Y1 -I1/1/37/10/19 20080418093700
FMT search:
> cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15 -m4.5/5.5/0.1 -Y1 -I10/10/37/10/19 20080418093700

The commands abouve find the best focal mechanism and moment magnitude of the
2008/4/18 Southern Illinois earthquake 20080418093700 using the central US
crustal velocity model cus with the earthquake at a depth of 15 km.
Here we assume that the Greens functions have already been computed and saved
in $green/cus/cus_15/.

DEPTH SEARCH
To find the best focal depth, repeat the inversion for different focal depths
either by using a for-loop or the '-A' flag
# for loop with depths 5 to 30 km at 5-km intervals
> for h in 5 10 15 20 25 30; do cap.pl -H0.2 -P1 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_$h/5.0 -E0 -K0 -Y2 20080418093700; done
# same as above but using the A flag in the format -Astart/end/incr
> cap.pl -H0.2 -P0.6 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15 -m4.5/5.5/0.1 -R0/0 -Y1 -I1/1/37/10/19 20080418093700 -A5/30/5

# to plot the depth test
> depth_test 20080418093700 cus
> gv dep_20080418093700.ps
-------------------------------------------------------------------------------------
[Header Info] 
The inversion results are saved in cus_15.out
Event 20080418093700 Model cus_015 FM  297 86.815262    0 Mw 5.20 rms 4.547e-05   112 ERR   0   0   0 CLVD -1.89 -nan ISO  10.212961 0.00 VR 80.3 data2 1.026e-04
# Hypocenter_sac_header elat 3.845000e+01 elon -8.789000e+01 edep 1.160000e+01
# tensor = 7.943e+23  0.9511 -0.5865 -0.0273 -0.6243  0.0474  0.1075
# norm L1    # Pwin 35 Swin 70    # N 8 Np 16 Ns 24

saying that the fault plane solution is strike 297, dip 87, rake 0, gamma (CLVD) -2, and delta (ISO) 10 degrees.
Also provided are misfit (rms), data norm (Data2) and Variance reduction (VR), and hypocenter location
tensor in Mxx Mxy Mxz Myy Myz Mzz (where x=North, y=East, z=Down).
norm info (L1 or L2); Duration of P (Pwin) and S windows; Number of stations (N); Number of P components (Np) and Number of S components (Ns)
  The rest of the files shows rms, cross-correlation coef., and time shift of individual waveforms.
  The waveform fits are plotted in file cus_15.ps in the event directory.
------------------------------------------------------------------------------------------------------

";

@ARGV > 0 || die $usage;

$ncom = 5;	# 5 segments to plot

open(INP,">$inp_cmd");
print INP "cap.pl ";
foreach $argnum (0 .. $#ARGV) {
    print INP "@ARGV[$argnum] ";}
close(INP);
system("chmod +x $inp_cmd");

#input options
foreach (grep(/^-/,@ARGV)) {
   $opt = substr($_,1,1);
   @value = split(/\//,substr($_,2));
   if ($opt eq "A") {
     $dep_min = $value[0];
     $dep_max = $value[1];
     $dep_inc = $value[2];
     if ($#value ==2){
       printf STDERR "Running cap for multiple depths: $dep_min to $dep_max at $dep_inc km increment\nWarning: overwriting -Mdepth\n";
     } else {
       $dep_inc=0;
       printf STDERR "Depth run flag -A not specified correctly\nUsing -Mdepth instead\n---------------------\n"; }
   } elsif ($opt eq "C") {
     ($f1_pnl, $f2_pnl, $f1_sw, $f2_sw) = @value;
   } elsif ($opt eq "D") {
     ($weight_of_pnl,$power_of_body,$power_of_surf)=@value;
   } elsif ($opt eq "F") {
     $fm_thr = $value[0] if $#value >= 0;
   } elsif ($opt eq "G") {
     $green = substr($_,2);
   } elsif ($opt eq "H") {
     $dt = $value[0];
   } elsif ($opt eq "I") {
       # Parameter for number of points or number of solutions
       # Two options.
       #    Length 1 = (RAND) number of solutions (nsol)
       #    Length 5 = (GRID) nv/nw/nk/nh/ns. nsol = nv*nw*nk*nh*ns
       #    
#      $deg = $value[0];
#      $dm = $value[1] if $#value > 0;
#      $dlune = $value[2] if $value[2]
       if ($#value==0) {
           $nI = 1;
           $nsol = $value[0];
           $nv = $nw = $nk = $nh = $ns = $nsol;
       } elsif ($#value==4) {
           $nI = 5;
           ($nv, $nw, $nk, $nh, $ns) = @value;
           $nsol = $nv * $nw * $nk * $nh * $ns;
           if ($nsol <= 0) {
               print STDERR "ERROR. Number of points per paramete should be >0.\n";
               exit(0);
           }
       }
#  elsif ($opt eq "J") {
#    $iso1   = $value[0] if $value[0];
#    $iso2  = $value[1] if $value[1];
#    $clvd1  = $value[2] if $value[2];
#    $clvd2 = $value[3] if $value[3];
#    $fmt_flag="true";     # used later for renaming output figures with fmt
#  }
   } elsif ($opt eq "K") {
     $isort = $value[0];
   } elsif ($opt eq "k") {
     $oldgrid = $value[0];
   } elsif ($opt eq "L") {
     $dura = join('/',@value);
   } elsif ($opt eq "m") {
       # Search parameter for magnitude
       # Two options.
       # 1. Length 1 = fixed magnitude 
       # 2. Length 3 = search over magnitude range 
       # 
       #($md_dep,$mg) = @value;
       #
       if ($#value==0) {
           ($mw1, $mw2, $dm) = ($value[0], $value[0], 0);
           $mg = $value[0];     # IS THIS NEEDED ANYMORE?
       } elsif ($#value==2) {
           ($mw1, $mw2, $dm) = @value;
           # nmw = number of magnitude points when running magnitude search.
           # note this allows to run at grid points dmw finer than 
           # 0.1 (eg 0.01, 0.001...), so run with care
           $nmw = sprintf("%.0f", ($mw2 - $mw1) / $dm +1 );
       }
   } elsif ($opt eq "M") {
       # ($md_dep,$mg) = @value;
       $md_dep = @value[0];
   } elsif ($opt eq "N") {
        $repeat = $value[0];
   } elsif ($opt eq "O") {
     $cmd = "cat";
   } elsif ($opt eq "P") {
     $plot = 1;
     $amplify = $value[0] if $#value >= 0;
     $spib = $value[1] if $value[1] > 0;
     $spis = $value[2] if $value[2] > 0;
     $keep = 1 if $#value > 2;
   } elsif ($opt eq "p") {
     $ampfact = $value[0];
   } elsif ($opt eq "Q") {
     $nof = $value[0];
   } elsif ($opt eq "R") {
       # Flag to set Ranges of parameters
       # Four options: 
       # 1. no flag  = FMT over full range
       # 2. Length 2 = fixed eigenvalue with grid search (v0/w0)
       # 3. Length 5 = fixed moment tensor (v0/w0/k0/h0/s0)
       # 4. Length 10 = subset case (v1/v2/w1/w2/k1/k2/h1/h2/s1/s2)
       if ($#value==1) {
           $nR = 2;
           ($v1, $v2) = @value;
           ($w1, $w2) = @value;
           $nv = $nw = 1;       # at least one lune point needs to exist
       } elsif ($#value==4) {
           $nR = 5;
           ($v1, $w1, $k1, $h1, $s1) = @value;
           ($v2, $w2, $k2, $h2, $s2) = @value;
           $h1 = $h2 = cos($value[3]*$deg2rad);  # cap expects h = cos(dip)
           #$h1 = $h2 = $value[3];         # cap expects dip
           $nsol = $nv = $nw = $nk = $nh = $ns = 1;
       } elsif ($#value==9) {
           $nR = 10;
           ($v1, $v2, $w1, $w2, $k1, $k2, $h1, $h2, $s1, $s2) = @value;
       }
   } elsif ($opt eq "S") {
     ($max_shft1, $max_shft2) = @value;
     $tie = $value[2] if $#value > 1;
   } elsif ($opt eq "T") {
     ($m1, $m2) = @value;
   } elsif ($opt eq "U") {
     ($rupDir) = @value;
     $pVel = 6.4;
     $sVel = 3.6;
     $rise = 0.4;
     $dirct = "_dir";
   } elsif ($opt eq "V") {
     ($vp, $love, $rayleigh) = @value;
   } elsif ($opt eq "W") {
     $disp = $value[0];
   } elsif ($opt eq "X") {
     $pol_wt = $value[0];
   } elsif ($opt eq "Y") {
     $norm = $value[0];
   } elsif ($opt eq "Z") {
     $weight = $value[0];
   } else {
     printf STDERR $usage;
     exit(0);
   }
}
@event = grep(!/^-/,@ARGV);

#-----------------------------------------------------------
#  Read parameter file instead
if (@ARGV == 1) {
   $parameter_file = $ARGV[0];
   sub_read_parameter_file($parameter_file)
}
#-----------------------------------------------------------

#-----------------------------------------------------------
# prepare additional parameters for CAP input
#-----------------------------------------------------------

# start to output some values
print STDERR "-------------------------------------------------------------\n";
print STDERR "cap.pl: input parameters for CAP:\n";

# convert strike/dip/rake to radians
if ($oldgrid == 0) {
    $k1 = $k1 * $deg2rad;
    $k2 = $k2 * $deg2rad;
    $s1 = $s1 * $deg2rad;
    $s2 = $s2 * $deg2rad;
}

unless ($dura) {
  $dura = int(10**(($mg-5)/2)+0.5);
  $dura = 1 if $dura < 1;
  $dura = 9 if $dura > 9;
}

# Flag to create regular grid as in Alvizuri & Tape (2016) and Silwal & Tape (2016).
# NOTE reproducibility may not be exact since grid spacing uses function gridvec (in the c code).
# Function gridvec does not implement the discretization of the previous version of 
# cap.c which uses rules to account for special grid points.
# Function gridvec also avoids endpoints in all parameters.
if( ($oldgrid == 1) && ($nI == 5)) {
    # Old grid
    print STDERR "\nWarning. Using flag -k1 for the lune grid. cap.c should\n";
    print STDERR "\tbe compiled with flag LUNE_GRID_INSTEAD_OF_UV = 1\n\n";
    # check that K flag works with flag R
    if (($nv == 1) && ($nw == 1) && ($nR == 0)) {
        # default is double couple if Range not specified and it's a single lune point
        ($dv, $dw) = (0, 0);
        ($v1, $v2) = (0, 0);
        ($w1, $w2) = (0, 0);
        ($k1, $k2) = (  0, 360);
        ($h1, $h2) = (  0, 90);
        ($s1, $s2) = (-90, 90);
    } elsif (($nv == 1) && ($nw == 1) && ($nR == 5)) {
        print STDERR "Fixed solution. Input values: $v1 $w1 $k1 $h1 $s1\n";
        # # if Range is set then its a subset
        # # lune points come from user input
        ($dv, $dw) = (0, 0);
    } elsif (($nv == 1) && ($nw == 1) && ($nR > 5)) {
        # If range=5 then it's a point solution.
        # Lune poionts come from user input
        ($dv, $dw) = (0, 0);
        ($k1, $k2) = (  0, 360);
        ($h1, $h2) = (  0, 90);
        ($s1, $s2) = (-90, 90);
    } else {
        # set the full range
        ($v1, $v2) = (-30, 30);
        ($w1, $w2) = (-90, 90);
        ($k1, $k2) = (  0, 360);
        ($h1, $h2) = (  0, 90);
        ($s1, $s2) = (-90, 90);
    }
    # NOTE values nX are grid spacings, NOT number of points
    # NOTE values dX should be integers (CAP expects integers)
    $dk = sprintf("%.0f", (($k2 - $k1) / $nk));  # strike -- do not include 360
    $dh = sprintf("%.0f", (($h2 - $h1) / $nh));  # dip    -- include 0 at start (though it will be offset later)
    $ds = sprintf("%.0f", (($s2 - $s1) / $ns));  # rake   -- include 0 point

    #print STDERR "Input parameters:\n$dv $dw $dk $dh $ds\n";  # for debugging
    $nsol = $nv * $nw * $nk * $nh * $ns;
} elsif( ($oldgrid == 1) && ($nI == 1)) {
    # random mode
    print STDERR "Warning. Using the old random setup.\n";
    print STDERR "cap.c should be compiled with flag LUNE_GRID_INSTEAD_OF_UV = 1\n";
    if (($nv == 1) && ($nw == 1) && ($nR == 0)) {
        # default is double couple if Range not specified and it's a single lune point
        ($dv, $dw) = (0, 0);
        ($v1, $v2) = (0, 0);
        ($w1, $w2) = (0, 0);
        ($k1, $k2) = (  0, 360);
        ($h1, $h2) = (  0, 90);
        ($s1, $s2) = (-90, 90);
        $nk = $nh = $ns = $nsol;
    } elsif (($nv == 1) && ($nw == 1) && ($nR == 5)) {
        print STDERR "Fixed solution. Input values: $v1 $w1 $k1 $h1 $s1\n";
        # # if Range is set then its a subset
        # # lune points come from user input
        ($dv, $dw) = (0, 0);
    } elsif (($nv == 1) && ($nw == 1) && ($nR > 0)) {
        # if Range is set then its a subset
        # lune points come from user input
        ($dv, $dw) = (0, 0);
        ($k1, $k2) = (  0, 360);
        ($h1, $h2) = (  0, 90);
        ($s1, $s2) = (-90, 90);
        $nk = $nh = $ns = $nsol;
    } else {
        # set the full range
        ($v1, $v2) = (-30, 30);
        ($w1, $w2) = (-90, 90);
        ($k1, $k2) = (  0, 360);
        ($h1, $h2) = (  0, 90);
        ($s1, $s2) = (-90, 90);
        $nv = $nw = $nk = $nh = $ns = $nsol;
    }
}

# plots for the DC don't have "fmt" in their filenames
if (($v1 == 0) && ($w1 == 0)) {
    $fmt_flag="false";   # double couple
    print STDERR "NOTE computing a double couple solution\n";
} else {
    $fmt_flag="true";
}

#-----------------------------------------------------------
# CHECK THAT USER INPUT MAKE SENSE
#-----------------------------------------------------------

# Set grid type to build (uniform rand, uniform grid) depending on entries in I (previously flag K).
# If using flag k3 then also need to set LUNE_GRID_INSTEAD_OF_UV=1 in cap.c
if ($nI == 1) {
    $grid_type = 2; # RAND
    $grid_type_label="random"
} elsif ($nI == 5) {
    $grid_type = 1; # GRID
    $grid_type_label="grid"
} else{
    print STDERR "STOP. Unable to set search type, check entries in flag R or I\n";
    exit(0);
}

#$grid_type = $type;

# Flag I: set defaults
if (($nI == 5) && ($nv == 1) && ($nw == 1) && ($oldgrid == 0)) {
    # default to the double couple
    ($v1, $v2) = (0, 0);
    ($w1, $w2) = (0, 0);
}

if ( -r $dura ) {	# use a sac file for source time function   
  $dt = 0;
  $riseTime = 1;
} else {
  $riseTime = $rise*$dura;
}

($model, $depth) = split('_', $md_dep);
unless ($depth) {
  $model = ".";
  $depth = 1;
}

if ($dep_inc==0) {
  $dep_min=$depth;
  $dep_max=$depth;
  $dep_inc=1;
}

# convert filter frequencies to periods
$Tf1 = 1/$f1_pnl;
$Tf2 = 1/$f2_pnl;
$Tf3 = 1/$f1_sw;
$Tf4 = 1/$f2_sw;
$filterBand = sprintf("Body:%.2f-%.2f. Surf:%.2f-%.2f",$Tf2,$Tf1,$Tf4,$Tf3);

#-----------------------------------------------------------
# END prepare additional parameters for CAP input
#-----------------------------------------------------------

for($dep=$dep_min;$dep<=$dep_max;$dep=$dep+$dep_inc) {
  foreach $eve (@event) {

    $md_dep = $model.'_'.$dep;
    next unless -d $eve;
    print STDERR "EVENT ID = $eve | EVENT DEPTH = $dep |  SOURCE DURATION = $dura\n";
    print STDERR "GRID TYPE = $grid_type ($grid_type_label search) | NSOL = $nsol\n";
    print STDERR "-------------------------------------------------------------\n\n";

# --------------------------
# Sort weight file by distance or azimuth
# Default is to use the weight file as it is
    $input_weight_file = "$eve/$weight";
    $clean_weight_file = "$eve/WEIGHT_CLEAN.dat";
    $station_info_file = "$eve/${eve}_station_list_ALL.dat";

    open(IN,$input_weight_file) || die "weight file not present\n";
    @weightlines = <IN>; $nsta = @weightlines;
    close(IN);

    if ($isort != 0){

        if ($isort == 1){
            print "Ouput file will be sorted by distance";
        }
        elsif ($isort == 2){
            print "Ouput file will be sorted by azimuth";
        }
        else {
            die "sorting can only be by distance (-K1) or azimuth (-K2)";
        }

        # Get data from station_list_ALL.dat
        # sta net   lat     lon         dist      azim
        # BKS BK 37.876200 -122.235600 518.062577 279.769407
        # CMB BK 38.034500 -120.386500 360.643229 285.609184
        # MHC BK 37.341600 -121.642600 462.460019 273.168208
        open(IN, $station_info_file) || die "STOP. {eid}/{eid}_station_list_ALL.dat file not present\n";
        @stnlines = <IN>;
        $nstn = @stnlines;
        close(IN);

        # Read specific columns from the weight file
        for ($i = 0; $i < $nsta; $i++){
            ($name,$dist,$pv,$pr,$sv,$sr,$st,$ptime,$plen,$stime,$slen,$shift)=split(" ",@weightlines[$i]);
            ($stnm,$pol) = split("/",$name);
            ($eve1,$net1,$name1,$loc1,$chan1) = split(/\./,$stnm);

            # sort by distance, azimuth, or leave as is
            for ($j = 0; $j < $nstn; $j++) {
                # read columns from station_list_ALL.dat
                ($name2, $net2, $lat2, $lon2, $dist2, $az2) = split(" ",@stnlines[$j]);
                if (($name1 eq $name2) && ($net1 eq $net2)) {
                    # by distance
                    if ($isort == 1) {
                        $col_to_sort[$i] = $dist2;
                    }
                    # by azimuth
                    elsif ($isort == 2) {
                        $col_to_sort[$i] = $az2;
                    }
                    # THIS CODE IS NEVER REACHED
                    # do nothing
                    else {
                        $col_to_sort[$i] = $i;
                    }
                }
            }
        }
        # GET INDICES OF SORTED ELEMENTS IN station_list_ALL.dat
        @sort_indx = sort{$col_to_sort[$a] <=> $col_to_sort[$b]} 0 .. $#col_to_sort;

        #-----------------------------------------------------------
        # Remove stations that have no information. 
        # DO THIS REGARDLESS. AT PRESENT THE SCRIPT SEEMS TO TAKE A STATION
        # INTO ACCOUNT WHETHER THE STATION IS USED OR NOT
        #-----------------------------------------------------------
        # NO polarity, NO P weight, NO S weight
        open(OUT,'>',$clean_weight_file);
        for ($i = 0; $i < $nsta; $i++) {
            $ipol = 1;
            ($name,$dist,$pv,$pr,$sv,$sr,$st,$ptime,$plen,$stime,$slen,$shift) = split(" ",@weightlines[$sort_indx[$i]]);
            ($stnm,$pol) = split("/",$name);

            # no polarity information
            if ($pol eq ''){
                $ipol = 0;
            }  
            # If pol or weight checks fail then keep the station
            # Skip IF and(pol checks) and(weight checks)
            if (($ipol==0 || $pol_wt==0 || $pol_wt==999) && $pv==0 && $pr==0 && $sv==0 && $sr==0 && $st==0) {
                print STDERR "$name \t $dist \t $pv \t $pr \t $sv \t $sr \t $st \t $ptime \t $plen \t $stime \t $slen \t $shift \n";
                next;
            } 
            # save in new weight file
            else {
                print OUT "$name \t $dist \t $pv \t $pr \t $sv \t $sr \t $st \t $ptime \t $plen \t $stime \t $slen \t $shift \n";
            }
        }
        close(OUT);

#        # Clean this weight file by removing stations that have no information (no polarity, no P weight, no S weight)
#        open(OUT,'>',$clean_weight_file);
#        @sort_indx = sort{$col_to_sort[$a] <=> $col_to_sort[$b]} 0 .. $#col_to_sort;
#        for ($i = 0; $i < $nsta; $i++){
#            $ipol = 1;
#            ($name,$dist,$pv,$pr,$sv,$sr,$st,$ptime,$plen,$stime,$slen,$shift)=split(" ",@weightlines[$sort_indx[$i]]);
#            ($stnm,$pol) = split("/",$name);
        #
#            if ($pol eq ''){$ipol = 0;}  # no polarity information
#            if (($ipol==0 || $pol_wt==0 || $pol_wt==999) && $pv==0 && $pr==0 && $sv==0 && $sr==0 && $st==0){
#                print STDERR "$name \t $dist \t $pv \t $pr \t $sv \t $sr \t $st \t $ptime \t $plen \t $stime \t $slen \t $shift \n";
#                next;} # No information available - skip this station
#            else {     # save in new weight file
#                print OUT "$name \t $dist \t $pv \t $pr \t $sv \t $sr \t $st \t $ptime \t $plen \t $stime \t $slen \t $shift \n";
#            }
#        }
#        close(OUT);
    }
    else {
        copy $input_weight_file, $clean_weight_file;
    }
# --------------------------

    open(WEI,$clean_weight_file);
    @wwf=<WEI>;
    close(WEI);
    $ncom = 2 if $wwf[0] =~ / -1 /;

    $cmd = "cap$dirct $eve $md_dep" unless $cmd eq "cat";    
    open(SRC, "| $cmd") || die "can not run $cmd\n";
    print SRC "$pVel $sVel $riseTime $dura $rupDir\n",$riseTime if $dirct eq "_dir";
    print SRC "$model $dep\n";          # first input in regular cap run
    print SRC "$m1 $m2 $max_shft1 $max_shft2 $repeat $fm_thr $tie\n";
    print SRC "@thrshd\n" if $repeat;   # no value in regular cap run
    print SRC "$vp $love $rayleigh\n";
    print SRC "$power_of_body $power_of_surf $weight_of_pnl $nof\n";
    print SRC "$plot\n";
    print SRC "$disp $pol_wt\n";
    print SRC "$green/$model/\n";
    print SRC "$grid_type\n";
    print SRC "$norm\n";
    print SRC "$dt $dura $riseTime\n";
    print SRC "$f1_pnl $f2_pnl $f1_sw $f2_sw\n";
    print SRC "$mw1 $mw2 $nmw $dm\n";
    print SRC "$v1 $v2 $nv $dv\n";
    print SRC "$w1 $w2 $nw $dw\n";
    print SRC "$k1 $k2 $nk $dk\n";
    print SRC "$h1 $h2 $nh $dh\n";
    print SRC "$s1 $s2 $ns $ds\n";
    print SRC "$nsol\n";
    printf SRC "%d\n",$#wwf + 1;
    print SRC @wwf;
    close(SRC);

    #-----save a copy of inpur command and weight file in the OUTPUT_DIR
    system("cp", $input_weight_file, './OUTPUT_DIR/weight.dat');
    system("cp", $inp_cmd, "./OUTPUT_DIR/${eve}_${model}_${dep}_caprun");
    system("git log | head -12 > ./OUTPUT_DIR/last_2git_commits.txt");

  plot:
    if ( $plot > 0 && ($? >> 8) == 0 ) {
      print STDERR "\n\ncap.pl: plot results ... \n";
      $odir = "./OUTPUT_DIR";
      chdir($odir);
      @dum = split('_', $md_dep);  # split mdl string
      $outfile = sprintf("%s_%s_%03d.out", @event, $model, int($dep));
      open(my $out,'>>',$outfile);
      say $out "INPUT_PAR $md_dep P_win $m1 S_win $m2 P $amplify p $ampfact NCOM $ncom spiB $spib spiS $spis $filterBand FMT $fmt_flag";

      &plot($md_dep, $m1, $m2, $amplify, $ampfact, $ncom, $spib, $spis, $filterBand, $fmt_flag, @event, $model, $dep, $dura, $riseTime, $pol_wt);
      unlink(<${md_dep}_*.?>) unless $keep;
      chdir("../");
      print STDERR "cap.pl: plotting finished.\n";
    } else {
      print STDERR "cap.pl: no plots generated.\n";
    }
  }
}
exit(0);
