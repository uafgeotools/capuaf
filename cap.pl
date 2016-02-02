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

# these are the only things one need to change based on the site installation
$home = $ENV{HOME};                        # my home directory
$caphome = $ENV{CAPHOME};                  # CAP home directory
$caprun = $ENV{CAPRUN};                    # run directory

require "$caphome/cap_plt.pl";             # include plot script

#================defaults======================================
$cmd = "cap";

$inp_cmd = "inp_cmd";

# green's function location
# within each of these directories are subdirectories of models (cus, scak, utuhalf, etc)
#$green = "$home/data/models/Glib";               # original
$green = "/store/wf/FK_synthetics";               # UAF linux network
#$green = "/import/c/d/ERTHQUAK/FK_synthetics ";  # UAF cluster
#$green = "$caprun/models";                       # user testing

$repeat = 0;
$bootstrap = 0;
$fm_thr = 0.01;
$dirct='';
$disp=0;
$mltp=0;
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
$parm = -1.0;
$type = -1.0;

# minimization (norm)
$norm = 2;

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

$nv = $dv = 10;    # gamma
$nw = $dw = 10;    # delta
$nk = $dk = 10;    # strike
$nh = $dh = 10;    # dip
$ns = $ds = 10;    # rake

# magnitude default is to run a single point
$nmw = 1;   
$dmw = 0;

# search type = RAND -- specify number of solutions to generate
$nsol    =  100000;     # full moment tensor
$nsol_fixlam =  100000;     # fixed lambda (includes DC)

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
  Usage: cap.pl -Mmodel_depth/mag [-A<dep_min/dep_max/dep_inc>] [-B] [-C<f1_pnl/f2_pnl/f1_sw/f2_sw>]
                  [-D<w1/p1/p2>] [-E<search>] [-F<thr>] [-Ggreen] [-Hdt] [-Idd[/dm]]
                  [-J[iso[/diso[/clvd[/dclvd]]]]] [-K<search_type>] [-L<tau>] [-M$model_$dep/$mw][-N<n>]
                  [-O] [-P[<Yscale[/Xscale_b[/Xscale_s[/k]]]]>] [-Qnof]
                  [-R<strike1/strike2/dip1/dip2/rake1/rake2>] [-S<s1/s2[/tie]>] [-T<m1/m2>]
                  [-Udirct] [-V<vp/vl/vr>] [-Wi] [-Xn] [-Y<norm>] [-Zstring] event_dirs

    -A  run cap for different depths. (dep_min/dep_max/dep_inc).
    -B  output misfit errors of all solutions for bootstrapping late ($bootrap).
    -C  filters for Pnl and surface waves, specified by the corner
	frequencies of the band-pass filter. ($f1_pnl/$f2_pnl/$f1_sw/$f2_sw).
    -D	weight for Pnl (w1) and distance scaling powers for Pnl (p1) and surface
   	waves (p2). If p1 or p2 is negative, all traces will be normalized. ($weight_of_pnl/$power_of_body/$power_of_surf).
    -E  Specify what kind of parameterization (0=Zhu; 1=Lune)
    -F	include first-motion data in the search. thr is the threshold ($fm_thr).
    	The first motion data are specified in $weight. The polarities
	can be specified using +-1 for P, +-2 for SV, and +-3 for SH after
	the station name, e.g. LHSA/+1/-3 means that P is up and SH is CCW.
	The Green functions need to have take-off angles stored in the SAC
	header.
        threshold should be NEGATIVE if polarities are allowed to conflict expected polarity (as mentioned in weight file).
    -G  Green's function library location ($green).
    -H  dt ($dt).
    -I  search interval in orientation (strike/dip/rake) and mag(dm)and non-DC(dlune) ($dorient/$dm/$dlune).
        If dm<0, the gain of each station will be determined by inversion. $dlune is ommited in case of search=0.
    -J  (if search=0)include isotropic and CLVD search using steps diso and dclvd (0/0/0/0).
        (if search=1,2) iso and clvd search range (lune parameterization)
        example -J10/10/5/5 will perform a direct search; and -J-10/10/-5/5 will perform a grid search (or the random search) wihtin the subset
    -K  Kind of search (0=line search; 1=Grid search; 2=random search)
    -L  source duration (estimate from mw, can put a sac file name here).
    -M	specify the model, source depth and initial magnitude.
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
    -R	grid-search range for strike/dip/rake (0/360/0/90/-90/90).
    -S	max. time shifts in sec for Pnl and surface waves ($max_shft1/$max_shift2) and
	tie between SH shift and SV shift:
	 tie=0 		shift SV and SH independently,
	 tie=0.5 	force the same shift for SH and SV ($tie).
    -T	max. time window lengths for Pnl and surface waves ($m1/$m2).
    -U  directivity, specify rupture direction on the fault plane (off).
    -V	apparent velocities for Pnl, Love, and Rayleigh waves (off).
    -W  use displacement for inversion; 1=> data in velocity; 2=> data in disp ($disp).
    -X  output other local minimums whose misfit-min<n*sigma ($mltp).
    -Y  specify norm (1 - L1 norm; 2 - L2 norm)
    -Z  specify a different weight file name ($weight).

=====================================================================================================
Examples:
> cap.pl -H0.2 -P0.3 -S2/5/0 -T35/70 -F -D2/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_15/5.0 20080418093700
  which finds the best focal mechanism and moment magnitude of tbe 2008/4/18 Southern Illinois earthquake
  20080418093700 using the central US crustal velocity model cus with the earthquake at a depth of 15 km.
  Here we assume that the Greens functions have already been computed and saved in $green/cus/cus_15/.
  The inversion results are saved in cus_15.out with the first line
Event 20080418093700 Model cus_15 FM 115 90  -2 Mw 5.19 rms 1.341e-02   110 ERR   1   3   4
  saying that the fault plane solution is strike 115, dip 90, and rake -2 degrees, with the
  axial lengths of the 1-sigma error ellipsoid of 1, 3, and 4 degrees.
  The rest of the files shows rms, cross-correlation coef., and time shift of individual waveforms.
  The waveform fits are plotted in file cus_15.ps in the event directory.
------------------------------------------------------------------------------------------------------
  To find the best focal depth, repeat the inversion for different focal depths:
> for h in 5 10 15 20 25 30; do cap.pl -H0.2 -P1 -S2/5/0 -T35/70 -F -D1/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_$h/5.0 -E0 -K0 -Y2 20080418093700; done
  and store all the results in a temporary file:
------------------------------------------------------------------------------------------------------
> grep -h Event 20080418093700/cus_*.out > junk1.out
> grep -h tensor 20080418093700/cus_*.out > junk2.out
  and then run
> ./depth.pl junk1.out junk2.out 20080418093700 > junk.ps
------------------------------------------------------------------------------------------------------
Instead of this you can run
> depth_test 20080418093700 cus
------------------------------------------------------------------------------------------------------
  The output from the above command
Event 20080418093700 Model cus_15 FM 115 90  -2 Mw 5.19 rms 1.341e-02   110 ERR   1   3   4 H  14.8 0.6
  shows that the best focal depth is 14.8 +/- 0.6 km.

  To include isotropic and CLVD in the inversion, use the -J option to specify the starting iso0, clvd0, and search steps. It requires
  that the Green's function library includes the explosion source components (.a, .b, .c).


";

@ARGV > 1 || die $usage;

$ncom = 5;	# 5 segemnts to plot

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
   } elsif($opt eq "B") {
     $bootstrap = 1;
   } elsif ($opt eq "C") {
     ($f1_pnl, $f2_pnl, $f1_sw, $f2_sw) = @value;
   } elsif ($opt eq "D") {
     ($weight_of_pnl,$power_of_body,$power_of_surf)=@value;
   } elsif ($opt eq "E") {
     $parm = $value[0];
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
           $nsol = $value[0];
           $nv = $nw = $nk = $nh = $ns = $nsol;
       } elsif ($#value==4) {
           ($nv, $nw, $nk, $nh, $ns) = @value;
           if ($nv == 0 || $nw == 0) {
               # NOTE! this only checks for (v,w) being zero
               $nsol = $nk * $nh * $ns;
           } else {
               $nsol = $nv * $nw * $nk * $nh * $ns;
           }

       }
   }
#  elsif ($opt eq "J") {
#    $iso1   = $value[0] if $value[0];
#    $iso2  = $value[1] if $value[1];
#    $clvd1  = $value[2] if $value[2];
#    $clvd2 = $value[3] if $value[3];
#    $fmt_flag="true";     # used later for renaming output figures with fmt
#  } 
   elsif ($opt eq "K") {
     $type = $value[0];
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
           ($v1, $v2) = @value;
           ($w1, $w2) = @value;
           $nv = $nw = 0;
           $fmt_flag="true";     # used later for renaming output figures with fmt
       } elsif ($#value==4) {
           ($v1, $w1, $k1, $h1, $s1) = @value;
           ($v2, $w2, $k2, $h2, $s2) = @value;
           $h1 = $h2 = cos($value[3]);  # cap expects h = cos(dip)
           $nsol = $nv = $nw = $nk = $nh = $ns = 1;
       } elsif ($#value==9) {
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
     $mltp = $value[0];
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
# prepare additional parameters for CAP input
#-----------------------------------------------------------

# convert strike/dip/rake to radians
$k1 = $k1 * $deg2rad;
$k2 = $k2 * $deg2rad;
$s1 = $s1 * $deg2rad;
$s2 = $s2 * $deg2rad;

unless ($dura) {
  $dura = int(10**(($mg-5)/2)+0.5);
  $dura = 1 if $dura < 1;
  $dura = 9 if $dura > 9;
}

# Additional settings for parameters when doing running regular grid.
# (dv, dw, dk, dh, ds) are grid SPACINGS for (gamma, delta, strike, dip, rake)
if( $type == 3 ) {
    use POSIX "fmod";   # may not be available outside unix/linux
    # full range for each parameter
    ($v1, $v2) = (-30, 30);
    ($w1, $w2) = (-90, 90);
    ($k1, $k2) = (  0, 360);
    ($h1, $h2) = (  0, 90);
    ($s1, $s2) = (-90, 90);

    # get total number of solutions
    # NOTE values nX are grid spacings, NOT number of points
    # NOTE values dX should be integers (CAP expects integers)
    ($dv, $dw, $dk, $dh, $ds) = ($nv, $nw, $nk, $nh, $ns);

    # get number of points per parameter
    if ($nv == 0 && $nw == 0) {
        print STDERR "double couple\n"; 
        $nv = 0;  # gamma
        $nw = 0;  # delta
    } else {
    # CHECK VERSION IN CAP AND +1, -1 VALUES
        $nv = (($v2 - $v1) / $dv) + 1;  # gamma -- include 0 point
        $nw = (($w2 - $w1) / $dw) + 1;  # delta -- include 0 point
    }
    # CHECK VERSION IN CAP AND +1, -1 VALUES
    $nk = (($k2 - $k1) / $dk) - 0;  # strike -- do not include 360
    $nh = (($h2 - $h1) / $dh);      # dip -- include 0 at start (though it will be offset later)
    $ns = (($s2 - $s1) / $ds) + 1;  # rake  -- include 0 point

    # check that spacings work
#   if (fmod($nv,2) || fmod($nw,2) || fmod($nk,2) || fmod($nh,2) || fmod($ns,2)) {
#       print STDERR "STOP. number of points not integer with current spacings and limits.\n";
#       print STDERR "$nv, $nw, $nk, $nh, $ns\n";
#       exit(0);
#   }
    $nsol = $nv * $nw * $nk * $nh * $ns;
}

$search = $type;

# CHECK THAT USER INPUT MAKE SENSE
# ONGOING
# I10000 (rand) goes with K2 else abort
# I10/10/10/10/10 goes with K1 else abort
#
# check parameterization type (zhu / lune) and search type (GRID/RANDOM) 
# NOTE parm is disabled so E is free
# if (($parm==0 ) && ($type==0)){
#     $search=0;
# } elsif (($parm==1 ) && ($type==1)){
#     $search=1;
# } elsif (($parm==1 ) && ($type==2)){
#     $search=2;
# } else {
#     printf STDERR "=====================================================
# This feature is not yet enabled. Try different combination of parameterization (-E) and search type (-K).
# =====================================================\n";
#     printf STDERR $usage;
#     exit(0);
# }


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

#----------------------------------------------------------
# delete 
#    exit(0);
#----------------------------------------------------------
 
for($dep=$dep_min;$dep<=$dep_max;$dep=$dep+$dep_inc) {
  foreach $eve (@event) {

    $md_dep = $model.'_'.$dep;
    next unless -d $eve;
    print STDERR "-------------------------------------------------------------\n";
    print STDERR "cap.pl: print some input parameters:\n";
    print STDERR "EVENT ID = $eve | EVENT DEPTH = $dep |  SOURCE DURATION = $dura\n";
    print STDERR "SEARCH TYPE  = $type | NSOL = $nsol\n";
    print STDERR "-------------------------------------------------------------\n\n";

    open(WEI, "$eve/$weight") || next;
    @wwf=<WEI>;
    close(WEI);
    $ncom = 2 if $wwf[0] =~ / -1 /;

    $cmd = "cap$dirct $eve $md_dep" unless $cmd eq "cat";    
    open(SRC, "| $cmd") || die "can not run $cmd\n";
    print SRC "$pVel $sVel $riseTime $dura $rupDir\n",$riseTime if $dirct eq "_dir";
    print SRC "$model $dep\n";          # first input in regular cap run
    print SRC "$m1 $m2 $max_shft1 $max_shft2 $repeat $bootstrap $fm_thr $tie\n";
    print SRC "@thrshd\n" if $repeat;   # no value in regular cap run
    print SRC "$vp $love $rayleigh\n";
    print SRC "$power_of_body $power_of_surf $weight_of_pnl $nof\n";
    print SRC "$plot\n";
    print SRC "$disp $mltp\n";
    print SRC "$green/$model/\n";
    print SRC "$search\n";
    print SRC "$norm\n";
    print SRC "$dt $dura $riseTime\n";
    print SRC "$f1_pnl $f2_pnl $f1_sw $f2_sw\n";
    #print SRC "$mg $dm $dlune\n";
    #print SRC "$iso1\n$iso2\n$clvd1\n$clvd2\n";
    #print SRC "$str1 $str2 $deg\n";
    #print SRC "$dip1 $dip2 $deg\n";
    #print SRC "$rak1 $rak2 $deg\n";
    #--- start parameters for search
    print SRC "$mw1 $mw2 $nmw $dm\n";
    print SRC "$v1 $v2 $nv $dv\n";
    print SRC "$w1 $w2 $nw $dw\n";
    print SRC "$k1 $k2 $nk $dk\n";
    print SRC "$h1 $h2 $nh $dh\n";
    print SRC "$s1 $s2 $ns $ds\n";
    print SRC "$nsol\n";
    #--- end parameters for search
    printf SRC "%d\n",$#wwf + 1;
    print SRC @wwf;
    close(SRC);

  plot:
    print STDERR "\n\ncap.pl: plot results ... \n";
    if ( $plot > 0 && ($? >> 8) == 0 ) {
      chdir($eve);
      @dum = split('_', $md_dep);  # split mdl string
      $outfile = sprintf("%s_%03d.out",@dum[0],int(@dum[1]));
      open(my $out,'>>',$outfile);
      say $out "INPUT_PAR $md_dep P_win $m1 S_win $m2 P $amplify p $ampfact NCOM $ncom spiB $spib spiS $spis $filterBand FMT $fmt_flag";

      #     &plot($md_dep, $m1, $m2, $amplify, $ncom, $sec_per_inch); # 20130102 calvizuri - original
      &plot($md_dep, $m1, $m2, $amplify, $ampfact, $ncom, $spib, $spis, $filterBand, $fmt_flag); # 20130102 calvizuri - added filter freq bands
      unlink(<${md_dep}_*.?>) unless $keep;
      chdir("../");
    }
  }
    print STDERR "cap.pl: plotting finished.\n";
}
exit(0);
