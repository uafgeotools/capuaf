# this subroutine plots waveform fits produced by source inversion srct
use List::Util qw[min max];

sub plot {

#  local($mdl, $t1, $t2 $am, $num_com, $spis) = @_; # original
  local($mdl, $t1, $t2, $am, $ampfact, $num_com, $spib, $spis, $filterBand, $fmt_flag, $evid, $model, $depth, $dura, $riseTime, $pol_wt) = @_;
  local($nn,$tt,$plt1,$plt2,$plt3,$plt4,$i,$nam,$com1,$com2,$j,$x,$y,@aa,$rslt,@name,@aztk);

# set this =1 if you want to plot time windows that have been excluded
  local $keepBad = 0;

# if you want to plot only polarities on the big beachball plot (No azimuth, station info or title)  - default 0
$only_pol = 0;
  
  @trace = ("1/255/255/255","3/0/0/0");       # plot data trace
  @name = ("P V","P R","Surf V"," Surf R","Surf T");

  $filterBand = "Filter periods (seconds): $filterBand";    # 20120719 - report filter bands
  $dura = sprintf("%.2f",$dura);
  $riseTime = sprintf("%.2f",$riseTime);
  $duration = "duration: $dura/$riseTime s";
  
#--------------------------
# plotting GMT file with synthetics and data
  
  # set GMT defaults
  # for all options (such as PAPER_MEDIA): http://gmt.soest.hawaii.edu/gmt/html/man/gmtdefaults.html
  # to check defaults, e.g.: gmtdefaults -L | grep MEASURE_UNIT
  @dum = split('_', $mdl);  # split mdl string
  $ftag=sprintf("%s_%s_%03d", $evid, $model, int($depth));
  $outfile = sprintf("%s.out", $ftag);

  # read in the output file results
#  open(FFF,"$mdl.out"); # original
  open(FFF,$outfile);    # 20130102 calvizuri - new file name
  @rslt = <FFF>;
  close(FFF);
  @meca = split('\s+',shift(@rslt));   # Line 1
  @hypo = split('\s+',shift(@rslt));   # Line 2
  @tensor = split('\s+',$rslt[0]);     # Line 3
  @others = grep(/^#/,@rslt);          # Line 4
  @ncomp = grep(/^#/,@rslt);           # Line 5
  @rslt=grep(!/^#/,@rslt);             # Remaing Lines
  $nrow = @rslt;

 # check if there are Input parameters in the last line
  @part = ();
  @last=split(' ',$rslt[$nrow-1]);
  if ($last[0] eq 'INPUT_PAR') {
      #$nrow = nrow-1;
      for $ii (0..$nrow-2) {
	   push @part, $rslt[$ii];}
       @rslt=@part;
       $nrow=$nrow-1;
   }

  # Page size
  $pheight_in = $nrow + 2;  # height of pape

  # positions of seismograms on the page
  # height of each seismogram
  $nn = int($pheight_in);
  $height = $pheight_in - 0.5;
  #($nn,$height) = (12,10.5);   # 10 rows of traces per 10.5 in.
  #print "\n$nn rows of traces per $height in";  
  print "\nseconds per inch = $spis";
  $sepb = 0.2*$spib;    # sec per inch (*1/2) bet body waves
  $seps = 0.2*$spis;    # separation bet surface waves

#  ($tt, $inc) = (2*$t1 + 3*$t2 + 4*$sepa, 1);
#  ($tt, $inc) = (2*$t1 + 3*$t2 + 2*$seps+2*$sepb, 1);
#  ($tt, $inc) = ($t1 + $t2 + $sepb, 4) if $num_com == 2;

  ($ttb, $inc) = (2*$t1 + 2*$sepb, 1);
  ($tts, $inc) = (3*$t2 + 2*$seps, 1);
  $tt = $ttb + $tts;

  ($tt, $inc) = ($t1 + $t2 + $sepb, 4) if $num_com == 2;
#  $width = 0.1*int(10*$tt/$sec_per_inch+0.5);
  $widthb = 0.1*int(10*$ttb/$spib+0.5);
  $widths = 0.1*int(10*$tts/$spis+0.5);
  $width = $widths+$widthb;
#  @x0 = ($t1+$sepa, $t1+$sepa, $t2+$sepa, $t2+$sepa, $t2);
  @x0 = ($t1+$sepb, $t1+$sepb, $t2+$seps, $t2+$seps, $t2);
  print "\n*** x0=@x0 *** \n";

  $pwidth_in = $width +1.5 ;  # width of paper    # orig 8.5
  print "\n$nrow rows to plot";
  print "\npaper is $pwidth_in inches wide and $pheight_in inches tall";
  system("gmtset BASEMAP_TYPE plain PAPER_MEDIA Custom_${pwidth_in}ix${pheight_in}i MEASURE_UNIT inch");

  # horizontal offset (why is it needed?) Because life is complicated bro!
  #  $xoffset="3.0";
  $xoffset=$widthb;

# pssac2 amplitude scaling:
# -M vertical scaling in sacfile_unit/MEASURE_UNIT = size<required> 
#           size: each trace will normalized to size (in MEASURE_UNIT)
#               scale =  (1/size) * [data(max) - data(min)]
#           size/alpha: plot absolute amplitude multiplied by (1/size)*r^alpha
#               where r is the distance range in km across surface
#               specifying alpha = 0.0 will give absolute amplitudes
#               scale = (1/size) * r^alpha
#           size/s: plot absolute amplitude multiplied by (1/size)*sqrt(sin(gcarc))
#               where gcarc is the distance in degrees.
#               scale = (1/size) * sqrt(sin(gcarc))

  # KEY: set amplitude scaling for seismograms
  if ($am>0.) {$stam = "$am/-1";} else {$stam=-$am;} # original line (with pssac, not pssac2)
  if ($am == 0x0){
      $amp = $am;}
  else{
      $amp = $am/$ampfact;}
  $stams = "$amp/0.";
  $stamb = "$am/0.";                                   # overwrite for absolute (to match default plotting)

  print "\namplitude scaling am = $am";
  print "\npssac2 amplitude scaling stam = $stam\n";
#  $outps = "$mdl.ps";   # original
  #$outps = sprintf("%s_%s_%03d.ps", $evid, $model, int($depth)); # reformatted filename
  #$outps = sprintf("%s_%s_%03d_fmt.ps", $evid, $model, ,int($depth)) if $fmt_flag eq "true";
  $outps = sprintf("${ftag}.ps");
  $outps = sprintf("${ftag}_fmt.ps") if $fmt_flag eq "true";

  # (1) plot cut seismograms with scaled amplitudes (first command: no -O appears)
  $tscale_x = 0.55;
  $tscale_y = $pheight_in - 2.0;
  $plt1b = "| pssac2 -JX${widthb}i/${height}i -L${spib} -l${tscale_x}/${tscale_y}/1/0.075/8 -R0/$ttb/0/$nn -Y0.2i -Ent-2 -M$stamb -K -P >> $outps";
  $plt1s = "| pssac2 -JX${widths}i/${height}i -L${spis} -l${tscale_x}/${tscale_y}/1/0.075/8 -R0/$tts/0/$nn -X${xoffset}i -Ent-2 -M$stams -O -K -P >> $outps";

  # (2) plot text labels
  $plt2_stn_info = "| pstext -JX -R -O -K -N -X-${xoffset}i >> $outps";
  $plt2_wf_info_b = "| pstext -JX${widthb}i/${height}i -R0/$ttb/0/$nn -O -K -N >> $outps";
  $plt2_wf_info_s = "| pstext -JX${widths}i/${height}i -R0/$tts/0/$nn -X${xoffset}i -O -K -N >> $outps";

  # (3) plot beachballs (solution, followed by possible local minima)
  $ballcolor = "150";
#  $dY = ${pheight_in} - 1.8;    # original
  $dY = ${pheight_in} - 1.6;
  $dX = -0.7-$xoffset;
#  $plt3 = "| psmeca -JX5i/1i -R-1/9/-1/1 -Sa5i -G$ballcolor -Y${dY}i -X-0.7i -O -K >> $outps";
  $plt3 = "| psmeca -JX5i/1i -R-1/9/-1/1 -Sa5i -G$ballcolor -Y${dY}i -X${dX}i -O -K >> $outps";
  $plt3 = "| psmeca -JX5i/1i -R-1/9/-1/1 -Sm8i -G$ballcolor -Y${dY}i -X${dX}i -O -K >> $outps" if $tensor[1] eq "tensor";

  # (4) plot markers on beachball
  # note: -JPa is a basemap for polar coordinates, clockwise from north

  # azimuths
  $plt4b = "| psxy -JPa1i -R0/360/0/1 -Sc0.02i -N -W0.5p,0/0/0 -G255 -O -K >> $outps";

  # supplemental: upper hemisphere piercing points on beachballs (o)
  #$plt4a = "| psxy -JPa1i -R0/360/0/1 -Sc0.08i -N -W0.5p,255/0/0 -O -K >> $outps";

  # default: lower hemisphere piercing points on beachballs (x) (last command: no -K appears)
  $plt4 = "| psxy -JPa1i -R0/360/0/1 -Sx0.10i -N -W0.5p,255/0/0 -G255 -O -K >> $outps";

#  $plt1=$plt2=$plt3="|cat";	# output GMT commands to command window for testing

  # (2.5) plot header information
  $dX = 0.8;
  $dY = 0.3;
  $plt4_5 = "| pstext -J -R -Y${dY}i -X${dX}i -O -N >> $outps";

#--------------------------

#  $outps2 = "${mdl}_beach.ps"; # original
  #$outps2 = sprintf("%s_%s_%03d_beach.ps", $evid, $model, int(depth));   # 20130102 calvizuri - revised filename
  #$outps2 = sprintf("%s_%s_%03d_beach_fmt.ps", $evid, $model, int($depth)) if $fmt_flag eq "true";
  $outps2 = sprintf("${ftag}_beach.ps");
  $outps2 = sprintf("${ftag}_beach_fmt.ps") if $fmt_flag eq "true";

  $fac = 6.5;
  $fac2 = 8.2*$fac;   # original: 5*$fac
  $JP = "-JPa${fac}i";

  # plot beachball
# $xplt3 = "| psmeca -JX${fac}i/${fac}i -R-1/1/-1/1 -N -G$ballcolor -W2p,0/0/0 -Sm${fac2}i -X1i -Y2i -K -P >> $outps2";
  $xplt3 = "| psmeca -JX${fac}i/${fac}i -R-1/1/-1/1 -N -G$ballcolor -W2p,0/0/0 -Sa${fac}i -X1i -Y2i -K -P >> $outps2";
  $xplt3 = "| psmeca -JX${fac}i/${fac}i -R-1/1/-1/1 -N -G$ballcolor -W2p,0/0/0 -Sm${fac2}i -X1i -Y2i -K -P >> $outps2" if $tensor[1] eq "tensor";

  # plot markers on beachball
  # note: -JPa is a basemap for polar coordinates, clockwise from north

  # azimuths
  $xplt4b = "| psxy $JP -R0/360/0/1 -Sc0.02i -N -W0.5p,0/0/0 -G255 -O -K >> $outps2";

  # supplemental: upper hemisphere piercing points on beachballs (o)
  $xplt4a = "| psxy $JP -R0/360/0/1 -Sc0.08i -N -W0.5p,255/0/0 -O -K >> $outps2";

  # default: lower hemisphere piercing points on beachballs (x)
  $xplt4 = "| psxy $JP -R0/360/0/1 -Sx0.10i -N -W0.5p,255/0/0 -G255 -O -K >> $outps2";
  $xplt4c = "| psxy $JP -R0/360/0/1 -St0.30i -N -W1p,0/255/0 -G255 -O -K >> $outps2";  # up polarity (green) - triangle
  $xplt4d = "| psxy $JP -R0/360/0/1 -Si0.30i -N -W1p,0/0/255 -G255 -O -K >> $outps2";  # down polarity (blue) - triangle
  $xplt4e = "| psxy $JP -R0/360/0/1 -St0.30i -N -W1p,255/0/0 -G255 -O -K >> $outps2";  # non-matching polarity red) - triangle
  $xplt4f = "| psxy $JP -R0/360/0/1 -Si0.30i -N -W1p,255/0/0 -G255 -O -K >> $outps2";  # non-matching polarity (red) - triangle

  # plot text labels
  $xplt5a = "| pstext $JP -R0/360/0/1 -N -O -K >> $outps2";
  $xplt5b = "| pstext $JP -R0/360/0/1 -N -O -K >> $outps2";
  $xplt5c = "| pstext $JP -R0/360/0/1 -N -O -K >> $outps2";

  # title (LAST COMMAND: no -K appears)
  $xplt6 = "| pstext -JX -R -N -O -Xa0 -Ya7.5 >> $outps2";

#--------------------------
# FIGURE 1: waveform fits with moment tensor
  print STDERR "cap_plt.pl: plotting summary for this solution ... \n";

  # uncomment output for debugging purposes
  print "\n-------------------";
  print "\nmeca:\n@meca";
  print "\ntensor:\n@tensor";
  print "\nothers:\n@others";
  print "\nrslt:\n@rslt";
  print "\n-------------------\n";

# get strike dip and rake
  $stk = @meca[0];
  $dip = @meca[1];
  $rak = @meca[2];

  # compute piercing points for beachballs
  $P_val=0; # maximum aplitude for pssac plotting (-P flag) - Body
  $S_val=0; # maximum aplitude for pssac plotting (-P flag) - Surface
  $i = 0; $j = 0; $i2 = 0; $j2 = 0;
  $pi = 3.14159265358979323846;
  @tklh=(); @tkuh=(); @staz=(); @az=(); @tklh_useweights=(); @staz_useweights=(); @tkuh_useweights=();
  foreach (@rslt) {
    @aa = split;
    if ($aa[7]>$P_val && $aa[2]!=0){$P_val=$aa[7];}   # maximum amplitude for pssac plotting (-P flag) [Maximum amplitude of vertical body wave]
    if ($aa[14]>$P_val && $aa[9]!=0){$P_val=$aa[14];} # Maximum amplitude of radial body wave
    if ($aa[21]>$S_val && $aa[16]!=0){$S_val=$aa[21];} # Maximum amplitude of vertical surface wave
    if ($aa[28]>$S_val && $aa[23]!=0){$S_val=$aa[28];} # Maximum amplitude of radial surface wave
    if ($aa[35]>$S_val && $aa[30]!=0){$S_val=$aa[35];} # maximum amplitude of love wave 
    $ifmp[$i] = $aa[37];  # first-motion polarity (input - data)
    $ifmpt[$i] = $aa[38];  # first-motion polarity (theoretical)
    $stnm = $aa[0];                              # station name
    #next if $aa[2] == 0;                        # skip if no body waves
    $x = `saclst az user1 f ${mdl}_$aa[0].0`;    # get the azimuth and P take-off angle
    @dd = @aa;
    @aa = split(' ', $x);       # outputs something like this: wes_1_HOYA.LL.TPH..LH.0 323.513 90.72
    @aa_pre = @aa;
    #print "\n--> saclst az user1 f ${mdl}_$aa[0].0";

    # compute polar coordinates azimuth and radius
    # NOTE this part outputs all azimuths in the weight file, even if the
    # station was not used in the inversion. Unless the input weigh files are
    # pre-sorted and clean.
    # WARNING if the weight files are not clean this line may cause a mismatch
    # between STNAME and AZIM 
    # CHECK
    $az[$i] = $aa[1];

    $azvec[$i] = sprintf("%s\n",$aa[1]);
    $staz[$i] = sprintf("%s %f %s\n",$aa[1],1.1,$stnm);      # station azimuth
    if ($aa[2]>90.) {                                         # upper hemisphere
       $rad = sqrt(2.)*cos($aa[2]*$pi/360);
       $tkuh[$j] = sprintf("%s %f %s\n",$aa[1],$rad,$stnm);
       $j++;
       # project piercing point to lower hemisphere
       $aa[1] += 180;
       $aa[2]=180-$aa[2];
    }
    $rad = sqrt(2.)*sin($aa[2]*$pi/360);
    $tklh[$i] = sprintf("%s %f %s\n",$aa[1],$rad,$stnm);        # lower hemisphere
    if (($dd[37] != 0  && $pol_wt != 0) || $dd[2]!=0 || $dd[9]!=0 || $dd[16]!=0 || $dd[23]!=0 || $dd[30]!=0 || $keepBad!=0){
	$tklh_useweights[$i2] = sprintf("%s %f %s\n",$aa[1],$rad,$stnm);
	$staz_useweights[$i2] = sprintf("%s %f %s\n",$aa_pre[1],1.1,$stnm);
	if ($aa_pre[2]>90.) {
	    $tkuh_useweights[$j2] = sprintf("%s %f %s\n",$aa_pre[1],$rad,$stnm);
	    $j2++;
	}
	$i2++;
    }
    $i++;
  }
#--------------------------compute pssac plotting info (scaling factor P_val)
  print "Maximum body wave amplitude = $P_val \n";
  print "Maximum surface wave amplitude = $S_val \n";
##---------------------------------------
# Three options for plotting (and scaling) the waveforms using -P flag (body waves) and -p flag (surface wave)
# Default for both -P and -p flag is 1 (i.e. option 1 in the following comments and using the scaling_factor=1)
# 1. Normalize by maximum body and surface amplitude separatetly, then apply a scaling factor
# 2. Normalized plotting -- data and synthetics have same maximum amplitude for all waveforms
# 3. The default plotting -- scale waveforms by given amplitude

  # -P flag for the body waves
  if ($am>=0.1){             # scale by the maximum body wave amplitude ($P_val) and then scale by $am factor (-P flag)
      $am1 = $P_val/$am;
    }
  elsif ($am==0){           # Normalized plotting (using the pssac2 bug) -- data and synthetics have same maximum amplitude
      $am1 = "0.5e+0.5";
    }
  else {                     # default plotting (FUTURE: find a better way to differentiate b/w exponents and rational number) -- scale by given amplitude $am  (-P flag)
      $am1 = $am;
  }
  # -p flag for the surface waves
  if ($ampfact>=0.1){        # scale by the maximum surface wave amplitude ($S_val) and then scale by $ampfact factor (-p flag) 
      $am2 = $S_val/$ampfact;
  }
  elsif ($ampfact==0){      # Normalized plotting (using the pssac2 bug) -- data and synthetics have same maximum amplitude
      $am2 = "0.5e+0.5";
  }
  else {                     # default plotting (FUTURE: find a better way to differentiate b/w exponents and rational number) -- scale by given amplitude $ampfact  (-P flag)
      $am2 = $ampfact;
  }
  $stamb = "$am1/0.";        # set the parameters
  $stams = "$am2/0.";

  print "pssac2 input for normalizing body waves = $am1 \n";
  print "pssac2 input for normalizing surface waves = $am2 \n";
#---------------------------------------
  # 20151025 cralvizuri - uncomment this command to normalize surf waves
  #                       This is for figures in Uturuncu FMT paper
  #$stams = $stamb;

  $plt1b = "| pssac2 -JX${widthb}i/${height}i -L${spib} -l${tscale_x}/${tscale_y}/1/0.075/8 -R0/$ttb/0/$nn -Y0.2i -Ent-2 -M$stamb -K -P >> $outps";
  $plt1s = "| pssac2 -JX${widths}i/${height}i -L${spis} -l${tscale_x}/${tscale_y}/1/0.075/8 -R0/$tts/0/$nn -X${xoffset}i -Ent-2 -M$stams -O -K -P >> $outps";

  # remove the file if it exists
  unlink($outps) if -e $outps;
  unlink($outps2) if -e $outps2;

  # save a copy for the second file
  @rslt0=@rslt;

  while (@rslt) {

#    # plot waveforms
#    open(PLT, $plt1);
#    $i = 0;
#    @aaaa = splice(@rslt,0,$nn-2);
#    foreach (@aaaa) {
#      @aa = split;
#      $nam = "${mdl}_$aa[0].";
#      $x=0;
#      for($j=0;$j<5;$j+=$inc) {
#        $com1=8-2*$j; $com2=$com1+1;
#	if ($aa[4*$j+2]>0) {
#	   printf PLT "%s %f %f 5/0/0/0\n",$nam.$com1,$x,$nn-$i-2;
#           printf PLT "%s %f %f 3/255/0/0\n",$nam.$com2,$x,$nn-$i-2;
#	} elsif ($keepBad) {
#	   printf PLT "%s %f %f 2/0/255/0\n",$nam.$com1,$x,$nn-$i-2;
#           printf PLT "%s %f %f 3/255/0/0\n",$nam.$com2,$x,$nn-$i-2;
#	}
#        $x = $x + $x0[$j];
#      }
#      $i++;
#    }
#    close(PLT);

    # plot waveforms body waves
    open(PLT, $plt1b);
    $i = 0;
    @aaaa = splice(@rslt,0,$nn-2);
    foreach (@aaaa) {   # go over each line in .out file
        @aa = split;
	if (($aa[37]!=0 && $pol_wt != 0) || ($aa[2]!=0 || $aa[9]!=0 || $aa[16]!=0 || $aa[23]!=0 || $aa[30]!=0 || $keepBad!=0)){
        $nam = "${mdl}_$aa[0].";
        $x=0;
        for($j=0;$j<2;$j+=$inc) {
            $com1=8-2*$j; $com2=$com1+1;    # seismogram extensions (.0, .1, .2...)
            if ($aa[7*$j+2]>0) {
#                printf "(j=$j) x=$x\t"; # debug
                printf PLT "%s %f %f 5/0/0/0\n",  $nam.$com1,$x+0,$nn-$i-2;     # data (black)
                printf PLT "%s %f %f 3/255/0/0\n",$nam.$com2,$x+0,$nn-$i-2;   # synthetic (red)
            } elsif ($keepBad) {
                printf PLT "%s %f %f 2/0/255/0\n",$nam.$com1,$x+0,$nn-$i-2;   # bad data (green)
                printf PLT "%s %f %f 3/255/0/0\n",$nam.$com2,$x+0,$nn-$i-2;   # synthetic (red)
            }
            $x = $x + $x0[$j];
        }
#                printf "\n"; # debug
        $i++;
    }
    }
    close(PLT);

    # plot waveforms surface waves
    open(PLT, $plt1s);
    $i = 0;
    foreach (@aaaa) {
        @aa = split;
	if (($aa[37]!=0 && $pol_wt != 0) || ($aa[2]!=0 || $aa[9]!=0 || $aa[16]!=0 || $aa[23]!=0 || $aa[30]!=0 || $keepBad!=0)){
        $nam = "${mdl}_$aa[0].";
#        $x=$x0[1];
        $x=0;
        for($j=2;$j<5;$j+=$inc) {
#                printf "(j=$j) x=$x\t"; # debug
            $com1=8-2*$j; $com2=$com1+1;
            if ($aa[7*$j+2]>0) {
                printf PLT "%s %f %f 5/0/0/0\n",  $nam.$com1,$x,$nn-$i-2;   # data (black)
                printf PLT "%s %f %f 3/255/0/0\n",$nam.$com2,$x,$nn-$i-2;   # synthetic (red)
            } elsif ($keepBad) {
                printf PLT "%s %f %f 2/0/255/0\n",$nam.$com1,$x,$nn-$i-2;   # bad data (green)
                printf PLT "%s %f %f 3/255/0/0\n",$nam.$com2,$x,$nn-$i-2;   # synthetic (red)
            }
            $x = $x + $x0[$j];
        }
#                printf "\n"; # debug
        $i++;
    }
    }
    close(PLT);

    
    # text labels
#    open(PLT, $plt2);
#    $y = $nn-2;
#    $i=0;
#    foreach (@aaaa) {
#      @aa = split;
#      $x = 0;
#      printf PLT "%f %f 10 0 0 1 $aa[0]\n",$x-0.8*$spis,$y;            # station label
#      printf PLT "%f %f 10 0 0 1 $aa[1]\n",$x-0.7*$spis,$y-0.2;        # distance_km/overal time shift
#      printf PLT "%f %f 10 0 0 1 %.1f\n",$x-0.7*$spis,$y-0.4,$az[$i];  # azimuth (see az above)
#      $i=$i+1;
#      for($j=0;$j<5;$j+=$inc) {
#	if ($aa[4*$j+2]>0 || $keepBad) {
#
#          printf "(j=$j) x=$x \t ";
#          printf PLT "%f %f 10 0 0 1 $aa[4*$j+5]\n",$x,$y-0.4;  # time shift each wave
#          printf PLT "%f %f 10 0 0 1 $aa[4*$j+4]\n",$x,$y-0.6;  # correl value
#	}
#        $x = $x + $x0[$j];
#      }
#      $y--;
#    }

    # plot station info
    
    open(PLT, $plt2_stn_info);
    $y = $nn-2;
    $i=0;
    foreach (@aaaa) {
# Ruler for reading CAP output
#    0          1       2     3   4   5      6      7       8     9    10  11   12    13      14       15   16   17  18   19    20    21        22    23   24  25   26    27    28     29 
#    |          |       |     |   |   |      |      |       |     |     |   |   |     |        |       |    |     |   |    |    |      |        |     |     |   |    |     |    |       |
# PLMK_XP    11.2/0.14  1   0.67 95 -0.08  0.64 8.19e-07 4.32e-07 1   0.79 80 -0.08 -0.09 8.05e-07 8.79e-07 1   3.48 79  1.89  0.84 9.47e-07 4.10e-07 1   4.47 75  1.89  1.24 1.03e-06 2.98e-07 1   3.68 80  0.23  1.61 7.56e-07 1.51e-07  1   0.45
#                                                                                                                                                                                               |     |   |   |      |    |         |      |    |
#                                                                                                                                                                                               30    31 32   33    34    35        36     37   38

        # variables from cap output
        @aa = split;
        @ab = split('/',$aa[1]);
        $dist_km = $ab[0];
        $tshift_all = $ab[1];
        $pol_syn = $aa[37];
        $pol_obs = $aa[38];

        # test if weight or polarity exists. if neither then print nothing and dont skip space
        if (($aa[37]!=0  && $pol_wt != 0) || ($aa[2]!=0 || $aa[9]!=0 || $aa[16]!=0 || $aa[23]!=0 || $aa[30]!=0 || $keepBad!=0)){

            $x = 0;
            @sensor_label = split('\.', $aa[0]);
            $inet = $sensor_label[1];
            $ista = $sensor_label[2];
            $iloc = $sensor_label[3];
            $icha = $sensor_label[4];
            # station label
            printf PLT "%f %f 10 0 0 1 $inet.$ista.$iloc.$icha\n", $x-0.8*$spis, $y;
            #printf PLT "%f %f 10 0 0 1 $sensor_label[1].$sensor_label[2].$sensor_label[3]\n", $x-0.8*$spis, $y;
            # station distance and overall time-shift
            if ($tshift_all==0.){
                printf PLT "%f %f 10 0 0 1 %d km\n", $x-0.8*$spis, $y-0.2, $dist_km;
                printf PLT "%f %f 10 0 0 1 %d\260 \n", $x-0.8*$spis, $y-0.4, $az[$i];  # azimuth (see az above)
            }
            else { 
                printf PLT "%f %f 10 0 0 1 %d km\n", $x-0.8*$spis, $y-0.2, $dist_km;
                printf PLT "%f %f 10 0 0 1 %d\260 \n", $x-0.8*$spis, $y-0.4, $az[$i];  # azimuth (see az above)
                printf PLT "%f %f 10 0 0 1 %.2f s\n", $x-0.8*$spis, $y-0.6, $tshift_all;  # tshift = Green_P_arrival - Input_P_arrival_weight_file
            }
            # azimuth
            # printf PLT "%f %f 10 0 0 1 %d\260 \n", $x-0.8*$spis, $y-0.6, $az[$i];  # azimuth (see az above)
            # polarities
            # NOTE if polarity is 0 or does not exist, then nothing is written
            if ($pol_syn || $keepBad==1) {
		if ($ab[1]==0.) {
		    printf PLT "%f %f 10 0 0 1 $pol_syn ($pol_obs)\n", $x-0.8*$spis, $y-0.6;
		}
		else {
		    printf PLT "%f %f 10 0 0 1 $pol_syn ($pol_obs)\n", $x-0.8*$spis, $y-0.8;
		}
            }
            $i=$i+1;
            $y--;
        } # end tests for weight and polarity
    }
    close(PLT);

    # plot data labels body waves
    open(PLT, $plt2_wf_info_b);
    $y = $nn-2;
    foreach (@aaaa) {
      @aa = split;
      if (($aa[37]!=0 && $pol_wt != 0) || ($aa[2]!=0 || $aa[9]!=0 || $aa[16]!=0 || $aa[23]!=0 || $aa[30]!=0 || $keepBad!=0)){
      $x = 0;
      for($j=0;$j<2;$j+=$inc) {
          if ($aa[7*$j+2]>0 || $keepBad) {
            # printf PLT "%f %f 10 0 0 1 $aa[4*$j+5]\n",$x,$y-0.4;  # time shift each wf
            # printf PLT "%f %f 10 0 0 1 $aa[4*$j+4]\n",$x,$y-0.6;  # correl value
              $fracmis=sprintf("%2.2f", $aa[7*$j+3]);
              $lamp=sprintf("%2.2f", $aa[7*$j+6]);
              printf PLT "%f %f 10 0 0 1 $aa[7*$j+5]\n", $x+0, $y-0.2;  # time shift each wf
              printf PLT "%f %f 10 0 0 1 $aa[7*$j+4]\n", $x+0, $y-0.4;  # correl value
              printf PLT "%f %f 10 0 0 1 $fracmis\n", $x+0, $y-0.6;     # fractional misfit
              printf PLT "%f %f 10 0 0 1 $lamp\n", $x+0, $y-0.8;        # log(max_amp_data/max_amp_syn)
          }
          $x = $x + $x0[$j];
      }
      $y--;
    }
    # plot labels PR and PV
    $x = 0.2*$spib;
    for($j=0;$j<2;$j+=$inc) {
      printf PLT "%f %f 12 0 0 1 $name[$j]\n",$x,$nn-1.5;
      $x = $x+$x0[$j];
    }
  }
    close(PLT);

    # plot data labels surface waves
    open(PLT, $plt2_wf_info_s);
    $y = $nn-2;
    foreach (@aaaa) {
      @aa = split;
      if (($aa[37]!=0 && $pol_wt != 0) || ($aa[2]!=0 || $aa[9]!=0 || $aa[16]!=0 || $aa[23]!=0 || $aa[30]!=0 || $keepBad!=0)){
#      $x = $x0[1];
      $x = 0;
      for($j=2;$j<5;$j+=$inc) {
          if ($aa[7*$j+2]>0 || $keepBad) {
              #printf PLT "%f %f 10 0 0 1 $aa[4*$j+5]\n",$x,$y-0.4;  # time shift each wave
              #printf PLT "%f %f 10 0 0 1 $aa[4*$j+4]\n",$x,$y-0.6;  # correl value
              $fracmis=sprintf("%2.2f", $aa[7*$j+3]);
              $lamp=sprintf("%2.2f", $aa[7*$j+6]);
              printf PLT "%f %f 10 0 0 1 $aa[7*$j+5]\n", $x, $y-0.2;  # time shift each wave
              printf PLT "%f %f 10 0 0 1 $aa[7*$j+4]\n", $x, $y-0.4;  # correl value
              printf PLT "%f %f 10 0 0 1 $fracmis\n", $x, $y-0.6;  # fractional misfit
              printf PLT "%f %f 10 0 0 1 $lamp\n", $x, $y-0.8;  # log(max_amp_data/max_amp_syn)
          }
          $x = $x + $x0[$j];     # original
      }
      $y--;
    }
    # -------------------- end plot data for each trace

    # plot labels PR PV SV SR SH for wave types

    # plot labels SV SR SH
    $x = 0.2*$spis;
#    $x = 0.2*$spis+$x0[2];
    for($j=2;$j<5;$j+=$inc) {
      printf PLT "%f %f 12 0 0 1 $name[$j]\n",$x,$nn-1.5;
      $x = $x+$x0[$j];
    }
  }
    close(PLT);

    
    # plot beachball
    # note: magnitude scale is "fixed" at 1e17 for psmeca -Sm and 1 for psmeca -Sa
    open(PLT, $plt3);
    if ($tensor[1] eq "tensor") {
        # moment tensor is converted from AkiRichads basis to GCMT basis, which is required for psmeca
        printf PLT "0 0 0 @tensor[9,4,7,6] %f %f 17\n",-$tensor[8],-$tensor[5];
    } else {
        # focal mechanism is plotted from the M0, strike/dip/rake values
        printf PLT "0 0 0 @meca[5,6,7] 1\n";  # 0.5*$spis,$nn-1;
    }
#    $x = 2;
#    foreach (@others) {
#       split;
#       printf PLT "%f -0.2 0 @_[1,2,3] 0.5 0 0 $_[6]\n",$x; $x+=1.5;
#    }
    close(PLT);

    # plot station azimuths beachballs (see staz above)
    #open(PLT, $plt4b);
    #foreach (@staz) {
    #  printf PLT;
    #}

    open(PLT, $plt4b);
    foreach (@staz_useweights) {
      printf PLT;
    }

    # Does this do anything??
    # plot station azimuths beachballs (see tkuh above)
    open(PLT, $plt4a);
    foreach (@tkuh) {
      printf PLT;
    }

    # plot piercing points on beachballs (see tklh above)
    #open(PLT, $plt4);
    #foreach (@tklh) {
    #  printf PLT;
    #}
    #close(PLT);

    open(PLT, $plt4);
    foreach (@tklh_useweights) {
      printf PLT;
    }
    close(PLT);

    # test start
#    $x = 0.5*$spis; 
    $x = 0; 
    $y = 0;
    $tgap=0.5;
    # plot four header labels (event type, focal mecha, var red, filters)
    # Event 19910914190000000 Model 19910914190000000_wes_001 FM  350 56.985645  -74 Mw 5.80 rms 2.673e-06     1 CLVD -4.08 ISO  -4.464618 VR 7.8 data2 2.783e-06
    #   0         1             2       3                     4    5      6        7  8   9   10      11      12  13    14  15       16    17  18   19     20
    open(PLT, $plt4_5);
    printf PLT "$x $y 12 0 0 0 Event $evid Model $model Depth $depth\n"; $y-=$tgap;
    printf PLT "$x $y 12 0 0 0 @meca[4] %d %d %d @meca[8,9] @~g@~ %3.0f @~d@~ %3.0f @meca[10,11] VR %3.1f pol_wt %0.2f \n",@meca[5], @meca[6], @meca[7], @meca[14],@meca[16],@meca[18],$pol_wt;$y-=$tgap;
    printf PLT "$x $y 12 0 0 0 $filterBand $duration\n" ; $y-=$tgap;  # 20120719 - filter bands
    printf PLT "$x $y 12 0 0 0 @ncomp[1]" ;
    close(PLT);

  print STDERR "cap_plt.pl: done. \n";
  }  # while (@rslt) {

#---------------------------------
# FIGURE 2: big moment tensor with station names at lower-hemisphere piercing points

  print STDERR "cap_plt.pl: plotting big beach ball ... \n";
  $pwidth_in = 8.5;  # width of paper
  $pheight_in = 11;  # height of paper
  system("gmtset BASEMAP_TYPE plain PAPER_MEDIA Custom_${pwidth_in}ix${pheight_in}i MEASURE_UNIT inch");

  # restore
  @rslt=@rslt0;

  while (@rslt) {
    @aaaa = splice(@rslt,0,$nn-2);

    # plot beachball (see notes above)
    open(XPLT, $xplt3);
    if ($tensor[1] eq "tensor") {
     printf XPLT "0 0 0 @tensor[9,4,7,6] %f %f 17\n",-$tensor[8],-$tensor[5];
    } else {
     printf XPLT "0 0 0 @meca[5,6,7] 1\n"; #0.5*$spis,$nn-1;
    }
    close(XPLT);

    # plot piercing points on beachballs (see tklh above)
    $i=0; $j=0; $k=0;
    open(XPLT, $xplt4);
    open(XPLTC, $xplt4c);
    open(XPLTD, $xplt4d);
    open(XPLTE, $xplt4e);
    open(XPLTF, $xplt4f);
    foreach (@tklh_useweights) {
	if ($ifmp[$i] * $ifmpt[$i] < 0) {     # mismatcing polarities
	    if ($ifmp[$i]>0){printf XPLTE;}   # input is UP (+1); theoretical is DOWN (-1)
	    else {printf XPLTF;}}             # input is DOWN (-1); theoretical is UP (+1)
	elsif ($ifmp[$i]>0){printf XPLTC;}    # both input and theoretical are UP (+1)
	elsif ($ifmp[$i]<0){printf XPLTD;}    # both input and theoretical are DOWN (-1)
        else {printf XPLT;}                   # no input polarity pick in the weight file
	$i=$i+1;
    }
    close(XPLT);
    close(XPLTC);
    close(XPLTD);
    close(XPLTE);
    close(XPLTF);

# Section for plotting azimuths and station name
if ($only_pol == 0) {
    # plot station azimuths beachballs (see staz above)
    open(XPLT, $xplt4b);
    foreach (@staz_useweights) {
      printf XPLT;
    }
    close(XPLT);

    # plot station azimuths beachballs (see tkuh above)
    open(XPLT, $xplt4a);
    foreach (@tkuh_useweights) {
      printf XPLT;
    }
    close(XPLT);
#------------

    open(XPLT, $xplt5a);
    foreach (@staz_useweights) {
      @aa = split;
      @aa_split = split('\.', $aa[2]);
      printf XPLT "%s %s 8 0 0 CB %s.%s.%s\n", 
          $aa[0], $aa[1], $aa_split[1], $aa_split[2], $aa_split[3];
    }
    close(XPLT);

#     open(XPLT, $xplt5b);
#     foreach (@tkuh) {
#       @aa = split;
#       printf XPLT "%s %s 8 0 0 CB (%s)\n",$aa[0],$aa[1],$aa[2]; 
#     }
#     close(XPLT);

    open(XPLT, $xplt5c);
    foreach (@tklh_useweights) {
      @aa = split;
      @aa_split = split('\.', $aa[2]);
      printf XPLT "%s %s 8 0 0 CB %s.%s.%s\n", 
          $aa[0], $aa[1], $aa_split[1], $aa_split[2], $aa_split[3];
    }
    close(XPLT);

    # TITLE
    $x = -1; 
    $y = 0;
    open(XPLT, $xplt6);
    printf XPLT "0 0 16 0 0 0 @meca[0..3]\n";
    # Event 19910914190000000 Model 19910914190000000_wes_001 FM  350 56.985645  -74 Mw 5.80 rms 2.673e-06     1 CLVD -4.08 ISO  -4.464618 VR 7.8 data2 2.783e-06
    #   0         1             2       3                     4    5      6        7  8   9   10      11      12  13    14  15       16    17  18   19     20
    printf XPLT "0 -0.05 16 0 0 0 @meca[4] %d %d %d @meca[8,9] @~g@~ %3.0f @~d@~ %3.0f @meca[10,11] VR %3.1f pol_wt %0.2f\n",@meca[5], @meca[6], @meca[7], @meca[14],@meca[16],@meca[18], $pol_wt;
    close(XPLT);

}
  print STDERR "cap_plt.pl: done.\n";
}
#---------------------------------

}
1;
