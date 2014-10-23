#!/usr/bin/env perl
#
# find best depth from CAP output and
# plot the misfit error and source mechanism
# as function of depth for each event
# Usage depth.pl result_file event_dir_names ...
#

$ballsize = 4;		    # controls default beachball size (psmeca)
$min0 = 1.0e+19;		# impossibly large misfit value
$Bscale = "-Ba5f1:\"\":/a20f5:\"\":";
$onlydc = 0;

#------------

#($rsl,@event) = @ARGV;     # original
($rsl,$rsl2,@event) = @ARGV; # input Event and Tensor lines from files cus_NNN.out, and event id

open(RSL,"$rsl") or die "couldn't open inputfile1 $rsl\n";
@aaa=<RSL>;
#print STDERR "debug: aaa=@aaa\n";  # 20130103 calvizuri - ok to delete
close(RSL);

#-----------
# GMT plot settings
$xoffset = "-X3.5c";
$yoffset = "-Y3.5c";

# default misfit
#$xtick1 = 50; $xtick2 = 5;
$xtick1 = 5; $xtick2 = 1;
$ytick1 = 20; $ytick2 = 5;
$B1 = "-Ba${xtick1}f${xtick2}:\"Depth, km\":/a${ytick1}f${ytick2}:\"Misfit relative to minimum\":";

# log (misfit/minimum misfit)
$xtick1 = 5; $xtick2 = 1;
$ytick1 = .01; $ytick2 = .005;
$B2 = "-Ba${xtick1}f${xtick2}:\" \":/a${ytick1}f${ytick2}:\"ln(ERR / ERR_min)\":nW";

# VR
$xtick1 = 5; $xtick2 = 1;
$ytick1 = 20; $ytick2 = 5;
$B3 = "-Ba${xtick1}f${xtick2}:\"Depth, km\":/a${ytick1}f${ytick2}:\"VR\":ES";

# size of plot
$zoom=3;
$magpsx=2.0*$zoom;
$magpsy=1.8*$zoom;
$J = "-JX$magpsx/$magpsy";

$pwidth_in = 8.5;		# width of paper
$pheight_in = 11;		# height of paper
$fsize = 16;
system("gmtset BASEMAP_TYPE plain PAPER_MEDIA Custom_${pwidth_in}ix${pheight_in}i MEASURE_UNIT inch LABEL_FONT_SIZE $fsize ANOT_FONT_SIZE $fsize");

#-----------
printf STDERR "============Starting depth search==========\n";
open(FH,"$rsl2") or die "couldn't open inputfile2 $rsl2\n";
@data_fmt=<FH>;
#print STDERR "debug: data_fmt=@data_fmt\n";
close(FH);

while (@event) {
  @aa = splice(@event,0,10); # 0=offset, 10=N elements (each element space separated)
  #print STDERR "debug. aa=@aa event=@event\n";
  $i=0;
  #$xx = "-K -Ba50f5/a20f5WSne";   # original <-- not enough ticks+labeling
  #$xx = "-K ${Bscale}WS -P";
  $xx = "-K -P $xoffset $yoffset";
  foreach $eve (@aa) {		# eve = current aa = event id
    $ii=1;
    $best=1;
    $min=$min0;
    foreach (grep(/$eve/,@aaa)) {
      chop;   
      # 20130103 calvizuri - example input to chop:
      # Event 20080418093700 Model cus_001 FM 291 51 -13 Mw 5.10 rms 3.748e-02   110 ERR   2   5   8 ISO 0.16 0.11 CLVD 0.14 0.08
      $line[$ii]=$_;
      @bb=split;
      ($aa,$dep[$ii])=split('_',$bb[3]); # split cus_001 to get depth (= 001km) -- $dep is used for y-axis range
          $strike[$ii]=$bb[5];    # not needed
          $dip[$ii]=$bb[6];       # not needed
          $rake[$ii]=$bb[7];      # not needed
      $mw[$ii]=$bb[9];
      #    printf STDERR "debug. mw[$ii]=%lf\n",$mw[$ii];
      $rms[$ii]=$bb[11];
      $vr[$ii]=$bb[24];
      if ($min>$rms[$ii]) {
	$best=$ii;$min=$rms[$ii];
	$depth = $dep[$ii];
      }
      $ii++;
  }
      # Find the ymax range for plotting log(misfit)
    #printf STDERR "%f %f \n",$rms[$ii-1],$min;
    $max = log($rms[$ii-1]/$min);
    if (log($rms[1]/$min) > log($rms[$ii-1]/$min)){
	$max = log($rms[1]/$min);}
    $max = sprintf("%1.2f",$max);   # suppress to 3 decimal places
    #printf STDERR "%f\n",$max;

    $jj=1;
    foreach (grep(/tensor/,@data_fmt)) {
      # We will go for consistency with cap_plt.pl, which has this line:
      # printf PLT "0 0 0 @tensor[9,4,7,6] %f %f 17\n",-$tensor[8],-$tensor[5];
      chop; # example input to chop: # tensor = 5.696e+23  0.838 -0.564 -0.335 -0.259  0.409 -0.185
      @kk=split;
      @dummy = split('\+',$kk[3]); # get only the exponent--needed for scaling beachballs on plot
      #$M0_floor="1E+$dummy[1]";
      #printf STDERR "debug. M0_floor=%e\n",2*$M0_floor;
      #$M0_exp[$jj] = $dummy[1];
      $M0_exp[$jj] = 17;	# see cap_plt.pl
      $mrr[$jj] = $kk[9];	# mrr=m33 -- eg box 8.3, p.351 (L&W95)
      $mtt[$jj] = $kk[4];	# m11
      $mff[$jj] = $kk[7];	# m22
      $mrt[$jj] = $kk[6];	# m13
      $mrf[$jj] =-$kk[8];	# -m23
      $mtf[$jj] =-$kk[5];	# -m12
          printf STDERR "debug. MT components: %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f M0_exp=%f\n",$mrr[$jj], $mtt[$jj], $mff[$jj], $mrt[$jj], $mrf[$jj], $mtf[$jj],$M0_exp[$jj];
      $jj++;
    }

    # compute a relative measure of error ($lerr = log(misfit/misfit_min))
    for ($jj=1;$jj<$ii;$jj=$jj+1) {
	$lerr[$jj] = $rms[$jj]/$rms[$best];
	$lerr[$jj] = log($lerr[$jj]);
	printf STDERR "%f %e %e %d %d %d\n",$lerr[$jj],$rms[$jj],$rms[$best],$dep[$best], $best,$ii;
    }

    # compute parabolic equation:  Y= aX^2 + bX + c (cleaner way to handle)
    $x1 = $dep[$best-1]; $y1 = $lerr[$best-1];   #(x1,y1)
    $x2 = $dep[$best]; $y2 = $lerr[$best];       #(x2,y2)
    $x3 = $dep[$best+1]; $y3 = $lerr[$best+1];   #(x3,y3)
    $denom = ($x1-$x2)*($x1-$x3)*($x2-$x3);
    $a = (($x3*($y2-$y1))+($x2 * ($y1 - $y3)) + ($x1 * ($y3 - $y2))) / $denom;
    $b = ($x3*$x3 * ($y1 - $y2) + $x2*$x2 * ($y3 - $y1) + $x1*$x1 * ($y2 - $y3)) / $denom;
    $c = ($x2 * $x3 * ($x2 - $x3) * $y1 + $x3 * $x1 * ($x3 - $x1) * $y2 + $x1 * $x2 * ($x1 - $x2) * $y3) / $denom;
    printf STDERR "point1 = (%f,%f); point2 = (%f,%f); point3 = (%f,%f)\n",$x1,$y1,$x2,$y2,$x3,$y3;
    printf STDERR "a %f, b = %f, c = %f\n", $a,$b,$c;

    # compute best fit parabola (NEW VERSION uses rms to compute parabola - complicated to debug
    #$dof = $bb[12];
    #$dof = 1;
    #next unless $ii>1;
    #if ($ii==2) {
    # $dep[0]=0.; $lerr[0]=$lerr[1];
    # } else {
    #   $dep[0] = 2*$dep[1]-$dep[2];		$lerr[0]=$lerr[2];
    # }
    # $dep[$ii] = 2*$dep[$ii-1]-$dep[$ii-2];
    # $lerr[$ii] = $lerr[$ii-2];
    # $best++ if $best==1 && $ii>2 && $min==$lerr[2];
    # $adj=0.; $adj=0.001*$lerr[$best] if $lerr[$best-1] eq $lerr[$best] and $lerr[$best+1] eq $lerr[$best];
    # $d1 = $dep[$best]-$dep[$best-1];
    # $d2 = $dep[$best+1]-$dep[$best];
    # #printf STDERR "%d %d %e %e %e %d\n",$d1,$d2,$lerr[$best-1],$lerr[$best],$lerr[$best+1], $best;
    # $sigma = $d1*$d2*($d1+$d2)/($d2*($lerr[$best-1]-$lerr[$best])+$d1*($lerr[$best+1]-$lerr[$best])+$adj*($d1+$d2));
    # $depth = 0.5*($lerr[$best+1]-$lerr[$best-1])*$sigma/($d1+$d2);
    # $min = $lerr[$best] - $depth*$depth/$sigma;
    # #printf STDERR "%s \n H %5.1f %5.1f %f\n", $line[$best],$depth,$sigma,$min;
    # $sigma = sqrt($sigma*$min/$dof);
    # $depth = $dep[$best] - $depth;
    # #printf STDERR "%s H %5.1f %5.1f\n", $line[$best],$depth,$sigma;

    # compute best fit parabola (OLD VERSION uses rms to compute parabola - complicated to debug)
    #$dof = $bb[12];
    #next unless $ii>1;
    #if ($ii==2) {
    #  $dep[0]=0.; $rms[0]=$rms[1];
    #} else {
    #  $dep[0]   = 2*$dep[1]-$dep[2];		$rms[0]=$rms[2];
    #}
    #$dep[$ii] = 2*$dep[$ii-1]-$dep[$ii-2];	$rms[$ii]=$rms[$ii-2];
    #$best++ if $best==1 && $ii>2 && $min==$rms[2];
    #$adj=0.; $adj=0.001*$rms[$best] if $rms[$best-1] eq $rms[$best] and $rms[$best+1] eq $rms[$best];
    #$d1 = $dep[$best]-$dep[$best-1];
    #$d2 = $dep[$best+1]-$dep[$best];
    #$sigma = $d1*$d2*($d1+$d2)/($d2*($rms[$best-1]-$rms[$best])+$d1*($rms[$best+1]-$rms[$best])+$adj*($d1+$d2));
    #$depth = 0.5*($rms[$best+1]-$rms[$best-1])*$sigma/($d1+$d2);
    #$min = $rms[$best] - $depth*$depth/$sigma;
    #$sigma = sqrt($sigma*$min/$dof);
    #$depth = $dep[$best] - $depth;
    #printf STDERR "%s H %5.1f %5.1f\n", $line[$best],$depth,$sigma;

    #-------------------------------
    $xmin = $dep[1]; $xmax = $dep[$ii-1];
    $ymin = -10; $ymax = 100;
    $R = "-R$xmin/$xmax/$ymin/$ymax";

    # define $R2
    $ymin2 = -.005; $ymax2 = $max;
    $R2 = "-R$xmin/$xmax/$ymin2/$ymax2";
    $xtick1 = 5; $xtick2 = 1;
    $ytick1 = $ymax2/5.; $ytick2 = $ymax2/10.;
    $B2 = "-Ba${xtick1}f${xtick2}:\" \":/a${ytick1}f${ytick2}:\"ln(ERR / ERR_min)\":nW";
    #-------------------------------
    
    $xinc = 0.5;
    $tmp = 1000000;   # temporary variable (start with very large misfit value to find the uncertainty)
    $err_ellipse = 0.01;   # to compute uncertainity (
    open(PLT, "| psxy $J $R2 $xx");     # key command that sets the scale (J) and region (R)
    printf STDERR "---min_dep=%f max_dep=%f----\n", $dep[1],$dep[$ii-1];
    for ($l=$dep[1];$l<$dep[$ii-1];$l+=$xinc) {
	$xcord = $l;
	$ycord = $a*$l*$l + $b*$l + $c;   # Y = ax^2 + bx + c
	#printf STDERR "%f %f\n", $xcord,$ycord;
	#$aa = ($l-$depth)/$sigma;
	printf PLT "%6.3f %6.3f\n",$xcord, $ycord; # plot parabola -- note: (depth) misfit function is not necessarily quadratic
	if (abs($ycord - $err_ellipse) < $tmp){
	    $tmp = abs($ycord - 0.01);
	    $unc = abs($l-$depth);
	}
    }
    close(PLT);

    # plot the misfit curve
    # $xinc = 1;
    # $tmp = 1000;   # temporary variable (start with very large misfit value to find the uncertainty)
    # open(PLT, "| psxy $J $R2 $xx");     # key command that sets the scale (J) and region (R)
    # for ($l=$dep[0];$l<$dep[$ii];$l+=$xinc) {
    #   $aa = ($l-$depth)/$sigma;
    #   printf PLT "%6.3f %6.3f\n",$l,$aa*$aa; # plot parabola -- note: (depth) misfit function is not necessarily quadratic
    #    if (abs($aa*$aa - 0.01) < $tmp){
    # 	  $tmp = abs($aa*$aa - 0.01);
    # 	  $unc = abs($l-$depth);
    #   }
    # }
    # close(PLT);
    # printf STDERR "unc = %d \n",$unc;

    # plot the log(err/min_err) with depth
    $xy = "-K $B2 -P -O";
    open(PLT, "| psxy $J $R2 $xy -Sc0.25c -Gred");     # key command that sets the scale (J) and region (R)
    for ($jj=1;$jj<$ii;$jj+=1) {
      $l=$dep[$jj];
      $aa = $lerr[$jj];
      printf PLT "%6.3f %6.3f\n",$l,$aa;
      #printf STDERR "%f %f %f\n",$l,$aa,$depth;
    }
    close(PLT);
    $xy = "-K $B2 -O";
    open(PLT, "| psxy $J $R2 $xy -Wthick,red");     # key command that sets the scale (J) and region (R)
    for ($jj=1;$jj<$ii;$jj+=1) {
      $l=$dep[$jj];
      $aa = $lerr[$jj];
      printf PLT "%6.3f %6.3f\n",$l,$aa;
    }
    close(PLT);

    # plot the VR (variance reduction) with depth
    $xy = "-K -O $B3";
    open(PLT, "| psxy $J $R $xy -Sc0.25c -Gblue");     # key command that sets the scale (J) and region (R)
    for ($jj=1;$jj<$ii;$jj+=1) {
      $l=$dep[$jj];
      $aa = $vr[$jj];
      printf PLT "%6.3f %6.3f\n",$l,$aa;
    }
    close(PLT);
    open(PLT, "| psxy $J $R $xy -Wthick,blue");     # key command that sets the scale (J) and region (R)
    for ($jj=1;$jj<$ii;$jj+=1) {
	$l=$dep[$jj];
      $aa = $vr[$jj];
      printf PLT "%6.3f %6.3f\n",$l,$aa;
    }
    close(PLT);

    # plot the beach balls
    #  open(PLT, "| psmeca -JX -R -O -K -Sa0.3");   # original
    open(PLT, "| psmeca $J $R2 -O -K -Sm${ballsize} -G100 -W0.5p,0 -P"); # plot Moment tensors
    if ($onlydc==1){
	open(PLT, "| psmeca $J $R2 -O -K -Sd${ballsize} -G100 -W0.5p,0 -P"); # plot Moment tensors
    }
    for ($l=1; $l<$ii; $l++) {
      $coordx=$dep[$l];		            # x-coord for moment tensor
      #$coordy=($rms[$l]-$min)/($min/$dof);  # y-coord for moment tensor
      #$coordy=($lerr[$l]-$min)/($min/$dof);  # y-coord for moment tensor
      #$coordy=log($rms[$l]/$rms[$best]);
      $coordy=$lerr[$l];
      #printf STDERR "%f %f \n",$coordx,$coordy;
      printf STDERR "%s \n", $line[$l];
      printf PLT "%6.1f %6.6f 0.0 %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %f 0.0 0.0 %4.2f\n",
        $coordx,$coordy,$mrr[$l], $mtt[$l], $mff[$l], $mrt[$l], $mrf[$l], $mtf[$l], $M0_exp[$l], $mw[$l];
      # output values to screen
      printf STDERR "%6.1f %6.1f 0.0 %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %f 0.0 0.0 %4.2f\n",
        $coordx,$coordy,$mrr[$l], $mtt[$l], $mff[$l], $mrt[$l], $mrf[$l], $mtf[$l], $M0_exp[$l], $mw[$l];
    }
    close(PLT);

    # plot event id, best depth, misfit value
    if ($i < $#aa) {
      open(PLT, "| pstext -JX $R -N -O -K");
    }			# when doing the depth test for multiple events (output for various depths should be precomputed for each event)
    else {
      open(PLT, "| pstext -JX $R -N -O -S2p,white -N");
    }

    # plot the title
    $xtitle = $dep[1];
    $ytitle = $ymax + ($ymax-$ymin)*0.05;
    $fsizet = $fsize+2;
    printf PLT "%f %f $fsizet 0 0 1 %s  h=%4.1f \261 %.1f km\n",$xtitle,$ytitle,$eve,$depth,$unc; # $sigma gives much larger estimation of uncertainties
    close(PLT);
    $xx = "-O -K -Y2 $B";           # shift up
    $i++;
    if ($i == 5) {
      #$xx = "-X3.5 -O -K -Y-8 -Ba50f5/a20f5WSne";
      #$xx = "-X3.5 -O -K -Y-8 ${Bscale}WSne";
      $xx = "-X3.5 -O -K -Y-8 $B";  # shift to a new column
    }

  }  # foreach $eve (@aa){
}    # while(@event) {
exit(0);
