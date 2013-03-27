#!/usr/bin/env perl
#
# find best depth from CAP output and
# plot the misfit error and source mechanism
# as function of depth for each event
# Usage depth.pl result_file event_dir_names ...
#

$ballsize = 4;		    # controls default beachball size (psmeca)
$min0 = 1.0e+19;		# impossibly large misfit value
$Bscale = "-Ba5f1:\"Depth, km\":/a20f5:\"Misfit\":";

#------------

#($rsl,@event) = @ARGV;     # original
($rsl,$rsl2,@event) = @ARGV; # input Event and Tensor lines from files cus_NNN.out, and event id

open(RSL,"$rsl") or die "couldn't open inputfile1 $rsl\n";
@aaa=<RSL>;
#print STDERR "debug: aaa=@aaa\n";  # 20130103 calvizuri - ok to delete
close(RSL);

#-----------
# GMT plot settings

#$xtick1 = 50; $xtick2 = 5;
$xtick1 = 4; $xtick2 = 1;
$ytick1 = 20; $ytick2 = 5;
$B = "-Ba${xtick1}f${xtick2}:\"Depth, km\":/a${ytick1}f${ytick2}:\"Misfit relative to minimum\":WSne";

$zoom=3;
$magpsx=2.0*$zoom;
$magpsy=1.8*$zoom;
$magpsme=1.0*$zoom;
$magpstx=1.0*$zoom;

$pwidth_in = 8.5;		# width of paper
$pheight_in = 11;		# height of paper
$fsize = 16;
system("gmtset BASEMAP_TYPE plain PAPER_MEDIA Custom_${pwidth_in}ix${pheight_in}i MEASURE_UNIT inch LABEL_FONT_SIZE $fsize ANOT_FONT_SIZE $fsize");

#-----------

open(FH,"$rsl2") or die "couldn't open inputfile2 $rsl2\n";
@data_fmt=<FH>;
#print STDERR "debug: data_fmt=@data_fmt\n";
close(FH);

while (@event) {
  @aa = splice(@event,0,10); # 0=offset, 10=N elements (each element space separated)
  #print STDERR "debug. aa=@aa event=@event\n";
  $i=0;
  #$xx = "-K -Ba50f5/a20f5WSne";   # original <-- not enough ticks+labeling
  $xx = "-K ${Bscale}WSne -P";
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
      #    $strike[$ii]=$bb[5];    # not needed
      #    $dip[$ii]=$bb[6];       # not needed
      #    $rake[$ii]=$bb[7];      # not needed
      $mw[$ii]=$bb[9];
      #    printf STDERR "debug. mw[$ii]=%lf\n",$mw[$ii];
      $rms[$ii]=$bb[11];
      if ($min>$rms[$ii]) {
	$best=$ii;$min=$rms[$ii];
      }
      $ii++;
    }
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
      #    printf STDERR "debug. MT components: %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f M0_exp=%f\n",
      #        $mrr[$jj], $mtt[$jj], $mff[$jj], $mrt[$jj], $mrf[$jj], $mtf[$jj],$M0_exp[$jj];
      $jj++;
    }

    # compute best fit parabola
    $dof = $bb[12];
    next unless $ii>1;
    if ($ii==2) {
      $dep[0]=0.; $rms[0]=$rms[1];
    } else {
      $dep[0]   = 2*$dep[1]-$dep[2];		$rms[0]=$rms[2];
    }
    $dep[$ii] = 2*$dep[$ii-1]-$dep[$ii-2];	$rms[$ii]=$rms[$ii-2];
    $best++ if $best==1 && $ii>2 && $min==$rms[2];
    $adj=0.; $adj=0.001*$rms[$best] if $rms[$best-1] eq $rms[$best] and $rms[$best+1] eq $rms[$best];
    $d1 = $dep[$best]-$dep[$best-1];
    $d2 = $dep[$best+1]-$dep[$best];
    $sigma = $d1*$d2*($d1+$d2)/($d2*($rms[$best-1]-$rms[$best])+$d1*($rms[$best+1]-$rms[$best])+$adj*($d1+$d2));
    $depth = 0.5*($rms[$best+1]-$rms[$best-1])*$sigma/($d1+$d2);
    $min = $rms[$best] - $depth*$depth/$sigma;
    $sigma = sqrt($sigma*$min/$dof);
    $depth = $dep[$best] - $depth;
    printf STDERR "%s H %5.1f %5.1f\n", $line[$best],$depth,$sigma;

#-------------------------------

  $xmin = $dep[0]; $xmax = $dep[$ii];
  $ymin = -10; $ymax = 100;
  $R = "-R$xmin/$xmax/$ymin/$ymax";

#-------------------------------

    # plot the misfit curve
    $xinc = 0.2;
    #  open(PLT, "| psxy -JX3/1.8 -R$dep[0]/$dep[$ii]/-10/100 $xx");    #original
    #open(PLT, "| psxy -JX16/16 $R $xx"); # portrait mode, larger+square plot
    open(PLT, "| psxy -JX$magpsx/$magpsy $R $xx");
    for ($l=$dep[0];$l<$dep[$ii];$l+=$xinc) {
      $aa = ($l-$depth)/$sigma;
      printf PLT "%6.3f %6.3f\n",$l,$aa*$aa; # plot parabola -- note: (depth) misfit function is not necessarily quadratic
    }
    close(PLT);
 
    # plot the beach balls
    #  open(PLT, "| psmeca -JX -R -O -K -Sa0.3");   # original
    open(PLT, "| psmeca -JX -R -O -K -Sm${ballsize} -G100 -W0.5p,0"); # plot Moment tensors
    for ($l=1;$l<$ii;$l++) {
      $coordx=$dep[$l];		# x-coord
      $coordy=($rms[$l]-$min)/($min/$dof); # y-coord for moment tensor
      #   printf PLT "%6.1f %6.1f 0 %s %s %s %s 0 0 %s\n",$dep[$l],($rms[$l]-$min)/($min/$dof),$strike[$l],$dip[$l],$rake[$l],$mw[$l],$mw[$l];   # original
      printf PLT "%6.1f %6.1f 0.0 %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %f 0.0 0.0 %4.2f\n",
        $coordx,$coordy,$mrr[$l], $mtt[$l], $mff[$l], $mrt[$l], $mrf[$l], $mtf[$l], $M0_exp[$l], $mw[$l];
      # output values to screen
      printf STDERR "%6.1f %6.1f 0.0 %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %f 0.0 0.0 %4.2f\n",
        $coordx,$coordy,$mrr[$l], $mtt[$l], $mff[$l], $mrt[$l], $mrf[$l], $mtf[$l], $M0_exp[$l], $mw[$l];
    }
    close(PLT);

    # plot event id, best depth, misfit value
    if ($i<$#aa) {
      open(PLT, "| pstext -N -JX -R -O -K");
    }			# when is this executed?
    else {
      open(PLT, "| pstext -JX -N -R -O -S2p,white -N");
    }

    # plot the title
    $xtitle = $dep[0];
    $ytitle = $ymax + ($ymax-$ymin)*0.05;
    $fsizet = $fsize+2;
    printf PLT "%f %f $fsizet 0 0 1 %s  h=%4.1f pm %.1f km\n",$xtitle,$ytitle,$eve,$depth,$sigma;
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
