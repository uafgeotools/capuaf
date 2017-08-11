#!/usr/bin/env perl
#
# find best depth from CAP output and
# plot the misfit error and source mechanism
# as function of depth for each event
# Usage depth.pl result_file event_dir_names ...
#

$ballsize = 2.5;		    # controls default beachball size (psmeca)
$min0 = 1.0e+19;	    # impossibly large misfit value
$Bscale = "-Ba5f1:\"\":/a20f5:\"\":";
$onlydc = 0; # only plots DC mechanism (Removes psmeca bugs that arises when plotting DC mechs)
$imodel = 1; # to draws bars at layer interfaces

#------------
# Structural model info : depth of layer interface
@tactmod = (3,11,24,31,76);
@scak = (4,9,14,19,24,33,49,66);
@cus =(1,10,20,30);
@wes = (2.5,32.5);
#------------

#($rsl,@event) = @ARGV;      # original
($rsl,$rsl2,@event) = @ARGV; # input Event and Tensor lines from files cus_NNN.out, and event id
#(@event,$smod)=@ARGV;


# name of junk files for storing header info while plotting
#$rsl = "OUTPUT_DIR/junk1.out";
#$rsl2 = "OUTPUT_DIR/junk2.out";

# system("grep -h Event OUTPUT_DIR/*$smod*.out > $rsl");
#system("grep -h tensor OUTPUT_DIR/*$smod*.out > $rsl2");

#------------
open(RSL,"$rsl") or die "couldn't open inputfile1 $rsl\n";
@aaa=<RSL>;
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

# VR
$xtick1 = 5; $xtick2 = 1;
$ytick1 = 20; $ytick2 = 5;
$B3 = "-Ba${xtick1}f${xtick2}:\"Depth, km\":/a${ytick1}f${ytick2}:\"VR (gray)\":ES";

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

#=====================READ output files=============================
while (@event)
{
    @aa = splice(@event,0,10); # 0=offset, 10=N elements (each element space separated)
    #$odir = "$aa[0]";  # If all the output files (*.out) are in the $eve directory
    $odir = "OUTPUT_DIR";  # If all the output files (*.out) are in the OUTPUT_DIR directory
    # print STDERR "debug. aa=@aa event=@event\n";
    $i=0;
    #$xx = "-K -Ba50f5/a20f5WSne";   # original <-- not enough ticks+labeling
    #$xx = "-K ${Bscale}WS -P";
    $xx = "-K -P $xoffset $yoffset";
    foreach $eve (@aa)
    {				# eve = current aa = event id
	$ii=1;
	$best=1;
	$min=$min0;
	$max=-100.;
	foreach (grep(/$eve/,@aaa))
    {
        chop;   
        # Split long names into 1, 2 or 3 strings
        # 20130103 calvizuri - example input to chop:
        # Event 20080418093700 Model cus_001 FM 291 51 -13 Mw 5.10 rms 3.748e-02   110 ERR   2   5   8 ISO 0.16 0.11 CLVD 0.14 0.08
        # New format:
        #   0           1          2               3                4   5    6         7  8   9   10   11        12    13   14  15     16      17  18   19    20
        # Event 20160124123742054 Model 20160124123742054_scak_100 FM   50 59.063911   46 Mw 4.40 rms 7.310e-01  3660 CLVD 0.00 ISO   0.000000 VR 46.6 data2 1.297e-06
        # OR
        # Event Little_Skull_Main Model Little_Skull_Main_wes_001 FM  243 67.006607  -21 Mw 5.10 rms 8.262e-01  2880 CLVD 0.00 ISO   0.000000 VR 31.7 data2 1.986e-05
        $line[$ii]=$_;
        @bb=split;
        $nwords = split('_',$bb[3]);
        if ($nwords == 3) {
            ($evname, $smodel, $dep[$ii]) = split('_',$bb[3]);
        }
        elsif ($nwords == 4) {
            ($evname1, $evname2, $smodel, $dep[$ii]) = split('_',$bb[3]); 
            $evname = join '', $evname1, " ", $evname2;
        }
        elsif ($nwords == 5) {
            ($evname1, $evname2, $evname3,$smodel,$dep[$ii]) = split('_',$bb[3]);
            $evname = join '', $evname1, " ", $evname2, " ", $evname3;
        }

        $strike[$ii]=$bb[5];		   # not needed
        $dip[$ii]=$bb[6];			   # not needed
        $rake[$ii]=$bb[7];			   # not needed
        $mw[$ii]=$bb[9];
        #    printf STDERR "debug. mw[$ii]=%lf\n",$mw[$ii];
        $rms[$ii]=$bb[11]; # this is the total misfit= polarity error + waveform error
        $vr[$ii]=$bb[18];

        if ($vr[$ii]>$max)
        { 
            $max=$vr[$ii]; $best=$ii;
        }
        $ii++;
    }
        # get the catalog depth from line #2 of the CAP output file
        # NEED A STATEMENT TO EXIT IF THE FILE DOES NOT EXIST
        $bfile = "./${odir}/${eve}_${smodel}_$dep[${best}].out";
	open(OUT,$bfile);
	@outfile=<OUT>;
       (undef,undef,undef,$elat,undef,$elon,undef,$edep)  =split(" ",$outfile[1]);
	printf STDERR "catalog depth (from sac header) is $edep\n";
    
	# ------------------read MT parameters --------------------------
	$jj=1;
	foreach (grep(/tensor/,@data_fmt))
	{
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
	    #printf STDERR "debug. MT components: %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f M0_exp=%f\n",$mrr[$jj], $mtt[$jj], $mff[$jj], $mrt[$jj], $mrf[$jj], $mtf[$jj],$M0_exp[$jj];
	    $jj++;
	}
    
	#===================== compute a relative measure of error ($lerr = log(misfit/misfit_min))
	$max=-1000;		# unrealisticly low maximum value
	for ($jj=1;$jj<$ii;$jj=$jj+1)
	{
	    $lerr[$jj] = $vr[$best]/$vr[$jj];
	    $lerr[$jj] = log($lerr[$jj]);
	    # printf STDERR "%f %e %e %d %d %d\n",$lerr[$jj],$rms[$jj],$rms[$best],$dep[$best], $best,$ii;
	    if ($lerr[$jj]>$max)
	    {
		$max=$lerr[$jj];
	    }
	}
	$max = sprintf("%1.3f",$max); # suppress to 3 decimal places

	#====================== compute parabolic equation:  Y= aX^2 + bX + c (cleaner way to handle) Using best 3 points only
	$x1 = $dep[$best-1]; $y1 = $lerr[$best-1]; #(x1,y1)
	$x2 = $dep[$best]; $y2 = $lerr[$best];	   #(x2,y2)
	$x3 = $dep[$best+1]; $y3 = $lerr[$best+1]; #(x3,y3)
	$denom = ($x1-$x2)*($x1-$x3)*($x2-$x3);
	# Below are the standard equation (check any maths textbook)
	$a = (($x3*($y2-$y1))+($x2 * ($y1 - $y3)) + ($x1 * ($y3 - $y2))) / $denom;
	$b = ($x3*$x3 * ($y1 - $y2) + $x2*$x2 * ($y3 - $y1) + $x1*$x1 * ($y2 - $y3)) / $denom;
	$c = ($x2 * $x3 * ($x2 - $x3) * $y1 + $x3 * $x1 * ($x3 - $x1) * $y2 + $x1 * $x2 * ($x1 - $x2) * $y3) / $denom;
	printf STDERR "point1 = (%f,%f); point2 = (%f,%f); point3 = (%f,%f)\n",$x1,$y1,$x2,$y2,$x3,$y3;
	printf STDERR "a %f, b = %f, c = %f\n", $a,$b,$c;
	# compute minima of parabola (at dY/dx = 0)
	$depth = -$b/(2*$a);	# minima occurs at this X
	$Ydepth = $a*$depth*$depth + $b*$depth +$c; # Yaxis value of parabola at minima

	#================GMT plotting ranges==========================
	$xmin = $dep[1]; $xmax = $dep[$ii-1];
	$ymin = -10; $ymax = 100;
	$R = "-R$xmin/$xmax/$ymin/$ymax";

    # define $R2
    $ymin2 = -$max/10; $ymax2 = $max;
    $R2 = "-R$xmin/$xmax/$ymin2/$ymax2";
    $xtick1 = 5; $xtick2 = 1;

    # show y-axis ticks up to 2 or 3 decimal places
    $ytick1 = $ymax2/5.;    # annotations
    $ytick2 = $ymax2/10.;   # locations
    # ticks for log err
    # y-axis ticks don't show for small log-err. This code handles those cases
    if ($ymax2 < 0.1) {
        $ytick1 = sprintf("%.3f", $ytick1);
        $ytick2 = sprintf("%.3f", $ytick2);
    }
    else {
        $ytick1 = sprintf("%.2f", $ytick1);
        $ytick2 = sprintf("%.2f", $ytick2); 
    }
    $B2 = "-Ba${xtick1}f${xtick2}:\" \":/a${ytick1}f${ytick2}:\"ln(VR_max / VR)\":nW";
    
	# Set the model for plotting layer interface
	if ($smodel eq "tactmod") {
	    @model=@tactmod;
	} elsif ($smodel eq "scak") {
	    @model=@scak;
	} elsif ($smodel eq "cus") {
	    @model=@cus;
	} elsif ($smodel eq "wes") {
	    @model=@wes;
	}

	#================== PLot misfit parabola
	$xinc = 0.1;
	$tmp = 1000000;	# temporary variable (start with very large misfit value to find the uncertainty)
	$err_cent = 0.01; # to compute uncertainity ($err_cent*100 percent confidence interval)
	open(PLT, "| psxy $J $R2 $xx -W1p,0/0/0,-"); # key command that sets the scale (J) and region (R)
	printf STDERR "---min_dep=%f max_dep=%f----\n", $dep[1],$dep[$ii-1];
	for ($l=$dep[1]; $l<$dep[$ii-1]; $l+=$xinc)
	{
	    $xcord = $l;
	    $ycord = $a*$l*$l + $b*$l + $c; # Y = ax^2 + bx + c  (value of Y at each $xcord using parabolic eq)
	    #printf STDERR "%f %f\n", $xcord,$ycord;
	    #$aa = ($l-$depth)/$sigma;
	    printf PLT "%f %f\n",$xcord, $ycord; # plot parabola -- note: (depth) misfit function is not necessarily quadratic
	    # compute the uncertainy based on $err_cent value 
	    if (abs($ycord - $err_cent) < $tmp)
	    {
		$tmp = abs($ycord - $err_cent);
		$unc = abs($l-$depth);
		$xcord2 = $l;      
		$ycord2 = $ycord;
	    }
	}
	close(PLT);

	# ==============Plot bars for uncertainty
	open(PLT, "| psxy $J $R2 -K -O -W2p,0/0/0");
	printf PLT "%f %f\n",$depth, $err_cent;
	printf PLT "%f %f\n",$xcord2, $err_cent;
	close(PLT);

	open(PLT, "| psxy $J $R2 -K -O -W1p,0/0/0,-");
	printf PLT "%f %f\n",$depth, $Ydepth;
	printf PLT "%f %f\n",$depth, $err_cent;
	close(PLT);
    
	# open(PLT, "| psxy $J $R2 -K -O -W0.5p,0/0/0,-");
	# printf PLT "%f %f\n",$xcord2, $ymin2;
	# printf PLT "%f %f\n",$xcord2, $err_cent;
	# close(PLT);

	# =============plot the log(err/min_err) with depth (perhaps don't need this)
	$xy = "-K $B2 -P -O";
	open(PLT, "| psxy $J $R2 $xy -Sc0.25c -Gred"); # key command that sets the scale (J) and region (R) - draws a circle
	for ($jj=1;$jj<$ii;$jj+=1)
	{
	    $l=$dep[$jj];
	    $aa = $lerr[$jj];
	    printf PLT "%f %f\n",$l,$aa;
	    #printf STDERR "%f %f %f\n",$l,$aa,$depth;
	}
	close(PLT);
	$xy = "-K $B2 -O";
	open(PLT, "| psxy $J $R2 $xy -Wthick,black"); # key command that sets the scale (J) and region (R) - drays a line
	for ($jj=1;$jj<$ii;$jj+=1)
	{
	    $l=$dep[$jj];
	    $aa = $lerr[$jj];
	    printf PLT "%f %f\n",$l,$aa;
	}
	close(PLT);

	#=================== plot the VR (variance reduction) with depth
	$xy = "-K -O $B3";
	open(PLT, "| psxy $J $R $xy -Sc0.25c -G150"); # key command that sets the scale (J) and region (R) - draws a circle
	for ($jj=1;$jj<$ii;$jj+=1)
	{
	    $l=$dep[$jj];
	    $aa = $vr[$jj];
	    printf PLT "%f %f\n",$l,$aa;
	}
	close(PLT);
	open(PLT, "| psxy $J $R $xy -Wthick,150"); # key command that sets the scale (J) and region (R)- drays a line
	for ($jj=1;$jj<$ii;$jj+=1)
	{
	    $l=$dep[$jj];
	    $aa = $vr[$jj];
	    printf PLT "%f %f\n",$l,$aa;
	}
	close(PLT);
    
	#================== Plot the bars where interface occurs in the structural model
	if ($imodel==1)
	{
	    for ($kk=0; $kk<=$#model; $kk++)
	    {
		open(PLT, "| psxy $J $R $xy -W2p,blue");
		printf PLT "%f %f\n",$model[$kk], -10;
		printf PLT "%f %f\n",$model[$kk], 0;
		close(PLT);
	    }
	}

	#================== Plot the depth(s) as inverted triangles
        
        # plot the catalog depth ($edep) as a RED inverted triangle
	open(PLT, "| psxy $J $R $xy -Si0.5c -G255/0/0 -W.05c");
	printf PLT "%f %f\n",$edep, -7.5;

    # Uncomment the following for depth test main event, FMT Uturuncu paper.
    # This is a tweak to plot red triangle at catalog depth.
    # NOTE Catalog depth for this event is 0.6 below sea level.
    # The depth test for the main event gives a best depth of 4.4km, then rounded to 4km.
    # But this is 4km from the surface, which is at elevation 4.6km. 
    # Therefore the inversion is with respect to the elevation:
    # 4.6 (elevation) - 4 (best depth) = 0.6 km above sea level.
    #printf PLT "%f %f\n",$edep + 4.6, -7.5;   # uncomment here for plotting catalog depth

	close(PLT);

        # plot the best-fitting depth ($depth) from the CAP grid search as a WHITE inverted triangle
        open(PLT, "| psxy $J $R $xy -Si0.5c -G255/255/255 -W.05c");
	printf PLT "%f %f\n",$depth, -7.5;
	close(PLT);

	#=============== plot the beach balls
	#  open(PLT, "| psmeca -JX -R -O -K -Sa0.3");   # original
	open(PLT, "| psmeca $J $R2 -O -K -Sm${ballsize} -G100 -W0.5p,0 -P"); # plot Moment tensors
	if ($onlydc==1)
	{
	    open(PLT, "| psmeca $J $R2 -O -K -Sd${ballsize} -G100 -W0.5p,0 -P"); # plot Moment tensors
	}
	for ($l=1; $l<$ii; $l++)
	{
	    $coordx=$dep[$l];	# x-coord for moment tensor
	    #$coordy=($rms[$l]-$min)/($min/$dof);  # y-coord for moment tensor
	    #$coordy=($lerr[$l]-$min)/($min/$dof);  # y-coord for moment tensor
	    #$coordy=log($rms[$l]/$rms[$best]);
	    $coordy=$lerr[$l];
	    #printf STDERR "%f %f \n",$coordx,$coordy;
	    printf STDERR "%s \n", $line[$l];
	    printf PLT "%f %f 0.0 %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %f 0.0 0.0 %4.2f\n",
	    $coordx,$coordy,$mrr[$l], $mtt[$l], $mff[$l], $mrt[$l], $mrf[$l], $mtf[$l], $M0_exp[$l], $mw[$l];
	    # output values to screen
	    printf STDERR "%6.1f %6.1f 0.0 %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %f 0.0 0.0 %4.2f\n",
	    $coordx,$coordy,$mrr[$l], $mtt[$l], $mff[$l], $mrt[$l], $mrf[$l], $mtf[$l], $M0_exp[$l], $mw[$l];
	}
	close(PLT);

	# ==============plot event id, best depth, misfit value
	if ($i < $#aa)
	{
	    open(PLT, "| pstext -JX $R -N -O -K");
	} # when doing the depth test for multiple events (output for various depths should be precomputed for each event)
	else
	{
	    open(PLT, "| pstext -JX $R -N -O -S2p,white -N");
	}

	# plot the title
	$xtitle = $dep[1];
	$ytitle = $ymax + ($ymax-$ymin)*0.05;
	$fsizet = $fsize+2;

    # $sigma gives much larger estimation of uncertainties
    # NOTE  for figure in Uturuncu FMT paper replace 'h' with 'depth', 'utuhalf' with more descriptive 'halfspace'
    #       Also follow comments #1 and #2 next
    printf PLT "%f %f $fsizet 0 0 1 %s | Model %s | Best depth %4.1f \261 %.1f km\n",$xtitle,$ytitle,$evname,$smodel,$depth,$unc;       # 1. comment for Uturuncu FMT paper
    #printf PLT "%f %f $fsizet 0 0 1 %s | halfspace | depth%4.1f \261 %.1f km\n",$xtitle,$ytitle, $evname,$depth,$unc; # 2. uncomment for Uturuncu FMT paper

	close(PLT);
	$xx = "-O -K -Y2 $B";	# shift up
	$i++;
	if ($i == 5)
	{
	    #$xx = "-X3.5 -O -K -Y-8 -Ba50f5/a20f5WSne";
	    #$xx = "-X3.5 -O -K -Y-8 ${Bscale}WSne";
	    $xx = "-X3.5 -O -K -Y-8 $B"; # shift to a new column
	}

    }				# foreach $eve (@aa){
}				# while(@event) {
exit(0);
