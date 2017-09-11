#!/usr/bin/env perl

=pod
Contents of an example parameter file (See: $caphome/EXAMPLES/flags.pf)
Usage:
     $ cap.pl flags.pf
XXX update paramter names to somthing more sensible

=cut

sub sub_read_parameter_file {

    my($parameter_file)=@_;

    open(IN,$parameter_file) || die "Could not read parameter file $parameter_file \n";
    printf STDERR "PARAMETER FILE: $parameter_file \n";
    @pflines = <IN>; $Npars = @pflines;
    close(IN);

    #  Read variables (these could be same as the input flags; but I think that will become confusing soon)
    for ($ii = 0; $ii < $Npars; $ii++) {
        ($opt,$val) = split("=",@pflines[$ii]);
        $opt =~ s/^\s+|\s+$//g;  # trim spaces
        $val=~ s/^\s+|\s+$//g;
        @value = split("/",$val);

        #printf STDERR "$opt @value \n";
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
	    # ($md_dep,$mg) = @value;
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
	} elsif ($opt eq "e") {
	    @event = grep(!/^-/,@value);
	} else {
	    printf STDERR $usage;
	    exit(0);
	}
    }
}

1;
