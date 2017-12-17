function[] = cap_fmt_misfit(datadir, evid, model, depth, nsol)
%
% Process waveform and polarity misfit. Output data files to show the misfit on
% the lune.
%
% USAGE
%   [gamma, delta, kappa, theta, sigma] = cap_fmt_misfit(datadir, evid, model, depth, nsol)
%
% EXAMPLES
%   cap_fmt_misfit('OUTPUT_DIR', '20100516063454464', 'utuhalf', 4, 22121190)
%   cap_fmt_misfit('OUTPUT_DIR', '20090407201255351', 'scak', 39, 22121190)
%   cap_fmt_misfit('OUTPUT_DIR', 'HOYA', 'wes', 1, 22121190)
%
% INPUT
%   Binary CAP files "grid_bb"  generated during inversion
%
% OUTPUT
%   Summary file for first motion polarities
%   Waveform misfit with psmeca files for plotting on the lune
%   The best solution in psmeca format
%
% WARNING
%   This function can read grid and random binary files.
%   This function creates subsets by polarities, and is fast when working with 
%   grid data. It's very slow with rand data.
%
% NOTE
%   Solutions are subset according to FMP==0 but this can be changed.
%   See section SUBSET FROM POLARITIES.
%   Magnitudes are assumed constant, so M0 = 22 (This is for plotting. See 'imag').
%   This code is based on the python scripts:
%       get_fmt_misfit_wf.py
%       get_fmt_misfit_fmp.py
% 
% Please cite this work as:
%
% C. Alvizuri and C. Tape. Full moment tensors for small events (Mw < 3) at
% Uturuncu volcano, Bolivia. Geophys. J. Int., 206(3):1761{1783, 2016. doi:
% 10.1093/gji/ggw247.
%
% 20160316 celso alvizuri - cralvizuri@alaska.edu 
%----------------------------------------------------------

close all;
TOL = 1e-10;
iplot = 0;

tic

%---------------------------------------
% filenames
%---------------------------------------
depth_str = sprintf('%03d', depth);
nsol0 = sprintf('%09d', nsol);

filename_key = [datadir, '/', evid,'_', model,'_', depth_str, '_'];
inputfile           = [filename_key, 'grid_bb', '_', nsol0, '.bin'];
out_misfit_wf       = [filename_key, 'misfit_wf_psmeca'];
out_misfit_fmp      = [filename_key, 'misfit_fmp'];
out_best_sol_psmeca = [filename_key, 'best_sol_psmeca'];

fprintf('\nRead cap binary output: %s \n', inputfile);
% revised function to read updated arrayMT
[v, w, kappa, theta, sigma, mag, misfit_wf, misfit_fmp, VR] = read_capbin_gd(inputfile, nsol);

%----------------------------------------------------------
% SUBSET FROM POLARITIES
%----------------------------------------------------------
fprintf('create subset of solutions ... ');
%ind = find(misfit_fmp ~= 0);    % all lune points except 0
%ind = find(misfit_fmp == 0);    % only lune points where FMP misfit==0.
ind = find(misfit_fmp >= 0);    % all lune points (new default)
% optionally use TOL 

nsol_pol = length(ind);
if(nsol_pol == 0)
    fprintf('\n\nWARNING. N pol solutions = %d. There are no solutions to process.\n\n', nsol_pol);
%    return;
end

% execute only when there are solutions to process
if(nsol_pol > 0)
    v_sub          = v(ind);
    w_sub          = w(ind);
    kappa_sub      = kappa(ind);
    theta_sub      = theta(ind);
    sigma_sub      = sigma(ind);
    mag_sub        = mag(ind);
    misfit_wf_sub  = misfit_wf(ind);
    misfit_fmp_sub = misfit_fmp(ind);
    VR_sub         = VR(ind);

    fprintf('N = %d \n', length(v_sub));

    % unique coordinates within set of solutions allowed by polarity
    fprintf('find unique coordinates in subset ... \n');
    %A = [v_sub, w_sub, kappa_sub, theta_sub, sigma_sub, mag_sub, misfit_wf_sub, misfit_fmp_sub];
    [UC, ia, ib] = unique([v_sub, w_sub], 'rows');
    vu = UC(:,1);
    wu = UC(:,2);

    %% max VR values
    nu = length(vu);  % should match nv*nu (see CAP output)
    ifinal0 = NaN(nu,1);

    fprintf('find best solutions at each coordinate ... \n');
    for ii=1:nu
        vtar = vu(ii);
        wtar = wu(ii);
        % indices in [v_sub, w_sub] that match unique coords UC
        imatch = find(and(vtar==v_sub, wtar==w_sub));
        [VRmax, iVRmax] = max(VR_sub(imatch));
        % iVRmax is index within imatch (a set), not within VR0
        ifinal0(ii) = imatch(iVRmax);
    end

    % get final data points
    % v, w, gamma, delta, kappa, theta, sigma, VR

    % v, w, VR
    vf  = v_sub(ifinal0);
    wf  = w_sub(ifinal0);
    VRf = VR_sub(ifinal0);

    [~, ibest_sol] = max(VRf)

    % gamma, delta, kappa, theta, sigma,
    rad2deg = 180.0 / pi;
    gammaf = v2gamma(v_sub(ifinal0));
    gammaf = gammaf * rad2deg;
    u = w_sub(ifinal0) + 3.0*pi/8.0;
    beta = u2beta(u);
    deltaf = beta - pi/2;
    deltaf = deltaf * rad2deg;
    kappaf = kappa_sub(ifinal0) * rad2deg;
    thetaf = theta_sub(ifinal0) * rad2deg;
    sigmaf = sigma_sub(ifinal0) * rad2deg;

    %% plot (for checking only)
    if iplot==1
        [vrmax, vrmax_ind] = max(VRf);
        [gamma_best, delta_best] = rect2lune(vf(vrmax_ind), wf(vrmax_ind));

        [~,imax] = max(VRmax);
        stlab = sprintf('max of %.2f at (%.1f, %.1f)',VRmax(imax),gamma_best,delta_best);
        figure; hold on;
        scatter(gammaf, deltaf, 8^2, VRf, 'filled');
        colormap(flipud(colormap)); caxis([4 20]);
        %title(stlab);
    end

    %----------------------------------------------------------
    % write results to files
    %----------------------------------------------------------
    fprintf('\nsave data to files\n');

    % output file EVID_MODEL_DEPTH_misfit_wf_psmeca
    % gamma, delta , VR  mrr, mtt, mff, mrt, mrf, mtf, exp

    fmtarray = TT2CMT(gammaf, deltaf, 1, kappaf, thetaf, sigmaf);
    ntensors = length(gammaf);
    imag = 22;  % beachball size (M0) for psmeca

    fid = fopen(out_misfit_wf, 'w');
    % fix issue 1. solutions with large -VR saturate lune plot
    % fix issue 2. psmeca doesn't plot values with -VR 
    for i=1:ntensors
        if (VRf(i) < 0)
            disp('WARNING, VR<0. Setting VR = 0');
            VRf(i) = 0;
        end
        fprintf(fid, '%14.7e %14.7e %14.7e  %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %4d\n',...
        gammaf(i), deltaf(i), VRf(i), fmtarray(:,i), imag);
    end
    fclose(fid);

    % write file for best solution
    %fid = fopen('EVID_MODEL_DEPTH_best_sol_psmeca', 'w');
    fid = fopen(out_best_sol_psmeca, 'w');
    fprintf(fid, '%14.7e %14.7e %14.7e  %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %4d\n',...
        gammaf(ibest_sol), deltaf(ibest_sol), VRf(ibest_sol), fmtarray(:,ibest_sol), imag);
    fclose(fid);
end % end check nsol==0

%----------------------------------------------------------
% output file EVID_MODEL_DEPTH_misfit_fmp
% gamma, delta, n(Lambda), N orientations, fraction orientations
%----------------------------------------------------------

% unique coordinates within entire set of solutions
% WARNING may be slow
nsol = length(v);
fprintf('\nfind unique coordinates in whole set (N = %d), may take a couple of minutes for large N ... \n', nsol);
[UC, ia, ib] = unique([v, w], 'rows');
vu = UC(:,1);
wu = UC(:,2);

% This is bruteforce but relatively quick ( <1 min for 2e7 on eagle). 
% Will leave as is for now.
nu = length(vu);  % should match nv*nu (see CAP output)

%fid = fopen('EVID_MODEL_DEPTH_misfit_fmp', 'w');
fid = fopen(out_misfit_fmp, 'w');
[gamma_fmp, delta_fmp] = rect2lune(vu, wu);
for ii=1:nu
    vtar = vu(ii);
    wtar = wu(ii);
    % indices in [v, w] that match unique coords UC
    imatch = find(and(vtar==v, wtar==w));
    [misfit_fmp_min, imisfit_fmp_min] = min(misfit_fmp(imatch));
    imatch2 = find(misfit_fmp(imatch) == misfit_fmp_min);
    norientations = length(misfit_fmp(imatch2));
%    [igamma_fmp, idelta_fmp] = rect2lune(vtar, wtar);
    fprintf(fid, '%11.6f %11.6f %d %8d %14.7e\n', gamma_fmp(ii), delta_fmp(ii), misfit_fmp_min, norientations, norientations/nsol);
%    fprintf(fid, '%11.6f %11.6f %d %8d %14.7e\n', vtar, wtar, misfit_fmp_min, norientations, norientations/nsol);
end
fclose(fid);

%% write results for testing
%dataf = [vf, wf, gammaf, deltaf,  VRf];
%fid = fopen('testdata.txt', 'w');
%for i=1:length(dataf)
%        fprintf(fid, '%11.6f %11.6f %11.6f %11.6f %11.6f\n', dataf(i, 1), dataf(i, 2), dataf(i, 3),  dataf(i, 4),  dataf(i, 5));
%end
%fclose(fid);

fprintf('Total number of solutions processed: %d\n', length(v));
fprintf('Max VR: %4.1f (should match CAP output if polarities are used)\n\n', max(VRf));
toc

%-----------------------------------------------------------
% dummy example
if 1==0
%    datadir = '~/REPOSITORIES/cap/EXAMPLES/OUTPUT_DIR/';
    evid = '19910914190000000';
    model = 'wes';
    depth = 1;
    nsol = 22121190
    cap_fmt_misfit(datadir, evid, model, depth, nsol);
end
%-----------------------------------------------------------
