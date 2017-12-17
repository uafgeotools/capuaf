function[] = stp(datadir, evid, model, depth, nsamples, Kfactor, nv, n)
% 
% Compute source type probability
%
% USAGE
%   stp(datadir, evid, model, depth, nsamples, Kfactor, nv, n);
% 
% EXAMPLES
%   stp('OUTPUT_DIR', '20100516063454464', 'utuhalf', 4, 22121190, 40, 11)
%   stp('OUTPUT_DIR', '20090407201255351', 'scak', 39, 22121190, 40, 11)
%   stp('OUTPUT_DIR', 'HOYA', 'wes', 1, 22121190, 40, 11)
% 
% INPUT 
%   datadir  - directory where CAP binary files are stored
%   evid     - event id
%   model    - model name
%   depth    - inversion depth
%   nsamples - number of solutions processed by CAP
%   Kfactor  - adjust prob density function
%   nv       - number of grid divisions (nw is derived from nv)
%   n        - (OPTIONAL) n-point smoothing
%
% OUTPUT
%   file 1  - prob(v, w)
%   file 2  - prob(gamma, delta)
%   file 3  - misfit(v, w) 
%   file 4  - 
%   file 5  - 
% 
% Please cite this work as:
%
%   Alvizuri, C., Silwal, V., Krischer, L., Tape, C. Estimation of full moment
%   tensors with uncertainties, for earthquakes, volcanic events, and nuclear
%   tests. Geophysics (in prep.).
%
% References:
%
%   Tape, W., and C. Tape, 2016, A confidence parameter for seismic moment
%   tensors: Geophys. J. Int., 205, 938–953.
%   Silwal, V., and C. Tape, 2016, Seismic moment tensors and estimated
%   uncertainties in southern Alaska: J. Geophys. Res. Solid Earth, 121,
%   2772–2797.
%
% 20160606 cralvizuri@alaska.edu
%-----------------------------------------------------------

close all;

ihist = 0;
ioutdata = 1;
icheckplot = 1;
iresidual = 0;
weight_pol = 0; 

% datadir, evid, model, depth, nsamples, Kfactor, nv, n
if nargin < 8
    ismooth = 0;
elseif nargin == 8
    ismooth = 1;
end

%---------------------------------------
% filenames
%---------------------------------------
nsamp_str = sprintf('%09d', nsamples);
depth_str = sprintf('%03d', depth);
model_depth = [model, '_', depth_str];
filename_key = [datadir, '/' ,evid, '_', model_depth];

% input
cap_out_file = [filename_key, '.out'];
capout_rand_bb = [filename_key, '_rand_bb_', nsamp_str, '.bin'];
capout_rand_mt = [filename_key, '_rand_mt_', nsamp_str, '.bin'];

% output
out_vw_p_density = [filename_key, '_', 'vw_p_density.txt'];
out_gd_p_density = [filename_key, '_', 'gd_p_density.txt'];
out_vw_misfit = [filename_key, '_', 'vw_misfit.txt'];

%---------------------------------------
% process data
%---------------------------------------
% read CAP binary output to get moment tensor and misfit data
fprintf('reading the data ...\n');
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,Nstn,~,~,~,~,~,~,~,~,~,~,~,~,~,pPol,ePol] = read_capout(cap_out_file);
[V, W, kappa, theta, sigma,  mag, misfit_wf, misfit_fmp, VR] = read_capbin_gd(capout_rand_bb, nsamples);
[mrr, mtt, mpp, mrt, mrp, mtp, mag, misfit_wf2, misfit_fmp2] = read_capbin_mt(capout_rand_mt, nsamples);
M = [mrr mtt mpp mrt mrp mtp]';

% compare best VR with VR from all solutions
if ihist==1
    figure;
    subplot(1,2,1);
    hist(VR, 1000);
    title('VR');
    subplot(1,2,2);
    hist(misfit_wf, 1000);
    title('misfit');
end

% combine misfit_wf and misfit_fmp
Np = length(find(pPol~=0)); % Number of station at which polarity is predicted

% TEST FUNCTION
%[V, W, misfit_wf, VR] =  get_gaussian3D(length(V));

fprintf('compute combined misfit ...\n')
[combined_misfit, misfit_wf, misfit_fmp] = total_misfit(misfit_wf, misfit_fmp, Np, weight_pol);

% Discretization for v and w
% NOTE  This value influences the location of best sol!
%       V is 9pi/8 times the size of W
nw = round(nv * 9.0 * pi / 8.0);
vmax = 1.0 / 3.0;
wmax = 3.0 * pi / 8.0;

% correct for numerical issue. apply offset 1e-8 at the boundaries. 
% note: offset >= 1e-8 
voffset = 1e-8;      
woffset = 1e-8;

v1 = -vmax - voffset; v2 = vmax + voffset; 
w1 = -wmax - woffset; w2 = wmax + woffset; 

% NOTE output center coordinates, not edges
dvcenter = ((vmax * 2) / (nv - 1)) * (1/2);
dwcenter = ((wmax * 2) / (nw - 1)) * (1/2);

vvec = linspace(v1, v2, nv);
wvec = linspace(w1, w2, nw);

[vgrid, wgrid] = meshgrid(vvec, wvec);
[vq, wq] = meshgrid(vvec, wvec);

fprintf('binning the data ...\n')
[count_in_ibinV, mapV2ibinV] = histc(V, vvec);
[count_in_ibinW, mapW2ibinW] = histc(W, wvec);

%-----------------------------------------------------------
% compute VR
%-----------------------------------------------------------
% find lowest misfit in each bin
vw_misfit = accumarray([mapW2ibinW mapV2ibinV], VR, [nw nv], @max);

% smoothing
if ismooth==1
    vw_misfit = conv2(vw_misfit, ones(n)/n.^2,'same');
end

% location of highest VR
% NOTE depends on bin size!
[VRmax, VRmaxindex] = max(vw_misfit(:));
[indWbestVR, indVbestVR] = ind2sub(size(vw_misfit), VRmaxindex);
vbestVR = vgrid(1, indVbestVR);
wbestVR = wgrid(indWbestVR, 1);
[gammaVR, deltaVR] = rect2lune(vbestVR, wbestVR);

%-----------------------------------------------------------
% compute probability
%-----------------------------------------------------------
misfit = combined_misfit * Kfactor;
p_unorm = exp(-misfit);      % unnormalized
P0 = sum(p_unorm);

% probability
vw_prob = accumarray([mapW2ibinW mapV2ibinV], p_unorm, [nw nv], @sum);
vw_prob = vw_prob / P0;
sum_pnorm = sum(sum(vw_prob));

% p density
n1 = (nv - 1) * (nw - 1);   % n1 cells in Q
qarea = pi / (2.0 * n1);
vw_p_density = vw_prob / qarea;
sum_p_density = sum(sum(vw_p_density)) * qarea;

% smoothing
if ismooth==1
    vw_p_density = conv2(vw_p_density, ones(n)/n.^2,'same');
end

% location of highest probability
% NOTE depends on bin size!
[Pmax, Pmaxindex] = max(vw_p_density(:));
[indWbestP, indVbestP] = ind2sub(size(vw_p_density), Pmaxindex);
vbestP = vgrid(1, indVbestP);
wbestP = wgrid(indWbestP, 1);
[gammaP, deltaP] = rect2lune(vbestP, wbestP);

%-----------------------------------------------------------
% compute residual (for checking only)
%-----------------------------------------------------------
% difference between misfit and probability
residual = (vw_misfit/max(max(vw_misfit))) - (vw_p_density/max(max(vw_p_density)));

if ismooth==1
    residual = conv2(residual, ones(n)/n.^2, 'same');
end

% DEBUG
%whos vgrid wgrid vw_p_density vw_misfit vw_misfit residual

%-----------------------------------------------------------
% save data to file
%-----------------------------------------------------------
if ioutdata == 1
    fprintf('writing data to files\n')
    % prob(v,w)
    fid = fopen(out_vw_p_density, 'w');
    for iv = 1:length(vvec)
        for iw = 1:length(wvec)
            fprintf(fid, '%11.7f %11.7f %11.7f\n', vvec(iv), wvec(iw), vw_p_density(iw, iv));
        end
    end
    fclose(fid);

    % prob(g,d)
    fid = fopen(out_gd_p_density, 'w');
    for iv = 1:length(vvec)
        for iw = 1:length(wvec)
            [igamma, idelta] = rect2lune(vvec(iv), wvec(iw));
            fprintf(fid, '%11.7f %11.7f %11.7f\n', igamma, idelta, vw_p_density(iw, iv));
        end
    end
    fclose(fid);

    % misfit(v,w)
    fid = fopen(out_vw_misfit, 'w');
    for iv = 1:length(vvec)
        for iw = 1:length(wvec)
            fprintf(fid, '%11.7f %11.7f %11.7f\n', vvec(iv), wvec(iw), vw_misfit(iw, iv));
        end
    end
    fclose(fid);
end

%-----------------------------------------------------------
% plot results (for checking only)
%-----------------------------------------------------------
if icheckplot == 1
    if iresidual == 1
        nfigs = 3;
    else
        nfigs = 2;
    end

    fprintf('plotting results (for checking only)\n');
    figure; 
    % VR or misfit on vw space
    subplot(1,nfigs,1)
    surfc(double(vgrid), double(wgrid), double(vw_misfit));
    hc1 = colorbar('southOutside');
    xlabel(hc1,'VR') 
    view(2)
    colormap(flipud(jet))
    grid on;
    axis equal;
    axis tight;
    label = sprintf('max(VR)\nv %6.2f w %6.2f\n gamma %5.1f delta %5.1f', vbestVR, wbestVR, gammaVR, deltaVR);
    title(label)
    hold on;
    plot3(vbestVR + dvcenter, wbestVR + dwcenter, VRmax, 'o', 'markerfacecolor', 'white')
    grid on;

    % probability on vw space
    subplot(1,nfigs,2)
    surfc(double(vgrid), double(wgrid), double(vw_p_density));
    view(2)
    hc2 = colorbar('southOutside');
    xlabel(hc2,'PD') 
    grid on;
    axis equal;
    axis tight;
    hold on;
    label = sprintf('K %d, int(p_d) %5.2f\nv %6.2f w %6.2f\n gamma %5.1f delta %5.1f',...
    Kfactor, sum_p_density, vbestP, wbestP, gammaP, deltaP);
    title(label)
    plot3(vbestP + dvcenter, wbestP + dwcenter, Pmax, 'o', 'markerfacecolor', 'white')

    if iresidual == 1
        subplot(1,nfigs,3)
        surfc(double(vgrid), double(wgrid), double(residual)); %
        shading interp;
        colorbar('southOutside');
        view(2)
        grid on;
        axis equal;
        axis tight;
        title('residual')
    end;
end; % icheckplot

fprintf('Done.\n');
