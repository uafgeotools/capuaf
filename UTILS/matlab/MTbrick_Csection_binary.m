function [M, mag, misfit_wf, misfit_fmp] = MTbrick_Csection_binary(eid,ddir,smodel,dep,nsamples,K,X,gmtdir)
% script to plot the cross section of strike-dip-rake space
% INPUT: log file for grid search (May be just use the brick elements)
% eid       event id 
% ddir      data directory containing the 2 grid binary files 
% gmtdir    gmt directory to save the cross-section (optional)
% Mref      Cross-sections at the moment tensor other than the minimum
%           misfit sample (optional)
%
% TASK THIS SCRIPT SHOULD PEFORM:
% (1) Cross-sections of strike-dip-rake brick (Input: kappa0,theta0,sigma0 (or Mref))
% (2) Contours of omega (from kappa0,theta0,sigma0 (or Mref)) on the cross-sections
% See run_CAP_unc.m for example
% Also see CAP_unc.m
%
% Vipul Silwal
% 27 July, 2015

%clear all, % close all

if nargin<8
    igmt = 0;
elseif nargin==8
    igmt = 1;
end

% Reference moment tensor (instead of the minimum misfit solution)
Mref = [];

% INPUT PAR 
% Number of samples
Nfname = sprintf('%09d',nsamples);

% model tag
dep = sprintf('%03d',dep);
model = strcat(smodel,'_',dep);

% filename types
ftype1 = '_grid_bb_';
ftype2 = '_grid_mt_';

%ddir = strcat(ddir,'/',eid,'/','OUTPUT_DIR','/');

% Read the binary grid file from the data directory
grid_log_filename1 = strcat(ddir,eid,'_',model,ftype1,Nfname,'.bin')
grid_log_filename2 = strcat(ddir,eid,'_',model,ftype2,Nfname,'.bin')

% read.out file
outfname = strcat(ddir,eid,'_',model,'.out');
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,Nstn,~,~,~,~,~,~,~,~,~,~,~,~,~,pPol,ePol] = read_capout(outfname);
%[otime,elat,elon,edep,strike,dip,rake,M,Mw,eid,capmod,capdep,rms,vr,Pwin,Swin,Nstn,Pstn,Sstn,...
%    stnm,lampPV,lampPR,lampSV,lampSR,lampSH,corrPV,corrPR,corrSV,corrSR,corrSH] = read_capout(outfname)

% Constants
deg = 180/pi;
% ---These are set in total_misfit.m ----
misfit_fact = K;         % misfit factor (this should be consistent with the scale factor in CAP_unc_binary.m)
%misfit_fact_wf = 40;         % misfit factor
%misfit_fact_fmp = 30;         % misfit factor (TBD)
misfit_pol_weight = X;   % weight for first motion  polarity misfit when adding to waveform misfit 
% --------------------------
%P = 0.8;        % Proability constants (See Walt's notes)
%ismooth = 0;    % to smooth the cdf
%imesa = 0;      % Use analytical formulation (1); estimate pdf from the sample distribution (0)
domega = 1/deg;
iTT2CMT = 0;
iminus = 0;

% Read log file
%[kappa,theta,sigma,F,Mw,delta,gamma,M0,M] = read_CAP_log(grid_log_filename);
[gamma, delta, kappa, theta, sigma, mag, misfit_wf, misfit_fmp,VR] = read_capbin_gd(grid_log_filename1,nsamples);
[mrr, mtt, mpp, mrt, mrp, mtp, mag, misfit_wf, misfit_fmp2] = read_capbin_mt(grid_log_filename2,nsamples);
M = [mrr mtt mpp mrt mrp mtp]';
kappa = kappa*deg;
theta = theta*deg;
sigma = sigma*deg;
%isequal(misfit_fmp,misfit_fmp2) % check

% Scale and combine misfit_wf and misfit_fmp (KEY)
%misfit_wf = misfit_fact_wf * misfit_wf;
%misfit_fmp = misfit_fact_fmp*misfit_fmp/Nstn;
%F = misfit_wf + misfit_fmp;
Np = length(find(pPol~=0)); % Number of station at which polarity is predicted
[F,misfit_wf,misfit_fmp] = total_misfit(misfit_wf,misfit_fmp,Np,misfit_pol_weight);
F = F*misfit_fact; % scale misfit

figure;
subplot(2,2,1)
hist(misfit_fmp,20)
subplot(2,2,2)
hist(misfit_wf,20)
subplot(2,2,3)
hist(F,20)

% Get the moment tensors (Another Option - takes longer time)
if iTT2CMT
    gamma = 0;
    delta = 0;
    M0 = 1/sqrt(2);     % such that |Lambda| = 1
    M = TT2CMT(gamma,delta,M0,kappa,theta,sigma);
end

% vectors of strike-dip-rak (edges of the brick)
kappa_vec = unique(kappa);
theta_vec = unique(theta);
sigma_vec = unique(sigma);
whos

% find global minimum (Mref)
[min_misfit,minindx] = min(F);
[min_misfit_wf,minindx_wf] = min(misfit_wf);
[min_misfit_fmp,minindx_fmp] = min(misfit_fmp);
% indx = 1; % For random checking

% cross-sections at these values of strike, dip and rake
if isempty(Mref) % Mref not specified
    kappa0 = kappa(minindx)
    sigma0 = sigma(minindx)
    theta0 = theta(minindx)
    Mref = M(:,minindx);
else
    [gamma0,delta0,M00,kappa0,theta0,sigma0] = CMT2TT(Mref);
    [ival,indx] = min(abs(kappa_vec-kappa0));
    kappa0 = kappa_vec(indx);
    [ival,indx] = min(abs(theta_vec-theta0));
    theta0 = theta_vec(indx);
    [ival,indx] = min(abs(sigma_vec-sigma0));
    sigma0 = sigma_vec(indx);
    Mref = TT2CMT(0,0,1, kappa0, theta0,sigma0);
end
% get the auxiliary fault plane
[kappa2,theta2,sigma2] = dcfaultpar2aux(kappa0,theta0,sigma0);

% Compute omega
omega = CMT2omega(Mref,M);

% scale misfit
%F = F * misfit_fact;

% Get sections at the minimum misfit solutions
ik = find(kappa==kappa0);
Z1 = [sigma(ik) theta(ik) F(ik)];
C1 = [sigma(ik) theta(ik) omega(ik)];
W1 = [sigma(ik) theta(ik) misfit_wf(ik)];
P1 = [sigma(ik) theta(ik) misfit_fmp(ik)];

it = find(theta == theta0);
Z2 = [sigma(it) kappa(it) F(it)];
C2 = [sigma(it) kappa(it) omega(it)];
W2 = [sigma(it) kappa(it) misfit_wf(it)];
P2 = [sigma(it) kappa(it) misfit_fmp(it)];

is = find(sigma==sigma0);
Z3 = [kappa(is) theta(is) F(is)];
C3 = [kappa(is) theta(is) omega(is)];
W3 = [kappa(is) theta(is) misfit_wf(is)];
P3 = [kappa(is) theta(is) misfit_fmp(is)];

% Plot waveform misfit sections
Msize = 30;
figure
subplot(2,2,1)
scatter(W1(:,1),W1(:,2),Msize,W1(:,3),'filled');
hold on; plot(sigma0,theta0,'p', 'MarkerSize',15,'MarkerFaceColor','w');
axis tight
subplot(2,2,2)
scatter(W2(:,1),W2(:,2),Msize,W2(:,3),'filled');
hold on; plot(sigma0,kappa0,'p', 'MarkerSize',15,'MarkerFaceColor','w');
axis tight
subplot(2,2,3)
scatter(W3(:,1),W3(:,2),Msize,W3(:,3),'filled');
hold on; plot(kappa0,theta0,'p', 'MarkerSize',15,'MarkerFaceColor','w');
axis tight
colorbar
subplot(2,2,4)
bb([kappa(minindx_wf) theta(minindx_wf) sigma(minindx_wf)]);
suptitle(sprintf('Waveform Misfit; kappa= %3.0f theta= %3.0f sigma=%3.0f',kappa0,theta0,sigma0));

if min(misfit_fmp)
    % Plot polarity misfit sections
    Ncols_colormap = Np;
    figure
    subplot(2,2,1)
    scatter(P1(:,1),P1(:,2),Msize,P1(:,3),'filled');
    h = colorbar;
    caxis([min(misfit_fmp) max(misfit_fmp)]);
    set(h,'ylim',[min(misfit_fmp) max(misfit_fmp)])
    
    %'ytick',linspace(min(misfit_fmp),max(misfit_fmp),Ncols_colormap));
    %# get current colormap
    map = colormap;
    %# adjust for number of colors you want
    rows = uint16(linspace(1, size(map,1), Ncols_colormap)) ;
    map = map(rows, :);
    %# and apply the new colormap
    colormap(map);
    hold on; %plot(sigma0,theta0,'p', 'MarkerSize',15,'MarkerFaceColor','w');
    axis tight
    %----------------------------
    subplot(2,2,2)
    scatter(P2(:,1),P2(:,2),Msize,P2(:,3),'filled');
    %colormap(map);
    h=colorbar;
    caxis([min(misfit_fmp) max(misfit_fmp)]);
    %set(h,'ylim',[0 0.5],'ytick',[0:0.01:0.5])
    hold on; %plot(sigma0,kappa0,'p', 'MarkerSize',15,'MarkerFaceColor','w');
    axis tight
    %----------------------------
    subplot(2,2,3)
    scatter(P3(:,1),P3(:,2),Msize,P3(:,3),'filled');
    hold on; %plot(kappa0,theta0,'p', 'MarkerSize',15,'MarkerFaceColor','w');
    axis tight
    %colormap(map);
    h = colorbar;
    caxis([min(misfit_fmp) max(misfit_fmp)]);
    %set(h,'ylim',[min(misfit_fmp) max(misfit_fmp)])
    %set(h,'ylim',[0 0.5],'ytick',[0:0.1:0.5])
    axis tight
    %----------------------------
    %subplot(2,2,4)
    %bb([kappa(minindx_fmp) theta(minindx_fmp) sigma(minindx_fmp)]);
    suptitle(sprintf('Polarity Misfit; kappa= %3.0f theta= %3.0f sigma=%3.0f',kappa0,theta0,sigma0));
end

% Plot total misfit sections
Msize = 30;
figure
subplot(2,2,1)
scatter(Z1(:,1),Z1(:,2),Msize,Z1(:,3),'filled');
hold on; plot(sigma0,theta0,'p', 'MarkerSize',15,'MarkerFaceColor','w');
axis tight
subplot(2,2,2)
scatter(Z2(:,1),Z2(:,2),Msize,Z2(:,3),'filled');
hold on; plot(sigma0,kappa0,'p', 'MarkerSize',15,'MarkerFaceColor','w');
axis tight
subplot(2,2,3)
scatter(Z3(:,1),Z3(:,2),Msize,Z3(:,3),'filled');
hold on; plot(kappa0,theta0,'p', 'MarkerSize',15,'MarkerFaceColor','w');
axis tight
colorbar
subplot(2,2,4)
bb([kappa0 theta0 sigma0]);
suptitle(sprintf('Total Misfit; kappa= %3.0f theta= %3.0f sigma=%3.0f',kappa0,theta0,sigma0));

% Plot contours of omega
figure
subplot(2,2,1)
scatter(C1(:,1),C1(:,2),Msize,C1(:,3),'filled');
axis tight
subplot(2,2,2)
scatter(C2(:,1),C2(:,2),Msize,C2(:,3),'filled');
axis tight
subplot(2,2,3)
scatter(C3(:,1),C3(:,2),Msize,C3(:,3),'filled');
axis tight
suptitle('Omega');

kappa0*deg
theta0*deg
sigma0*deg

%% script to save (Z's and C's) input file for gmt

if and(igmt,~isempty(gmtdir))
    %plotdir=strcat(gmtdir,eid,'/',norm,'/',wts{ii},'/');
    %plotdir = '~/rough/';
    % gmtdir = '/home/vipul/gmt/data/cap/test/';
    %plotdir = '~/rough/';
    %eid='1';
    %plotdir=strcat(gmtdir,eid,'/',norm,'/',wts{ii},'/');
    mkdir(gmtdir);
    
    %copyfile(strcat(datadir,'FIGS/',eid,'/',norm,'/',wts{ii},'/',strcat(model,'_',edep,'.out')),strcat(plotdir,eid,'.out'));
    %copyfile(strcat(evtdir,eid,'/',strcat(eid,'_event.dat')),strcat(plotdir));
    %copyfile(strcat(evtdir,eid,'/',strcat(eid,'_station.dat')),strcat(plotdir));
    
    % copy CAP output to GMT plotting directory
    %outfile = dir(strcat(ddir,'out'));
    copyfile(outfname,strcat(gmtdir,eid,'.out'));
    % copy event and station location file to GMT plotting directory
    %copyfile(strcat(ddir,'..',eid,'/',eid,'_event.dat'),gmtdir);
    copyfile(strcat(ddir,'../',eid,'/',eid,'_station_list_ALL.dat'),gmtdir);
    
    
    % XYZ files
    fid = fopen(strcat(gmtdir,eid,'_stkfile.dat'),'w');
    fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',Z1');
    fclose(fid);
    fid = fopen(strcat(gmtdir,eid,'_dipfile.dat'),'w');
    fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',Z2');
    fclose(fid);
    fid = fopen(strcat(gmtdir,eid,'_rakfile.dat'),'w');
    fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',Z3');
    fclose(fid);
    
    % XYZ for omega contour
    fid = fopen(strcat(gmtdir,eid,'_stkcont.dat'),'w');
    fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',C1');
    fclose(fid);
    fid = fopen(strcat(gmtdir,eid,'_dipcont.dat'),'w');
    fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',C2');
    fclose(fid);
    fid = fopen(strcat(gmtdir,eid,'_rakcont.dat'),'w');
    fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',C3');
    fclose(fid);
    
    % Write the Mref
    fid = fopen(strcat(gmtdir,eid,'_Mref.dat'),'w');
    fprintf(fid,'%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f',[kappa0;theta0;sigma0;kappa2;theta2;sigma2]);
    fclose(fid);
end

%% EXAMPLES
if 0
    % Update example for IRIS data
    eid = '20090407201253480';  % IRIS catalog id
    ddir = '/home/vipul/REPOSITORIES/cap';
    dep = 41;
    nsamples = 50616;     % Number of samples (grid)
    smodel = 'scak';        % velocity model
    gmtdir = '/home/vipul/REPOSITORIES/cap/20090407201253480_gmt_data/';
    %gmtdir = []; 
    MTbrick_Csection_binary(eid,ddir,dep,gmtdir);
    
    eid = '20160124123742054';  % IRIS catalog id
    ddir = '/home/vipul/PROJECTS/josh/';
    dep =107;
    nsamples = 49248;     % Number of samples (grid)
    smodel = 'scak';        % velocity model
    gmtdir = '/home/vipul/PROJECTS/josh/20160124123742054/gmt_data/';
    MTbrick_Csection_binary(eid,ddir,dep,gmtdir);
    
    % Default example 
    eid='20090407201255351';
    ddir = '/home/vipul/CAP/inv/scak/';
    dep = 41;
    gmtdir = '~/CAP/inv/scak/20090407201255351_gmt_data/';
    MTbrick_Csection_binary(eid,ddir,dep);

    %------NEHRP exmaple
    eid = '20140124120703813';
    ddir = '/home/vipul/CAP/inv/scak/NEHRP/20140124120703813/';
    MTbrick_Csection_binary(eid,ddir);
    
    % Kaltag example
    eid='19990517173210860';
    ddir = '/home/vipul/CAP/inv/tact/kaltag/19990517173210860/';
    MTbrick_Csection_binary(eid,ddir);
end

