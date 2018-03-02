function [Mref, post_samples_info,unc_result]= CAP_unc_binary(eid,ddir,smodel,dep,nsamples,K,X,gmtdir)
% script for doing posterior analysis and uncertainty estimate on MT
% Input : Log files from CAP (contains MT info and misfit values)
%   eid     - event id
%   ddir    - directory where the log files are stored
%   smodel  - velocity model (eg: scak)
%   dep     - depth
%   nsmaples- Number of sample points (binary files need this as input)
%   K       - Misfit scale factor
%   X       - polarity misfit weight factor
%   gmtdir  - (optional) directory where you want to save the file fot GMT plotting
%   
%
% Also see:
%   MTbrick_Csection_binary.m   : For getting misfit cross-secitons of
%                                 strike-dip-rake grid
%   CAP_summary_plot.pl         : For making summary plots using GMT
%  
% Vipul Silwal
% 20 July, 2015

% clear all,% close all

if nargin<8
    igmt = 0;
elseif nargin==8
    igmt = 1;
end

% Reference moment tensor (instead of the minimum misfit solution)
Mref = [];

% Output directory
% odir = 'output_dir/';
% ddir = strcat(ddir,'/',eid,'/',odir); % Example: /home/vipul/CAP/inv/scak/MOOS/20090407201255351/output_dir/

% File tag
Nfname = sprintf('%09d',nsamples);
dep = sprintf('%03d',dep);
model = strcat(smodel,'_',dep);

% filename types
ftype1 = '_grid_bb_';
ftype2 = '_grid_mt_';

%ddir = strcat(ddir,'/',eid,'/','OUTPUT_DIR','/');

% filename types
ftype1 = '_rand_bb_';
ftype2 = '_rand_mt_';

% Read the binary random file from the data directory
% XXX: Check what takes less time to read (binary or mat file)
rand_log_filename1 = strcat(ddir,eid,'_',model,ftype1,Nfname,'.bin')
rand_log_filename2 = strcat(ddir,eid,'_',model,ftype2,Nfname,'.bin');

% read.out file
outfname = strcat(ddir,eid,'_',model,'.out')
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,Nstn,~,~,~,~,~,~,~,~,~,~,~,~,~,pPol,ePol] = read_capout(outfname);

% Constants and Input parameters
deg = 180/pi;
%K = 5;         % misfit factor
P = 0.8;        % Probability constants (See Walt's notes)
misfit_pol_weight = X;   % weight for first motion  polarity misfit when adding to waveform misfit 
ismooth = 1;    % to smooth the cdf
imesa = 1;      % (1) Use analytical formulation; (0)estimate pdf from the sample distribution
idc = 1;        % (1) Double couple only;  (0) Full moment tensor
domega = 1/(deg);
iTT2CMT = 0;
Nsamp_reject = 10000;
iplot = 1;
ipost_samples_file = 0;
ipost_beach_file = 0;

% % Read files and Get the moment tensors
% [kappa,theta,sigma,F,Mw,delta,gamma,M0,M] = read_CAP_log(rand_log_filename);

% [gamma, delta, kappa, theta, sigma, mag, misfit_wf, misfit_fmp] = read_capbin_gd(rand_log_filename1,nsamples);
% [M(1,:),M(2,:),M(3,:),M(4,:),M(5,:),M(6,:), mag, misfit_wf, misfit_fmp] = read_capbin_mt(rand_log_filename2,nsamples);
% n = length(kappa);
[V, W, kappa, theta, sigma, mag, misfit_wf, misfit_fmp,VR] = read_capbin_gd(rand_log_filename1,nsamples);
[mrr, mtt, mpp, mrt, mrp, mtp, mag, misfit_wf, misfit_fmp] = read_capbin_mt(rand_log_filename2,nsamples);
M = [mrr mtt mpp mrt mrp mtp]';
kappa = kappa*deg;
theta = theta*deg;
sigma = sigma*deg;

% combine misfit_wf and misfit_fmp (KEY)
Np = length(find(pPol~=0)); % Number of station at which polarity is predicted
[F,misfit_wf,misfit_fmp] = total_misfit(misfit_wf,misfit_fmp,Np,misfit_pol_weight);

if iTT2CMT      % Another Option (takes longer time)
    V = 0;
    W = 0;
    M0 = 1/sqrt(2);     % such that |Lambda| = 1
    M = TT2CMT(V,W,M0,kappa,theta,sigma);
end

% find global minimum (Mref)
[ival,indx] = min(F)
% indx = 1;
kappa0 = kappa(indx);
sigma0 = sigma(indx);
theta0 = theta(indx);
Mmin = M(:,indx);

% You can choose any Mref (reference moment tensor) !
if isempty(Mref)
    Mref = Mmin; % default option
    iMref = 0; % to change the position of Mref (default: Mref = Mmin)
else
    iMref = 1; % to change the position of Mref (default: Mref = Mmin)
end

% omega angle from Mref to all M
omega_deg = CMT2omega(M,Mref);
omega = omega_deg/deg;

%% probability density function
% Kvec = 10:10:40; % test for multiple scale factors (these are now set in
% total_misfit.m)
Kvec = K;
figure;
for ii=1:length(Kvec)
    K = Kvec(ii);
    misfit = F*K;
    p = exp(-misfit);      % Get probability (not normalized yet) from scaled misfit
    % Test to see when probability is same throughout
    % p = 0.5*ones(length(F),1);
    % Monte carlo integration: integrate 'p' over random 'omega' in 'domega' increment
    [omega_mid, prob1, cum_prob, N] = MC_int(omega,p,domega);
    % normalized cdf (MC_int don't do that)
    prob1 = prob1/(sum(prob1)*domega);
    cum_prob = cumsum(prob1)*domega;
    post_pdf = [(omega_mid*deg)'  smooth(prob1,20)]';
    post_cdf = [(omega_mid*deg)'  cum_prob]';
    % After this: sum(prob1)*domega = 1(i.e. normalized = area under the
    % curve is 1)
    % cum_prob(end) = 1

    % Analytical mesa
    if imesa == 1
        if idc ==1  % for double couple
            mesa = omegadcpdf(omega_mid,false);
        else        % for full moment tensor
            mesa = omegapdf(omega_mid,false);
        end
    else % compute the mesa function from the sample histogram
        [N,mesa,~] = plot_histo(omega,omega_mid(1)-domega/2:domega:omega_mid(end)+domega/2,3);
        %omega_mid=omega_mid';
        mesa = mesa';
        %whos mesa omega_mid
    end
    
    mesa = mesa/(sum(mesa)*domega);
    homo_prob = cumsum(mesa)*domega;
    mesa_pdf = [(omega_mid*deg)' mesa']';
    mesa_cdf = [(omega_mid*deg)' homo_prob']';
    whos mesa_pdf mesa_cdf mesa 
    
    % Compute volume v(P) for Mref (M0 or best solution)
    [ival,indx] = min(abs(cum_prob - P));
    OMEGA = omega_mid(indx);
    vP = homo_prob(indx);
    %unc_result = [P vP OMEGA*deg]';
    
    % Get the the volume v(P) for range of P
    Pvec = [0:0.01:1];
    [v,vindx] = volumetric_unc(cum_prob,homo_prob,Pvec);
    vP_vec = [Pvec' v']';
    end_points = [1 0;0 0];
    vP_vec_avg = [vP_vec end_points];
    
    % Plotting
    % Data from these figures will go to subfigure(b) of the Summary plots
    % figure(10);
    %----------------------------------------------------
    subplot(2,2,1)
    if ismooth
        prob1 = smooth(prob1,20);
    end
    plot(omega_mid,prob1,'g','LineWidth',2);
    hold on
    plot(omega_mid,mesa,'b','LineWidth',2);
    xlim([omega_mid(1)-domega/2 omega_mid(end)+domega/2]);
    axis equal
    axis([0 pi 0 2])
    title(sprintf('PDF, p = exp(-K x misfit); K = %3.0f; max(misfit) = %3.3f', K, max(misfit)))
    %----------------------------------------------------
    subplot(2,2,2)
    if ismooth
        prob1 = smooth(prob1,10);
    end
    plot(omega_mid,smooth(prob1,20),'g--','LineWidth',2);
    hold on
    plot(omega_mid,mesa,'b','LineWidth',2);
    xlim([omega_mid(1)-domega/2 omega_mid(end)+domega/2]);
    axis equal
    axis([0 pi 0 2])
    title(sprintf('smoothed PDF, p = exp(-K x misfit); K = %3.0f', K))
    %----------------------------------------------------
    subplot(2,2,3)
    plot(omega_mid,cum_prob,'g','LineWidth',2);
    hold on
    plot(omega_mid,homo_prob,'b','LineWidth',2);
    plot([0 omega_mid(indx)],[P P],'k--');
    plot([omega_mid(indx) omega_mid(indx)],[vP cum_prob(indx)],'k--');
    plot([0 omega_mid(indx)],[vP vP],'k--');
    %title(sprintf('Cumulative PDF, P = %2.2f, v(P) = %2.3f, Omega = %3.0f',P,vP,OMEGA*deg))
    title(sprintf('Cumulative PDF'));
    xlabel('omega');
    xlim([omega_mid(1)-domega/2 omega_mid(end)+domega/2])
    axis equal
    axis([0 pi 0 2]);
    %----------------------------------------------------
    subplot(2,2,4)
    plot(Pvec,v,'r-','LineWidth',2);
    hold on
    %plot([P P],[0 v(Pvec==P)],'k--');
    %plot([0 P],[v(Pvec==P) v(Pvec==P)],'k--');
    plot([0 1],[0 1],'k--')
    axis([0 1 0 1]);
    xlabel('V');
    ylabel('P(V)');
    %title(sprintf('P = %2.2f, v(P) = %1.3f',P,v(Pvec==P)));
    title(sprintf('P_{av} = %1.3f',sum(v)/length(v)));
    unc_result=sum(v)/length(v);
    axis
end

%% Check with the Rejection method (always keep number of posterior samples
% = 2000)
% (UPDATE)
%
chance = rand(nsamples,1);
ikeep = [];
while (length(ikeep) < Nsamp_reject)
    ik = find((p/max(p))>chance);
    ikeep = [ikeep ik];
end
ikeep = ikeep(1:Nsamp_reject);
Mpost = M(:,ikeep);

figure
plot_histo(omega(ikeep),0:2*domega:pi,3);
hold on
plot(omega_mid,prob1,'r','LineWidth',2);
xlabel('omega');
title('cumulative probability method (Red); Rejection method (Green)');
axis equal, axis tight;

if iplot
    % HISTOGRAM for the posterior samples in strike (kappa), dip (theta) and
    % rake (sigma)
    % Data from this figure will go to subfigure(e) of the Summary plots
    % TO DO: Add gamma delta histograms for FMT analysis
    figure
    subplot(3,2,1)
    plot_histo(kappa(ikeep),0:5:360);
    hold on
    plot(kappa0,0,'vr','MarkerSize',5)
    xlabel('Strike')
    %----------------------------------------------------
    subplot(3,2,2)
    plot_histo(theta(ikeep),0:5:90);
    hold on
    plot(theta0,0,'vr','MarkerSize',5)
    xlabel('Dip')
    %----------------------------------------------------
    subplot(3,2,3)
    plot_histo(sigma(ikeep),-90:5:90);
    hold on
    plot(sigma0,0,'vr','MarkerSize',5)
    xlabel('Rake')
    %----------------------------------------------------
    subplot(3,2,4)
    plot_histo(V(ikeep),-1/3:.005:1/3);
    hold on
    plot(sigma0,0,'vr','MarkerSize',5)
    xlabel('v')
    %----------------------------------------------------
    subplot(3,2,5)
    plot_histo(W(ikeep),-3*pi/8:.05:3*pi/8);
    hold on
    plot(sigma0,0,'vr','MarkerSize',5)
    xlabel('w')
end

% Auxiliary fault plane (for the same set of posterior moment tensor samples)
kappa_samples = kappa(ikeep);
theta_samples = theta(ikeep);
sigma_samples = sigma(ikeep);
[kappa_samples2,theta_samples2,sigma_samples2]=dcfaultpar2aux(kappa_samples,theta_samples,sigma_samples);

%========== get outline points for omega1 and rnd_err=========
nbin=60;
[omega_outline,misfit_outline,M_outline] = cappts2outline(omega_deg,misfit,M,nbin);

%% Pick few samples from posterior space for plotting
% For 400 randomly picked samples from posterior
Npost = 200;
ikeep_subset  = randperm(Nsamp_reject,Npost);
Mpost_subset = Mpost(:,ikeep_subset);
omega_subset = omega_deg(ikeep(ikeep_subset));
misfit_subset = misfit(ikeep(ikeep_subset));

% if you want to keep subset posterior samples same for multiple figures
% ipost_samples_file = 1;
if ipost_samples_file
    % read misfit_subset and Mpost_subset here
    % compute omega_subset = CMTomega(Mref,Mpost_subset)
    % Add file here
    fname = '/home/vipul/gmt/data/cap/NW/20140418184418152/L1/M111/20140418184418152_Mpost_file.dat';
    fname = '/home/vipul/gmt/data/cap/MOOS/20090407201255351/L1/M111/20090407201255351_Mpost_file.dat';
    file = dlmread(fname);
    Mpost_subset = file(:,1:6)';
    misfit_subset = file(:,7);
    omega_subset = CMT2omega(Mref,Mpost_subset);
end
Mpost_file = [Mpost_subset; misfit_subset'];

% compute P-T axis points
[lam,U] = CMTdecom(Mpost_subset);
Usamp = convertv(1,5,U);  % convert GCMT to south-east-up
Upasamp = U2pa(Usamp,1,0);
Upasamp(:,3:4) = [];     % cut the nodal axis
%Upasamp(end+1,:) = [0 0 0 90];
Tp = Upasamp(:,1);
Ta = Upasamp(:,2);
Pp = Upasamp(:,3);
Pa = Upasamp(:,4);
[Txs,Tys] = pa2xy(Tp,Ta);
[Pxs,Pys] = pa2xy(Pp,Pa);

% omega and misfit for Npost sample points (green sample points in GMT fig)
post_samples = [omega_subset misfit_subset Pxs Pys Txs Tys]';
% PT axis for Npost sample points
%PTaxis = [Pxs Pys Txs Tys]';
post_samples_info = [Mpost_subset' omega_subset misfit_subset];
% figure;plot(omega_subset,misfit_subset,'.')

% 20 random samples from the posterior space
Nbeach = 20;
ikeep_beach  = randperm(Nsamp_reject,Nbeach);
Mbeach = Mpost(:,ikeep_beach);
omega_beach = omega_deg(ikeep(ikeep_beach));
misfit_beach = misfit(ikeep(ikeep_beach));
kappa_beach = kappa(ikeep(ikeep_beach));
theta_beach = theta(ikeep(ikeep_beach));
sigma_beach = sigma(ikeep(ikeep_beach));

% if you want to keep same 20 beachballs subset
% ipost_beach_file = 1;
if ipost_beach_file
    % read misfit_beach and Mbeach here
    % compute omega_beach = CMTomega(Mref,Mbeach)
    % Add file here
    fname = '/home/vipul/gmt/data/cap/NW/20140418184418152/L1/M111/20140418184418152_Mbeach_file.dat';
    fname = '/home/vipul/gmt/data/cap/MOOS/20090407201255351/L1/M111/20090407201255351_Mbeach_file.dat';
    file = dlmread(fname);
    Mbeach = file(:,1:6)';
    misfit_beach = file(:,7);
    omega_beach = CMT2omega(Mref,Mbeach);
    [~,~,~,kappa_beach,theta_beach, sigma_beach] = CMT2TT(Mbeach)
end
whos Mbeach misfit_beach
Mbeach_file = [Mbeach; misfit_beach'];

whos omega_beach misfit_beach kappa_beach theta_beach sigma_beach
post_beach = [omega_beach misfit_beach kappa_beach theta_beach sigma_beach];
post_beach = [sortrows(post_beach,2)]';  % sort by misfit for plotting beachballs

% Figure (c) from the summary plot
if iplot
    figure
    plot(omega_outline,misfit_outline)
    hold on
    plot(omega_subset,misfit_subset,'.g')
    xlabel('Omega');
    ylabel('Misfit');
end

%% Output for GMT

if igmt
    %gmtdir = '/home/vipul/gmt/data/cap/test/';
    %gmtdir = '/home/vipul/gmt/data/cap/test/';
    %plotdir = '~/rough/';
    %eid='20140416202423770';
    %plotdir=strcat(gmtdir,eid,'/',norm,'/',wts{ii},'/');
    mkdir(gmtdir);
    
    % Data for plotting histograms of strik-dip-rake
    data_gmt = [kappa_samples';theta_samples';sigma_samples'];
    fid = fopen(strcat(gmtdir,eid,'_histo.dat'),'w');
    fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',data_gmt);
    fclose(fid);
    % Data for plotting histograms of auxiliary planes of strik-dip-rake
    data_gmt = [kappa_samples2';theta_samples2';sigma_samples2'];
    fid = fopen(strcat(gmtdir,eid,'_histo2.dat'),'w');
    fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',data_gmt);
    fclose(fid);
    
    % omega-misfit outline (instead of using all 100000 sample points)
    fid = fopen(strcat(gmtdir,eid,'_omega_misfit_outline.dat'),'w');
    fprintf(fid,'%f\t%f\n',[omega_outline'; misfit_outline']);
    % fprintf(fid,'%f\t%f\n',[omega_outline'; ptry_outline']);
    fclose(fid);
    
    % file1[omega, misfit] and PTaxis for 400 posterior samples
    fid = fopen(strcat(gmtdir,eid,'_omega_err_post.dat'),'w');
    fprintf(fid,'%3.2f\t%3.3f\t%2.4f\t%2.4f\t%2.4f\t%2.4f\n',post_samples);
    fclose(fid);
    % save the 400 subsetted MTs and misfit
    fid = fopen(strcat(gmtdir,eid,'_Mpost_file.dat'),'w');
    fprintf(fid,'%2.4f\t%2.4f\t%2.4f\t%2.4f\t%2.4f\t%2.4f\t%3.4f\n',Mpost_file);
    fclose(fid);
    
    % [omega misfit strike dip rake omega] and misfit for 20 samples
    fid = fopen(strcat(gmtdir,eid,'_samples.dat'),'w');
    fprintf(fid,'%3.2f\t%3.3f\t%3.3f\t%3.2f\t%3.2f\n',post_beach);
    fclose(fid);
    % save the 20 subsetted MTs and misfit for plotting beachballs
    fid = fopen(strcat(gmtdir,eid,'_Mbeach_file.dat'),'w');
    fprintf(fid,'%2.4f\t%2.4f\t%2.4f\t%2.4f\t%2.4f\t%2.4f\t%3.4f\n',Mbeach_file);
    fclose(fid);
    
    % posterior pdf (not cdf)
    fid = fopen(strcat(gmtdir,eid,'_post_pdf.dat'),'w');
    fprintf(fid,'%3.4f\t%3.4f\n',post_pdf);
    fclose(fid);
    
    % plot mesa  pdf (not cdf)
    fid = fopen(strcat(gmtdir,eid,'_mesa_pdf.dat'),'w');
    fprintf(fid,'%3.4f\t%3.4f\n',mesa_pdf);
    fclose(fid);
    
    % posterior cdf
    fid = fopen(strcat(gmtdir,eid,'_post_cdf.dat'),'w');
    fprintf(fid,'%3.4f\t%3.4f\n',post_cdf);
    fclose(fid);
    
    % plot mesa cdf
    fid = fopen(strcat(gmtdir,eid,'_mesa_cdf.dat'),'w');
    fprintf(fid,'%3.4f\t%3.4f\n',mesa_cdf);
    fclose(fid);
    
    % plot uncertainty curve P(V)
    fid = fopen(strcat(gmtdir,eid,'_vP.dat'),'w');
    fprintf(fid,'%3.4f\t%3.4f\n',vP_vec);
    fclose(fid);
    fid = fopen(strcat(gmtdir,eid,'_vP_avg.dat'),'w');
    fprintf(fid,'%3.4f\t%3.4f\n',vP_vec_avg);
    fclose(fid);
    
    % save the main result (P, v(P), OMEGA)
    fid = fopen(strcat(gmtdir,eid,'_result.dat'),'w');
    fprintf(fid,'%1.3f\n',unc_result);
    fclose(fid);
    
    % Ouput the Mref points (for Noatak event only)
    if (iMref)
        Mref_set = [num2cell(Mref_set);Mtag];
        fid = fopen(strcat(gmtdir,eid,'_Mref_pts.dat'),'w');
        fprintf(fid,'%f\t%f\t%s\n',Mref_set{:});
        fclose(fid);
    end
end

    function [PV,indx] = volumetric_unc(Phat,Vhat,V_vec)
        % Returns PV on basis of cumulative probability Phat (cdf), homogeneous
        % sample distribution, Vhat. For each 'prob' in P, it returns the
        % corresponding volume (i.e. homogeneous sample distribution) of V
        % Phat = cumulative cdf (blue curve in the original doc)
        % Vhat = Homogeneous cdf (red curve in the original doc)
        % V_vec = V (see figure 4(b) of SilwalTape2016)
     
        for ii = 1:length(V_vec)
            Vi = V_vec(ii);
            % Find the Vhat which is closest to the Vi
            % V_vec and Vhat may not the same point
            [ival,indx(ii)] = min(abs(Vhat - Vi));
            PV(ii) = Phat(indx(ii));     % Key: Find PV
        end
        
        % Plotting
        if 0
            figure
            plot(prob,v);
            hold on
            plot(prob,v,'ro');
                xlabel('V');
                ylabel('PV');
        end
    

%% EXAMPLES
if 0    
    % SilwalTape2016 example
    % check 'smodel' and 'nsamples' above 
    %eid = '20090407201255351'; % AEC catalog id
    eid = '20090407201253480';  % IRIS catalog id
    ddir = '/home/vipul/CAP/inv/scak/20090407201253480_save/OUTPUT_DIR_1/';    
    dep = 41;
    smodel = 'scak';      % velocity model
    %gmtdir = '/home/vipul/CAP/inv/scak/20090407201253480_save/20090407201253480_gmt_data/';
    gmtdir = []; 
    K = 40/1.8003;
    X = 0;   % no polarity
    nsamples = 50000;     % Number of samples (random)
    CAP_unc_binary(eid,ddir,smodel,dep,nsamples,K,X,gmtdir);
    % Now get cross-section at minimum misfit (this is only for GMT plots)
    % NOT for uncertainty estimate
    nsamples = 49248;     % Number of samples (grid)
    MTbrick_Csection_binary(eid,ddir,smodel,dep,nsamples,K,X,gmtdir)
 
    % Kodiak offshore earthquake Aftershock
    eid = '20180131200145000';
    ddir = '/home/vipul/CAP/inv/scak/Kodiak_offshore/aftershocks/20180131200145000/OUTPUT_DIR/';    
    dep = 10;
    nsamples = 49248;     % Number of samples (random)
    smodel = 'scak';      % velocity model
    gmtdir = '/home/vipul/CAP/inv/scak/Kodiak_offshore/aftershocks/20180131200145000/20180131200145000_gmt_data/';
    %gmtdir = [];
    K = 20;
    X = 0;
    nsamples = 50000;
    CAP_unc_binary(eid,ddir,smodel,dep,nsamples,K,X,gmtdir);
    % Now get cross-section at minimum misfit (this is only for GMT plots)
    % NOT for uncertainty estimate
    nsamples = 49248;     % Number of samples (grid)
    MTbrick_Csection_binary(eid,ddir,smodel,dep,nsamples,K,X,gmtdir);

    
    % Iniskin Aftershock
    eid='20160124123742054';
    nsamples = 50000;     % Number of samples
    smodel = 'scak';     % velocity model
    ddir = '/home/vipul/PROJECTS/josh/20160124123742054/OUTPUT_DIR/';
    dep = 107;
    K = 40/1.8003;
    X = 0.2;
    gmtdir = '/home/vipul/PROJECTS/josh/20160124123742054/gmt_data/';
    CAP_unc_binary(eid,ddir,smodel,dep,nsamples,K,X,gmtdir);
    
    %% OBSOLETE EXAMPLES
    % 20090407201255351 for Celso's thesis
    % NOTE: Change the gmtdir
    eid = '20090407201255351';  % IRIS catalog id
    %ddir = '/home/vipul/CAP/inv/scak/';
    ddir = '/home/vipul/celso_thesis/';
    dep = 39;
    nsamples = 50000;     % Number of samples (grid)
    smodel = 'scak';        % velocity model
    K = 1;
    gmtdir = strcat(ddir,eid,'/gmt_data/');
    gmtdir = []; 
    CAP_unc_binary(eid,ddir,smodel,dep,nsamples,gmtdir,K);
    nsamples = 49248;
    MTbrick_Csection_binary(eid,ddir,smodel,dep,nsamples,gmtdir,K);
    
    % 20090407201255351 for testing new CAP
    % NOTE: Change the gmtdir
    eid = '20070911234634153';  % IRIS catalog id
    ddir = '/home/vipul/CAP/inv/scak/MOOS_updated/waveforms_alaska';
    dep = 94;
    nsamples = 50000;     % Number of samples (grid)
    smodel = 'scak';        % velocity model
    K = 1;
    gmtdir = strcat(ddir,'/',eid,'/gmt_data/');
    gmtdir = []; 
    CAP_unc_binary(eid,ddir,smodel,dep,nsamples,gmtdir,K);
    nsamples = 49248;
    MTbrick_Csection_binary(eid,ddir,smodel,dep,nsamples,gmtdir,K);
    
    %====================================================================
    % Beluga events (NEHRP)
    eid = '20080126042942584'; dep = 11; X = 0.3; K = 100;
    eid = '20120629110739385'; dep = 8; X = 0.5; K = 40;
    ddir = strcat('/home/vipul/CAP/inv/scak/NEHRP/beluga/',eid,'/OUTPUT_DIR/');
    eid = '20080418041458669'; dep = 3; X = 0.1; K = 10;
    ddir = strcat('/home/vipul/CAP/inv/scak/NEHRP/north_susitna/M_greater_than_3/',eid,'/OUTPUT_DIR/');
    nsamples = 50000;     % Number of samples (grid)
    smodel = 'scak';        % velocity model
    gmtdir = strcat(ddir,'../gmt_data/');
    %gmtdir = []; 
    CAP_unc_binary(eid,ddir,smodel,dep,nsamples,K,X,gmtdir);
    nsamples = 51840;
    MTbrick_Csection_binary(eid,ddir,smodel,dep,nsamples,K,X,gmtdir);
    %====================================================================
    
    % Little_Skull_Main for Celso's thesis
    % NOTE: Change the gmtdir
    eid = 'Little_Skull_Main';
    ddir = '/home/vipul/celso_thesis/';
    dep = 10;
    nsamples = 50000;     % Number of samples (random)
    smodel = 'wes';        % velocity model
    K = 20;
    gmtdir = strcat(ddir,eid,'/gmt_data/');
    %gmtdir = []; 
    CAP_unc_binary(eid,ddir,smodel,dep,nsamples,gmtdir,K);
    nsamples = 48618;     % Number of samples (grid)
    MTbrick_Csection_binary(eid,ddir,smodel,dep,nsamples,gmtdir,K);

    % AlvizuriTape2016 example
    % check 'smodel' and 'nsamples' above 
    eid='20100516063454464';
    nsamples = 1311975;     % Number of samples
    smodel = 'utuhalf';     % velocity model
    ddir = '/home/vipul/CAP/inv/utu/';
    dep = 4;
    CAP_unc_binary(eid,ddir,smodel,dep,nsamples);

    % OBSOLETE EXAMPLES
    % CAP log_random file should be present
    eid = '20090407201255351';
    %ddir = '~/CAP/inv/scak/';
    %ddir = strcat(ddir,eid,'_output_dir/');
    ddir = '~/CAP/inv/scak/20090407201255351_output_dir/'
    %gmtdir = '/home/vipul/gmt/data/cap/testbinary/';
    gmtdir = '~/CAP/inv/scak/20090407201255351_gmt_data/';
    CAP_unc_binary(eid,ddir,gmtdir);
    % Just input the log file
    %     CAP_unc('/home/vipul/CAP/inv/scak/MOOS/RESULTS/20090407201255351/L1/M111/log_039_rand');
    %     CAP_unc('/home/vipul/CAP/inv/scak/MOOS/20090407201255351/log_039_rand_ext');
    

end

