function [Mref, post_samples_info,Pav,PV,PV_diff]= CAP_unc(eid,ddir,gmtdir,Mref)
% script for doing posterior analysis and uncertainty estimate on MT
% Input : Log files from CAP (contains MT info and misfit values)
%   eid  - event id
%   ddir - directory where the log files are stored
%   Mref -  Reference moment tensor
%   gmtdir - directory where you want to save the file fot GMT plotting
%
% Vipul Silwal
% 20 July, 2015

% clear all,% close all

% FILENAME SHOULD BE AN INPUT PARAMETER

%mat_filename = '/home/vipul/CAP/inv/scak/MOOS/20070911234634153/log_100_rand.mat';
%mat_filename = '/home/vipul/CAP/inv/scak/MISC/20140416202423770/log_072_rand.mat';

% dont create files for GMT plotting if gmt data dir is not provided
% if nargin~=4
%     igmt = 0;
% else
%     igmt = 1;
% end
%
% if nargin==1
%     Mref=[];
% end

%
% TO DO:
% Move the p(m0) block from mtvipul to here

if nargin<3
    igmt = 0;
    Mref=[];
elseif nargin==3
    igmt = 1;
    Mref=[];
end

% If there is only one input parameter it must be the log file
if (nargin==1)
    rand_log_filename = eid;
else
    % read the rand_log file from the data directory and change into mat file
    log_file = dir(strcat(ddir,'*rand'));
    if isempty(log_file)
        warning('CAP random search log file in %s do not exist',ddir);
        return
    end
    rand_log_filename = strcat(ddir,log_file.name);
end

% Create .mat if it does not exist already
if ~exist(strcat(rand_log_filename,'.mat'),'file')
    rand_log_filename = ascii2mat(rand_log_filename);
else
    rand_log_filename = strcat(rand_log_filename,'.mat');
end

% Constants (or Inputs varianbles)
deg = 180/pi;
K = 40;         % misfit factor
P = 0.8;        % Probability constants (See Walt's notes)
ismooth = 0;    % to smooth the cdf
imesa = 1;      % (1) Use analytical formulation; (0)estimate pdf from the sample distribution
idc = 1;        % use pdf for omega double couple or full moment tensor
domega = 1/(deg);
iTT2CMT = 0;
Nsamp_reject = 2000;
iplot = 0;
ipost_samples_file = 0;
ipost_beach_file = 0;
Vmax = 1; % parameter to govern what fraction of V(w) we want to consider for Pav estimation
%Vmax = 0.1;
%Vmax = 0.05;

% Read files and Get the moment tensors
[kappa,theta,sigma,F,Mw,delta,gamma,M0,M] = read_CAP_log(rand_log_filename);
n = length(kappa);

if iTT2CMT      % Another Option (takes longer time)
    gamma = 0;
    delta = 0;
    M0 = 1/sqrt(2);     % such that |Lambda| = 1
    M = TT2CMT(gamma,delta,M0,kappa,theta,sigma);
end

% find global minimum (Mref)
[ival,indx] = min(F);
% indx = 1;
kappa0 = kappa(indx);
sigma0 = sigma(indx);
theta0 = theta(indx);
Mmin = M(:,indx)

% You can choose any Mref (reference moment tensor) !
if isempty(Mref)
    Mref = Mmin; % default option
    iMref = 0; % to change the position of Mref (default: Mref = Mmin)
else
    iMref = 1; % to change the position of Mref (default: Mref = Mmin)
end

% This section is for NOATAK event only
if (iMref)
    % might want to keep beachballs and posterior samples same
    % ipost_samples_file = 1;
    % ipost_beach_file = 1;
    M1 = Mmin;
    M2 = M(:,3); % randomly chosen - because M is random anyways         % Omega 80
    %%M3 = [-0.1592   -0.1531    0.3123   -0.7805    0.5409   -0.1585]'; % Omega 90
    %M4 = [-0.2618   .5315    -0.2697 0.6372 -0.6145 -0.0671]';          % Omega 72 (local minima)
    %M5 = [0.0250   -0.6454    0.6205    0.0437    0.1540    0.7571]';   % Omega 175
    %M6 = [-4.620 3.740 0.876 -0.880 0.982 -1.730 ]'; % GCMT solution    % Omega 62
    
    Mref = M2; % CHANGE THIS: M1 M2 M3 M4 M5
    %Mtag = {'M1','M2','M3','M4','M5'};
    Mtag = {'M1','M2'};
    
    % Find these moment tensors in the full MT set (to get their omega and
    % misfit values - used for plotting red points in CAP_summary_plot.pl
    [Mset(:,1),~,set_indx(:,1)] = omega2Mset(0,M1,M);
    [Mset(:,2),~,set_indx(:,2)] = omega2Mset(0,M2,M);
    %     [Mset(:,3),~,set_indx(:,3)] = omega2Mset(0,M3,M);
    %     [Mset(:,4),~,set_indx(:,4)] = omega2Mset(0,M4,M);
    %     [Mset(:,5),~,set_indx(:,5)] = omega2Mset(0,M5,M);
    
    omega_set = CMT2omega(Mref,Mset);
    misfit_set = F(set_indx)*K;
    Mref_set = [omega_set misfit_set]';
end

% omega angle from Mref to all M
%Mref = [-0.9765    0.1361    0.8404    0.1423    0.1425    0.3468]';
omega_deg = CMT2omega(M,Mref);
omega = omega_deg/deg;

%
iMark = 1;
if iMark
    Mark = [-0.9765    0.1361    0.8404    0.1423    0.1425    0.3468]';
    omega_mark = CMT2omega(Mark,Mref);
    omega_mark = omega_mark/deg
end

%% probability density function
%Kvec = 40:10:40; % test for multiple scale factors
Kvec = 40;
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
        [N,mesa,omega_mid] = plot_histo(omega,omega_mid(1)-domega/2:domega:omega_mid(end)+domega/2,3);
        omega_mid=omega_mid';
    end
    
    %     if idc ==1  % for double couple
    %         mesa_eq = omegadcpdf(omega_mid,false);
    %         else        % for full moment tensor
    %         mesa_eq = omegapdf(omega_mid,false);
    %     end
    mesa = mesa/(sum(mesa)*domega);
    homo_prob = cumsum(mesa)*domega;
    mesa_pdf = [(omega_mid*deg)' mesa']';
    mesa_cdf = [(omega_mid*deg)' homo_prob']';
    
    % Compute volume v(P) for Mref (M0 or best solution)
    [ival,indx] = min(abs(cum_prob - P));
    OMEGA = omega_mid(indx);
    vP = homo_prob(indx);
    %unc_result = [P vP OMEGA*deg]';
    
    % Get the PV curve (see SilwalTape2016 Figure 4)
    % TO DO: PV_vec and PV_vec_end cann be reduced to one
    Vvec = [0:0.01:1];
    [PV,vindx] = volumetric_unc(cum_prob,homo_prob,Vvec);
    PV_vec = [Vvec' PV']';
    end_points = [1 0;0 0];
    PV_vec_end = [PV_vec end_points]; % with end points (why is it needed?)
    PV_diff = diff(PV)./diff(Vvec); % derivative of PV
    cum_prob_diff = diff(cum_prob);
    homo_prob_diff = diff(homo_prob);
    PV_diff_theory = cum_prob_diff'./homo_prob_diff; % Sanity check: This should be same as PV_diff
    whos cum_prob_diff homo_prob_diff PV_diff omega_mid PV_diff_theory
    
    % Compute uncertainty Pav
    % cum_probArguement: Pav is too conservative (Vmax is the knob to adjust the
    % samples taken for estiamting Pav)
    V_total = find(Vvec==Vmax);
    Pav = sum(PV(1:V_total))/V_total;
    
    % Plotting
    % Data from these figures will go to subfigure(b) of the Summary plots
    % figure(10);
    %----------------------------------------------------
    % TO DO We perhaps wants to add option to save as pdf (input flag isave)
    subplot(2,2,1)
    if ismooth
        prob1 = smooth(prob1,20);
    end
    plot(omega_mid,prob1,'b','LineWidth',2);
    hold on
    plot(omega_mid,mesa,'r','LineWidth',2);
    %plot(omega_mid,mesa_eq,'-k','LineWidth',2);
    if iMark
        plot(omega_mark,0.1,'v','MarkerFaceColor','m','Markersize',7);
    end
    xlim([omega_mid(1)-domega/2 omega_mid(end)+domega/2]);
    axis equal
    axis([0 pi 0 2])
    title(sprintf('PDF, p = exp(-K x misfit); K = %3.0f', K))
    %----------------------------------------------------
    subplot(2,2,2)
    if ismooth
        prob1 = smooth(prob1,10);
    end
    plot(omega_mid,smooth(prob1,20),'b--','LineWidth',2);
    hold on
    plot(omega_mid,mesa,'r','LineWidth',2);
    xlim([omega_mid(1)-domega/2 omega_mid(end)+domega/2]);
    axis equal
    axis([0 pi 0 2])
    title(sprintf('smoothed PDF, p = exp(-K x misfit); K = %3.0f', K))
    %----------------------------------------------------
    subplot(2,2,3)
    plot(omega_mid,cum_prob,'b','LineWidth',2);
    hold on
    plot(omega_mid,homo_prob,'r','LineWidth',2);
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
    plot(Vvec,PV,'k-','LineWidth',2);
    hold on
    %plot([P P],[0 v(Pvec==P)],'k--');
    %plot([0 P],[v(Pvec==P) v(Pvec==P)],'k--');
    plot([0 1],[0 1],'k--');
    %     PV_x = linspace(0,1,length(PV_diff_theory));
    %     plot(PV_x,PV_diff_theory,'r--');
    %PV_diff = 20*(PV_diff/max(PV_diff));
    PV_x = linspace(0,1,length(PV_diff));
    PV_diff = smooth(PV_diff,5);
    plot(PV_x,smooth(PV_diff,5),'g-','LineWidth',2);
    %axis([0 1 0 1]);
    xlabel('V');
    ylabel('P(V)');
    %title(sprintf('P = %2.2f, v(P) = %1.3f',P,v(Pvec==P)));
    title(sprintf('P_{av} = %1.3f',Pav));
    axis
end

if (0)
    % Option: Another intutive way for estimating uncertainity
    % Omega at which the difference between the cdf (blue) and the mesa (red) is the maximum
    figure
    plot(omega_mid,cum_prob-homo_prob)
    xlabel('omega');
    title('cdf - mesa');
end

%% Check with the Rejection method (always keep number of posterior samples
% = 2000)
% (UPDATE)
%
chance = rand(n,1);
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
    figure
    subplot(2,2,1)
    plot_histo(kappa(ikeep),0:5:360);
    hold on
    plot(kappa0,0,'vr','MarkerSize',5)
    xlabel('Strike')
    %----------------------------------------------------
    subplot(2,2,2)
    plot_histo(theta(ikeep),0:5:90);
    hold on
    plot(theta0,0,'vr','MarkerSize',5)
    xlabel('Dip')
    %----------------------------------------------------
    subplot(2,2,3)
    plot_histo(sigma(ikeep),-90:5:90);
    hold on
    plot(sigma0,0,'vr','MarkerSize',5)
    xlabel('Rake')
end

% Auxiliary fault plane (for the same set of posterior moment tensor samples)
kappa_samples = kappa(ikeep);
theta_samples = theta(ikeep);
sigma_samples = sigma(ikeep);
[kappa_samples2,theta_samples2,sigma_samples2]=dcfaultpar2aux(kappa_samples,theta_samples,sigma_samples);

%========== get outline points for omega1 and rnd_err=========
nbin=60;
[omega_outline,misfit_outline,M_outline] = cappts2outline(omega_deg,misfit,M,nbin);
if iplot
    figure
    subplot(2,1,1)
    plot(omega_outline,misfit_outline)
    xlabel('Omega');
    ylabel('Misfit');
end

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
    
    % posterior pdf (not cdf)
    fid = fopen(strcat(gmtdir,eid,'_post_cdf.dat'),'w');
    fprintf(fid,'%3.4f\t%3.4f\n',post_cdf);
    fclose(fid);
    
    % plot mesa  pdf (not cdf)
    fid = fopen(strcat(gmtdir,eid,'_mesa_cdf.dat'),'w');
    fprintf(fid,'%3.4f\t%3.4f\n',mesa_cdf);
    fclose(fid);
    
    % plot mesa  pdf (not cdf)
    fid = fopen(strcat(gmtdir,eid,'_vP.dat'),'w');
    fprintf(fid,'%3.4f\t%3.4f\n',PV_vec);
    fclose(fid);
    fid = fopen(strcat(gmtdir,eid,'_vP_avg.dat'),'w');
    fprintf(fid,'%3.4f\t%3.4f\n',PV_vec_end);
    fclose(fid);
    
    % save the main result (P, v(P), OMEGA)
    fid = fopen(strcat(gmtdir,eid,'_result.dat'),'w');
    fprintf(fid,'%1.3f\n',Pav);
    fclose(fid);
    
    % Ouput the Mref points (for Noatak event only)
    if (iMref)
        Mref_set = [num2cell(Mref_set);Mtag];
        fid = fopen(strcat(gmtdir,eid,'_Mref_pts.dat'),'w');
        fprintf(fid,'%f\t%f\t%s\n',Mref_set{:});
        fclose(fid);
    end
end

    function [v,indx] = volumetric_unc(P,V,prob)
        % Returns v(P) on bases of in cumulative probability P (pdf), homogeneous
        % sample distribution, V. For each 'prob' in P, it returns the
        % corresponding volume (i.e. homogeneous sample distribution) of V
        % P = cumulative pdf (blue curve in the original doc)
        % V = Homogeneous pdf (red curve in the original doc)
        %
        
        PV = 1;
        
        if PV
            temp = P;
            P = V;
            V = temp;
        end
        
        for ii = 1:length(prob)
            p = prob(ii);
            [ival,indx(ii)] = min(abs(P - p));
            v(ii) = V(indx(ii));
        end
        if 0
            figure
            plot(prob,v);
            hold on
            plot(prob,v,'ro');
            if PV
                xlabel('V');
                ylabel('PV');
            else
                xlabel('P');
                ylabel('v(P)');
            end
        end
    end

%% EXAMPLES
if 0
    % CAP log_random file should be present
    eid = '20070911234634153';
    ddir = '/home/vipul/CAP/inv/scak/MOOS/RESULTS/20070911234634153/L1/M111/';
    gmtdir = '/home/vipul/gmt/data/cap/test/';
    CAP_unc(eid,ddir,gmtdir);
    %--------------------------------------
    eid = '20090407201255351';
    ddir = '/home/vipul/CAP/inv/scak/MOOS/20090407201255351/stn1/';
    gmtdir = '/home/vipul/gmt/data/cap/test_stn/';
    CAP_unc(eid,ddir,gmtdir);
    
    %--------------------------------------
    % Just input the log file
    CAP_unc('/home/vipul/CAP/inv/scak/MOOS/RESULTS/20090407201255351/L1/M111/log_039_rand');
    CAP_unc('/home/vipul/CAP/inv/scak/MOOS/20090407201255351/log_039_rand_ext');
    %=========================================================
end
end
