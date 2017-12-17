function [] = MTbrick_Csection(eid,ddir,gmtdir,Mref)
% script to plot the cross section of strike-dip-rake space
% INPUT: log file for grid search (May be just use the brick elements)
%
% TASK THIS SCRIPT SHOULD PEFORM:
% (1) Cross-sections of strike-dip-rake brick (Input: kappa0,theta0,sigma0 (or Mi))
% (2) Contours of omega (from kappa0,theta0,sigma0) on the cross-sections
% See run_CAP_unc.m for example
% Also see CAP_unc.m
%
% Vipul Silwal
% 27 July, 2015

%clear all, % close all

% MAKE FILENAME AS INPUT PARAMETER
%mat_filename = '/home/vipul/CAP/inv/scak/MOOS/20070911234634153/log_100_grid.mat';
%mat_filename = '/home/vipul/CAP/inv/scak/MISC/20140416202423770/log_072_grid.mat';

% dont create files for GMT plotting if gmt data dir is not provided
if nargin<3
    igmt = 0;
    Mref=[];
elseif nargin==3
    igmt = 1;
    Mref=[];
end

% Get the grid_log file from the data directory and change into mat file
if (nargin==1)
    grid_log_filename = eid;
else
    log_file = dir(strcat(ddir,'*grid'));
    if isempty(log_file)
        warning('CAP grid search log file in %s do not exist',ddir);
        return
    end
    grid_log_filename = strcat(ddir,log_file.name);
end

% Create .mat if it does not exist already
if ~exist(strcat(grid_log_filename,'.mat'),'file')
    grid_log_filename = ascii2mat(grid_log_filename);
else
    grid_log_filename = strcat(grid_log_filename,'.mat');
end

% Constants
deg = 180/pi;
misfit_fact = 40;         % misfit factor
%P = 0.8;        % Proability constants (See Walt's notes)
%ismooth = 0;    % to smooth the cdf
%imesa = 0;      % Use analytical formulation (1); estimate pdf from the sample distribution (0)
domega = 1/deg;
iTT2CMT = 0;
iminus = 0;

% Read log file
[kappa,theta,sigma,F,Mw,delta,gamma,M0,M] = read_CAP_log(grid_log_filename);
n = length(kappa);

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
whos kappa_vec theta_vec sigma_vec
% find global minimum (Mref)
[min_misfit,minindx] = min(F);
% indx = 1; % For random checking

% cross-sections at these values of strike, dip and rake
if isempty(Mref) % Mref not specified
    kappa0 = kappa(minindx);
    sigma0 = sigma(minindx);
    theta0 = theta(minindx);
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

% change misfit (F) from vector to mesh format
[ERR] = vec2mesh(F,kappa_vec,theta_vec,sigma_vec);
ERR = ERR*misfit_fact;   % scale by misfit factor
% Get the probability
prob = exp(-ERR);
% Omega in mesh fromat
[OMEGA] = vec2mesh(omega,kappa_vec,theta_vec,sigma_vec);

% Plot Cross-Sections and Omega contours
if iminus     % Plot cross-sections at completely opposite solution
    [kappa1,theta1,sigma1] =dcfaultpar2minusM(kappa0,theta0,sigma0); % 180 (opposite fault plane)
end
icont = 0; % plot contour
Csection_variable = ERR;
[Z1,Z2,Z3] = plot_MTbrick_Csection(Csection_variable,kappa_vec,theta_vec,sigma_vec,kappa0,theta0,sigma0,icont);

Csection_variable = OMEGA;
[C1,C2,C3] = plot_MTbrick_Csection(Csection_variable,kappa_vec,theta_vec,sigma_vec,kappa0,theta0,sigma0,1);

%% script to save (Z's and C's) input file for gmt

if igmt
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
    %copyfile(strcat(ddir,outfile.name),strcat(gmtdir,eid,'.out'));
    copyfile(strcat(ddir,eid,'.out'),gmtdir);
    % copy event and station location file to GMT plotting directory
    copyfile(strcat(ddir,eid,'_event.dat'),gmtdir);
    copyfile(strcat(ddir,eid,'_station.dat'),gmtdir);
    
    
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

%----------------------------------------------------
    function [C] = vec2mesh(c,x,y,z)
        % SCRIPT SHOULD
        % Input: c vector with x,y,z edges of brick
        % Output C in mesh format
        %
        % Vipul Silwal
        % July 30, 2015
        %
        % Future Work: WHAT IF, length(c) ~= length(x)*length(y)*length(z)
        
        % get length of edges
        x_len = length(x);
        y_len = length(y);
        z_len = length(z);
        
        % precreate the mesh
        C = zeros(y_len,x_len,z_len);
        % storing elements in mesh format
        for ii=1:z_len
            for jj=1:y_len
                for kk=1:x_len
                    C(jj,kk,ii)=c(x_len*y_len*(ii-1)+x_len*(jj-1)+kk);
                end
            end
        end
    end

    function  [Z1,Z2,Z3] = plot_MTbrick_Csection(plot_variable,kappa_vec,theta_vec,sigma_vec,kappa0,theta0,sigma0,icontour)
        %
        % Given the brick (kappa-theta-sigma) elements, plots the cross-sections
        % SCRIPT SHOULD
        % (1) Input: kappa-sigma-theta brick in mesh format;
        %            [kappa0, theta0, sigma0] where cross-sectiosn needs to be made
        %            [kappa_vec,theta_vec,sigma_vec] edges of brick
        % (2) Output: Plots of cross-sections, contours of omega on the Csections
        % (3) Return: 2D cross-sections files
        %
        % Vipul Silwal
        % July 30, 2015
        %
        % Future Work: Fix the plotting axis
        
        iplot  = 1;
        if nargin == 7
            icontour = 0;
        end
        % get the cross-sections
        stk_sec=squeeze(plot_variable(:,(kappa_vec==kappa0),:));
        dip_sec=squeeze(plot_variable((theta_vec==theta0),:,:));
        rak_sec=squeeze(plot_variable(:,:,(sigma_vec==sigma0)));
        
        if iplot
            figure
            subplot(2,2,1)
            if icontour
                contour(stk_sec)
            else
                imagesc(stk_sec)
            end
            subplot(2,2,2)
            if icontour
                contour(dip_sec)
            else
                imagesc(dip_sec)
            end
            subplot(2,2,3)
            if icontour
                contour(rak_sec)
            else
                imagesc(rak_sec)
            end
        end
        
        % return in vector for gmt input
        [x1,y1] = meshgrid(sigma_vec,theta_vec);
        Z1 = [x1(:) y1(:) stk_sec(:)];
        [x2,y2] = meshgrid(sigma_vec,kappa_vec);
        Z2 = [x2(:) y2(:) dip_sec(:)];
        [x3,y3] = meshgrid(kappa_vec,theta_vec);
        Z3 = [x3(:) y3(:) rak_sec(:)];
    end

%% EXAMPLES

if 0
    MTbrick_Csection('/home/vipul/CAP/inv/scak/MOOS/RESULTS/20090407201255351/L1/M111/log_039_grid');
    %--------------------------------------
    eid = '20090407201255351';
    ddir = '/home/vipul/CAP/inv/scak/MOOS/20090407201255351/stn1/';
    gmtdir = '/home/vipul/gmt/data/cap/test_stn/';
    MTbrick_Csection(eid,ddir,gmtdir);
end
end
