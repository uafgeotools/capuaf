function [w,otime,slon,slat,depc,Mw,eidlist] = get_rs(oran_vec,iregion,itest,irs)
% Contains important parameter for extracting waveform. However most of
% them stay unused of the time. 4 most frquently encountered options
% included
% oran      =   origin time, (single element or a vector)
% iregion   =   1 (SCAK); 2(MFSZ); 3 (interior); 4 (NEHRP) - OPTIONAL
% itest     =   1 test mode (only Z component of database 1) - OPTIONAL
% irs       =   1 plotting record section - OPTIONAL
%
% See get_rs_short.m for more examples
% 
% Vipul Silwal
% 3/27/2014

%clear all
%clc
%close all

spdy = 86400;
w= []; otime=[]; slon=[]; slat = []; dep=[];Mw = []; eidlist = [];

% If only eid is given - Search all AEC eq catalog and let it be a test run
if nargin==1
    iregion=0; itest=1; irs=1;
end


% CHOOSE region (or Add more) 
% This is to avoid the case in which multiple events occur at a same time
% at different place
switch iregion
    case 0
        % All Alaska (most of it!)
        ax3 = [-170 -130 50 70 -10 700];
        Mwran = [2 10];
    case 1
        % SCAK
        axmoos = [-154 -146 58 62.5];
        ax3 = [axmoos -10 700];
        Mwran = [2 10];
    case 2
        % MFSZ
        axmfz = [-151 -147.5 63.5 65.5];
        ax3 = [axmfz -10 700];
        Mwran = [2 10];
    case 3
        % interior (qyu)
        axint = [-153.5 -145.5 62.2 66.2];
        ax3 = [axint -10 700];
        Mwran = [2 10];
    case 4
        % NEHRP
        ax3 = [-153.5 -149 59 62 -10 200];
        Mwran = [2 10];
    case 5 % AKTOMO 
        ax3 = [-154 -139 59 66 -10 200];
        Mwran = [4 10];
end

%========================================================================
% Choose Stations (epicentral range and channels)
if itest
        % Smaller subset for quick checking
    stasub = [0 200]; idatabase = [1]; % Change as needed
    %chan = {'SHZ','SHE','SHN','SH1','SH2','SH*'};
    %chan = {'BHZ','EHZ','SHZ'};
    chan = {'BHZ'};
    iprocess = 1; sacdir = './'; 
else    % Save waveforms
    stasub = [0 300];
    idatabase = [1 2 3 4];      % UAF(1) BEAAR(2) ARCTIC(3) MOOSE(4) YAHTSE(5)
    %chan = {'BHZ','BHE','BHN','BH1','BH2','SHZ','SHE','SHN','SH1','SH2','SH*','EHZ','EHE','EHN','EH1','EH2'};
    %chan = {'BH*','EH*','SH*'};
    chan = {'BH*'};
    iprocess = 2;                % iprocess = 2 to deconvolve
    sacdir = './';
end

% Set the constants and plotting parameters here
% Usually you will have to change it only once for a project (for
% consistency)
% See getwaveform_input.m for info about the parameters 
isort = 2;          % =1 by azimuth, =2 by distance
iabs = 0;
T1 = [];
T2 = [];
trshift = 0;
tmark = [];
pmax = 50;
iintp = 0;
inorm = 1;
nfac = 1;
azcen = [];
iunit = 1;
tlims = [-40 200];       	% time limits for plotting
imap = 1;

% signal processing parameters
samplerate = 50;
cutoff = [];
iint = 0;                   % 1 to integrate waveform to displacement
duration_s = 300;
oshift = 100;

%========================================================================
% In case you want to extract multiple waveforms 
% Event will be extracted based on origin time vector and chosen region
for kk=1:length(oran_vec)
    oran=oran_vec(kk);
    % get earthquake information from the database
    oran = [oran 1] % Find events within +-1 sec window of oran_vec(kk)
    [otime,slon,slat,depc,Mw,eidlist] = read_eq_AEC(oran,ax3,Mwran)
    eidlist
    
    % This may seem like a redundant loop - just in case there are multiple
    % events in 1 sec window (see read_eq_AEC.m) 
    for ii=1:length(otime)
        close all; whos
        originTime = otime(ii);
        elat = slat(ii);
        elon = slon(ii);
        %dep = sdep(ii);
        mag = Mw(ii);
        eid = eidlist(ii);
        edep_km = depc(ii);
        
        % Set time limits for data extraction
        tshift = oshift + trshift;
        startTime = originTime - oshift/spdy;
        endTime   = originTime + duration_s/spdy;
        dur_dy = endTime-startTime;
        fprintf('origin time is %s\n',datestr(originTime,'yyyy-mm-dd HH:MM:SS.FFF'))
        fprintf('startTime is %s\n',datestr(startTime,31));
        fprintf('total length of time requested: %.2f s (= %.2f min = %.2f hours)\n',...
            dur_dy*spdy,dur_dy*3600,dur_dy*24);
        
        
        % DATA EXTRACTION
        tic
        [w,s,site,sitechan] = getwaveform(idatabase,startTime,endTime,chan,iint,...
            iprocess,cutoff,samplerate,stasub,sacdir,originTime,elat,elon,edep_km,mag,eid);
        toc
        disp(sprintf('%.1f s to execute getwaveform.m from getwaveform_cap.m',toc));
        
        % PLOT record section
        if and(irs==1,~isempty(w))
            % assume all waveforms have the same event ID
            keid = get(w(1),'KEVNM');
            
            % plot record section
            plotw_rs(w,isort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azcen,iunit,imap);
            % note: set imap = 0 in plotw_rs.m to NOT plot the map (and print the
            %       record section with the following command)
            iprint = 0;
            if iprint==1, orient tall; print(gcf,'-dpsc',sprintf('%s_%i_%i_isort%i_rs',get(w(1),'KEVNM'),round(T1),round(T2),isort)); end
            %orient tall; print(gcf,'-dpsc','sumatra2012hf_fullA');
        end
        
        % Post processing (Rotate to RTZ)
        % CAUTION: only one of these commands at a time
        %          Running with different channels overwrites older weight
        %          files
        if 0
            rt_cap('./ProcessedData2',eid{1},{'BHZ'});
            %rt_cap('./ProcessedData2',eid{1},{'BHZ','SHZ','EHZ'});
        end
        
        % write kmlfiles
        iwritekml=0;
        if iwritekml
            kmlfname=strcat(eid{1},'_kmlfile');
            stnm=get(w,'station');
            stla=get(w,'STLA');
            stlo=get(w,'STLO');
            kmlwrite(kmlfname,stla,stlo,'Name',stnm);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE
if 1==0
    iregion = 1;    % 1 = SACK; 2 = MFSZ
    itest  = 1;     % in test mode it takes only one component and one database
    irs = 1;        % plot record section
    oran = [datenum(2007,9,19,11,22,26)];
       
    [w,otime,slon,slat,depc,Mw,eidlist] = get_rs(oran,iregion,itest,irs);
    %=============================================
    iregion = 5;    % 1 = SACK; 2 = MFSZ
    itest  = 1;     % in test mode it takes only one component and one database
    irs = 1;        % plot record section
    oran = [datenum(2007,09,11,23,46,34)];
       
    [w,otime,slon,slat,depc,Mw,eidlist] = get_rs(oran,iregion,itest,irs);
    %=============================================
    iregion = 4;    % 1 = SACK; 2 = MFSZ
    itest  = 0;     % in test mode it takes only one component and one database
    irs = 1;        % plot record section
    oran = [datenum(2012,03,08,10,57,43)];
    oran = [datenum(2014,01,24,12,07,03)];  
    oran = [datenum(2012,03,06,06,12,58)];
       
    [w,otime,slon,slat,depc,Mw,eidlist] = get_rs(oran,iregion,itest,irs);
    
    oran = [datenum(2008,02,05,03,51,42)];
    [w,otime,slon,slat,depc,Mw,eidlist] = get_rs(oran)
end
%======================================================================
