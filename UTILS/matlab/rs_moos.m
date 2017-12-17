%
% rs_moos.m
% generate record section plots of MOOS events
%
% To make composite PDF:
% for file in `ls *_best_rs.ps` ; do ps2pdf $file ; done ; pdcat -r *_rs.pdf rs_moos_best.pdf
% for file in `ls *_small_rs_isort1.ps` ; do ps2pdf $file ; done ; pdcat -r *_rs_isort1.pdf rs_moos_small_isort1.pdf
% for file in `ls *_small_rs_isort2.ps` ; do ps2pdf $file ; done ; pdcat -r *_rs_isort2.pdf rs_moos_small_isort2.pdf
% for file in `ls 319798*_rs.ps` ; do ps2pdf $file ; done ; pdcat -r 319798*_rs.pdf rs_319798.pdf
% 
%

clear
close all
clc

iprint = 0;
irs = 1;

ibest = 0;      % PROBABLY THE CODE FOR THIS NEEDS TO BE UPDATED
ismall = 1;

spdy = 86400;

%-----------------------

if ibest==1

% load AEIC MT catalog
[otime0,elat0,elon0,edep0,M,M0,Mw0,eid0,eidc,depc] = read_mech_AEIC;
[~, Mw] = CMT2all(M);

T1 = 10; T2 = 50; tshift = 140; duration = 140;
T1 = 5; T2 = 10; tshift = 30; duration = 140;
T1 = 2; T2 = 30; tshift = 30; duration = 140;

% SUBSET OF EVENTS
%ax0 = [-151.5 -147.5 59.5 63.0];
%ax0 = [-150.5461 -148.4908   59.9452   62.2171];
ax0 = [-156 -144 58 62.5];
ax3 = [ax0 -10 700];
oran = [datenum(2007,8,15) datenum(2009,8,15)];
Mwran = [0 10];

[otime0,elat0,elon0,edep0,Mw0,isubset] = get_seis_subset(otime0,elat0,elon0,edep0,Mw0,oran,ax3,Mwran);
eid0 = eid0(isubset);
eidc = eidc(isubset);
depc = depc(isubset);

[~,isort1] = sort(depc,'ascend');
display_eq_list(isort1,otime0,elon0,elat0,edep0,Mw0,eid0,depc,eidc);

% PLOT
rfile = '/home/carltape/gmt/stations/seismic/Matlab_output/STATIONS_ALASKA_moos_specfem';
read_station_specfem(rfile);
hold on; plot(elon0,elat0,'rp','markersize',12);

% bounding box for stations -- southern Alaska
%axb = [-156 -144 56 64];
%axb = [-154.5 -147.5 58.5 64.5];

% display_eq_list.m: listing 21 events
%   1   4 otime 2007-10-10 (2007283) 18:03:26 lon -147.41 lat 59.96 dep  15.00 ( 11.80) km Mw 4.05 280436 (280253)
%   2   2 otime 2007-09-19 (2007262) 11:22:26 lon -146.11 lat 61.38 dep  35.00 ( 30.82) km Mw 4.45 279084 (278884)
%   3  15 otime 2009-04-07 (2009097) 20:12:55 lon -149.74 lat 61.45 dep  45.00 ( 33.03) km Mw 4.59 319798 (319605)
%   4  12 otime 2009-02-15 (2009046) 19:35:00 lon -146.33 lat 61.60 dep  35.00 ( 37.24) km Mw 4.49 316444 (316248)
%   5   8 otime 2008-08-28 (2008241) 23:14:18 lon -149.60 lat 62.12 dep  50.00 ( 42.95) km Mw 4.11 305302 (305104)
%   6  21 otime 2009-07-30 (2009211) 22:39:10 lon -151.09 lat 59.93 dep  30.00 ( 44.08) km Mw 4.38 326766 (326568)
%   7   3 otime 2007-10-03 (2007276) 14:06:12 lon -151.29 lat 58.28 dep  40.00 ( 45.46) km Mw 5.17 279857 (279709)
%   8  17 otime 2009-04-30 (2009120) 04:54:57 lon -151.31 lat 58.99 dep  40.00 ( 52.73) km Mw 4.88 321298 (321087)
%   9  20 otime 2009-06-26 (2009177) 16:48:20 lon -150.64 lat 61.91 dep  60.00 ( 59.48) km Mw 4.25 324973 (324781)
%  10  19 otime 2009-06-22 (2009173) 19:28:05 lon -150.70 lat 61.94 dep  80.00 ( 64.59) km Mw 5.51 324684 (324487)
%  11   7 otime 2008-03-27 (2008087) 23:07:45 lon -152.17 lat 59.01 dep  75.00 ( 68.53) km Mw 5.26 290389 (290230)
%  12   5 otime 2007-11-28 (2007332) 23:57:03 lon -151.13 lat 61.91 dep  75.00 ( 69.61) km Mw 4.82 283325 (283131)
%  13  13 otime 2009-02-23 (2009054) 00:04:27 lon -153.63 lat 58.92 dep  85.00 ( 87.75) km Mw 4.86 316902 (316685)
%  14  10 otime 2008-12-28 (2008363) 07:13:10 lon -151.05 lat 62.35 dep  80.00 ( 89.31) km Mw 4.42 313450 (313240)
%  15  14 otime 2009-03-17 (2009076) 01:13:33 lon -152.15 lat 60.24 dep  90.00 ( 90.14) km Mw 4.22 318329 (318113)
%  16   9 otime 2008-09-18 (2008262) 19:43:53 lon -152.79 lat 59.50 dep  85.00 ( 90.15) km Mw 4.53 306659 (306480)
%  17  11 otime 2009-01-24 (2009024) 18:09:50 lon -152.89 lat 59.43 dep  95.00 ( 97.87) km Mw 5.72 314984 (314787)
%  18   1 otime 2007-09-11 (2007254) 23:46:34 lon -151.53 lat 61.53 dep 100.00 (100.88) km Mw 4.40 278678 (278485)
%  19  16 otime 2009-04-14 (2009104) 17:14:27 lon -153.06 lat 60.16 dep 105.00 (117.78) km Mw 4.27 320317 (320117)
%  20  18 otime 2009-05-24 (2009144) 09:40:04 lon -153.25 lat 59.78 dep 125.00 (125.48) km Mw 4.60 322995 (322785)
%  21   6 otime 2008-03-14 (2008074) 09:38:21 lon -152.64 lat 61.07 dep 165.00 (143.65) km Mw 4.97 289517 (289317)

% moos
%teid = eid0;
%teid = eidc;
%teid = {'278678' '279084' '280436' '283325' '289517' '290389' '305302' '306659' '313450' '314984' '316444' '318329' '319798' '320317' '322995' '324684' '324973' '326766'};
%teid = {'319798'};
%db = '/home/admin/databases/MOOS/wf/moos'; ds = datasource('antelope',db);

% AEIC
%teid = {'349057'}; ds = datasource('uaf_continuous'); db = '/aerun/sum/params/Stations/master_stations';

neid = length(eid0);
%neid = 1;  % testing

%neid = 1;
k1 = 1; k2 = neid;
%k1 = 4; k2 = k1;

for xx = k1:k2
   
    kk = isort1(xx);
    eid = eidc{kk};
    %eid = eid0{kk};
    %imatch = find(strcmp(teid{isort1(kk)},eid0)==1);
    elon = elon0(kk);
    elat = elat0(kk);
    %edep_km = edep0(kk);
    edep_km = depc(kk);
    mag = Mw(kk);
    otime = otime0(kk);
    
    originTime = otime;
    startTime = originTime - tshift/spdy;
    endTime   = originTime + duration/spdy;

    % bounding box containing stations
    %dlon = 6;
    %dlat = 4;
    %axb = [elon-dlon elon+dlon elat-dlat elat+dlat];
    %axb = [-156.1 -144 56 64];  % extend to include TTA
    axb = [0 500];
    
    chan= {'BHZ'};
    iraw = 0;
    iint = 0;
    sacdir = [];
    cutoff = [];
    samplerate = [];
    [w,s,site,sitechan] = getwaveform(0,startTime,endTime,chan,iint,iraw,cutoff,samplerate,axb,sacdir,originTime,elat,elon,edep_km,eid);

    % plot record section
    if irs==1
        nw = length(w);
        for ii = 1:nw
            rlat(ii) = site(ii).lat;
            rlon(ii) = site(ii).lon;
            relev_km(ii) = 0;
        end 

        % add fields
        w = wset(w,rlat,rlon,elat,elon,eid,edep_km,relev_km,mag,tshift,get(w,'KNETWK'));

        % INPUT PARAMETERS FOR RECORD SECTION
        pmax = 1000;
        isort = 2;  % =1 by azimuth, =2 by distance
        iabs = 0;
        iint = 0;
        inorm = 1;
        nfac = 1.5;
        azcen = [];
        iunit = 1;
        plotw_rs(w,isort,iabs,T1,T2,pmax,iint,inorm,[],nfac,azcen,iunit);
        if iprint==1, orient tall; print(gcf,'-dpsc',sprintf('%2.2i_%s_%i_%i_isort%i_best_rs',xx,eid,round(T1),round(T2),isort)); end
    end
end

% for kk=1:neid
% 
%     imatch = find(strcmp(teid{kk},eid0)==1);
%     elon = elon0(imatch);
%     elat = elat0(imatch);
%     edep_km = edep0(imatch);
%     mag = Mw(imatch);
%     otime = otime0(imatch);
% 
%     tstart = otime - tshift/spdy;
%     tend = otime + duration/spdy;
%     scnl=scnlobject('*','BHZ','','');
%     w = waveform(ds,scnl,tstart,tend);
% 
%     % collect event and station info
%     sta = get(w,'station');
%     nsta = length(sta);
%     for ii=1:nsta
%         site0 = db_get_site_info(sta(ii),get(w(ii),'start'),db);
%         rlat(ii) = site0.lat;
%         rlon(ii) = site0.lon;
%         relev_km(ii) = site0.elev;
%     end
% 
%     % add fields to waveform object
%     w = wset(w,rlat,rlon,elat,elon,teid{kk},edep_km,relev_km,mag,tshift);
% 
%     % KEY parameters: record sections
%     isort = 2;
%     inorm = 1; nfac = 1.0;
%     %inorm = 0; nfac = 1.5;
%     T1 = 2; T2 = 30;
%     %T1 = []; T2 = [];
%     pmax = 40;
%     iint = 0;
%     azcen = 0;
%     iunit = 1;
%     plotw_rs(w,isort,T1,T2,pmax,iint,inorm,[],nfac,azcen,iunit);
%     
%     if iprint==1, orient tall; print(gcf,'-dpsc',sprintf('%s_best_rs',teid{kk})); end
% end
   
end

%--------------------------------------------------------------------------

if ismall==1
    % open MOOS array to see bounding box for target events
    rfile = '/home/carltape/gmt/stations/seismic/Matlab_output/STATIONS_ALASKA_moos_specfem';
    [rlon,rlat,relev,rburial,stnm,netwk] = read_station_specfem(rfile);
    
    % SUBSET OF EVENTS
    %ax0 = [-151.5 -147.5 59.5 63.0];
    ax0 = [-150.5461 -148.4908   59.9452   62.2171];
    ax3 = [ax0 -10 30];
    oran = [datenum(2007,8,15) datenum(2009,8,15)];
    magran = [2.5 2.9];

    [otime0,elon0,elat0,edep0,mag0,eid0] = read_eq_AEC(oran,ax3,magran);
    [~,isort1] = sort(mag0,'descend');
    display_eq_list(isort1,otime0,elon0,elat0,edep0,mag0,eid0);
    
    % subset of events
    if 1==1
        %ax0 = [-150.5461 -148.4908   59.9452   62.2171];
        %ax3 = [ax0 -10 30];
        %oran = [datenum(2007,8,15) datenum(2009,8,15)];
        %magran = [2.5 2.9];
        %ipick = [32 10 5 20 9 28 17];
        ipick = [32 35 38 14 12 15 16];   % 2nd index in list
        
        otime0 = otime0(ipick);
        elon0 = elon0(ipick);
        elat0 = elat0(ipick);
        edep0 = edep0(ipick);
        eid0 = eid0(ipick);
        mag0 = mag0(ipick);
        isort1 = [1:length(ipick)]';
    end
    n = length(otime0);
    hold on;
    plot(elon0,elat0,'rp','markersize',12);
    text(elon0,elat0,eid0,'fontsize',8);
    
    display_eq_list(isort1,otime0,elon0,elat0,edep0,mag0,eid0);
    
%     % eids
%     sday = 1 + floor(otime0 - datenum(year(otime0),1,1));
%     eid0 = repmat(cellstr(''),1,n);
%     for ii=1:n
%         x = otime0(ii);
%         eid0(ii) = cellstr(sprintf('%4i%3.3i%2.2i%2.2i%2.2i',year(x),sday(ii),hour(x),minute(x),round(second(x))));
%     end
    
    tshift = 0; duration = 100;
    %tshift = 0; duration = 35;

    % moos
    db = '/home/admin/databases/MOOS/wf/moos';
    ds = datasource('antelope',db);

    %n = 4;
    
    keid = [];
    for jj=1:n
    %for jj=1:1
        %kk = isort1(jj);
        disp('------------------------------');
        kk = jj;
        eid = eid0{kk}
        elon = elon0(kk);
        elat = elat0(kk);
        edep_km = edep0(kk);
        mag = mag0(kk);
        otime = otime0(kk);

        % bounding box containing stations
        %dlon = 3;
        %dlat = 2;
        %axb = [elon-dlon elon+dlon elat-dlat elat+dlat]
        stasub = [0 300];
        
        tstart = otime - tshift/spdy;
        tend = otime + duration/spdy;
        
        chan = {'BHZ'};
        samplerate = [];
        sacdir = [];
        iint = 0;
        iprocess = 1;   % =2 for deconvolution
        cutoff = [1/100 8];
        
        %continue    % to skip the rest of the loop
        
        [w,s,site,sitechan] = getwaveform([1 4],tstart,tend,chan,iint,iprocess,...
            cutoff,samplerate,stasub,sacdir,otime,elat,elon,edep_km,mag,eid);
        
        if irs==1
            % save only if there is at least one station within DMIN km of event
            DMIN = 50;
            dists = get(w,'DIST');
            if min(dists) < DMIN
                disp(sprintf('%3i %3i %s',jj,kk,eid));
                % KEY parameters: record sections
                isort = 2;  % DISTANCE (=2) or AZIMUTH (=1)
                iabs = 0;
                tmark = [];
                inorm = 1; nfac = 1.5;
                %inorm = 0; nfac = 1.5;
                T1 = 0.5;
                T2 = 4;
                pmax = 80;
                iintp = 0;
                tlims = [];
                azcen = 0;
                iunit = 1;
                imap = 0;

                plotw_rs(w,isort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azcen,iunit,imap);

                if iprint==1, orient tall; print(gcf,'-dpsc',sprintf('%s_small_rs_isort%i',eid,isort)); end
            end
        end
    end

end

%==========================================================================

if 0==1
    % 316248
    startTime = 7.338198156261343e+05;
    endTime = 7.338198175937268e+05;
    datestr(startTime)
    scnl = scnlobject({'COLA','COR','KDAK'},'BHZ','','');
    ds = datasource('uaf_continuous');
    w = waveform(ds,scnl,startTime,endTime);
end

%==========================================================================