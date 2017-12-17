function rt_cap(data_dir,eid,chanl)
%RT_CAP adapts waveforms from getwaveform.m for CAP moment tensor inversion
%
% Operations done:
%	Changing right handed system to left handed (by multiplying t-axis by-1)
% 	Rotation to RTZ (left-handed system)
%	Changing amplitude from nm/s to cm/s (*e-7)
%	Adding headers fields: CMPAZ (changed after rotation)
%
% Input:
%	data_dir  Folder containg the raw waveforms (e.g., 'ProcessedData2/')
%   eid       event ID (string)
%   chanl     Z channel code (e.g., 'BHZ', 'HHZ', ...)
%
% Output (written as sac files):
%	station_name.r	station_name.t	station_name.z
%
% Example:
%   rt_cap('ddir_name','event_name',{'chans'})
%   rt_cap('./ProcessedData1/','20090407201255351',{'BHZ'})
%   rt_cap('/home/carltape/PROJECTS/nenana/','20120411091837_HL',{'BHZ'})
%
% calls rt_rot.m
%

% format for weight files
%stfmt = '%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n';
stfmt = '%12s %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n';

% if nargin == 2
%     cha = 'B';
% end

if isnumeric(eid), eid = num2str(eid); end
evtdir = strcat(data_dir,'/',eid,'/');
if ~exist(evtdir,'dir'), error([evtdir ' does not exist']); end

if nargin==2
    disp('input channel is not specified -- set cha to wildcard');
    cha = '*';
end

mkdir(strcat(evtdir,'/weight/'));
evtdir1 = strcat(evtdir,'/weight/');

% generate template weight files for CAP (st_name  dist w1 w2 w3 w4 w5 t1 t2)
% Example : BMR 282 0 0 1 1 1 41.0 0
% CAP WEIGHT FILES WITH STATIONS SORTED BY AZIMUTH
fid1 = fopen(strcat(evtdir1,'w1_az.dat'),'w');
fid2 = fopen(strcat(evtdir1,'weight_full_az.dat'),'w');
fid3 = fopen(strcat(evtdir1,'weight_body_az.dat'),'w');
fid4 = fopen(strcat(evtdir1,'weight_surf_az.dat'),'w');
fid5 = fopen(strcat(evtdir1,'weight_near_az.dat'),'w');
fid6 = fopen(strcat(evtdir1,'weight_far_az.dat'),'w');
fid7 = fopen(strcat(evtdir,eid,'_','station.dat'),'w');
fid8 = fopen(strcat(evtdir,eid,'_','event.dat'),'w');

count = 0;

dist=[];stla=[];stlo=[];evla=[];evlo=[];az=[];evdp=[];kevnm=[];
IN=[]; IE=[];

ii_carryon = 0;     % to keep track of different channels
tmp = 1;            % to keep track of stations with same name (ex: POKR_TA and POKR_TA1)
for kk=1:length(chanl)
    Zcha = chanl{kk};
    chatag = Zcha(1:2);
    % Get all vertical one - example = *.BHZ.sac
    sacfiles = dir(strcat(evtdir,'*','.','*',Zcha,'.sac'))
    % Number of *Z.sac files. This way we might be missing files which only
    % has N or E (rare but possible - example if sitechan info is not
    % available for Z)
    nseis = length(sacfiles);
    in = zeros(nseis,1);
    ie = zeros(nseis,1);
    
    % Loop over the number of stations (one can do this with using the input
    % file for stations and use dir command instead
    for ii_stn = 1:nseis
        %disp(sprintf('%s', s{ii}));
        fname=textscan(sacfiles(ii_stn).name,'%s','delimiter','.');
        a=fname{1};
        ii = ii_stn + ii_carryon;    % ii_carryon is for number of stations of another sensor location
        ntwk{ii} = a{2};    % network id...... these variables are used for to picking correct set of records
        s{ii} = a{3};       % station name
        loc{ii} = a{4};     % sensor location tag

        % New fromat for naming sacfiles
        ftag{ii} = strcat(eid,'.',ntwk{ii},'.',s{ii},'.',loc{ii},'.',chatag);
        fname_start = strcat(evtdir,ftag{ii})
        fname_end = 'sac';
        
        %===================LOAD WAVEFORMS=================================
        % north
        % there are multiple ways in which station components are described
        % examples: BHN, BH1, BHN_01. If you encounter other ways modify
        % underneath
        
        sacevt = dir(strcat(fname_start,'N.',fname_end));
        sacevt1 = dir(strcat(fname_start,'1.',fname_end));  % these are just 2 ways
        if (isempty(sacevt) && isempty(sacevt1))
            warning('could not find N component file! - creating null N file ');
            % if only Z component is present create N and E
            sacevt = dir(strcat(fname_start,'Z.',fname_end));
            in(ii_stn) = 1;
        elseif isempty(sacevt)
            sacevt = sacevt1;
        end
        filename = strcat(evtdir,sacevt(1).name);
        n = rsac(filename);
        % Change the headers so that N and E can pass the QC before rotation (in
        % rt_rot)
        if in(ii_stn)
            n = ch(n,'CMPINC',90);
            n = ch(n,'CMPAZ',0);
            n(:,2) = n(:,2)*0;      % Null North comp file; amplitde = 0
        end
        
        % east
        sacevt = dir(strcat(fname_start,'E.',fname_end));
        sacevt1 = dir(strcat(fname_start,'2.',fname_end));
        if (isempty(sacevt) && isempty(sacevt1))
            warning('could not find E component file! - creating null E file ');
            % if only Z component is present create N and E
            sacevt = dir(strcat(fname_start,'Z.',fname_end));
            ie(ii_stn) = 1;
        elseif isempty(sacevt)
            sacevt = sacevt1;
        end
        filename = strcat(evtdir,sacevt(1).name);
        e = rsac(filename);
        % Change the headers so that N and E can pass the QC before rotation (in
        % rt_rot)
        if ie(ii_stn)
            e = ch(e,'CMPINC',90);
            e = ch(e,'CMPAZ',90);
            e(:,2) = e(:,2)*0;      % Null East comp file; amplitde = 0
        end
        
        % vertical
        sacevt = dir(strcat(fname_start,'Z.',fname_end));
        filename = strcat(evtdir,sacevt(1).name);
        z = rsac(filename);
        
        %========================ROTATION==================================
        
        disp(sprintf('Before Rotation		%f(N)	%f(E)	%f(BAZ)',lh(n,'CMPAZ'),lh(e,'CMPAZ'),lh(n,'BAZ')));
        
        % get sample rate from channel
        chan = lh(z,'KCMPNM');
        rcha = chan(1);
        
        % rotation (KEY COMMAND)
        [r,t] = rt_rot(n,e,rcha);
        if r==0
            continue
        end
        
        disp(sprintf('After Rotation		%f(R)	%f(T)   %f(AZ)',lh(r,'CMPAZ'),lh(t,'CMPAZ'),lh(r,'AZ')));
        
        %=======================CHANGE HEADERS=============================
        % modification for CAP (nm/s to cm/sec)
        r(:,2) = r(:,2)*1e-7;
        t(:,2) = t(:,2)*1e-7;
        z(:,2) = z(:,2)*1e-7;
        % change header from VELOCITY (NM/SEC) to UNKNOWN to avoid confusion.
        % Note that 5=UNKNOWN but this is only viewable using lh from inside
        % sac itself; you cannot see it from 'saclst', for example.
        r = ch(r,'IDEP',5);
        t = ch(t,'IDEP',5);
        z = ch(z,'IDEP',5);
        
        % changing headers
        %r = ch(r,'LPSPOL',1,'CMPINC',90);
        %t = ch(t,'LPSPOL',1,'CMPINC',90);
        %z = ch(z,'LPSPOL',1,'CMPINC',0,'CMPAZ',0);
        
        % getting header info
        % note: here we get the headers from the radial component file (r)
        knet(ii,:) = lh(r,'KNETWK');
        dist(ii,:) = lh(r,'DIST');
        stla  = lh(r,'STLA');
        stlo  = lh(r,'STLO');
        evla  = lh(r,'EVLA');
        evlo  = lh(r,'EVLO');
        az(ii,:)    = lh(r,'AZ');
        evdp  = lh(r,'EVDP');
        kevnm = lh(r,'KEVNM');
        
        %=======================SAVING RTZ=================================
        % saving RTZ sac files
        %wsac(strcat(evtdir,s{ii},'_',knet(ii,:),'.z'),z);
        %wsac(strcat(evtdir,s{ii},'_',knet(ii,:),'.r'),r);
        %wsac(strcat(evtdir,s{ii},'_',knet(ii,:),'.t'),t);
        wsac(strcat(fname_start,'.z'),z);
        wsac(strcat(fname_start,'.r'),r);
        wsac(strcat(fname_start,'.t'),t);
       
        % making staions and event file
        fprintf(fid7, '%s\t%s\t%f\t%f\t%f\t%f\n',s{ii},knet(ii,:),stla,stlo,dist(ii),az(ii));
        
        count = count + 1;
        
        disp('------------------------------------------------------------------');
    end
    ii_carryon = ii;
    IN = [IN;in];
    IE = [IE;ie];
end

fprintf(fid8, '%s\t%f\t%f\t%f\n',kevnm,evla,evlo,evdp);

%=====================CREATING WEIGHT FILES========================
% FUTURE: MAKE SCRIPT FOR THIS
% creating weight file for CAP
% Azimuth sorted weight files
whos
length(az)
[azs, iaz] = sort(az);

for jj = 1:length(azs)
    if azs(jj)~=0
        %stnm = strcat(s{iaz(jj)},'_',knet(iaz(jj),:));
        stnm = ftag{iaz(jj)};
        fprintf(fid1,stfmt,stnm,0,0,0,0,0,0,0,0,0,0,0);
        %-------------------------------------------------------------
        if ~and(IN(iaz(jj)),IE(iaz(jj)))
            fprintf(fid2,stfmt,stnm,0,1,1,1,1,1,0,0,0,0,0);
        else
            fprintf(fid2,stfmt,stnm,0,1,0,1,0,0,0,0,0,0,0);
        end
        %-------------------------------------------------------------
        if ~and(IN(iaz(jj)),IE(iaz(jj)))
            fprintf(fid3,stfmt,stnm,0,1,1,0,0,0,0,0,0,0,0);
        else
            fprintf(fid3,stfmt,stnm,0,1,0,0,0,0,0,0,0,0,0);
        end
        %-------------------------------------------------------------
        if ~and(IN(iaz(jj)),IE(iaz(jj)))
            fprintf(fid4,stfmt,stnm,0,0,0,1,1,1,0,0,0,0,0);
        else
            fprintf(fid4,stfmt,stnm,0,0,0,1,0,0,0,0,0,0,0);
        end
        %-------------------------------------------------------------
        if dist(jj) < 200       % Near (stations less than 200km away) -only body
            if ~and(IN(iaz(jj)),IE(iaz(jj)))
                fprintf(fid5,stfmt,stnm,0,1,1,0,0,0,0,0,0,0,0);
            else
                fprintf(fid5,stfmt,stnm,0,1,0,0,0,0,0,0,0,0,0);
            end
        else                % Far distance 
            if ~and(IN(iaz(jj)),IE(iaz(jj)))
                fprintf(fid6,stfmt,stnm,0,1,1,1,1,1,0,0,0,0,0);
            else
                fprintf(fid6,stfmt,stnm,0,1,0,1,0,0,0,0,0,0,0);
            end
        end
    end
end

count = count + 1;
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
fclose(fid6);
fclose(fid7);
fclose(fid8);

% ==========sort by dist================
% CAP WEIGHT FILES WITH STATIONS SORTED BY DISTANCE
fid1 = fopen(strcat(evtdir1,'w1_dist.dat'),'w');
fid2 = fopen(strcat(evtdir1,'weight_full_dist.dat'),'w');
fid3 = fopen(strcat(evtdir1,'weight_body_dist.dat'),'w');
fid4 = fopen(strcat(evtdir1,'weight_surf_dist.dat'),'w');
fid5 = fopen(strcat(evtdir1,'weight_near_dist.dat'),'w');
fid6 = fopen(strcat(evtdir1,'weight_far_dist.dat'),'w');

length(dist)
[dists, idist] = sort(dist);

if isempty(dist)
    return
    disp('No weight files generate');
end

for jj = 1:length(dists)
    if dists(jj)~=0
        %stnm = strcat(s{idist(jj)},'_',knet(idist(jj),:));
        stnm = ftag{idist(jj)};
        fprintf(fid1,stfmt,stnm,0,0,0,0,0,0,0,0,0,0,0);
        %-------------------------------------------------------------
        if ~and(IN(idist(jj)),IE(idist(jj)))
            fprintf(fid2,stfmt,stnm,dists(jj),1,1,1,1,1,0,0,0,0,0);
        else
            fprintf(fid2,stfmt,stnm,dists(jj),1,0,1,0,0,0,0,0,0,0); % Only Vertical component
        end
        %-------------------------------------------------------------
        if ~and(IN(idist(jj)),IE(idist(jj)))
            fprintf(fid3,stfmt,stnm,0,1,1,0,0,0,0,0,0,0,0);
        else
            fprintf(fid3,stfmt,stnm,0,1,0,0,0,0,0,0,0,0,0); % Only Vertical component
        end
        %-------------------------------------------------------------
        if ~and(IN(idist(jj)),IE(idist(jj)))
            fprintf(fid4,stfmt,stnm,0,0,0,1,1,1,0,0,0,0,0);
        else
            fprintf(fid4,stfmt,stnm,0,0,0,1,0,0,0,0,0,0,0); % only Surf V
        end
        %-------------------------------------------------------------
        if dist(jj) < 200       % Near (stations less than 200km away) -only body
            if ~and(IN(idist(jj)),IE(idist(jj)))
                fprintf(fid5,stfmt,stnm,0,1,1,0,0,0,0,0,0,0,0);
            else
                fprintf(fid5,stfmt,stnm,0,1,0,0,0,0,0,0,0,0,0);
            end
        else                    % Stations farther than 200km
            if ~and(IN(idist(jj)),IE(idist(jj)))
                fprintf(fid6,stfmt,stnm,0,1,1,1,1,1,0,0,0,0,0);
            else
                fprintf(fid6,stfmt,stnm,0,1,0,1,0,0,0,0,0,0,0);
            end
        end
        %-------------------------------------------------------------
    end
end

fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
fclose(fid6);

disp(sprintf('\n Total %d x3 files generated',count-1));
name = strcat(evtdir,s{ii},'.z');

% Creates weight files with first motion information (first motion and P arrival time)
% stored in AEIC database
idb = 1;                % set database id
iUseReviewed = 1;       % use reviewed polarity info
picks2weight(data_dir,eid,idb,iUseReviewed)
end