% This script performs 2 tasks:
% 1. Gets events from the catalog based on the selection criteria
% 2. Extracts the selected events for CAP inversion 
% 
% This script is built upon getwaveform.m (See getwaveform.m;
% getwaveform_input.m)
%
% README:
% Select a region ('iregion' - SCAK or MFFZ) and selection critera
% Perform a test run (itest) first for a couple of events. Select couple of
% events and extract the waveform. 
% This script will get_rs.m. get_rs.m has been adapted to generate waveforms 
% required for CAP inversion :
% 1. Rotate to RTZ (radial-transverse-vertical component) from NEZ (north-east-vertical)
% 2. resample to 50sps
% 3. create weight files needed for running CAP

% Vipul Silwal 
% May 30, 2015

%%%%% KEY parameter %%%%%%
iregion = 2;    % 1 = SCAK (south-secntal Alaska) ; 2 = MFSZ (interior Alaska)
itest  = 1;     % for test (extract raw waveforms; see get_rs.m); Chose 1 for waveforms for CAP
irs = 0;        % plot record sections or not
imoos = 0;      % 21 main events from SilwalTape2016
iscak = 0;      % 85 smaller events from SilwalTape2016
iint = 1;       % for running in interior alaska (qyu and Josh)
iBeluga = 1;

% Region bounds
axmoos = [-154 -146 58 62.5];
axint = [-153.5 -145.5 62.2 66.2]; % qyu and Josh
axmfz = [-151 -147.5 63.5 65.5];   % vipul

%====Step 1: Get events that you want to extract============
switch iregion
    case 1
        % MOOSE
        oran = [datenum(2007,09,19,11,22,26)];
        %oran = [datenum(2007,09,11,23,46,34) 1];
        
        % GCMT
        % oran = [datenum(2009,01,02,14,15,36) 1];
        
        % SCAK region
        oran = [datenum(2008,06,02,17,27,40)];
        oran = [datenum(2009,06,17,03,48,49)];
        oran = [datenum(2014,05,05,04,59,4)];
        oran = [datenum(2008,11,18,19,56,51)];
        oran = [datenum(2014,04,16,20,24,24)];
        oran = [datenum(2013,08,01,21,32,47)];
        
        % subset example (MOOS) - 21 events in SilwalTape2016 (Moment
        % tensor solution from AEC)
        if imoos
            otime = [datenum(2007,8,15) datenum(2009,8,15)];
            ax3 = [axmoos -10 700];
            Mwran = [0 10];
            [oran,slat,slon,sdep,M,M0,Mw,eid] = read_mech_AEC(otime,ax3,Mwran);
            [~,isort] = sort(Mw,'ascend');
            %display_eq_list(isort,oran,slon,slat,sdep,Mw,eid);
        end
        
        % 85 events in SilwalTape2016 (fault-plane solution from AEC)
        if iscak
            otime = [datenum(2007,8,15) datenum(2009,8,15)];
            ax3 = [axmoos -10 700];
            Mwran = [3.5 10];
            [oran,slat,slon,sdep,M,M0,Mw,eid] = read_mech_AECfp(otime,ax3,Mwran);
            %[~,isort] = sort(eid);
            [~,isort] = sort(sdep,'ascend');
        end
        
        % small events from studying seismic clustering
        % events obtained from get_events_NEHRP.m)
        if iBeluga 
            oran = [datenum(2008,02,05,03,51,42),datenum(2010,03,28,16,05,36),datenum(2009,05,16,01,51,04),datenum(2012,06,29,11,07,39)...
                ,datenum(2014,01,24,12,07,03),datenum(2008,01,26,04,29,42),datenum(2014,07,14,06,04,10),datenum(2012,03,06,06,12,58)];
        end
        
    case 2
        % MFSZ region
        % EAST
        oran = [datenum(2001,03,25,11,34,50)];
        oran = [datenum(2002,12,29,20,38,30)];    % only 7 stations
        oran = [datenum(2002,08,13,06,23,52)];    % only 7 stations
        oran = [datenum(2013,03,05,21,55,58)];
        oran = [datenum(2001,05,02,23,53,15)];
        % WEST
        oran = [datenum(2000,02,07,04,17,53)];    % Seems like clustered events       
        oran = [datenum(2008,07,16,10,12,00)];
        oran = [datenum(2012,04,11,09,21,57)];
        oran = [datenum(2009,07,28,12,13,15)];
        % NORTH
        oran = [datenum(2013,07,12,07,59,17)];
        oran = [datenum(2008,07,20,08,09,37)];
        % SOUTH
        oran = [datenum(2000,11,29,10,35,47)];
        oran = [datenum(2000,12,06,18,40,26)];
        % SOUTHWEST
        oran = [datenum(2013,06,05,18,58,23)];
        oran = [datenum(2010,08,29,18,08,06)];
        if 0
            % USED in MFFZ paper
            oran = [datenum(2000,11,29,10,35,47) datenum(2000,12,06,18,40,26) datenum(2001,03,25,11,34,50)...
                datenum(2001,06,30,09,41,42)...
                datenum(2008,07,16,10,12,00) datenum(2009,07,28,12,13,15) datenum(2012,04,11,09,21,57)...
                datenum(2013,03,05,21,55,58) datenum(2013,06,05,18,58,23) datenum(2013,07,12,07,59,17)...
                datenum(2014,12,13,15,47,31)...
                datenum(2014,08,31,03,06,57) datenum(2014,08,31,12,24,58) datenum(2014,10,21,00,36,58) datenum(2014,10,23,16,30,24)
                ];
        end
        if iint
            if 0            % Events for Qingping's project
                otime = [datenum(2009,1,1) datenum(2011,06,1)]; % Qingping
                ax3 = [axint -10 800];
                Mwran = [0 10];
                [oran,slat,slon,sdep,M,M0,Mw,eid] = read_mech_AEC(otime,ax3,Mwran);
                [~,isort] = sort(oran,'ascend');
                display_eq_list(isort,oran,slon,slat,sdep,Mw,eid);
            end
            
            % Events for Josh's project
            % Option 1: Interior Alaska event
            otime = [datenum(2011,06,1) datenum(2015,06,1)];
            ax3 = [axint -10 800];
            Mwran = [4 10];     % events between magnitude 4 and 10
            [oran,slon,slat,sdep,Mw,eid,etype] = read_eq_AEC(otime,ax3,Mwran);
            [~,isort] = sort(oran,'ascend');
            display_eq_list(isort,oran,slon,slat,sdep,Mw,eid);
        end
end

%====Step 2: Extract selected eevnts============
% See get_rs.m
oran = datenum(2009,4,7,20,12,55);
[w,otime,slon,slat,depc,Mw,eidlist] = get_rs(oran,iregion,itest,irs);