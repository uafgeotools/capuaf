function picks2weight(evtdir,eid,idb,iUseReviewed)
% Function to insert the P arrival times and first motion into the weight
% file for CAP inversion. 
% Input 
%   evtdir      parent directory containing the directory of event data
%   eid         event id
%   idb         database index
% Output directories
%   weight_pick     include P arrival time
%   weight_fm       include first motion information
%   weight_pick_fm  include all P wave information (both first motion and
%                   arrival time 
%  
% Example:
% evtdir = '/home/vipul/GEOTOOLS/matlab_util/util_data_syn/ProcessedData2/'
% eid = '20070911234634000';
% idb = 1;
% picks2weight(evtdir,eid,idb)
% 
% CAUTION !  get_picks may give multiple first motions for different
% channels. Generally all of them are same. If they are not same any of
% them will be chosen at random. We need a better way of deciding! 
% 
%
% Vipul Silwal
% July 17, 2015

spdy = 86400;
iPonly = 1;
iolddb = 0;     % use old database paths (They still work - as of 03/30/2016)

% call get_picks
otar = eid2otime(eid);

if iolddb
    [sta,chan,timepick,fm,ifm,delta_deg,rlon,rlat,relev,elon,elat,edep,eml,originTime,dist_hyp] = get_picks_OLD(otar,idb,iPonly,iUseReviewed);
else
    [sta,chan,timepick,iphase,fm,ifm,delta_deg,rlon,rlat,relev,elon,elat,edep,eml,originTime,dist_hyp] = get_picks(otar,idb,iPonly,iUseReviewed);
end


% if using plutons database then also read from external database of
% polarities and arrival times
if idb==6
    %hwin_sec = 1    % half window of time interval to search
    %get_picks([otar hwin_sec], idb, 1)
    [x_eid,x_otime,x_evid,x_orid,x_elon,x_elat,x_edep_km,x_emag,x_sta,x_stlo,x_stla,x_stel,x_staz,x_tdist,x_vdist,x_hdist,x_Ppol,x_Ptime] = read_utu_fmp;
    imatch = find(strcmp(eid,x_eid)==1);
    if ~isempty(imatch)
        eid_xdb = x_eid(imatch);
        stn_xdb = x_sta(imatch);
        pol_xdb = x_Ppol(imatch);
        Ptime_xdb = x_Ptime(imatch);
        Ptime_xdb(isnan(Ptime_xdb)) = 0;
    else
       disp('WARNING: no eid match in PLUTONS fmp database'); 
    end
end

% KEY: this time shift must be relative to the target origin time, which
% might differ from the origin time that is in the database; the reason is
% that our sac files have the origin time set as otar
p_time = (timepick-originTime)*spdy;

% directory containing weight files
wtdir='weight';
evtdir1 = strcat(evtdir,'/',eid,'/',wtdir,'/');
filename_vec={'w1_az.dat';'weight_body_az.dat';    'weight_far_az.dat';   'weight_full_az.dat';    'weight_near_az.dat';    'weight_surf_az.dat';...
        'w1_dist.dat';  'weight_body_dist.dat'; 'weight_far_dist.dat';  'weight_full_dist.dat' ; 'weight_near_dist.dat'; 'weight_surf_dist.dat'};
    
% output directory
wtdir_new = {'weight_pick' ;'weight_fm';'weight_pick_fm'};

% note end of line
%stfmt = '%s\t%d\t%d\t%d\t%d\t%d\t%d\t%4.2f\t%d\t%d\t%d\t%d\n';
stfmt = '%12s %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %6.2f %4.0f %4.0f %4.0f %4.0f\n';

for kk=1:length(wtdir_new)
    
    wtdir2 = wtdir_new{kk};
    evtdir2 = strcat(evtdir,'/',eid,'/',wtdir2,'/');
    mkdir(evtdir2);
    %length(filename_vec)
    
    %filename = 'weight_full_dist.dat';
    for jj=1:length(filename_vec)
        % changes for each jj
        filename=filename_vec{jj};
        
        % get weight file for reading
        fid = fopen(strcat(evtdir1,filename));
        
        % open file for writing
        fid1 = fopen(strcat(evtdir2,filename),'w');
        
        %S = textscan(fid,'%s%s','delimiter','_');               % read station names from file
        S = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%f%f');
        fclose(fid);
        stnm_net = S{1};
        
        % probably this could be more efficient
        % get name of each station
        stnm = stnm_net;    % initialize
        for xx=1:length(stnm_net)
            stemp = stnm_net{xx};
            stnm_wt=textscan(stemp,'%s','delimiter','.');
            tmp=stnm_wt{1};
            ntwk = tmp{2};    % network id...... these variables are used for to picking correct set of records
            stnm(xx) = {tmp{3}};       % station name
            loc = tmp{4};     % sensor location tag
            chan_wt = strcat(tmp{5},'Z');
        end
        
        [istnm,ista]=ismember(stnm,sta);                        % compare station names from weight file with those from get_picks
        % istnm (0,1)
        % ista - index in sta of common elements
        
        fid = fopen(strcat(evtdir1,filename));
        W = textscan(fid,'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n');
        fclose(fid);
        dist=W{2};
        w1=W{3};
        w2=W{4};
        w3=W{5};
        w4=W{6};
        w5=W{7};
        w6=W{8};    % Ptime
        w7=W{9};
        w8=W{10};
        w9=W{11};
        w10=W{12};

        P_pick = zeros(length(W{1}),1);      % w6 = W{8}
        P_fm = zeros(length(W{1}),1);
        for ii=1:length(stnm)       % loop over stations
            
            % if the station is listed in the arrival table (get_picks.m)
            % (note: this does NOT mean that the ptime and fm fields are present)
            if istnm(ii)
                P_pick(ii) = p_time(ista(ii));                 % assign P_pick
                P_fm(ii) = ifm(ista(ii));                      % assign first motion
            end
            
            if idb==6
                % note: this would over-write P_fm from get_picks.m (if fm is present in Total)
                imatch = strcmp(stnm(ii),stn_xdb);
                if ~isempty(imatch)
                    P_fm(ii) = pol_xdb(imatch);
                    % if the station is not listed, then use our own pick
                    %if and(~istnm(ii),~P_fm(ii))
                    if ~istnm(ii)
                        P_pick(ii) = Ptime_xdb(imatch);
                        disp('WARNING: custom P time is being used');
                    end    
                end
            end
            
            stag = stnm_net{ii};
            if kk==1                                % weight_pick 
                fprintf(fid1,stfmt,stag,...
                    dist(ii),w1(ii),w2(ii),w3(ii),w4(ii),w5(ii),P_pick(ii),w7(ii),w8(ii),w9(ii),w10(ii));

            elseif kk==2                            % weight_fm
                if P_fm(ii)==0
                    w1(ii)=0; w2(ii)=0; w3(ii)=0; w4(ii)=0; w5(ii)=0; % don't use stations without the first-motion info
                else
                    stag = strcat(stnm_net{ii},'/',num2str(P_fm(ii)));
                end
                fprintf(fid1,stfmt,stag,...
                    dist(ii),w1(ii),w2(ii),w3(ii),w4(ii),w5(ii),w6(ii),w7(ii),w8(ii),w9(ii),w10(ii));
                
            else                                    % weight_pick_fm
                if P_fm(ii)==0
                    % w1(ii)=0; w2(ii)=0; w3(ii)=0; w4(ii)=0; w5(ii)=0; % don't use stations without the first-motion info
                else
                    stag = strcat(stnm_net{ii},'/',num2str(P_fm(ii)));
                end   
                fprintf(fid1,stfmt,stag,...
                    dist(ii),w1(ii),w2(ii),w3(ii),w4(ii),w5(ii),P_pick(ii),w7(ii),w8(ii),w9(ii),w10(ii));
            end
        end
        fclose(fid1);

    end  % end loop over files (jj)
end  % end loop over sets of files (kk)

