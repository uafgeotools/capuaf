function [] = add_picks2weight(weight_input,weight_output,idb,iUseReviewed,HWIN_SEC)
% Adds first-motion polarity and P arrival times to the weight file

%idb = 1;            % Use AEC db
iPonly = 1;         % Get only P arrival time
%HWIN_SEC = 2;       % this is really important because otime is different for AEC and IRIS
%iUseReviewed = 0;   % Analyst reviewed or unreviewed polarity
spdy = 86400;

% Read the weight file
weight_input
if ~exist(weight_input,'file'), error('file does not exist'); end
[stn_tag, pol, edist, PV_wt, PR_wt, SV_wt, SR_wt, ST_wt, P_arrival, P_len, S_arrival, S_len, waveform_shft] = read_cap_weight(weight_input);

% Get the origin time and add the time range (HWIN_SEC) to match with the AEC otime
tmp = strsplit(stn_tag{1},'.');
eid = tmp{1};
otime = eid2otime(eid);
otar = [otime HWIN_SEC];

% Get polarity and time picks
[sta_aec,chan,timepick,iphase,fm,ifm,delta_deg,rlon,rlat,relev,elon,elat,edep,eml,originTime,dist_hyp] = get_picks(otar,idb,iPonly,iUseReviewed);

% KEY: this time shift must be relative to the target origin time, which
% might differ from the origin time that is in the database; the reason is
% that our sac files have the origin time set as otar
p_time_aec = (timepick-originTime)*spdy;

% Difference between AEC origin time and IRIS origin time
odiff = (originTime-otime)*spdy;

% add the difference between otimes' to the P arrival time w.r.t.  AEC
% origin time
p_time_iris = p_time_aec + odiff;

% Match station name and add polarity (ifm) and P arrival time (timepick)
for ii = 1:length(stn_tag)
    tmp = strsplit(stn_tag{ii},'.');
    stnm_iris(ii) = tmp(3);
    for jj = 1:length(sta_aec)
        if strcmp(stnm_iris{ii},sta_aec{jj})
            pol_aec(ii) = ifm(jj);
            P_arrival_aec(ii) = p_time_iris(jj);
        end
    end
end

% write the weight file
write_cap_weight(weight_output, stn_tag, pol_aec, edist, PV_wt,PR_wt, SV_wt, SR_wt, ST_wt, P_arrival_aec, P_len, S_arrival, S_len, waveform_shft)

% Examples

if 0
    weight_input = '/home/vipul/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/20090407201253480/weight.dat';
    weight_output = '/home/vipul/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/20090407201253480/weight_test.dat';
    idb = 1;            % Use AEC db
    HWIN_SEC = 2;       % this is really important because otime is different for AEC and IRIS
    iUseReviewed = 0;   % Analyst reviewed or unreviewed polarity
    add_picks2weight(weight_input,weight_output,idb,iUseReviewed,HWIN_SEC)
    %-----------------------------------------------------------------
%     weight_input = '/home/ksmith/REPOSITORIES/capuaf/Nov20MT/20151120105348168/weight_body.dat';
%     weight_output = '/home/ksmith/REPOSITORIES/capuaf/Nov20MT/20151120105348168/weight_body_picks2.dat';
%     idb = 1;            % Use AEC db
%     HWIN_SEC = 2;       % this is really important because otime is different for AEC and IRIS
%     iUseReviewed = 1;   % Analyst reviewed or unreviewed polarity
%     add_picks2weight(weight_input,weight_output,idb,iUseReviewed,HWIN_SEC)
end
