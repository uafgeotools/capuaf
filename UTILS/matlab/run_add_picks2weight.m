clear all, close all, clc

% Nenana basin study
% evid = '20181113152641907';
% evid = '20180916191154565';
% evid = '20190925134513809';
% evid = '20171230114316278';
% evid = '20160711200557702';
% evid = '20170613073936181';
% evid = '20190117121355727';
% evid = '20190624090423301';
% evid = '20190624090423301';
% evid = '20170527163305640';
% evid = '20180825181551595';
% evid = '20170429111548898';
% evid = '20190326212718519';
% evid = '20190306213313991';
evid = '20190113164555437';

ftag0 = 'weight';
% ftag0 = 'weight_orig';
%ftag0 = 'WEIGHT_CLEAN';

%'/home/ksmith/REPOSITORIES/capuaf/MTs/Apr29MT/V1/';
cap_path = '/home/ksmith/REPOSITORIES/capuaf/';
%cap_path = '/home/ksmith/REPOSITORIES/capuaf/Mar06MT/V1/';

case_no = 2;
switch case_no 
    case 1, ftag = '';
    case 2, ftag = '_body';
    case 3, ftag = '_surf';
    case 4, ftag = '_body_nobasin';
end

weight_input = sprintf('%s%s/%s%s.dat',cap_path,evid,ftag0,ftag)
weight_output = sprintf('%s%s/%s%s_picks_and_pol.dat',cap_path,evid,ftag0,ftag)
% testing
weight_output = sprintf('~/%s%s_picks_and_pol.dat',ftag0,ftag)

% testing
%[stn_tag, pol, edist, PV_wt, PR_wt, SV_wt, SR_wt, ST_wt, P_arrival, P_len, S_arrival, S_len, waveform_shft] = read_cap_weight(weight_input);

idb = 1;            % Use AEC db
HWIN_SEC = 2;       % this is really important because otime is different for AEC and IRIS
iUseReviewed = 0;   % Analyst reviewed or unreviewed polarity
add_picks2weight(weight_input,weight_output,idb,iUseReviewed,HWIN_SEC);

%==========================================================================
