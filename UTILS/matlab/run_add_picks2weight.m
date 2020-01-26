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
evid = '20170429111548898';
cap_path = '/home/ksmith/REPOSITORIES/capuaf/MTs/Apr29MT/V1/';

case_no = 4;
switch case_no 
    case 1
    	weight_input = [cap_path evid '/weight.dat'];
    	weight_output = [cap_path evid '/weight_wpicks_and_pol.dat'];
        weight_output = '~/weight_wpicks_and_pol.dat';  % testing
    case 2
    	weight_input = [cap_path evid '/weight_body.dat'];
    	weight_output = [cap_path evid '/weight_body_wpicks_and_pol.dat']
    case 3 
    	weight_input = [cap_path evid '/weight_surf.dat'];
    	weight_output = [cap_path evid '/weight_surf_wpicks_and_pol.dat']
    case 4 
    	weight_input  = [cap_path evid '/weight_body_nobasin.dat'];
        weight_output = [cap_path evid '/weight_body_nobasin_wpicks_and_pol.dat'];
        % TESTING
        weight_output = '~/weight_body_nobasin_wpicks_and_pol.dat'
end

idb = 1;            % Use AEC db
HWIN_SEC = 2;       % this is really important because otime is different for AEC and IRIS
iUseReviewed = 1;   % Analyst reviewed or unreviewed polarity
add_picks2weight(weight_input,weight_output,idb,iUseReviewed,HWIN_SEC)
