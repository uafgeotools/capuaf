clear all
close all
clc

% Nenana basin study
weight_input = '/home/ksmith/REPOSITORIES/capuaf/20181113152641907/weight_body.dat'
weight_output = '/home/ksmith/REPOSITORIES/capuaf/20181113152641907/weight_body_wpicks.dat'
idb = 1;            % Use AEC db
HWIN_SEC = 2;       % this is really important because otime is different for AEC and IRIS
iUseReviewed = 1;   % Analyst reviewed or unreviewed polarity
add_picks2weight(weight_input,weight_output,idb,iUseReviewed,HWIN_SEC)
