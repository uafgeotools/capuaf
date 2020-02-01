function [stnm, pol, edist, PV_wt, PR_wt, SV_wt, SR_wt, ST_wt, P_arrival, P_len, S_arrival, S_len, waveform_shft] = read_cap_weight(weight_filename)

% get number of columns
% THIS STEP OUGHT TO BE UNNECESSARY IF WE CAN USE fscanf OR textscan PROPERLY
[~,ncol] = system(sprintf("head -1 %s | tr '|' ' ' | wc -w",weight_filename));
ncol = str2num(ncol);

fid = fopen(weight_filename);
if ncol==12
    C = textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f');
elseif ncol==13
    C = textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f %f');
else
    ncol
    weight_filename
    error('can only have 12 or 13 columns in cap weight file');
end
fclose(fid);

stnm_pol = C{1};
% Get polarities
for ii = 1:length(stnm_pol)
    tmp = strsplit(stnm_pol{ii},'/');
    stnm{ii} = tmp{1};
    if length(tmp)==2
        pol(ii) = str2double(tmp{2});
    else
        pol(ii) = 0;
    end
end
edist = C{2};
PV_wt = C{3};
PR_wt = C{4};
SV_wt = C{5};
SR_wt = C{6};
ST_wt = C{7};
P_arrival = C{8};
P_len = C{9};
S_arrival = C{10};
S_len = C{11};
waveform_shft = C{12};

% EXamples
if 0
    weight_filename = '/home/vipul/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/20090407201253480/weight.dat';
    [stnm, pol, edist, PV_wt, PR_wt, SV_wt, SR_wt, P_arrival, P_len, S_arrival, S_len, waveform_shft] = read_cap_weight(weight_filename);
    weight_filename = '/home/vipul/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/20090407201253480/weight_test.dat';
    write_cap_weight(weight_filename, stnm, pol, edist, PV_wt,PR_wt, SV_wt, SR_wt, ST_wt, P_arrival, P_len, S_arrival, S_len, waveform_shft)
end
