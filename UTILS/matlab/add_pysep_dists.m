clear all
close all
clc

% obtaining MT distances retrospectively and writing to new obspy like file but specifies max stn dist grabs
aggwtfile = '/home/ksmith/REPOSITORIES/capuaf/nenanabasin_MTs/weight_files/old_wts/tail_file.txt';

[stnm, pol, edist, PV_wt, PR_wt, SV_wt, SR_wt, ST_wt, P_arrival, P_len, S_arrival, S_len, waveform_shft] = read_cap_weight(aggwtfile);
eidlength = 4+2+2+2+2+2+3; % corresponds to yyyy mm dd ...
eiddaylength = 4+2+2;
for s_i = 1:length(stnm)
    old_eidset{s_i} = stnm{s_i}(1:eidlength);
    old_eiddayset{s_i} = stnm{s_i}(1:eiddaylength);
end
edistmax = ceil(edist/25)*25;

obspy_evs = '/home/ksmith/REPOSITORIES/manuscripts/kyle/papers/nenanabasin/data/nenanabasin';

[~,eid,~,lon,lat,dep_m,mag] = textread([obspy_evs '_catalogMdep_obspy.txt'],'%f%s%s%f%f%f%f');
dep = dep_m / 1000;
otime = eid2otime(eid);

%[otime,lon,lat,dep,mag,eid] = read_seis_obspy(obspy_evs);

for e_i = 1:length(eid) 
    eid_day{e_i} = eid{e_i}(1:eiddaylength); 
end
% start writing new file
n = length(lon);

%filename = [obspy_evs '_obspy_pysep_dwnld.txt'];
filename = [obspy_evs '_catalogMdep_obspy_pysep_dwnld.txt'];

disp(['writing ' filename]);

stfmt = '%3i  %s  %s %12.6f %12.6f %8.1f %6.2f %1.0f';

fid = fopen(filename,'w');
for ii = 1:n
    nd_i = find(strcmp(eid_day{ii},old_eiddayset));
    if ~isempty(nd_i)
    	newedist = edistmax(nd_i);
    else
    	newedist = 999;
    end 
    % note: depth is in meters
    fprintf(fid,[stfmt '\n'],...
        ii,eid{ii},datestr(otime(ii),'yyyy-mm-ddTHH:MM:SS.FFF'),...
        lon(ii),lat(ii),dep(ii)*1000,mag(ii),newedist);
end
fclose(fid);
