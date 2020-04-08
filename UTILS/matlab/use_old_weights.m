clear all
close all
clc

% by Kyle Smith
% created on 2020/03/26
% this code will copy weights and polarity from one weight file into another

for evnum = 8:21
    
    old_and_new_wts
    
    oldfile = [kylecapdir 'nenanabasin_MTs/weight_files/old_wts/' neweid '_weight.dat'];
    std_wt = [kylecapdir 'MTs/' datestr(evdate,'mmmdd') 'MT/' Vtype '/' neweid '/weight.dat'];
    newfile = [kylecapdir 'MTs/' datestr(evdate,'mmmdd') 'MT/' Vtype '/' neweid '/' wtfilename];
    
    [stnm_o, pol_o, edist_o, PV_wt_o, PR_wt_o, SV_wt_o, SR_wt_o, ST_wt_o, P_arrival_o, P_len_o, S_arrival_o, S_len_o, waveform_shft_o] = read_cap_weight(oldfile);
    [stnm, pol, edist, PV_wt, PR_wt, SV_wt, SR_wt, ST_wt, P_arrival, P_len, S_arrival, S_len, waveform_shft] = read_cap_weight(std_wt);
    
    for kk = 1:length(stnm_o)
        stnm_o{kk}(1:length(neweid)) = neweid;
    end
    
    missing_stns = setdiff(stnm_o,stnm)
    added_stns = setdiff(stnm,stnm_o)
    
    fid = fopen([kylecapdir 'MTs/' datestr(evdate,'mmmdd') 'MT/' Vtype '/station_diff.txt' ],'w');
    fprintf(fid,'missing stations: %s\n',strjoin(missing_stns))
    fprintf(fid,'added stations: %s\n',strjoin(added_stns))
    fclose(fid)
    
    % copy same information for "colocated stations"
    DEstns = {[neweid '.DE.UAF01..HH'],[neweid '.DE.UAF02..HH']};
    nearDEstns = {[neweid '.XV.FAPT..HH'],[neweid '.AK.NEA2..BH']};
    if ~isempty(intersect(added_stns,DEstns))
        [~,Ade,~] = intersect(stnm,DEstns,'stable');
        [~,Ande,~] = intersect(stnm,nearDEstns,'stable');
        pol(Ade) = pol(Ade);
        PV_wt(Ade) = PV_wt(Ande);
        PR_wt(Ade) = PR_wt(Ande);
        SV_wt(Ade) = SV_wt(Ande);
        SR_wt(Ade) = SR_wt(Ande);
        ST_wt(Ade) = ST_wt(Ande);
        %for de_i = 1:length(DEstns)
        %	[~,Ade,Bde] = intersect(stnm,DEstns);
        %end
	%error('stopping')
    end
    
    [C,IA,IB] = intersect(stnm,stnm_o,'stable');
    
    pol(IA) = pol_o(IB);
    PV_wt(IA) = PV_wt_o(IB); PR_wt(IA) = PR_wt_o(IB); SV_wt(IA) = SV_wt_o(IB); SR_wt(IA) = SR_wt_o(IB); ST_wt(IA) = ST_wt_o(IB);
    write_cap_weight(newfile, stnm, pol, edist, PV_wt,PR_wt, SV_wt, SR_wt, ST_wt, P_arrival, P_len, S_arrival, S_len, waveform_shft);
end
