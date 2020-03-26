% this code might not work for stations with polarity informtion
oldfile = '/home/ksmith/REPOSITORIES/capuaf/MTs/Aug28MT/V1/20180828151844000/weight.dat';

neweid = '20180828151843464';
newfile = ['/home/ksmith/REPOSITORIES/capuaf/MTs/Aug28MT/V2/' neweid '/weight_picks_and_pol.dat'];
[stnm_o, pol_o, edist_o, PV_wt_o, PR_wt_o, SV_wt_o, SR_wt_o, ST_wt_o, P_arrival_o, P_len_o, S_arrival_o, S_len_o, waveform_shft_o] = read_cap_weight(oldfile);
[stnm, pol, edist, PV_wt, PR_wt, SV_wt, SR_wt, ST_wt, P_arrival, P_len, S_arrival, S_len, waveform_shft] = read_cap_weight(newfile);

for kk = 1:length(stnm_o)
    stnm_o{kk}(1:length(neweid)) = neweid;
end

missing_stns = setdiff(stnm_o,stnm)
added_stns = setdiff(stnm,stnm_o)

DEstns = {[neweid '.DE.UAF01..HH'],[neweid '.DE.UAF02..HH']};
nearDEstns = {[neweid '.XV.FAPT..HH'],[neweid '.AK.NEA2..BH']};
if ~isempty(intersect(added_stns,DEstns))
    [~,Ade,~] = intersect(stnm,DEstns,'stable');
    [~,Ande,~] = intersect(stnm,nearDEstns,'stable');
    PV_wt(Ade) = PV_wt(Ande); 
    PR_wt(Ade) = PR_wt(Ande); 
    SV_wt(Ade) = SV_wt(Ande); 
    SR_wt(Ade) = SR_wt(Ande); 
    ST_wt(Ade) = ST_wt(Ande);
    %for de_i = 1:length(DEstns)
%	[~,Ade,Bde] = intersect(stnm,DEstns);
    %end
end

[C,IA,IB] = intersect(stnm,stnm_o,'stable');

PV_wt(IA) = PV_wt_o(IB); PR_wt(IA) = PR_wt_o(IB); SV_wt(IA) = SV_wt_o(IB); SR_wt(IA) = SR_wt_o(IB); ST_wt(IA) = ST_wt_o(IB);

write_cap_weight(newfile, stnm, pol, edist, PV_wt,PR_wt, SV_wt, SR_wt, ST_wt, P_arrival, P_len, S_arrival, S_len, waveform_shft);
