function [] = write_cap_weight(weight_filename, stnm, pol, edist, PV_wt,...
    PR_wt, SV_wt, SR_wt, ST_wt, P_arrival, P_len, S_arrival, S_len, waveform_shft)

% Combine station name and polarity into one string
for ii = 1:length(pol)
    if pol(ii)==0
        stnm_pol{ii} = stnm{ii};
    else
        stnm_pol{ii} = strcat(stnm{ii},'/',num2str(pol(ii)));
    end
end

% Write everything to a weight file
fid = fopen(weight_filename,'w');
stfmt = '%34s %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.4f %4.0f %4.0f %4.0f %4.0f \n';

for ii = 1:length(stnm_pol)
    fprintf(fid, stfmt, stnm_pol{ii}, edist(ii), PV_wt(ii),...
        PR_wt(ii), SV_wt(ii), SR_wt(ii), ST_wt(ii), P_arrival(ii), P_len(ii),...
        S_arrival(ii), S_len(ii), waveform_shft(ii));
end
fclose(fid);