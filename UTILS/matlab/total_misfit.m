function [F,misfit_wf,misfit_fmp] = total_misfit(misfit_wf, misfit_fmp, Np, misfit_pol_weight)
% combine various misfit measures (Add more)
% See : MTbrick_Csection_binary.m for example
% Final misfit is multiplied by the same scale factor used in SilwalTape2016

% to prevent division by 0 while normalizing (line 19)
if Np==0;
    Np=1;
end

% weight for each misfit component
%misfit_pol_weight = 0.25; % belongs [0,1] 
misfit_wf_weight = (1-misfit_pol_weight);

%Max_pol = max(misfit_pol);
%Max_wf = max(misfit_wf);

misfit_wf = misfit_wf_weight * misfit_wf;
misfit_fmp = misfit_pol_weight * misfit_fmp/Np;

% Compute total misfit (NOTE: we might want to normalize each misfit component)
F = misfit_wf + misfit_fmp;
