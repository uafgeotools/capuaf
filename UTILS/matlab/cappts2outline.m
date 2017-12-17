function [omega_outline,misfit_outline,M_outline] = cappts2outline(omega,misfit,M,nbin)

iplot = 0;
[omega,iind] = sort(omega);
misfit = misfit(iind);
M = M(:,iind);

edges = linspace(0,180,nbin+1);
nedges = length(edges);

misfit_min = NaN(nedges,1);
misfit_max = NaN(nedges,1);
omega_vec_min = NaN(nedges,1);
omega_vec_max = NaN(nedges,1);

[omega_max,imax] = max(omega);
misfit_omega_max = misfit(imax);
disp(sprintf('max omega is %.2f with misfit %.3f',omega_max,misfit_omega_max));

imatch_inc=0;
for ii=1:nbin
   
    imatch = find(and(omega >= edges(ii), omega < edges(ii+1))==1);
    
    if isempty(imatch)
        warning('no omega points in this bin (assigning minmax values from previous bin)');
        misfit_min(ii) = misfit_min(ii-1);
        misfit_max(ii) = misfit_max(ii-1);
    else
        [misfit_min(ii),iind] = min( misfit(imatch) );
        omega_vec_min(ii) = omega(imatch_inc+iind);
        M_min(:,ii) = M(:,imatch_inc+iind);
        [misfit_max(ii),iind] = max( misfit(imatch) );
        omega_vec_max(ii) = omega(imatch_inc+iind);
        M_max(:,ii) = M(:,imatch_inc+iind);
    end
    imatch_inc=imatch_inc+length(imatch);
end

% final point
% BEST OPTION: take misfit of the largest omega
% ALTERNATIVE: take misfit of the mean of the left bin edge (which could
%               have been extrapolated from previous bins)
%misfit_min(ii+1) = mean([ misfit_min(ii-1) misfit_max(ii-1)]);
misfit_min(ii+1) = misfit_omega_max;
misfit_max(ii+1) = misfit_min(ii+1);
%omega_vec_min(ii+1) = omega_max;
%omega_vec_max(ii+1) = omega_max;

% reorder points in anticlockwise sense
%omega_outline = [omega_vec_min ;flip(omega_vec_max)];
misfit_outline = [misfit_min; flip(misfit_max)];
omega_outline = [edges flip(edges)]';
M_outline = [M_min flip(M_max')'];


% plot
if iplot
    plot(omega_outline,misfit_outline)
end