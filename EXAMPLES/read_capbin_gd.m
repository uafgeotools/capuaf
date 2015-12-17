function [gamma, delta, kappa, theta, sigma, misfit_wf, mw, misfit_pol] = read_capbin_gd(inputfile)
% script to read cap output capout_gd.bin
% data is put into columns
%
% usage:
% [gamma, delta, kappa, theta, sigma, misfit_wf, mw, misfit_fmp] = read_capbin_gd('capout_gd.bin') ;
%
% 20151216 celso alvizuri - cralvizuri@alaska.edu 
%----------------------------------------------------------

% there is a way to not have to specify total number of rows.
% [ongoing]
nrows = 102960;

% read data
m = memmapfile(inputfile, 'Format', {'single', [8 nrows] 'floats'});
m.data(1).floats(1:16) ;
a = m.data(1).floats()';
gamma = a(:,1);
delta = a(:,2);
kappa = a(:,3);
theta = a(:,4);
sigma = a(:,5);
mw         = a(:,6);
misfit_wf  = a(:,7);
misfit_pol = a(:,8);

