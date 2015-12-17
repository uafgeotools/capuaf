function[mrr, mtt, mpp, mrt, mrp, mtp, mag, misfit_wf, misfit_fmp] = read_capbin_mt(inputfile);
% script to read cap output capout_gd.bin
% data is put into columns
%
% usage:
% [mrr, mtt, mpp, mrt, mrp, mtp, misfit_wf, mag, misfit_fmp] = read_capbin_mt('capout_mt.bin');
%
% 20151216 celso alvizuri - cralvizuri@alaska.edu 
%----------------------------------------------------------

% there is a way to not have to specify total number of rows.
% [ongoing]
nrows = 102960;

% read data
m = memmapfile(inputfile, 'Format', {'single', [9 nrows] 'floats'});
m.data(1).floats(1:16) ;
a = m.data(1).floats()';
mrr = a(:,1);
mtt = a(:,2);
mpp = a(:,3);
mrt = a(:,4);
mrp = a(:,5);
mtp = a(:,6);
mag = a(:,7);
misfit_wf = a(:,8);
misfit_fmp = a(:,9);

