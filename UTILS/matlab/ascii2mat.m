function [mat_filename]  = ascii2mat(ascii_filename, mat_filename)
% simple script to convert file from ascii to mat format

% same name with '.mat' appended if no mat_filename is provided
if nargin==1
    mat_filename = strcat(ascii_filename,'.mat');
end
% read file
file = dlmread(ascii_filename);
% save with the variable name 'file'
save(mat_filename,'file');

if 0
    ascii_filename = '/home/vipul/CAP/inv/scak/NW/20140418184418152/L1/M111/log_020_grid';
    mat_filename = '/home/vipul/CAP/inv/scak/NW/20140418184418152/L1/M111/log_020_grid.mat';
    
    ascii_filename = '/home/vipul/CAP/inv/scak/MISC/20140416202423770/L1/M111/log_072_grid';
    mat_filename = '/home/vipul/CAP/inv/scak/MISC/20140416202423770/L1/M111/log_072_grid.mat'
    
    ascii_filename = '/home/vipul/CAP/inv/scak/MISC/20140416202423770/L1/M111/log_072_rand';
    mat_filename = '/home/vipul/CAP/inv/scak/MISC/20140416202423770/L1/M111/log_072_rand.mat'
    ascii2mat(ascii_filename, mat_filename)
    
    ascii_filename = '/home/vipul/CAP/inv/scak/MOOS/20070911234634153/log_100_grid';
    ascii2mat(ascii_filename)
    
    ascii_filename = '/home/vipul/CAP/inv/scak/MISC/20140416202423770/log_072_grid';
    ascii2mat(ascii_filename)
end
end


    