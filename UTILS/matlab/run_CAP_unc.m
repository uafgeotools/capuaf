% SCRIPT SHOULD
% Input: CAP datadir info
% call run_dcMT_unc - for random log files
% call MT brick_section - for grid files
%
% Vipul Silwal
% 3 Aug, 2015
% REMEMBER to to use 'close all' inside the loop when running for all
% events

clear all; close all

% Input parameters
itest = 1;  % to test the script
igmt = 1;

% Get event info 
CAPdir = '/home/vipul/CAP/inv/scak/MOOS/RESULTS/'; % event dir(s) is here
%CAPdir = '/home/vipul/CAP/inv/scak/NW/';
%CAPdir = '/home/vipul/CAP/inv/scak/MISC/';

eids = dir(strcat(CAPdir,'20*')); %(criteria for events in 21st century)

% main GMT dir (all events directories are inside it)
GMTdir = '/home/vipul/gmt/data/cap/MOOS/'; % main GMT dir (all events 
%GMTdir = '/home/vipul/gmt/data/cap/NW/';
%GMTdir = '/home/vipul/gmt/data/cap/MISC/';
GMTdir = '/home/vipul/gmt/data/cap/test/';

% choose the norm and weight files
norms = {'L1' 'L2'};
wts = {'M110' 'M111' 'M011' 'M112' 'M012' 'M101'};
norms = {'L1'};
wts = {'M111'};

% get all the event eids 
for ii = 1:length(eids)
    eid(ii) = {eids(ii).name};
end

% run only for one event (for testing)
if itest
    %eid = {eids(1).name}; % just run for one event for testing
    %eid = {'20081228071310738'};
    %eid = {'20140418184418152'};
    eid = {'20090407201255351'};
    %eid = {'20070911234634153'};
    %eid = {'20140416202423770'};
    norms = {'L1'};
    wts = {'M111'};
end

for ii = 1:length(norms)
    nrm=norms{ii};
    for jj = 1:length(wts)
        weight = wts{jj};
        for kk = 1:length(eid)
            % directory containing CAP log files (random and grid search)
            ddir = strcat(CAPdir,eid{kk},'/',nrm,'/',weight,'/');
            % directory for saving files for GMT plotting
            gmtdir = strcat(GMTdir,eid{kk},'/',nrm,'/',weight,'/');
            if igmt
                % Estimate the uncertainty in solution
                [Mref, post_samples_info] = CAP_unc(eid{kk},ddir,gmtdir);
                % Generate cross-section in strike-dip-rake (kappa-theta-sigma
                % space)
                MTbrick_Csection(eid{kk},ddir,gmtdir,Mref);
                %close all
            else % No need to do cross-sections nor save files
                [Mref, post_samples_info]= CAP_unc(eid{kk},ddir);
                [Mref, post_samples_info]= sortrows(post_samples_info,7); % sort by omega
            end
            close all
        end
    end
end

