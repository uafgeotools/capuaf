function [otime,elat,elon,edep,strike,dip,rake,M,Mw,eid,capmod,capdep,rms,vr,Pwin,Swin,Nstn,Pstn,Sstn,...
    stnm,lampPV,lampPR,lampSV,lampSR,lampSH,corrPV,corrPR,corrSV,corrSR,corrSH,pPol,ePol] = read_capout_old(filename)
% OLDER VERSION of read_capout.m (git commit: commit 3066f1a)
% 
% READ_CAPOUT read the output file produced by CAP
%
% EXAMPLE:
%   filename = '/usr/local/GEOTOOLS_copy/tomo_util/cap/20080418093700_check/cus_015.out';
%   filename = '/home/vipul/GEOTOOLS/tomo_util/cap/20090407201255351/scak_041.out';
%   filename = '/home/vipul/CAP/inv/scak/MFSZ/RESULTS/20001129103547240.out';
%   filename = [getenv('CAPRUN') '/inv/scak/20090407201255351/scak_041.out'];
%   [otime,elat,elon,edep,M,eid,rms,vr] = read_capout(filename);
%   figure; plot_beachballs(M,elon,elat,1);
%
% FUTURE WORK:
%   - add a line with the origin time as # origin_time_sac_header yyyy mm dd etc
%        (in case the eid is NOT the origin time)
%   - add output for stations
%   - possibly adapt this to read a concatenated set of output file
%        (or loop over the output files you want to read)
%   - fix otime2eid to check the number of characters in eid
%

NHEADER_LINES = 4;
bplot=0;            % multuiple plots
write_meca = 0;     % write and psmeca file for GMT plot    
idep = 1;           % 0=use AEC dep; 1=use CAP dep

% check if file exists
if ~exist(filename,'file'), error(sprintf('%s does not exist',filename)); end

% read in concatenated CMTSOLUTION files
lines = textread(filename,'%s','delimiter','\n');
nlines = length(lines);

nstation = nlines - NHEADER_LINES;

%disp(['File : ' filename ]);
%disp(['Number of stations : ' num2str(nstation) ]);

% extract the header information
[~,eid,~,ctag,~,strike,dip,rake,~,Mw,~,rms,X,~,strike_sig,dip_sig,rake_sig,~,delta,~,~,gamma,~,~,vr,~,dat2] = ...
    strread(lines{1},'%s%s%s%s%s%f%f%f%s%f%s%f%f%s%f%f%f%s%f%f%s%f%f%s%f%s%f');
[~,~,~,elat,~,elon,~,edep] = strread(lines{2},'%s%s%s%f%s%f%s%f');
[~,~,~,M0_dyne_cm,Mxx,Mxy,Mxz,Myy,Myz,Mzz] = strread(lines{3},'%s%s%s%f%f%f%f%f%f%f');
[~,~,norm,~,~,Pwin,~,Swin,~,~,Nstn,~,Pstn,~,Sstn] = strread(lines{4},'%s%s%s%s%s%f%s%f%s%s%f%s%f%s%f');

% Removing trailing '/' from eid (if present)
eid = eid{1};               % change cell to char
if strcmp(eid(end),'/')
    eid(end)=[];
    disp('hi');
end
eid = cellstr(eid);         % back to cell

% WARNING: this assumes that the eid has the format yyyymmddHHMMSSFFF
otime = eid2otime(eid);

whos eid
% depth used in CAP inversion
ctemp = char(ctag);
capdep = str2num(ctemp(end-2:end));
% model used in CAP inversion
capmod = cellstr(ctemp(1:end-4));

M0 = M0_dyne_cm * 1e-7;     % N-m
% note: Mw is calculated in CAP here: amp=pow(10.,1.5*mt[0].par+16.1-20);

% moment tensor tensor computed from TT parameters
% note: angles are rounded to integers, so this will not be as accurate as
%       the MT entries listed in the output file
Mcheck = TT2CMT(gamma,delta,M0,strike,dip,rake);

% moment tensor listed in CAP
Maki = M0 * [Mxx Myy Mzz Mxy Mxz Myz]';
[Mcap,T] = convert_MT(2,1,Maki);

disp('comparing MT computed from TT and MT listed in output');
disp('both are in the convention of GCMT');
disp([Mcheck Mcap]);

M = Mcap;

errsum=0;
% READ station info
for ii=1:Nstn
    iline = ii+4;
    [stnm(ii),~,wtPV(ii),errPV(ii),corrPV(ii),shiftPV(ii),lampPV(ii),DampPV(ii),SampPV(ii),...
        wtPV(ii),errPR(ii),corrPR(ii),shiftPR(ii),lampPR(ii),DampPR(ii),SampPR(ii),...
        wtSV(ii),errSV(ii),corrSV(ii),shiftSV(ii),lampSV(ii),DampSV(ii),SampSV(ii),...
        wtSR(ii),errSR(ii),corrSR(ii),shiftSR(ii),lampSR(ii),DampSR(ii),SampSR(ii),...
        wtSH(ii),errSH(ii),corrSH(ii),shiftSH(ii),lampSH(ii),DampSH(ii),SampSH(ii),...
        pPol(ii),ePol(ii)] = strread(lines{iline},'%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
    errsum =errsum+ wtPV(ii)*errPV(ii)+wtPV(ii)*errPR(ii)+ wtSV(ii)*errSV(ii)+wtSR(ii)*errSR(ii)+wtSH(ii)*errSH(ii);
end
disp(sprintf('Sum of error: %f',errsum));

% plotting various measures of misfit: err, corr, shift, amp
if bplot
    % NOTE: example for PR time window only
    figure;
    subplot(3,3,1); plot(errPR,corrPR,'.');     xlabel('err'); ylabel('corr');
    subplot(3,3,2); plot(errPR,shiftPR,'.');    xlabel('err'); ylabel('shift');
    subplot(3,3,3); plot(errPR,lampPR,'.');     xlabel('err'); ylabel('amp');
    subplot(3,3,4); plot(lampPR,corrPR,'.');    xlabel('amp'); ylabel('corr');
    subplot(3,3,5); plot(shiftPR,corrPR,'.');   xlabel('shift'); ylabel('corr');
    subplot(3,3,6); plot(shiftPR,lampPR,'.');   xlabel('shift'); ylabel('amp');
    subplot(3,3,7); plot(errPR,lampPR,'.');     xlabel('err'); ylabel('amp');
    subplot(3,3,8); plot(errPR,corrPR,'.');     xlabel('err'); ylabel('corr');
end

%--------------------------------------------------------------------------
if write_meca
    filename = eid{1};
    if idep==0
        write_psmeca(filename,otime,elat,elon,edep,M,eid);
    else
        write_psmeca(filename,otime,elat,elon,capdep,M,eid);
    end
end
%==========================================================================
