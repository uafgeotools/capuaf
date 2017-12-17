function [otim,lat,lon,dep,stk,dp,rak,Mwi] = read_mech_CAP(ed)
% Reads CAP MT solutions and outputs as psmeca, table and latex format
%
% Vipul Silwal
% 06/04/2013

% Set input catalog direcotry and output dire/file here
if 0
    ddir = '/home/vipul/CAP/inv/scak/MOOS/RESULTS/out/';    % 21 main events
    gmtdir = '/home/vipul/gmt/data/cmt/CAP/';
    filename = strcat(gmtdir,'CAP_21'); % AECMT
end

if 1
    ddir = '/home/vipul/CAP/inv/scak/SCAK/RESULTS/out_full/';  % 106 main events
    gmtdir = '/home/vipul/gmt/data/cmt/CAP/';
    gmtdir = './';
    filename = strcat(gmtdir,'CAP_106'); % AECMT + AECfp
end

if 0
    ddir = '/home/vipul/CAP/inv/scak/SCAK/RESULTS/out/';   % 88 other events
    gmtdir = '/home/vipul/gmt/data/cmt/CAP/';
    filename = strcat(gmtdir,'CAP_88'); % AECfpn
end

% Set the index here (Caution: Files may get overwritten)
igmt = 0;
itable = 1;
ilatex = 0;

% read cap*.out
q = dir(strcat(ddir,'*out'));
for ii=1:length(q)
    [otime(ii,1),elat(ii,1),elon(ii,1),edep(ii,1),strike(ii,1),dip(ii,1),rake(ii,1),M(:,ii),Mw(ii,1),...
        eid(ii,1),capmod,capdep(ii,1),rms,vr,Pwin,Swin,Nstn(ii,1),Pstn,Sstn]=...
        read_capout(strcat(ddir,q(ii).name));
end

% find the match
% [C,imatch,IB] = intersect(eid,ed);
% otim = otime(imatch);
% lat = elat(imatch);
% lon = elon(imatch);
% dep = capdep(imatch);
% Mwi = Mw(imatch);
% stk = strike(imatch);
% dp = dip(imatch);
% rak = rake(imatch);


% psmeca output for GMT
if igmt
    write_psmeca(filename,otime,elat,elon,capdep,M,eid);
end

% mecha table output
if itable
    write_mech_table(filename,otime,elat,elon,capdep,M,eid);
end
% output for latex table
if ilatex 
    write_latex_table(elat,elon,capdep,M,eid,Nstn);
end

%%%%%%%%%%%%%%%%%%%%%%% OTHER EXAMPLES %%%%%%%%%%%%%%%%%%%%%%%%
if 0
    gmtdir = '/home/vipul/gmt/data/cmt/';
    igmt = 0;
    itable = 0;
    ilatex = 1;
    axmoos = [-154 -146 58 62.5];
    ax3 = [axmoos -10 200];
    oran = [datenum(2007,8,15) datenum(2009,8,15)];
    Mwran = [3.5 10];
    
    % ====== Events for CAP that are in AEC fault plane catalog
    [otime,slat,slon,sdep,M,M0,Mw,eid2] = read_mech_AECfp(oran,ax3,Mwran);
    filename = strcat(gmtdir,'AEC/AECfm');
    if igmt % psmeca output for GMT
        write_psmeca(filename,otime,slat,slon,sdep,M,eid2);
    end
    [~,isort] = sort(sdep,'ascend');
    %[~,isort] = sort(Mw,'descend');
    display_eq_list(isort,otime,slon,slat,sdep,Mw,eid2);
    if itable % mech table output
        write_mech_table(filename,otime,slat,slon,sdep,M,eid2);
    end
    if ilatex % latex
        write_latex_table(slat,slon,sdep,M,eid2)
    end
    
    
    % ===== Events for CAP that are in AEC MT catalog
    [otime,slat,slon,sdep,M2,M0,Mw,eid,depc] = read_mech_AEC(oran,ax3,Mwran);
    filename = strcat(gmtdir,'AEC/Maec');
    if igmt
        write_psmeca(filename,otime,slat,slon,sdep,M2,eid);
    end
    [~,isort] = sort(sdep,'ascend');
    %[~,isort] = sort(Mw,'descend');
    display_eq_list(isort,otime,slon,slat,sdep,Mw,eid);
    if itable
        write_mech_table(filename,otime,slat,slon,sdep,M2,eid);
    end
    if ilatex
        write_latex_table(slat,slon,sdep,M2,eid)
    end
    
    % ====== use intersect.m instead
    for ii=1:length(eid)
        for jj=1:length(eid2)
            if strcmp(eid(ii),eid2(jj))
                Mi(:,ii) = M(:,jj);
                omega(ii) = CMT2omegadc_xi0(M2(:,ii),M(:,jj),0,0)
            end
        end
    end
    filename = strcat(gmtdir,'AEC/Mfm');
    if igmt
        write_psmeca(filename,otime,slat,slon,sdep,Mi,eid);
    end
    if itable
        write_mech_table(filename,otime,slat,slon,sdep,Mi,eid);
    end
end