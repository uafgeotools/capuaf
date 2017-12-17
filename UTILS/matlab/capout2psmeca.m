function [otime,elon,elat,edep,M,M0,Mw,eid] = capout2psmeca(ddir)
% INPUT
%   ddir      Data directory containing CAP output files 
% OUTPUT can be used to create psmeca files, catalog.txt or CMTSOLTIONS
% (See EXAMPLES)

q = dir(strcat(ddir,'*out'));

% read from the CAP.out files
for ii=1:length(q)
    [otime(ii),elat(ii),elon(ii),edep(ii),strike,dip,rake,M(:,ii),Mw,eid(ii),capmod,capdep(ii),rms,vr,Pwin,Swin,Nstn(ii),Pstn,Sstn,...
    stnm,lampPV,lampPR,lampSV,lampSR,lampSH,corrPV,corrPR,corrSV,corrSR,corrSH,pPol,ePol] = read_capout(strcat(ddir,q(ii).name));
end
[Mw,M0] = CMT2mw(M);

% write latex table
write_latex_table(elat,elon,capdep,M,eid,Nstn);


% ====================================
% EXAMPLES
if 0
    % Cook Inlet events
    ddir = '/home/vipul/CAP/inv/scak/NEHRP/ci/out_files/';
    filename = '/home/vipul/gmt/data/NEHRP/ci/ci';
    [otime,elon,elat,edep,M,M0,Mw,eid] = capout2psmeca(ddir);
    % Beluga events 
    ddir = '/home/vipul/CAP/inv/scak/NEHRP/beluga/out_files/';
    filename = '/home/vipul/gmt/data/NEHRP/beluga/beluga';
    [otime,elon,elat,edep,M,M0,Mw,eid] = capout2psmeca(ddir);  
    % north susitna events 
    ddir = '/home/vipul/CAP/inv/scak/NEHRP/north_susitna/M_greater_than_3/out_files/';
    filename = '/home/vipul/gmt/data/NEHRP/ns/ns';
    [otime,elon,elat,edep,M,M0,Mw,eid] = capout2psmeca(ddir)
    % All NEHRP events
    ddir = '/home/vipul/CAP/inv/scak/NEHRP/out_files/';
    filename = '/home/vipul/gmt/data/NEHRP/cap/nehrp';
    [otime,elon,elat,edep,M,M0,Mw,eid] = capout2psmeca(ddir)
    % Josh SCAK catalog; (2009-2016)
    ddir = '/home/jcpurba/PROJECTS/CAP/inv/scak/out_scak_all/';
    filename = '/home/jcpurba/PROJECTS/CAP/inv/scak/gmt/scak_all';
    fname = '/home/vipul/data/josh/Josh'
    [otime,elon,elat,edep,M,M0,Mw,eid] = capout2psmeca(ddir)
    % northern thrust and fold event
    ddir = '/home/vipul/CAP/inv/tact/20070217145140160/out_files/';
    filename = [];
    [otime,elon,elat,edep,M,M0,Mw,eid] = capout2psmeca(ddir);
    %--------------------------------
    
    ipsmeca = 0;
    iwritemech = 0;
    iwriteCMTSOLUTION = 0;
    
    if ipsmeca
        write_psmeca(filename,otime,elat,elon,edep,M,eid);
    end
    
    if iwritemech
        %fname = '/home/vipul/REPOSITORIES/manuscripts/vipul/papers/2016nehrp/data/NEHRP';
        write_mech_table(fname,otime,elat,elon,edep,M,eid)
    end
    
    if iwriteCMTSOLUTION
        tshift = zeros(length(otime),1);
        hdur = tshift;
        write_CMTSOLUTION('./',0,otime,tshift,hdur,elat,elon,edep,M,eid,[],[],[]);
    end
end