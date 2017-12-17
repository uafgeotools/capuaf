% Updates weight files for CAP
% Uses:
% 1. update (or scale) weights for all stations (no sorting)
% 2. sort by any of the sac header (currenty azimuth only). Build upon if
% you need any other sorting (example distance)
%
% Usage: enter event directory, event id and weight file name
% new weight file is created in the same directory as old one
%
% Vipul Silwal
% 06/09/2014

clear all
isort=0;    % if 1 then sort (else use script for scaling weights)
sortpar=1;  % 1 for azimuthal sorting; 2 for distance sorting


% EXAMPLES
% eid = '20090407201255351';
% evtdir = '/home/vipul/GEOTOOLS/tomo_util/cap/';

%  eid = '20090730223910267';
%  evtdir = strcat('/home/vipul/CAP/inv/scak/MOOS/',eid,'/');

% eid = '20140831030657325';
% eid = '20140831122458286';
% eid = '20080720080937779';
% eid = '20140831030657325';
% eid = '20141021003658579';
% eid = '20141023163024175';

eidlist={'20070911234634153' '20070919112226549' '20071003140612444' '20071010180326301' '20071128235703849'...
    '20080314093821771' '20080327230745201' '20080828231418631' '20080918194353069' '20081228071310738'...
    '20090124180950811' '20090215193500098' '20090223000427175' '20090317011333066' '20090407201255351'...
    '20090414171427415' '20090430045457938' '20090524094004552' '20090622192805162' '20090626164820729' '20090730223910267'};

%eidlist={'20150131173903205'};
%eidlist={'20080501234643121'};
eidlist={'20070919112226549'};
eidlist={'20141023163023968'};

%edir = '/home/vipul/GEOTOOLS/tomo_util/cap/';
edir = '/home/vipul/CAP/inv/scak/MOOS/';
%edir = '/home/vipul/CAP/inv/scak/bering/';
%edir = '/home/vipul/CAP/inv/scak/SCAK/';
edir = '/home/vipul/CAP/inv/scak/MFSZ/inv2/';

for  id=1:length(eidlist)
    eid=eidlist{id}
    evtdir = strcat(edir,eid,'/',eid,'/');
    
    % wt_filename = 'weight111.dat';
    wt_filename = 'weight111.dat';
    % new_wt_filename = strcat('UPT_',wt_filename);
    new_wt_filename = 'weight011.dat';
    %new_wt_filename = 'weight101.dat';
    
    % read weight file
    fid=fopen(strcat(evtdir,wt_filename));
    wt = textscan(fid,'%s %f %d %d %d %d %d %f %f %f %f %f');
    fclose(fid);
    
    % set the scaling factor
    if strcmp(new_wt_filename,'weight011.dat')
        PV_fact = 0;
        PR_fact = 0;
        SV_fact = 1;
        SR_fact = 1;
        SH_fact = 1;
    elseif strcmp(new_wt_filename,'weight101.dat')
        PV_fact = 1;
        PR_fact = 1;
        SV_fact = 0;
        SR_fact = 0;
        SH_fact = 0;
    end

    
    % update weights
    %[stnm,pol]=regexp(wt{1},'\.','split')
    stnm = wt{1};
    dist = wt{2};
    w1 = wt{3}*PV_fact;
    w2 = wt{4}*PR_fact;
    w3 = wt{5}*SV_fact;
    w4 = wt{6}*SR_fact;
    w5 = wt{7}*SH_fact;
    tp = wt{8};
    ts = wt{9};
    
    % choice of sorting
    if sortpar==1
        head = 'AZ';
    elseif sortpar==2
        head = 'DIS';
    end
    
    % sort by SAC header information (example azimuthal sorting)
    if isort
        for ii=1:length(wt{1})
            %[stnm(ii),pol]=regexp(stnm2{ii},'\.','split')
            zfile=strcat(evtdir,stnm{ii},'.z')
            zcomp = rsac(zfile);
            sort_list(ii) = lh(zcomp,head);
            %az(ii) = lh(zcomp,'DIST');
        end
        [~,ind]=sort(sort_list);
    else
        ind=1:length(wt{1});
    end
    
    % write new weight file
    fid=fopen(strcat(evtdir,new_wt_filename),'w');
    for jj=1:length(wt{1})
        ii=ind(jj);
        fprintf(fid,'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%3.2f\t%3.2f\t%d\t%d\t%d\n',...
            stnm{ii},dist(ii),w1(ii),w2(ii),w3(ii),w4(ii),w5(ii),tp(ii),ts(ii),0,0,0);
    end
    fclose(fid);
end