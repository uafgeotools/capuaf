clear all
close all
clc

misfit_factor = 15;  % works only if greater than 1
icov = 0;       % will calculate the covariance of posterior (not needed!)
err_norm=1;
igmt =0;
iminus = 0;
ifig = 0;
label = 'MOOS'; % OR 'MFSZ'
model = 'scak'; % OR 'tactmod'
itest = 1;

% choose the norm and weight files
norms = {'L1' 'L2'};
wts = {'M110' 'M111' 'M011' 'M112' 'M012' 'M101'};
if itest
    norms = {'L1'};
    wts = {'M111'};
end

% direcoty for input files for gmt plotting.
gmtdir = '/home/vipul/gmt/data/cap/';
gmtdir = strcat(gmtdir,label,'/');

% data directory (with log files and output files from CAP)
datadir = '/home/vipul/CAP/inv/scak/MOOS/';
% datadir = '/home/vipul/CAP/inv/scak/MFSZ/north/';

% NEW id's
eidlist={'20070911234634153' '20070919112226549' '20071003140612444' '20071010180326301' '20071128235703849'...
    '20080314093821771' '20080327230745201' '20080828231418631' '20080918194353069' '20081228071310738'...
    '20090124180950811' '20090215193500098' '20090223000427175' '20090317011333066' '20090407201255351'...
    '20090414171427415' '20090430045457938' '20090524094004552' '20090622192805162' '20090626164820729' '20090730223910267'};
if itest 
    %eidlist={'20071010180326301'};
    % eidlist={'20070911234634153'};
    % eidlist={'20090407201255351'}; % EXAMPLE EVENT PAPER
end

% log file format
grd_form = 'log*grid';
rnd_form = 'log*rand';

%---test events-------------

% datadir = '/home/vipul/GEOTOOLS/tomo_util/cap/';
% eidlist = {'20080418093700'};
%
% datadir = '/home/vipul/CAP/inv/scak/MISC/';
% eidlist = {'20140416202424814'};
%
% datadir = '/home/vipul/CAP/inv/scak/MISC/';
% eidlist = {'20140418184418317'};

%% Block for getting the misfit for all 21 events
for in=1:length(norms)
    norm=norms{in};
    for il=1:length(eidlist)
        eid = eidlist{il};
        
        evtdir = strcat(datadir,eid,'/');
        
        % filename = strcat('/home/vipul/CAP/inv/scak/319605/log/log_grid_fine');
        
        for ii=1:length(wts)
            close all
            fname = dir(strcat(datadir,eid,'/',norm,'/',wts{ii},'/',grd_form));
            filename = fname.name;
            fid = fopen(strcat(datadir,eid,'/',norm,'/',wts{ii},'/',filename));
            s = textscan(filename,'%s','delimiter','_');
            s = s{1};
            edep = s{2};
            
            wt=textscan(fid,'%f %f %f %f %f %f');       % read the log file
            fclose(fid);
            
            stk=wt{1};
            dip=wt{2};
            rak=wt{3};
            err=wt{4};
            mw=wt{5};
            misfit(:,il) = err;
            %---------------------------
            
            [ival,ix]=min(err);
            stk_min = stk(ix);
            dip_min = dip(ix);
            rak_min = rak(ix);
            
            if il==1  % get omega for only the first event (only for seeing distribution of misfit with omega across all events: saves time)
            M_min0 = dcfaultpar2CMT(stk_min,dip_min,rak_min,0);
        
            % Distane from (M)
            M_min = repmat(M_min0,1,length(stk));
            [Mi,k1,d1,n1,p1,p2,p3] = dcfaultpar2CMT(stk,dip,rak,0);
            [omega1,xi1] = CMT2omegadc_xi0(Mi,M_min,0,0);
            end

            % scaling error
            continue

        end
    end
end

%% Block for (1) scaling, (2) getting probability, (3) creating samples
% 
nbin=60;
iscale = 1;
iopts = 4;
 %close all
clear omega_outline scaled_misfit samples_outline
for ii = 1:3
    clear ikeep 
    hold on
    [omega_outline,misfit_outline(:,ii)] = cappts2outline(omega1,misfit(:,ii),nbin);
    
    max_misfit = max(misfit(:,ii));
    min_misfit = min(misfit(:,ii));
    if iopts==1
        % scaling option : 1
        % No sclaing of the misfit function
        scaled_misfit = misfit;
    elseif iopts==2
        % scale the misfit functions 
        % scaling option : 2
        % This is currently used in the paper
        scaled_misfit(:,ii) = 15*(log(misfit(:,ii)/min_misfit)/log(max_misfit/min_misfit));
        [~,misfit_outline(:,ii)] = cappts2outline(omega1,scaled_misfit(:,ii),nbin);
    elseif iopts==3
        % scaling option : 3
        % Simplification of what is currently used in paper (iopts 3)
        scaled_misfit(:,ii) = ((misfit(:,ii) - min_misfit)/(max_misfit - min_misfit)).*(max_misfit./misfit(:,ii));
        [~,misfit_outline(:,ii)] = cappts2outline(omega1,scaled_misfit(:,ii),nbin);
    elseif iopts==4
        % scaling option : 4
        % Simply scaling the misfit. This is same as using different sigma
        % for all of them. 
        scaled_misfit(:,ii) = 15*misfit(:,ii);
        [~,misfit_outline(:,ii)] = cappts2outline(omega1,scaled_misfit(:,ii),nbin);
    end
    
    % get probability
    ptry = exp(-scaled_misfit(:,ii));
    ptrysum = sum(ptry);
    ptry = ptry/ptrysum;
    [~,ptry_outline(:,ii)] = cappts2outline(omega1,ptry,nbin);
    
    % get posterior samples
    % max_ptry = max(ptry);
    % min_ptry = min(ptry);
    % chance = min_ptry + (max_ptry - min_ptry)*rand(length(stk),1);
    % ikeep = find(chance>ptry);
    % [samples_omega,samples_outline(:,ii)] = cappts2outline(omega1(ikeep),ptry(ikeep),nbin);
end

figure
plot(omega_outline,misfit_outline)
xlabel('omega')
ylabel('misfit');
figure
plot(omega_outline,ptry_outline)
%hold on
%plot(samples_omega,samples_outline,'-s')
xlabel('omega')
ylabel('probability');


% Add sampling algorithm 
