% Intended to perform posterior analaysis by reading in the summary files
%
% Vipul Silwal
% 9 Sept, 2015

clear all
close all
iwrite=0;
ddir = '/home/vipul/CAP/inv/scak/MOOS/RESULTS/';
dirs = dir(strcat(ddir,'2*'));
iAEC = 1;
inorm2 = 0;
if iAEC
    nrow=8;
else
    nrow=6;
end
wt_order = [6 5 1 4 2 3 7 8];
K  = 40; % misfit factor
N = 1e5;  % Number of random samples
dV = (4*pi*pi)/N;       % Volume of the brick

% See write_psmeca_MOOS_catalog
%AECMTfile = '/home/vipul/gmt/data/cmt/AEC/Maec_psmeca';
%AECfpfile = '/home/vipul/gmt/data/cmt/AEC/Mfm_psmeca';

%[lat1,lon1,dep1,M1,px1,py1,slabel1] = read_psmeca(AECMTfile);
%[lat2,lon2,dep2,M2,px2,py2,slabel2] = read_psmeca(AECfpfile);

axmoos = [-154 -146 58 62.5];
ax3 = [axmoos -10 200];
oran = [datenum(2007,8,15) datenum(2009,8,15)];
Mwran = [3.5 10];

[otime1,lat1,lon1,dep1,M1,M0,Mw1,eid1,depc] = read_mech_AEC(oran,ax3,Mwran);
% Events for CAP (call read_eq_AEC) - AECfp
[otime2,lat,lon,dep,M,M0,Mw,eid] = read_mech_AECfp(oran,ax3,Mwran);

for ii=1:length(eid1)
    for jj=1:length(eid)
        if strcmp(eid1(ii),eid(jj))
            Mw2(ii) = Mw(jj);
            dep2(ii) = dep(jj);
            M2(:,ii) = M(:,jj);
            omega(ii) = CMT2omegadc_xi0(M1(:,ii),M2(:,ii),0,0);
        end
    end
end
%DEP = [dep1 dep2' depc];

datadir = '/home/vipul/gmt/data/cmt/cmt_analysis/';
% MT tags for findings the files
Mlist = {'M011','M112','M012','M101','M110','M111','Maec','Mfm'};
% MT tags for latex printing
Mlist2 = {'M_{011}','M_{112}','M_{012}','M_{101}','M_{110}','M_{111}','M_{AEC}','M_{FM}'};
Mlist2 = {'\mb','\md','\me','\mc','\mf','\ma','\maec','\mfm'};
omegadir = '/home/vipul/gmt/data/cap/MOOS/';

nele = nrow*length(dirs); 

for ii=1:length(dirs)
    eid{ii} = dirs(ii).name;
    ffname = strcat(ddir,eid{ii},'/summary.dat')
    fid = fopen(ffname);
    w = textscan(fid,'%f %f %f %f %f %f %f %f');
    fclose(fid);
    norm = w{1};
    wt = w{2};
    stk = w{3};
    dip = w{4};
    rak = w{5};
    mag = w{6};
    dep = w{7};
    vr = w{8};
    Mii = dcfaultpar2CMT(stk,dip,rak,0);
    %fname = omega_dif
    k1 = 1;
    k2 = 1;
    lat_vec(1:nrow,ii) = lat1(ii);
    lon_vec(1:nrow,ii) = lon1(ii);
    otime(1:nrow,ii) = eid2otime(eid{ii});
    for jj = 1:length(norm)
        if norm(jj) == 1   % L1
            dep_vec(k1,ii) = dep(jj);  % AEC MT depth is different from fp depth
            mag_vec(k1,ii) = mag(jj);
            vr_vec1(k1,ii) = vr(jj);
            M_vec(:,k1) = Mii(:,jj);
            k1 = k1+1;
        end
        if inorm2
            if norm(jj) == 2   % L2
                dep_vec2(k2,ii) = dep(jj); % AEC MT depth is different from fp depth
                mag_vec2(k2,ii) = mag(jj);
                vr_vec2(k2,ii) = vr(jj);
                M_vec2(:,k2) = Mii(:,jj);
                k2 = k2+1;
            end
        end
    end
    
    if iAEC
        % Add AEC MT info
        dep_vec(k1,ii) = dep1(ii); % L1
        dep_vec2(k2,ii) = dep1(ii); % L2
        mag_vec(k1,ii) = Mw1(ii);
        mag_vec2(k1,ii) = Mw1(ii);
        M_vec(:,k1) = M1(:,ii);
        M_vec2(:,k2) = M_vec(:,k1);
        vr_vec1(k1,ii) = 0;
        vr_vec2(k1,ii) = 0;
        % Add AEC fp info
        dep_vec(k1+1,ii) = dep2(ii); % L1
        dep_vec2(k2+1,ii) = dep2(ii); % L2
        mag_vec(k1+1,ii) = Mw2(ii);
        mag_vec2(k1+1,ii) = Mw2(ii);
        M_vec(:,k1+1) = M2(:,ii);
        M_vec2(:,k2+1) = M_vec(:,k1+1);
        vr_vec1(k1+1,ii) = 0;
        vr_vec2(k1+1,ii) = 0;
    end
    
    % Compute the p(M0) - probability at the minimum
    for jj=1:6
        LogDir = strcat(ddir,eid{ii},'/L1/',Mlist{jj},'/');
        f = dir(strcat(LogDir,'*rand.mat'));
        fname =strcat(LogDir,f.name);
        [kappa,theta,sigma,F,Mw,delta,gamma,M0,M] = read_CAP_log(fname);
        Fsum = sum(exp(-F*K));
        [fmin,indx1] = min(F);
        Mmin = M(:,indx1);
        
        pM0 = exp(-fmin*K);     % unnormalized probability at the minimum
        kfac = 1/(dV*Fsum);     % normalization factor (Walt2015b notes)
        
        PM0(jj,ii) = kfac*pM0;         % normalized probability at the minimum
        
    end
        
    % Get the uncertainty/confidence measure values (get Pav - precomputed
    % and stored)
    dirname = strcat(omegadir,eid{ii},'/');
    for jj=1:6
        fname = strcat(dirname,'L1/',Mlist{jj},'/',eid{ii},'_result.dat');
        fid = fopen(fname);
        w = dlmread(fname) % get the confidence measure value
        fclose(fid);
        OMEGA1(jj,ii) = w;
        % for norm 2
        if 0
            fname = strcat(dirname,'L2/',Mlist{jj},'/',eid{ii},'_result.dat');
            fid = fopen(fname);
            w = dlmread(fname); % get the confidence measure value
            fclose(fid);
            OMEGA2(jj,ii) = w{1};
        end
    end
    
    if iAEC
        OMEGA1(jj+1,:) = 0;
        OMEGA1(jj+2,:) = 0;
        OMEGA2(jj+1,:) = 0;
        OMEGA2(jj+2,:) = 0;
        PM0(jj+1,:) = 0;
        PM0(jj+2,:) = 0;
    end
    
    M_vec
    if iwrite==1
        write_psmeca(strcat(datadir,eid{ii},'_L1'),otime(:,ii),lat_vec(:,ii),lon_vec(:,ii),dep_vec(:,ii),M_vec)
        write_psmeca(strcat(datadir,eid{ii},'_L2'),otime(:,ii),lat_vec(:,ii),lon_vec(:,ii),dep_vec(:,ii),M_vec)
    end
    
    [omegadc1(:,ii)] = CMT2omegadc_xi0(repmat(M_vec(:,6),1,nrow),M_vec,0,0);
    if inorm2
        [omegadc2(:,ii)] = CMT2omegadc_xi0(repmat(M_vec2(:,6),1,nrow),M_vec2,0,0);
    end
end

% Latex table
for ii=1:length(dirs)
    fprintf('\\multicolumn{4}{l}{Event %s (%d)} \\\\ \\hline \n',eid{ii},ii);
    for kk=1:nrow
        jj=wt_order(kk);
        fprintf('L1 & %s & %3.0f & %1.1f & %3.0f & %3.2f & %1.2f  \\\\ \n',Mlist2{jj},dep_vec(jj,ii),mag_vec(jj,ii),omegadc1(jj,ii),PM0(jj,ii), OMEGA1(jj,ii));
    end
    fprintf('\\hline \n');
    if 0
        for kk=1:nrow
            jj=wt_order(kk);
            fprintf('L2 & %s & %3.0f & %1.1f & %3.1f & %3.0f & %1.2f \\\\ \n',Mlist2{jj},dep_vec2(jj,ii),mag_vec2(jj,ii),vr_vec2(jj,ii),omegadc2(jj,ii),OMEGA2(jj,ii));
        end
        fprintf('\\hline \n');
    end
end

if iwrite==1
    fid = fopen(strcat(datadir,'Omega_L1'),'w');
    fprintf(fid,'%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\n',omegadc1)
    fclose(fid);
    if inorm2
        fid = fopen(strcat(datadir,'Omega_L2'),'w');
        fprintf(fid,'%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\n',omegadc2)
        fclose(fid);
    end
end

% RESHAPE
vrl1=reshape(vr_vec1,nele,1);
ol1=reshape(OMEGA1,nele,1);
if inorm2
    vrl2=reshape(vr_vec2,nele,1);
    ol2=reshape(OMEGA2,nele,1);
end

% =========
iplot=0;
if iplot
    figure
    subplot(2,1,1)
    plot(dep_vec')
    title('depth');
    subplot(2,1,2)
    plot(mag_vec')
    title('mag');
    if inorm2
        figure
        subplot(2,1,1)
        plot(dep_vec2')
        title('depth');
        subplot(2,1,2)
        plot(mag_vec2')
        title('mag');
    end
    figure(100)
    set(gca,'FontSize',22)
    plot(vrl1,ol1,'.','Markersize',16)
    xlabel('VR')
    ylabel('P_{av}')
    title('L1')
    if inorm2
        figure(101)
        set(gca,'FontSize',22)
        plot(vrl2,ol2,'.','Markersize',15)
        xlabel('VR')
        ylabel('\Omega')
        title('L2')
    end
end

if 0
    figure
    plot(vr_vec1(1,:), OMEGA1(1,:),'.r','MarkerSize',15)    % M011
    hold on
    plot(vr_vec1(2,:), OMEGA1(2,:),'.g','MarkerSize',15)   % M112
    plot(vr_vec1(3,:), OMEGA1(3,:),'.b','MarkerSize',15)    % M012
    plot(vr_vec1(4,:), OMEGA1(4,:),'.c','MarkerSize',15)   % M101
    plot(vr_vec1(5,:), OMEGA1(5,:),'.m','MarkerSize',15)    % M110
    plot(vr_vec1(6,:), OMEGA1(6,:),'.k','MarkerSize',15)    % M111
end
%=====regression
if 0
    N=100;
    xmin = 1;
    xmax = 100;
    xvec = linspace(xmin,xmax,N);
    
    % y = a*exp(-b*x)
    %close all
    amin = 350;
    amax = 350;
    ainc = 10;
    bmin = .12; % .12 (L2)  .08(L1)
    bmax = .12; % .12 (L2)  .08(L1)
    binc = .01;
    for a=amin:ainc:amax
        for b = bmin:binc:bmax
            yvec = (a*exp(-b*xvec))+13;
            %plot(xvec,yvec,'r','Linewidth',2)
            hold on
        end
    end
    
    ylim([0 180])
end

% For making a table of [pmin , domega and conf]
% Latex table
if 0
    for jj=1:21
        fprintf('%d & ',jj)
        for kk=1:6
            ii=wt_order(kk);
            fprintf('%3.2f & ',PM0(ii,jj))
        end
        for kk=1:6
            ii=wt_order(kk);
            fprintf('%3.0f & ',omegadc1(ii,jj))
        end
        for kk=1:6
            ii=wt_order(kk);
            if kk == 6
                fprintf('%2.2f  ',OMEGA1(ii,jj))
            else
            fprintf('%2.2f & ',OMEGA1(ii,jj))
            end
        end
        fprintf('\\\\ \n')
    end
    corrPM0 = corr(log(PM0(1:6,:)'));
    corrOMEGA1 = corr(log(OMEGA1(1:6,:)'));
    corromegadc1 = corr(log(omegadc1(1:6,:)'));
    
    % log difference in probability of different moment tensor subsets
    figure;
    ymax = 4.5; xmin = -3; xmax = 3; ytick = 0:1:5; fontsize = 20;
    
    subplot(3,2,1); 
    %figure;
    set(gca,'FontSize',fontsize)
    plot_histo(log(PM0(6,:))-log(PM0(5,:)),[xmin:0.1:xmax],1)
    hold on; plot([0 0],[0 ymax],'--r','Linewidth',2); box on;
    %title('M_{111} vs M_{110}'); 
    title('(a) Log(pM_{111}/pM_{110})');
    %xlabel('Log(pM_{111}/pM_{110})');
    ylim([0 ymax]); ylabel('Events (N=21)'); 
    set(gca,'yTick',ytick)
    
    subplot(3,2,2); 
    %figure;
    set(gca,'FontSize',fontsize)
    plot_histo(log(PM0(6,:))-log(PM0(1,:)),[xmin:0.1:xmax],1)
    hold on; plot([0 0],[0 ymax],'--r','Linewidth',2); box on;
    %title('M_{111} vs M_{011}'); 
    title('(b) Log(pM_{111}/pM_{011})');
    %xlabel('Log(pM_{111}/pM_{011})');
    ylim([0 ymax]);ylabel('Events (N=21)');
    set(gca,'yTick',ytick)
    
    subplot(3,2,3); 
    %figure;
    set(gca,'FontSize',fontsize)
    plot_histo(log(PM0(6,:))-log(PM0(2,:)),[xmin:0.1:xmax],1)
    hold on; plot([0 0],[0 ymax],'--r','Linewidth',2); box on;
    %title('M_{111} vs M_{112}'); 
    title('(c) Log(pM_{111}/pM_{112})');
    %xlabel('Log(pM_{111}/pM_{112})');
    ylim([0 ymax]); ylabel('Events (N=21)');
    set(gca,'yTick',ytick)
    
    subplot(3,2,4); 
    %figure;
    set(gca,'FontSize',fontsize)
    plot_histo(log(PM0(4,:))-log(PM0(1,:)),[xmin:0.1:xmax],1)
    hold on; plot([0 0],[0 ymax],'--r','Linewidth',2); box on;
    %title('M_{101} vs M_{011}'); 
    title('(d) Log(pM_{101}/pM_{011})');
    %xlabel('Log(pM_{101}/pM_{011})');
    ylim([0 ymax]); ylabel('Events (N=21)');
    set(gca,'yTick',ytick)
    
    subplot(3,2,5); 
    %figure;
    set(gca,'FontSize',fontsize)
    plot_histo(log(PM0(3,:))-log(PM0(1,:)),[xmin:0.1:xmax],1)
    hold on; plot([0 0],[0 ymax],'--r','Linewidth',2); box on;
    %title('M_{012} vs M_{011}'); 
    title('(e) Log(pM_{012}/pM_{011})')
    %xlabel('Log(pM_{012}/pM_{011})');
    ylim([0 ymax]); ylabel('Events (N=21)');
    set(gca,'yTick',ytick)
end


% For making a table of [omegadc]
% Updated latex table
if 0
    for jj=1:21
        fprintf('%s (%d) & ',eid1{jj},jj)
        for kk=1:8
            ii=wt_order(kk);
            fprintf('%3.0f & ',omegadc1(ii,jj))
        end
        fprintf('\\\\ \n')
    end
end












