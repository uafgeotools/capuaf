% This script reads summary files and plots variation of uncertainty with
% VR
%
% Perhaps put it in mtvipul.m (for southern Alaska moment tensor paper)
%
% Vipul Silwal
% Sept 7, 2015

clear all
close all
ddir = '/home/vipul/CAP/inv/scak/MOOS/RESULTS/';
events = dir(strcat(ddir,'2*'));
norm = {'L1' 'L2'};
norm = {'L1'};
    
Mlist = {'011','112','012','101','110','111'};
cindx = {'r','g','b','c','m','k'};
%Mlist = {'111'}
gmtdir = '/home/vipul/gmt/data/cap/MOOS/';
K = 40;
N =1e5;
dV = (4*pi*pi)/N;       % Volume of the brick

VR = [];
conf = [];
PMall = [];
for ii=1:length(events)
    eid = events(ii).name;
    fname = strcat(ddir,eid,'/summary.dat');
    fid = fopen(fname);
    w = textscan(fid,'%d %d %d %f %d %f %d %f');
    fclose(fid);
    nor = w{1};
    wt = w{2};
    stk = w{3};
    dip = w{4};
    rak = w{5};
    mag = w{6};
    dep = w{7};
    vr = w{8};
    %vr = [vr; w{8}];
    for jj = 1:length(norm)
        for kk = 1:length(Mlist)
            VR = [VR; vr(2*kk-jj)];
            LogDir = strcat(ddir,eid,'/L1/M',Mlist{kk},'/');
            f = dir(strcat(LogDir,'*rand.mat'));
            fname =strcat(LogDir,f.name);
            [kappa,theta,sigma,F,Mw,delta,gamma,M0,M] = read_CAP_log(fname);
            
            Fsum = sum(exp(-F*K));
            [fmin,indx1] = min(F);
            Mmin = M(:,indx1);
            
            pM0 = exp(-fmin*K);     % unnormalized probability at the minimum
            kfac = 1/(dV*Fsum);     % normalization factor (Walt2015b notes)
            
            PM0 = kfac*pM0;         % normalized probability at the minimum
            PMall = [PMall; PM0];
            
            gdir = strcat(gmtdir,eid,'/',norm{jj},'/M',Mlist{kk},'/');
            omega_file = strcat(gdir,eid,'_result.dat');
            fid = fopen(omega_file);
            o = textscan(fid,'%f');
            conf = [conf; o{1}];
            
            figure(1);
            semilogx(PM0,vr(2*kk-jj),strcat('.',cindx{kk}),'MarkerSize',15)
            %plot(PM0,vr(2*kk-jj),strcat('.',cindx{kk}),'MarkerSize',15)
            hold on
            
            figure(2);
            semilogx(PM0,o{1},strcat('.',cindx{kk}),'MarkerSize',15);
            %plot(PM0,o{1},strcat('.',cindx{kk}),'MarkerSize',15);
            hold on
        end
    end
end

figure(1);
hold on
semilogx([12.29],[58.9],'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',15); % Denali
semilogx([1.77],[16.8],'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',15); % Noatak
%plot([12.29],[58.9],'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',15); % Denali
%plot([1.77],[16.8],'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',15); % Noatak
set(gca,'FontSize',18)
xlabel('p(M_0)','FontSize',18);
ylabel('VR','FontSize',18);

figure(2);
semilogx(12.29,0.98,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',15); % Denali
text(20,0.95,'Denali')
semilogx(1.77,0.61,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',15); % Noatak
text(1.77,0.61,'Noatak')
%plot(12.29,0.98,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',15); % Denali
%plot(1.77,0.61,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',15); % Noatak
set(gca,'FontSize',18)
xlabel('p(M_0)','FontSize',18);
ylabel('P_{av}','FontSize',18);


break
% Plot omega vs VR
set(gca,'FontSize',18)
plot(VR,conf,'.','MarkerSize',15)
xlim([0 100]);
ylim([0.5 1]);
xlabel('VR','FontSize',18);
ylabel('P_{av}','FontSize',18);