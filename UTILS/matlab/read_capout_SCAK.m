close all
clear all
cap_dir = '/home/vipul/CAP/inv/scak/MOOS/RESULTS/';
%cap_dir = '/home/vipul/CAP/inv/scak/SCAK/';
events=dir(strcat(cap_dir,'200*'));
MTtype = 'M110';

stnm_all=[];
lampPV_all=[];
lampPR_all=[];
lampSV_all=[];
lampSR_all=[];
lampSH_all=[];
corrPV_all=[];
corrPR_all=[];
corrSV_all=[];
corrSR_all=[];
corrSH_all=[];
for ii=1:length(events)
    edir = events(ii).name;
    outdir = strcat(cap_dir,edir,'/L1/',MTtype,'/');
    outfilename = dir(strcat(outdir,'2*.out'));
    outfile = strcat(outdir,outfilename.name);
    [otime(ii),elat(ii),elon(ii),edep(ii),strike(ii),dip(ii),rake(ii),M,Mw(ii),eid(ii),...
        capmod,capdep(ii),rms(ii),vr(ii),Pwin(ii),Swin(ii),Nstn(ii),Pstn(ii),Sstn(ii),...
        stnm,lampPV,lampPR,lampSV,lampSR,lampSH,corrPV,corrPR,corrSV,corrSR,corrSH]= read_capout_old(outfile);
    stnm_all=horzcat(stnm_all,stnm);
    lampPV_all=horzcat(lampPV_all,lampPV);
    lampPR_all=horzcat(lampPR_all,lampPR);
    lampSV_all=horzcat(lampSV_all,lampSV);
    lampSR_all=horzcat(lampSR_all,lampSR);
    lampSH_all=horzcat(lampSH_all,lampSH);
    corrPV_all=horzcat(corrPV_all,corrPV);
    corrPR_all=horzcat(corrPR_all,corrPR);
    corrSV_all=horzcat(corrSV_all,corrSV);
    corrSR_all=horzcat(corrSR_all,corrSR);
    corrSH_all=horzcat(corrSH_all,corrSH);
    %[MDC,strike1(ii),dip1(ii),rake1(ii),strike2(ii),dip2(ii),rake2(ii),k1,d1,n1,k2,d2,n2,U,lams]= CMT2dcfaultpar(M,0);
end

figure
bin=.5;
subplot(3,2,1)
plot_histo(lampPV_all,-10:bin:10,2,1)
title('PV')
subplot(3,2,2)
plot_histo(lampPR_all,-10:bin:10,2,1)
title('PR')
subplot(3,2,3)
plot_histo(lampSV_all,-10:bin:10,2,1)
title('SR')
subplot(3,2,4)
plot_histo(lampSR_all,-10:bin:10,2,1)
title('SV')
subplot(3,2,5)
plot_histo(lampSH_all,-10:bin:10,2,1)
title('SH')
suptitle(sprintf('%s - Best 21 event; Ln(data/syn)',MTtype));

figure
bin=5;
subplot(3,2,1)
plot_histo(corrPV_all,-0:bin:100,2,1)
ylim([0 0.15])
title('PV')
subplot(3,2,2)
plot_histo(corrPR_all,0:bin:100,2,1)
ylim([0 0.15])
title('PR')
subplot(3,2,3)
plot_histo(corrSV_all,0:bin:100,2,1)
ylim([0 0.3])
title('SR')
subplot(3,2,4)
plot_histo(corrSR_all,0:bin:100,2,1)
ylim([0 0.3])
title('SV')
subplot(3,2,5)
plot_histo(corrSH_all,0:bin:100,2,1)
ylim([0 0.3])
title('SH')
suptitle(sprintf('%s - Best 21 event; cross-corr',MTtype));

%=================================================================
[C,iA,iC]=unique(stnm_all);
figure(100);
hold on
[N,Nplot,centers]=plot_histo((iC),1:1:length(C)+1,1,1);

lampPVmat = zeros(length(C),21);
lampPRmat = zeros(length(C),21);
lampSVmat = zeros(length(C),21);
lampSRmat = zeros(length(C),21);
lampSHmat = zeros(length(C),21);
corrPVmat = zeros(length(C),21);
corrPRmat = zeros(length(C),21);
corrSVmat = zeros(length(C),21);
corrSRmat = zeros(length(C),21);
corrSHmat = zeros(length(C),21);

indx = [1;cumsum(N)];
for ii=1:length(C)
    k=1;
    for jj=1:length(stnm_all)
        if strcmp(C(ii),stnm_all(jj))
           lampPVmat(ii,k)=lampPV_all(jj);
           lampPRmat(ii,k)=lampPR_all(jj);
           lampSVmat(ii,k)=lampSV_all(jj);
           lampSRmat(ii,k)=lampSR_all(jj);
           lampSHmat(ii,k)=lampSH_all(jj);
           corrPVmat(ii,k)=corrPV_all(jj);
           corrPRmat(ii,k)=corrPR_all(jj);
           corrSVmat(ii,k)=corrSV_all(jj);
           corrSRmat(ii,k)=corrSR_all(jj);
           corrSHmat(ii,k)=corrSH_all(jj);
           k = k+1;
        end
    end
    lampPV_avg(ii,1) = ii;
    lampPV_avg(ii,2) = sum(lampPVmat(ii,:));
    lampPR_avg(ii,1) = ii;
    lampPR_avg(ii,2) = sum(lampPRmat(ii,:));
    lampSV_avg(ii,1) = ii;
    lampSV_avg(ii,2) = sum(lampSVmat(ii,:));
    lampSR_avg(ii,1) = ii;
    lampSR_avg(ii,2) = sum(lampSRmat(ii,:));
    lampSH_avg(ii,1) = ii;
    lampSH_avg(ii,2) = sum(lampSHmat(ii,:));
    corrPV_avg(ii,1) = sum(corrPVmat(ii,:));
    corrPR_avg(ii,1) = sum(corrPRmat(ii,:));
    corrSV_avg(ii,1) = sum(corrSVmat(ii,:));
    corrSR_avg(ii,1) = sum(corrSRmat(ii,:));
    corrSH_avg(ii,1) = sum(corrSHmat(ii,:));
end

lampPV_avg(:,2) = lampPV_avg(:,2)./N;
lampPR_avg(:,2) = lampPR_avg(:,2)./N;
lampSV_avg(:,2) = lampSV_avg(:,2)./N;
lampSR_avg(:,2) = lampSR_avg(:,2)./N;
lampSH_avg(:,2) = lampSH_avg(:,2)./N;

figure;
subplot(3,2,1)
plot(lampPV_avg(:,2),'.-')
xlim([0 105])
grid on
subplot(3,2,2)
plot(lampPR_avg(:,2),'.-')
xlim([0 105])
grid on
subplot(3,2,3)
plot(lampSV_avg(:,2),'.-')
xlim([0 105])
grid on
subplot(3,2,4)
plot(lampSR_avg(:,2),'.-')
xlim([0 105])
grid on
subplot(3,2,5)
plot(lampSH_avg(:,2),'.-')
xlim([0 105])
grid on

figure
subplot(3,2,1)
plot(corrPV_avg./N,'.-')
xlim([0 105])
subplot(3,2,2)
plot(corrPR_avg./N,'.-')
xlim([0 105])
subplot(3,2,3)
plot(corrSV_avg./N,'.-')
xlim([0 105])
subplot(3,2,4)
plot(corrSR_avg./N,'.-')
xlim([0 105])
subplot(3,2,5)
plot(corrSH_avg./N,'.-')
xlim([0 105])


% [~,isort] = sort(capdep,'ascend');
% % Output in different formats
% for jj=1:length(events)
%     ii=isort(jj);
%     % Output for Latex
%     %fprintf('%3.0f & %s & %3.2f & %3.2f & %3.0f & %3.0f & %3.0f & %2.1f & %3.0f & %3.0f & %3.0f & %3.0f\\\\ \n',...
%     %    jj,eid{ii},elat(ii),elon(ii),strike1(ii),dip1(ii),rake1(ii),Mw(ii),capdep(ii),Pstn(ii),Sstn(ii),Nstn(ii))
%     %fprintf(' & & & & %3.0f & %3.0f & %3.0f &  &  &  &  & \\\\ \n',strike2(ii),dip2(ii),rake2(ii))
%     % Output for GMT psmeca plotting (use -Sa flag)
%     fprintf('%3.2f \t %3.2f \t %3.0f \t %3.0f \t %3.0f \t %3.0f \t %2.1f \t %3.2f \t %3.2f \n',...
%         elon(ii),elat(ii),capdep(ii),strike1(ii),dip1(ii),rake1(ii),Mw(ii),elon(ii),elat(ii))
% end


if 1==0
    plotdir = '/home/vipul/gmt/data/lamp/';
    fid = fopen(strcat(plotdir,'lampPV.dat'),'w');
    fprintf(fid,'%d %3.3f\n',lampPV_avg');
    fclose(fid);
    fid = fopen(strcat(plotdir,'lampPR.dat'),'w');
    fprintf(fid,'%d %3.3f\n',lampPR_avg');
    fclose(fid);
    fid = fopen(strcat(plotdir,'lampSV.dat'),'w');
    fprintf(fid,'%d %3.3f\n',lampSV_avg');
    fclose(fid);
    fid = fopen(strcat(plotdir,'lampSR.dat'),'w');
    fprintf(fid,'%d %3.3f\n',lampSR_avg');
    fclose(fid);
    fid = fopen(strcat(plotdir,'lampSH.dat'),'w');
    fprintf(fid,'%d %3.3f\n',lampSH_avg');
    fclose(fid);
end
