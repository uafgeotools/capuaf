clear all
close all

% Get CAP results for:
%---------------------------------------------------------------------
% 21 AEC moment tensor events - (2007-2009)
%ddir = '/home/vipul/CAP/inv/scak/MOOS/RESULTS/L1_M111_out/'; 
% 106 SilwalTape2016 events (85 AECfp + 21 AECMT) - (2007-2009)
ddir = '/home/vipul/CAP/inv/scak/SCAK/RESULTS/out_full/';
ax2 = [wrapTo360(-154) wrapTo360(-146) 58 62.5];   % SilwalTape2016
oran = [datenum(2007,8,1) datenum(2009,8,15)];
ioutformat = 0;
%---------------------------------------------------------------------
% 72 AEC moment tensors (2007-2014) - for aktomo
ddir = '/home/vipul/CAP/inv/aktomo/scak/RESULTS/out/';
ioutformat = 0; 

%---------------------------------------------------------------------
% 90 Josh's events here (2 events are not in AECfp catalog - Iniskin aftershocks)
ddir = '/home/jcpurba/PROJECTS/CAP/inv/scak/out_scak_all/';
ax2 = [wrapTo360(-154) wrapTo360(-139) 59 62.5];   % aktomo
oran = [datenum(2009,8,12) datenum(2016,6,1)];
ioutformat = 1;

q = dir(strcat(ddir,'*out'));
for ii=1:length(q)
    if ioutformat == 0
        [otime2(ii,1),elat(ii,1),elon(ii,1),edep(ii,1),strike(ii,1),dip(ii,1),rake(ii,1),M1(:,ii),Mw1(ii,1),...
            eid1(ii,1),capmod,capdep(ii,1),rms,vr,Pwin,Swin,Nstn(ii,1),Pstn,Sstn]=...
            read_capout_old(strcat(ddir,q(ii).name));
    else
        [otime2(ii,1),elat(ii,1),elon(ii,1),edep(ii,1),strike(ii,1),dip(ii,1),rake(ii,1),M1(:,ii),Mw1(ii,1),...
            eid1(ii,1),capmod,capdep(ii,1),rms,vr,Pwin,Swin,Nstn(ii,1),Pstn,Sstn,...
            stnm,lampPV,lampPR,lampSV,lampSR,lampSH,corrPV,corrPR,corrSV,corrSR,corrSH,pPol,ePol,shiftSavg,shiftPavg] = read_capout(strcat(ddir,q(ii).name));
    end
end


% Get AECfp solutions
ax3 = [ax2 -10 700];
Mwran = [3.5 10];
[otime,slat,slon,sdep,M2,M0,Mw2,eid2] = read_mech_AECfp(oran,ax3,Mwran);
axis(ax2)

% Get AECMT solutions
if 1==1 % set equal to 0 if (you want to compute omega between fault plane and SCAK events (cap )
    Mwran = [0 10];
    [otime,slat,slon,sdep,M3,M0,Mw3,eid3,depc] = read_mech_AEC(oran,ax3,Mwran);
end

% Compute Omega distance between any 2 from above
% difference between AECfp and AEC moment tensor solutions
[c,A,B]=intersect(eid2,eid3);
Mi1 = M2(:,A);
Mi2 = M3(:,B);
omega1 = CMT2omega(Mi1,Mi2);

% difference between CAP and AEC moment tensor solutions
[c,A,B]=intersect(eid1,eid3);
Mi1 = M1(:,A);
Mi2 = M3(:,B);
omega2 = CMT2omega(Mi1,Mi2);

% difference between CAP and AECfp solutions
[c,A,B]=intersect(eid1,eid2);
Mi1 = M1(:,A);
Mi2 = M2(:,B);
omega3 = CMT2omega(Mi1,Mi2);

% Same as above figure (generate by CMT2omega) but wiht slight
% modifications (figures for papers or posters)
if 1
    figure
    set(gca,'FontSize',22)
    plot_histo(omega1,0:10:180,1,1)
    title('AEC fault-plane vs AEC moment tensor solutions','FontSize', 15)
    xlabel('\omega');
    
    figure
    set(gca,'FontSize',22)
    plot_histo(omega2,0:10:180,1,1)
    title('CAP vs AEC moment tensor solutions','FontSize', 15)
    xlabel('\omega');
    
    figure
    set(gca,'FontSize',22)
    plot_histo(omega3,0:10:180,1,1)
    title('CAP vs AEC fault-plane solutions','FontSize', 15)
    xlabel('\omega');
end

% average difference in mag
if  1
    [c,A,B]=intersect(eid1,eid2);
    Mi1 = Mw1(A);
    Mi2 = Mw2(B);
    dMw = Mi1 - Mi2;
    ddMw = std(dMw);
    figure
    box on
    hold on
    set(gca,'FontSize',22)
    plot_histo(dMw,-0.5:.10:1,1,1)
    %plot([mean(dMw) mean(dMw)],[0 35],'--r','Linewidth',3)
    plot([0 0],[0 35],'--r','Linewidth',3)
    xlabel('m_{CAP} - m_{AEC}');
    % set(gca,'FontSize',22)
    % hold on
    % box on
    % plot(dMw,'.','Markersize',14);
    % plot([0 100],[mean(dMw) mean(dMw)],'--','Linewidth',3)
    % plot([0 100],[mean(dMw)+ddMw mean(dMw)+ddMw],'--r','Linewidth',3);
    % plot([0 100],[mean(dMw)-ddMw mean(dMw)-ddMw],'--r','Linewidth',3);
    % xlim([0 90]);
    % xlabel('Events');
    % ylabel('\Delta Mw');
end

if 0
    [a,b] = sort(omega);
    eid1(b)'
    a
    Mw(b)'
end