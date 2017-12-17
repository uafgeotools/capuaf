% Figure omega_dc vs Delta(VR) figure (main paper)
if (1)
    % Read CAP solutions
    clear all
    ddir = '/home/vipul/CAP/inv/scak/MOOS/RESULTS/L1_M111_out/';    % 21 main events (CAP solutions)
    ddir = '/home/vipul/CAP/inv/scak/MOOS/RESULTS04/out/';
    q = dir(strcat(ddir,'*out'));
    for ii=1:length(q)
        [otime(ii,1),elat(ii,1),elon(ii,1),edep(ii,1),strike(ii,1),dip(ii,1),rake(ii,1),M1(:,ii),Mw(ii,1),...
            eid(ii,1),capmod,capdep(ii,1),rms1(ii,1),vr1(ii),Pwin,Swin,Nstn(ii,1),Pstn,Sstn]=...
            read_capout_old(strcat(ddir,q(ii).name));
    end
    
    % Read AEC MT solutions
    ddir = '/home/vipul/CAP/inv/scak/MOOS/RESULTS_AEC/AECMT/out/';  % Solutions at AEC MT and searching over Mw
    q = dir(strcat(ddir,'*out'));
    for ii=1:length(q)
        [otime(ii,1),elat(ii,1),elon(ii,1),edep(ii,1),strike(ii,1),dip(ii,1),rake(ii,1),M2(:,ii),Mw(ii,1),...
            eid(ii,1),capmod,capdep(ii,1),rms2(ii,1),vr2(ii),Pwin,Swin,Nstn(ii,1),Pstn,Sstn]=...
            read_capout_old(strcat(ddir,q(ii).name));
    end
    
    ddir = '/home/vipul/CAP/inv/scak/MOOS/RESULTS_AEC/AECMT_Mw/out/';  % Solutions at AEC MT and fixed Mw
    q = dir(strcat(ddir,'*out'));
    for ii=1:length(q)
        [otime(ii,1),elat(ii,1),elon(ii,1),edep(ii,1),strike(ii,1),dip(ii,1),rake(ii,1),M2(:,ii),Mw(ii,1),...
            eid(ii,1),capmod,capdep(ii,1),rms3(ii,1),vr3(ii),Pwin,Swin,Nstn(ii,1),Pstn,Sstn]=...
            read_capout_old(strcat(ddir,q(ii).name));
    end
    
    omega = CMT2omega(M1,M2);
    dvr1 = vr1 - vr2;
    dvr2 = vr1 - vr3;
    dvr1 = log(rms2./rms1);
    dvr2 = log(rms3./rms1);
    
    figure(1)
    hold on
    box on
    set(gca,'FontSize',22)
    %plot(dvr1,omega,'.','Markersize',18);
    %plot(dvr2,omega,'.g','Markersize',18);
    xlabel('Ln(pCAP/pAEC)')
    ylabel('\omega')
    %xlim([0 32])
    
    figure(2)
    hold on
    box on
    set(gca,'FontSize',22)
    %plot(dvr1,omega,'.','Markersize',18);
    %plot(dvr2,omega,'.g','Markersize',18);
    xlabel('Ln(pCAP/pAEC)')
    ylabel('\omega')
    %xlim([0 32])

    %[xf,yf,a,b] = linfit(dvr2',omega);
    for ii=1:length(vr1)
        if ii ==3
            figure(1)
            plot([dvr1(ii) dvr2(ii)],[omega(ii) omega(ii)],'k--','Linewidth',2)
            figure(2)
            semilogx([dvr1(ii) dvr2(ii)],[omega(ii) omega(ii)],'k--','Linewidth',2)
        else
            figure(1)
            plot([dvr1(ii) dvr2(ii)],[omega(ii) omega(ii)],'k-','Linewidth',2)
            figure(2)
            semilogx([dvr1(ii) dvr2(ii)],[omega(ii) omega(ii)],'k-','Linewidth',2)
        end
    end
    figure(1)
    plot(dvr2,omega,'.r','Markersize',18);
    plot(dvr1,omega,'.','Markersize',18);
    %plot(dvr2,omega,'.g','Markersize',18);
    %plot(xf,yf,'k','Linewidth',2)
end

%% Pav vs omega figures for all events - in mtvipul/extra.tex
if (0)
    eidlist={'20070911234634153' '20070919112226549' '20071003140612444' '20071010180326301' '20071128235703849'...
        '20080314093821771' '20080327230745201' '20080828231418631' '20080918194353069' '20081228071310738'...
        '20090124180950811' '20090215193500098' '20090223000427175' '20090317011333066'...
        '20090414171427415' '20090430045457938' '20090524094004552' '20090622192805162' '20090626164820729' '20090730223910267'};
    for jj = 1:length(eidlist)
        eid = eidlist{jj};
        ddir = strcat('/home/vipul/CAP/inv/scak/MOOS/RESULTS/',eid,'/L1/M111/');
        %ddir = strcat('/home/vipul/CAP/inv/scak/NW/',eid,'/L1/M111/');
        f = dir(strcat(ddir,'*rand.mat'));
        fname =strcat(ddir,f.name);
        [kappa,theta,sigma,F,Mw,delta,gamma,M0,M] = read_CAP_log(fname); % read rand log files for misfit values
        f = dir(strcat(ddir,'*rand'));
        fname =strcat(ddir,f.name);
        
        [fmin,indx1] = min(F);
        Mmin = M(:,indx1);
        [fmax,indx2] = max(F);
        Mmax = M(:,indx2);
        Mref = Mmin;
        
        Nref = 600;
        ivec = 1:Nref;
        Mvec = M(:,ivec);
        Fref = F(ivec);
        om_ref = CMT2omega(Mvec,Mref);
        F_samp= F(ivec);
        VR = 100*(1-F.*F);
        VR_samp = VR(ivec);
        
        % For running the uncertainty analysis (estimating Pav) AND saving
        % the results in *.mat file
        if 0
            for jj = 1:length(Mvec)
                [~,~, Pav(jj,1)]= CAP_unc(fname,Mvec(:,jj));
                %close all
            end

            % to save the omegadc and Pav values (CAUTIOUS!! will overwrite)
            save(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmega.mat'),'om_ref','Pav')
            % save(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmegaVR.mat'),'om_ref','Pav','VR_samp')
            % save(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmegaVR3.mat'),'om_ref','Pav','VR_samp')
        end
    end 
        
        % Now plot all of them
        % For figures in extra (20 events (one is in supplement). 
        close all
        set(gca,'FontSize',22)
        for ii = 1:length(eidlist)
            eid = eidlist{ii};
            PO = load(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmega.mat'));
            PavMax = dlmread(strcat('/home/vipul/gmt/data/cap/MOOS/L1_M111_unc/',eid,'_result.dat'));
            Nsubplots = 10;
            ifig = fix(ii/Nsubplots)+1;
            isubplot = mod(ii,Nsubplots);
            if isubplot==0;
                ifig = ifig-1;
                isubplot=Nsubplots;
            end
            figure(ifig)
            subplot(5,2,isubplot)
            plot(PO.om_ref,PO.Pav,'.','Markersize',15);
            grid on
            xtickValues = 0:90:180;
            ytickValues = 0:0.5:1;
            set(gca,'XTick',xtickValues)
            set(gca,'YTick',ytickValues)
            %subplot(4,2,isubplot)
            %plot(PO.om_ref,PO.Pav,'.','Markersize',15);
            xlim([0 180]); ylim([0 1]);
            xlabel('\omega');
            ylabel('P_{av}');
            title(sprintf('%s; max(P_{av}) = %1.2f',eid,PavMax));
        end
end
%% Pav vs omegadc with changing Mref (omega = [0 11 50 175]) - in mtvipul/extra.tex
if (0)
    eidlist = {'20090407201255351'};
    for jj = 1:length(eidlist)
        eid = eidlist{jj};
        ddir = strcat('/home/vipul/CAP/inv/scak/MOOS/RESULTS/',eid,'/L1/M111/');
        %ddir = strcat('/home/vipul/CAP/inv/scak/NW/',eid,'/L1/M111/');
        f = dir(strcat(ddir,'*rand.mat'));
        fname =strcat(ddir,f.name);
        % read rand log files for misfit values
        [kappa,theta,sigma,F,Mw,delta,gamma,M0,M] = read_CAP_log(fname); 
        f = dir(strcat(ddir,'*rand'));
        fname =strcat(ddir,f.name);
        
        [fmin,indx1] = min(F);
        Mmin = M(:,indx1);
        [fmax,indx2] = max(F);
        Mmax = M(:,indx2);
        indx3 = 3831; % -M0 (completely opposite to minimum misfit solution)
        Mref = Mmin;
        isamp_vec = [indx1 indx2 indx3];
        Msamp = M(:,isamp_vec); % for samples on top
        
        
        Nsampl = 600;  % Number of sample points (kind of meaningless since we are using presaved files)
        ivec = 1:Nsampl;
        Mvec = M(:,ivec);
        Fref = F(ivec);
        om_ref = CMT2omega(Mvec,Mref);
        F_samp= F(ivec);
        VR = 100*(1-F.*F);
        VR_samp = VR(ivec);
        % for samples on top
        om_samp = CMT2omega(Msamp,Mref);
        
        % For running the uncertainty analysis (estimating Pav) AND saving
        % the results in *.mat file
        % BUT SEE ABOVE (ALREADY PRESAVED FOR ALL)
        if 0
            for jj = 1:length(Mvec)
                [~,~, Pav(jj,1)]= CAP_unc(fname,Mvec(:,jj));
                %close all
            end
            
            % to save the omegadc and Pav values (CAUTIOUS!! will overwrite)
            save(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmega.mat'),'om_ref','Pav')
            % save(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmegaVR.mat'),'om_ref','Pav','VR_samp')
            % save(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmegaVR3.mat'),'om_ref','Pav','VR_samp')
        end
        
         if 1
            % get the 3 points to plot on top
            a = load(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmegaVR.mat'));
            b = load(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmegaVR3.mat'));
            
            om_all = a.om_all; VR_all = a.VR_all; Pav_all = a.Pav_all;
            om3 = b.om3; VR3 = b.VR3; Pav3 = b.Pav3;
            om_all = om_ref;
            om3 = om_samp;
            
            figure(1); 
            plot(om_all(1:Nsampl),Pav_all(1:Nsampl),'.','Markersize',15)
            set(gca,'FontSize',22); xlabel('\omega'); ylabel('P_{av}');
            hold on
            plot(om3(1),Pav3(1),'o','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'Markersize',11);
            plot(om3(2),Pav3(2),'o','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'Markersize',11);
            plot(om3(3),Pav3(3),'o','MarkerEdgeColor','k','MarkerFaceColor',[.5 0.5 .5],'Markersize',11);
            xlim([0 180]); ylim([0 1]);
            xtickValues = 0:30:180;
            ytickValues = 0:0.2:1;
            set(gca,'XTick',xtickValues)
            set(gca,'YTick',ytickValues)
            
            ref_indx = 208;       % found by sort(omega(M,Mmin))  
            Mref = M(:,ref_indx);
            om_all = CMT2omega(Mvec,Mref);
            om3 = CMT2omega(Msamp,Mref);
            figure(2); 
            plot(om_all(1:Nsampl),Pav_all(1:Nsampl),'.','Markersize',15)
            set(gca,'FontSize',22); xlabel('\omega'); ylabel('P_{av}');
            hold on
            plot(om3(1),Pav3(1),'o','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'Markersize',11);
            plot(om3(2),Pav3(2),'o','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'Markersize',11);
            plot(om3(3),Pav3(3),'o','MarkerEdgeColor','k','MarkerFaceColor',[.5 0.5 .5],'Markersize',11);
            xlim([0 180]); ylim([0 1]);
            xtickValues = 0:30:180;
            ytickValues = 0:0.2:1;
            set(gca,'XTick',xtickValues)
            set(gca,'YTick',ytickValues)
            
            ref_indx = 321;         % found by sort(omega(M,Mmin))
            Mref = M(:,ref_indx);
            om_all = CMT2omega(Mvec,Mref);
            om3 = CMT2omega(Msamp,Mref);
            figure(3); 
            plot(om_all(1:Nsampl),Pav_all(1:Nsampl),'.','Markersize',15)
            set(gca,'FontSize',22); xlabel('\omega'); ylabel('P_{av}');
            hold on
            plot(om3(1),Pav3(1),'o','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'Markersize',11);
            plot(om3(2),Pav3(2),'o','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'Markersize',11);
            plot(om3(3),Pav3(3),'o','MarkerEdgeColor','k','MarkerFaceColor',[.5 0.5 .5],'Markersize',11);
            xlim([0 180]); ylim([0 1]);
            xtickValues = 0:30:180;
            ytickValues = 0:0.2:1;
            set(gca,'XTick',xtickValues)
            set(gca,'YTick',ytickValues)
            
            ref_indx = 3831;         % found by sort(omega(M,Mmin))
            Mref = M(:,ref_indx);
            om_all = CMT2omega(Mvec,Mref);
            om3 = CMT2omega(Msamp,Mref);
            figure(4); 
            plot(om_all(1:Nsampl),Pav_all(1:Nsampl),'.','Markersize',15)
            set(gca,'FontSize',22); xlabel('\omega'); ylabel('P_{av}');
            hold on
            plot(om3(1),Pav3(1),'o','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'Markersize',11);
            plot(om3(2),Pav3(2),'o','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'Markersize',11);
            plot(om3(3),Pav3(3),'o','MarkerEdgeColor','k','MarkerFaceColor',[.5 0.5 .5],'Markersize',11);
            xlim([0 180]); ylim([0 1]);
            xtickValues = 0:30:180;
            ytickValues = 0:0.2:1;
            set(gca,'XTick',xtickValues)
            set(gca,'YTick',ytickValues)
        end
    end
end

%% misfit vs omegadc; VR vs omegadc; Pav vs omegadc; Pav vs VR figures - in mtvipul/supp.tex
if (0)
    eidlist = {'20090407201255351'};
    
    for jj = 1:length(eidlist)
        eid = eidlist{jj};
        ddir = strcat('/home/vipul/CAP/inv/scak/MOOS/RESULTS/',eid,'/L1/M111/');
        %ddir = strcat('/home/vipul/CAP/inv/scak/NW/',eid,'/L1/M111/');
        f = dir(strcat(ddir,'*rand.mat'));
        fname =strcat(ddir,f.name);
        % read rand log files for misfit values
        [kappa,theta,sigma,F,Mw,delta,gamma,M0,M] = read_CAP_log(fname); 
        f = dir(strcat(ddir,'*rand'));
        fname =strcat(ddir,f.name);
        
        [fmin,indx1] = min(F);
        Mmin = M(:,indx1);
        [fmax,indx2] = max(F);
        Mmax = M(:,indx2);
        Mref = Mmin;
        indx3 = 3831;
        
        Nsampl = 600;  % Number of sample points (kind of meaningless since we are using presaved files)
        ivec = 1:Nsampl;
        Mvec = M(:,ivec);
        Fref = F(ivec);
        om_ref = CMT2omega(Mvec,Mref);
        F_samp= F(ivec);
        VR = 100*(1-F.*F);
        VR_samp = VR(ivec);
        
        % For running the uncertainty analysis (estimating Pav) AND saving
        % the results in *.mat file
        % BUT SEE ABOVE (ALREADY PRESAVED FOR ALL)
        if 0
            for jj = 1:length(Mvec)
                [~,~, Pav(jj,1)]= CAP_unc(fname,Mvec(:,jj));
                %close all
            end
            
            % to save the omegadc and Pav values (CAUTIOUS!! will overwrite)
            save(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmega.mat'),'om_ref','Pav')
            % save(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmegaVR.mat'),'om_ref','Pav','VR_samp')
            % save(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmegaVR3.mat'),'om_ref','Pav','VR_samp')
        end
        
         if 1
            % get the 3 points to plot on top
            a = load(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmegaVR.mat'));
            b = load(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmegaVR3.mat'));
            
            om_all = a.om_all; VR_all = a.VR_all; Pav_all = a.Pav_all;
            om3 = b.om3; VR3 = b.VR3; Pav3 = b.Pav3;
            
            figure(1); 
            plot(om_all(1:Nsampl),VR_all(1:Nsampl),'.','Markersize',15)
            set(gca,'FontSize',22); xlabel('\omega'); ylabel('VR');
            hold on
            plot(om3(1),VR3(1),'o','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'Markersize',11);
            plot(om3(2),VR3(2),'o','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'Markersize',11);
            plot(om3(3),VR3(3),'o','MarkerEdgeColor','k','MarkerFaceColor',[.5 0.5 .5],'Markersize',11);
            %xlim([0 180]); ylim([-20 80]);
            xtickValues = 0:30:180;
            ytickValues = -20:20:80;
            set(gca,'XTick',xtickValues);xlim([0 180]);
            set(gca,'YTick',ytickValues);ylim([-20 80]);
            
            figure(2); 
            plot(om_all(1:Nsampl),Pav_all(1:Nsampl),'.','Markersize',15)
            set(gca,'FontSize',22); xlabel('\omega'); ylabel('P_{av}');
            hold on
            plot(om3(1),Pav3(1),'o','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'Markersize',11);
            plot(om3(2),Pav3(2),'o','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'Markersize',11);
            plot(om3(3),Pav3(3),'o','MarkerEdgeColor','k','MarkerFaceColor',[.5 0.5 .5],'Markersize',11);
            xlim([0 180]); ylim([0 1]);
            xtickValues = 0:30:180;
            ytickValues = 0:0.2:1;
            set(gca,'XTick',xtickValues);xlim([0 180]);
            set(gca,'YTick',ytickValues)
            
            figure(3); 
            plot(VR_all(1:Nsampl),Pav_all(1:Nsampl),'.','Markersize',15)
            set(gca,'FontSize',22); xlabel('VR'); ylabel('P_{av}');
            hold on
            plot(VR3(1),Pav3(1),'o','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'Markersize',11);
            plot(VR3(2),Pav3(2),'o','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'Markersize',11);
            plot(VR3(3),Pav3(3),'o','MarkerEdgeColor','k','MarkerFaceColor',[.5 0.5 .5],'Markersize',11);
            xlim([-20 80]); ylim([0 1]);
            xtickValues = -20:20:80;
            ytickValues = 0:0.2:1;
            set(gca,'XTick',xtickValues);
            set(gca,'YTick',ytickValues)
            
            figure(4); 
            isamp_vec = [indx1 indx2 indx3];
            F3 =40*F(isamp_vec);
            plot(om_all(1:Nsampl),40*F(1:Nsampl),'.','Markersize',15)
            set(gca,'FontSize',22); xlabel('\omega'); ylabel('misfit');
            hold on
            plot(om3(1),F3(1),'o','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'Markersize',11);
            plot(om3(2),F3(2),'o','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'Markersize',11);
            plot(om3(3),F3(3),'o','MarkerEdgeColor','k','MarkerFaceColor',[.5 0.5 .5],'Markersize',11);
            %xlim([0 180]); ylim([0 1]);
            xtickValues = 0:30:180;
            %ytickValues = 0:0.2:1;
            set(gca,'XTick',xtickValues);xlim([0 180]);
            %set(gca,'YTick',ytickValues)
        end
    end
end
%% Probability vs omega and Probability vs Pav figures (equal and logy axis) - used in mtvipul/extra.tex
if (0)
    eidlist = {'20090407201255351'};
    if 0
    eidlist = {'20070911234634153'};
    eidlist = {'20071003140612444'};
    eidlist = {'20071010180326301'};
    eidlist = {'20140418184418152'}; % NW
    eidlist={'20070911234634153' '20070919112226549' '20071003140612444' '20071010180326301' '20071128235703849'...
        '20080314093821771' '20080327230745201' '20080828231418631' '20080918194353069' '20081228071310738'...
        '20090124180950811' '20090215193500098' '20090223000427175' '20090317011333066' '20090407201255351'...
        '20090414171427415' '20090430045457938' '20090524094004552' '20090622192805162' '20090626164820729' '20090730223910267'};
    end
    for ii = 1:length(eidlist)
        eid = eidlist{ii};
        ddir = strcat('/home/vipul/CAP/inv/scak/MOOS/RESULTS/',eid,'/L1/M111/');
        %ddir = strcat('/home/vipul/CAP/inv/scak/NW/',eid,'/L1/M111/');
        f = dir(strcat(ddir,'*rand.mat'));
        fname =strcat(ddir,f.name);
        [kappa,theta,sigma,F,Mw,delta,gamma,M0,M] = read_CAP_log(fname);
        f = dir(strcat(ddir,'*rand'));
        fname =strcat(ddir,f.name);
        
        [fmin,indx1] = min(F);
        Mmin = M(:,indx1);
        [fmax,indx2] = max(F);
        Mmax = M(:,indx2);
        indx3 = 3831; % for opposite solution (-Mref) (obtianed by [a,b] = sort(om_ref))
        %indx3 = 208; 
        %indx3= 321;   % 
        Mref = Mmin;
        %Mref=M(:,indx3); %
        
        Nref = 600;
        ivec = 1:Nref;
        Mvec = M(:,ivec);
        Fref = F(ivec);
        om_ref = CMT2omega(Mvec,Mref);
        fal = F(ivec);
        prob = exp(-40*fal); % probabbility of all samples 
        
        if 1 % uncomment for getting (min, max and -M)
            isamp_vec = [indx1 indx2 indx3];
        end
        F3 = F(isamp_vec);
        P3 = exp(-40*F3); 
        om3 = CMT2omega(M(:,isamp_vec),Mref);
        b = load(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmegaVR3.mat'));
        Pav3 = b.Pav3;
        
        if 1
            figure(1); set(gca,'FontSize',22)
            plot(om_ref,prob,'.','Markersize',15)
            hold on
            xlabel('\omega');
            xtickValues = 0:30:180;
            set(gca,'XTick',xtickValues);xlim([0 180]);
            ylabel('p(M)');
            plot(om3(1),P3(1),'o','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'Markersize',11);
            plot(om3(2),P3(2),'o','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'Markersize',11);
            plot(om3(3),P3(3),'o','MarkerEdgeColor','k','MarkerFaceColor',[.5 0.5 .5],'Markersize',11);
            
            figure(2); set(gca,'FontSize',22)
            semilogy(om_ref,prob,'.','Markersize',15)
            hold on
            xlabel('\omega');
            xtickValues = 0:30:180;
            set(gca,'XTick',xtickValues);xlim([0 180]);
            ylabel('p(M)');
            semilogy(om3(1),P3(1),'o','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'Markersize',11);
            semilogy(om3(2),P3(2),'o','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'Markersize',11);
            semilogy(om3(3),P3(3),'o','MarkerEdgeColor','k','MarkerFaceColor',[.5 0.5 .5],'Markersize',11);
            
            a = load(strcat('/home/vipul/manuscripts/2014mt/data/',eid,'_PavOmegaVR.mat'));
            pall = a.Pav_all;
            figure(3); set(gca,'FontSize',22)
            plot(pall,prob,'.','Markersize',15)
            hold on
            xlabel('P_{av}');
            xtickValues = 0:0.2:1;
            set(gca,'XTick',xtickValues);xlim([0 1]);
            ylabel('p(M)');
            plot(Pav3(1),P3(1),'o','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'Markersize',11);
            plot(Pav3(2),P3(2),'o','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'Markersize',11);
            plot(Pav3(3),P3(3),'o','MarkerEdgeColor','k','MarkerFaceColor',[.5 0.5 .5],'Markersize',11);
            
            
            figure(4); set(gca,'FontSize',22)
            semilogy(pall,prob,'.','Markersize',15)
            hold on
            xlabel('P_{av}');
            xtickValues = 0:0.2:1;
            set(gca,'XTick',xtickValues);xlim([0 1]);
            ylabel('p(M)');
            semilogy(Pav3(1),P3(1),'o','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'Markersize',11);
            semilogy(Pav3(2),P3(2),'o','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'Markersize',11);
            semilogy(Pav3(3),P3(3),'o','MarkerEdgeColor','k','MarkerFaceColor',[.5 0.5 .5],'Markersize',11);
        end
    end
end


%% =====================================================================
% generate uncertainty summary files for -Mref (and the maximum misfit solution)

if (0)
    fname ='/home/vipul/CAP/inv/scak/MOOS/RESULTS/20090407201255351/L1/M111/log_039_rand.mat';
    [kappa,theta,sigma,F,Mw,delta,gamma,M0,M] = read_CAP_log(fname);
    
    eid = '20090407201255351';
    ddir =strcat('/home/vipul/CAP/inv/scak/MOOS/RESULTS/',eid,'/L1/M111/');
    GMTdir = '/home/vipul/gmt/data/cap/test/';
    
    [fmin,indx] = min(F);
    Mmin = M(:,indx);
    [fmax,indx] = max(F);
    Mmax = M(:,indx);
    
    %Mref = M(:,3831);   % for opposite solution (-Mref) (obtianed by [a,b] = sort(om_ref))
    %Mref = Mmax;  % for maximum misfit solution
    Mref = Mmin;
    
    [~,~, Pav]= CAP_unc(eid,Mref,ddir,GMTdir);
    MTbrick_Csection(eid,ddir,GMTdir,Mref);
end

%% =====================================================================
% Compute p(M0) for all the 21+2 events

if (0)
    clear all
    
    K  = 40; % misfit factor
    N = 100000;  %
    mtag = 'MOOS';
    eidlist={'20070911234634153' '20070919112226549' '20071003140612444' '20071010180326301' '20071128235703849'...
        '20080314093821771' '20080327230745201' '20080828231418631' '20080918194353069' '20081228071310738'...
        '20090124180950811' '20090215193500098' '20090223000427175' '20090317011333066' '20090407201255351'...
        '20090414171427415' '20090430045457938' '20090524094004552' '20090622192805162' '20090626164820729' '20090730223910267'};
    % eidlist={'20090317011333066'};
    % eidlist = {'20140418184418152'}; mtag = 'NW';     % NW
    % eidlist = {'20140416202423770'}; mtag = 'MISC';   % Denali
    
    for ii = 1:length(eidlist)
        eid = eidlist{ii};
        ddir = strcat('/home/vipul/CAP/inv/scak/MOOS/RESULTS/',eid,'/L1/M111/');
        %ddir = strcat('/home/vipul/CAP/inv/scak/NW/',eid,'/L1/M111/');
        %ddir = strcat('/home/vipul/CAP/inv/scak/MISC/',eid,'/L1/M111/');
        % Read the *.out file
        [otime,elat,elon,edep,strike,dip,rake,M,Mw,eid,capmod,capdep,rms,vr(ii),Pwin,Swin,Nstn,Pstn,Sstn,...
            stnm,lampPV,lampPR,lampSV,lampSR,lampSH,corrPV,corrPR,corrSV,corrSR,corrSH] = read_capout_old(strcat(ddir,eid,'.out'));
        % Read the log file
        f = dir(strcat(ddir,'*rand.mat'));
        fname =strcat(ddir,f.name);
        [kappa,theta,sigma,F,Mw,delta,gamma,M0,M] = read_CAP_log(fname);
        
        Fsum = sum(exp(-F*K));
        [fmin(ii),indx1] = min(F);
        Mmin = M(:,indx1);
        
        pM0(ii) = exp(-fmin(ii)*K); % unnormalized probability at the minimum
        dV = (4*pi*pi)/N;       % Volume of the brick
        kfac = 1/(dV*Fsum);     % normalization factor (Walt2015b notes)
        
        PM0(ii) = kfac*pM0(ii); % normalized probability at the minimum
        
        SlopeM0(ii) = 4*pi*pi*PM0(ii);          % Slope of the confidence curve at V=0
        Ang(ii) = atan(SlopeM0(ii))*180/pi;     % Angle of the confidence curve at V=0
    end
    
    for ii = 1:length(eidlist)
        fprintf('%d & %s & %f & %2.2f & %4.2f & %3.2f & %2.1f \\\\ \n',ii, eidlist{ii}, fmin(ii), PM0(ii), SlopeM0(ii), Ang(ii), vr(ii));
    end
end
%% =====================================================================

    % For running uncertainty test with increasing stations (both body
    % and/0r surface waves were used for all 10 stations)
    datadir = '/home/vipul/CAP/inv/scak/MOOS/20090407201255351/rand_log/';
    outputfilename = strcat(datadir,'log_rand_1');
    [Mref, post_samples_info,Pav,PV,PV_diff]= CAP_unc(outputfilename);
    outputfilename = strcat(datadir,'log_rand_2');
    [Mref, post_samples_info,Pav,PV,PV_diff]= CAP_unc(outputfilename);
    outputfilename = strcat(datadir,'log_rand_3');
    [Mref, post_samples_info,Pav,PV,PV_diff]= CAP_unc(outputfilename);
    outputfilename = strcat(datadir,'log_rand_4');
    [Mref, post_samples_info,Pav,PV,PV_diff]= CAP_unc(outputfilename);
    outputfilename = strcat(datadir,'log_rand_5');
    [Mref, post_samples_info,Pav,PV,PV_diff]= CAP_unc(outputfilename);
    outputfilename = strcat(datadir,'log_rand_6');
    [Mref, post_samples_info,Pav,PV,PV_diff]= CAP_unc(outputfilename);
    outputfilename = strcat(datadir,'log_rand_7');
    [Mref, post_samples_info,Pav,PV,PV_diff]= CAP_unc(outputfilename);
    outputfilename = strcat(datadir,'log_rand_8');
    [Mref, post_samples_info,Pav,PV,PV_diff]= CAP_unc(outputfilename);
    outputfilename = strcat(datadir,'log_rand_9');
    [Mref, post_samples_info,Pav,PV,PV_diff]= CAP_unc(outputfilename);
    outputfilename = strcat(datadir,'log_rand_10');
    [Mref, post_samples_info,Pav,PV,PV_diff]= CAP_unc(outputfilename);
%% =====================================================================    
    datadir = '/home/vipul/CAP/inv/scak/MOOS/20090407201255351/rand_log/';
    
    Nst=10;
    for ii=1:Nst
        outputfilename = strcat(datadir,'log_rand_',num2str(ii));
        [Mref, post_samples_info,Pav,PV,PV_diff(:,ii)]= CAP_unc(outputfilename);
    end
    
    Vvec = 0:0.01:1;
    Vfrac = 0.1;
    for ii=1:Nst
        iindx = find(Vvec==Vfrac);
        P_sum(ii) = sum(PV_diff(1:iindx,ii));
    end
    
    Vvec = linspace(0,1,100);
    figure;
    subplot(2,2,1)
    plot(Vvec,PV_diff)
    xlabel('V'); ylabel('P''(V)');
    subplot(2,2,2)
    plot(sum(PV_diff),'o-')
    xlabel('Nst'); ylabel('sum(P''(V))');
    subplot(2,2,3)
    plot(PV_diff(1,:),'o-')
    

%% =====================================================================
% Prepare data for all 10 multiple station suns for GMT plotting
    eid = '20090407201255351';
    Nst=10;
    for ii=1:Nst
        ddir = strcat('/home/vipul/CAP/inv/scak/MOOS/20090407201255351/stn',num2str(ii),'/');
        gmtdir = strcat('/home/vipul/gmt/data/cap/test_stn',num2str(ii),'/');
        CAP_unc(eid,ddir,gmtdir);
        MTbrick_Csection(eid,ddir,gmtdir);
    end 
 %% ===============
 % For test effect of adding SH component with only one station
    eid = '20090407201255351';
    Nst=1;
    for ii=1:Nst
        ddir = strcat('/home/vipul/CAP/inv/scak/MOOS/20090407201255351/stn',num2str(ii),'_SH/');
        gmtdir = strcat('/home/vipul/gmt/data/cap/test_stn',num2str(ii),'_SH/');
        CAP_unc(eid,ddir,gmtdir);
        MTbrick_Csection(eid,ddir,gmtdir);
    end 
    
 %% ============
 % For checking uncertainty when piercing points are at P-T axis
 eid = '20090407201255351';
    Nst=2;
    for ii=2:Nst
        ddir = strcat('/home/vipul/CAP/inv/scak/MOOS/20090407201255351/stn',num2str(ii),'_PT/');
        gmtdir = strcat('/home/vipul/gmt/data/cap/test_stn',num2str(ii),'_PT/');
        CAP_unc(eid,ddir,gmtdir);
        MTbrick_Csection(eid,ddir,gmtdir);
    end 
    
%% ============
 %  For checking uncertainty when inverting using No Shift
 eid = '20090407201255351';
 ddir = '/home/vipul/CAP/inv/scak/MOOS/20090407201255351/log_NoShift/stn1C/';
 gmtdir = '/home/vipul/gmt/data/cap/test_NoShift/stn1C/';
 CAP_unc(eid,ddir,gmtdir);
 MTbrick_Csection(eid,ddir,gmtdir);
 