clear all
close all
clc

Nfiles = 1;
multi_run=0;   % run for more than one log files
misfit_factor = 15;  % works only if greater than 1
true_per=1;         % works if misfit_factor==1; if true_per is 1 mean of ptryN will be used
log_len = 0:Nfiles-1;
iomega = 0;     % plot error vs omega (and omega histo) for both full grid and random samples
icov = 0;       % will calculate the covariance of posterior (not needed!)
err_norm=1;
prob_norm=0;
igmt =0;
iminus = 0;
ifig = 0;
irun =1;
cap_aeic= 0;
%misfit_fact = [5 10 15 20 25 30];

% direcoty for input files for gmt plotting
gmtdir = '/home/vipul/gmt/data/';

% data directory (with log files and output files from CAP)
datadir = '/home/vipul/CAP/inv/scak/MOOS/';
% NEW id's

eidlist = {'20090407201255','20090215193500','20090430045457','20090730223910',...
    '20080828231418','20090626164820','20090622192805','20080327230745','20090223000427',...
    '20081228071310','20090317011333','20090124180950','20090414171427','20090524094004'...
    '20070911234634','20080314093821','20071003140612','20071128235703','20070919112226'...
    '20080918194353'};

dep = {'042','056','042','051','048','056','068','078','081','080','095','095','104','116',...
    '095','142','044','070','036','068'};

eidlist = {'20090407201255','20090215193500','20090430045457','20090730223910',...
    '20090626164820','20090622192805','20080327230745','20090223000427',...
    '20081228071310','20090317011333','20090124180950','20090414171427','20090524094004'...
    '20080918194353'};

dep = {'042','056','042','051','056','068','081','080','095','095','104','116',...
    '068'};

%---test events-------------
% eidlist = {'20090407201255'}; dep={'042'};
% eidlist = {'20090215193500'}; dep={'056'};
% eidlist = {'20090430045457'}; dep={'042'};
% eidlist = {'20090730223910'}; dep={'051'};
% eidlist = {'20080828231418'}; dep={'051'};
% eidlist = {'20090626164820'}; dep={'057'};
% eidlist = {'20090622192805'}; dep={'069'};
% eidlist = {'20080327230745'}; dep={'078'};
% eidlist = {'20090223000427'}; dep={'085'};
% eidlist = {'20081228071310'}; dep={'088'};
% eidlist = {'20090317011333'}; dep={'095'};
% eidlist = {'20090124180950'}; dep={'100'};
% eidlist = {'20090414171427'}; dep={'104'};
% eidlist = {'20090524094004'}; dep={'111'};
  eidlist = {'20070911234634'}; dep={'095'};
% eidlist = {'20080314093821'}; dep={'150'};
% eidlist = {'20071003140612'}; dep={'041'};
% eidlist = {'20071128235703'}; dep={'074'};
% eidlist = {'20070919112226'}; dep={'042'};
% eidlist = {'20080918194353'}; dep={'068'};
% eidlist = {'20071010180326'}; dep={'027'};
% eidlist = {'20090102141536'}; dep={'052'}; % GCMT
% %--------------------------------

 datadir = '/home/vipul/GEOTOOLS/tomo_util/cap/';
 eidlist = {'20080418093700'};
 
datadir = '/home/vipul/CAP/inv/scak/MISC/';
 eidlist = {'20140416202424814'};
 
 datadir = '/home/vipul/CAP/inv/scak/MISC/';
 eidlist = {'20140418184418317'};
%%
for il=1:length(eidlist)
        
        close all
        %clc
        eid = eidlist{il}
        edep = dep{il};

        evtdir = strcat(datadir,eid,'/');
        
        % filename = strcat('/home/vipul/CAP/inv/scak/319605/log/log_grid_fine');
        
        for ii=1:1:length(log_len)
            if multi_run==1
                if ii<=10
                    % filename = strcat('/home/vipul/GEOTOOLS/tomo_util/cap/data/data2/log_00',num2str(log_len(ii)));
                    filename = strcat(evtdir,'log/log_',edep,'_00',num2str(log_len(ii)));
                    % 319605
                    %filename = strcat('/home/vipul/CAP/inv/scak/319605/log/log_00',num2str(log_len(ii)));
                else
                    filename = strcat(evtdir,'log/log_',edep,'_0',num2str(log_len(ii)));
                end
            else
                % filename = '/home/vipul/GEOTOOLS/tomo_util/cap/data/log_grid_fine';
                % filename = '/home/vipul/GEOTOOLS/tomo_util/cap/data/log_grid';
                %filename = '/home/vipul/GEOTOOLS/tomo_util/cap/data/log_grid_coarse';
                %        filename = '/home/vipul/GEOTOOLS/tomo_util/cap/log_000';
                % SCAK
                filename = strcat(evtdir,'log/log_',edep,'_grid')
                if cap_aeic
                    switch jj
                        case 1
                            filename = strcat(evtdir,'cap_aeic_log/M11_grid');
                        case 2
                            filename = strcat(evtdir,'cap_aeic_log/M12_grid');
                        case 3
                            filename = strcat(evtdir,'cap_aeic_log/M21_grid');
                        case 4
                            filename = strcat(evtdir,'cap_aeic_log/M22_grid');
                    end
                end
                %filename = strcat(evtdir,'../log/log_015_grid')
                % filename = strcat(evtdir,'/log/log_077_grid')
                filename = strcat(evtdir,'/log/log_023_grid')
                
            end
            fid = fopen(filename);
            
            wt=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f');       % read the log file
            fclose(fid);
            
            stk=wt{1};
            dip=wt{2};
            rak=wt{3};
            err=wt{4};
            mw=wt{5};
            vared = err;
            %---------------------------
            figure
            subplot(2,2,1)
            plot(stk,err,'.')
            subplot(2,2,2)
            plot(dip,err,'.')
            subplot(2,2,3)
            plot(rak,err,'.')
            subplot(2,2,4)
            plot_histo(err,linspace(min(err),max(err),20));
            %---------------------------
            
            [ival(ii),ix]=min(err);
            stk_min = stk(ix);
            dip_min = dip(ix);
            rak_min = rak(ix);
            M_min0 = dcfaultpar2CMT(stk_min,dip_min,rak_min,0);
            if iminus ==1
                [kappa1,theta1,sigma1] =faultpar2minusM(stk_min,dip_min,rak_min); % 180 (opposite fault plane)
                sigma1= real(sigma1);
            end
            % Distane from (M)
            M_min = repmat(M_min0,1,length(stk));
            [Mi,k1,d1,p1,p1,p2,p3] = dcfaultpar2CMT(stk,dip,rak,0);
            [omega1,xi1] = CMT2omega_xi0(Mi,M_min,0,0);
            
            % Distance from (-M)
            if iminus ==1
                M_min1 = faultpar2CMT(kappa1,theta1,sigma1,0);
                M_min = repmat(M_min1,1,length(stk));
                [omega2,xi2] = CMT2omega_xi0(Mi,M_min,0,0);
            end
            
            %-----------------------------------------
            % PLOT in moment tensor space and other vectors (see
            % faultpar2CMT.m)
            %         figure
            %         subplot(3,2,1)
            %         plot(Mi(1,:),err,'.')
            %         subplot(3,2,2)
            %         plot(Mi(2,:),err,'.')
            %         subplot(3,2,3)
            %         plot(Mi(3,:),err,'.')
            %         subplot(3,2,4)
            %         plot(Mi(4,:),err,'.')
            %         subplot(3,2,5)
            %         plot(Mi(5,:),err,'.')
            %         subplot(3,2,6)
            %         plot(Mi(6,:),err,'.')
            %         figure
            %         subplot(2,2,1)
            %         scatter3(p1(1,:),p1(2,:),p1(3,:),50,err,'.')
            %         subplot(2,2,2)
            %         scatter3(p2(1,:),p2(2,:),p2(3,:),50,err,'.')
            %         subplot(2,2,3)
            %         scatter3(p3(1,:),p3(2,:),p3(3,:),50,err,'.')
            
            %-----------------------------------------
            % scaling error
            if err_norm ==1
                %err = err-ival(ii);
                %%%err = ((err-ival(ii))/min(err))*misfit_factor;
                
                err = log((err)/min(err))/log(max(err)/min(err));
                err = err*misfit_factor;
                %err = err*misfit_factor;
                %err = err - min(err);
            end
            %
            %err = ((err)*facto);
            %----------------------------------------
            stk_len  = length(unique(stk));
            dip_len  = length(unique(dip));
            rak_len  = length(unique(rak));
            
            % create vectors
            %dip1 = linspace(0,1,dip_len);
            %dip2 = acosd(dip1);
            stk1 = linspace(min(stk),max(stk),stk_len);
            rak1 = linspace(min(rak),max(rak),rak_len);
            dip1 = unique(dip);
            
            % find minimum
            istk = find(stk1==stk_min);
            idip = find(dip1==dip_min);
            irak = find(rak1==rak_min);
            
            ddip=1/(dip_len-1);
            
            % increment
            dstk=360/(stk_len-1);
            drak=180/(rak_len-1);
            
            % create mesh
            [STK,DIP,RAK]=meshgrid(stk1,dip1,rak1);
            ERR(1:dip_len,1:stk_len,1:rak_len)=0;
            omega0 = zeros(dip_len,stk_len,rak_len);
            omega0_opp = zeros(dip_len,stk_len,rak_len);
            xi0 = zeros(dip_len,stk_len,rak_len);
            xi0_opp = zeros(dip_len,stk_len,rak_len);
            
            % organize ERR into a Mesh format
            K=0;
            K1=0;
            for i=1:rak_len
                for j=1:dip_len
                    for k=1:stk_len
                        ERR(j,k,i)=err(stk_len*dip_len*(i-1)+stk_len*(j-1)+k);
                        sig(j,k,i)=exp((-ERR(j,k,i)));
                        K= K+ exp(-ERR(j,k,i))*(drak*ddip*dstk*((pi/180)^2));
                        K1 = K1 + (drak*ddip*dstk*((pi/180)^2));
                        %[STK(j,k,i) DIP(j,k,i) RAK(j,k,i)]
                        omega0(j,k,i) = omega1(stk_len*dip_len*(i-1)+stk_len*(j-1)+k);
                        xi0(j,k,i) = xi1(stk_len*dip_len*(i-1)+stk_len*(j-1)+k);
                        if iminus==1
                            omega0_opp(j,k,i) = omega2(stk_len*dip_len*(i-1)+stk_len*(j-1)+k);
                            xi0_opp(j,k,i) = xi2(stk_len*dip_len*(i-1)+stk_len*(j-1)+k);
                        end
                        % plot(omega_orthoB(j,k,i),ERR(j,k,i),'.')
                        % temp(stk_len*(i-1)+dip_len*(j-1)+k)=omega_orthoB(j,k,i);
                        % plot(xi_orthoB(j,k,i),ERR(j,k,i),'.')
                    end
                end
            end
            
            %--------------------plotting misift-----------
            
            if ifig==1
                figure;
                subplot(3,1,1)
                hold on
                plot(stk1,ERR(idip,:,irak),'-x')
                text(stk1(end),ERR(idip,stk_len,irak),num2str(mw(1)));
                xlabel('STK');
                ylabel('ERR');
                title(sprintf('Dip=%3.0f Rak=%d',dip1(idip),rak1(irak)));
                subplot(3,1,2)
                hold on
                plot(dip1,ERR(:,istk,irak),'-x')
                text(dip1(end),ERR(dip_len,istk,irak),num2str(mw(1)));
                xlabel('DIP');
                ylabel('ERR');
                title(sprintf('Stk=%d Rak=%d',stk1(istk),rak1(irak)));
                subplot(3,1,3)
                hold on
                plot(rak1,squeeze(ERR(idip,istk,:)),'-x')
                text(rak1(end),ERR(idip,istk,rak_len),num2str(mw(1)));
                xlabel('RAK');
                ylabel('ERR');
                
                title(sprintf('Stk=%d Dip=%3.0f',stk1(istk),dip1(idip)));
                
                %figure
                %plot_histo(ERR(:),linspace(min(ERR(:)),max(ERR(:)),20));
            end
            
            % --------------------- misfit and pdf mesh--------------------
            if 1==1
                
                % Sliced plots for misfit and pdf
                PROB = exp(-ERR)/K;
                figure
                subplot(2,1,1);
                slice(STK,DIP,RAK,ERR,stk_min,dip_min,rak_min)
                xlabel('STK');
                ylabel('DIP');
                zlabel('RAK');
                title('Misfit Error');
                shading flat
                axis equal
                colorbar
                subplot(2,1,2);
                slice(STK,DIP,RAK,PROB/20,stk_min,dip_min,rak_min)
                xlabel('STK');
                ylabel('DIP');
                zlabel('RAK');
                title('pdf');
                shading flat
                axis equal
                colorbar
                max_sigma = max(max(max(PROB)));
                min_sigma = min(min(min(PROB)));
                fprintf('stk_inc=%f \t cos(dip_inc)=%f \t rak_inc=%f\n',dstk,ddip,drak);
                fprintf('K=%f \t max_sigma=%f\n',K,max_sigma);
                axis equal
                
                % Slice misift plot through (-M) instead
                %             figure
                %             subplot(2,1,2);
                %             slice(STK,DIP,RAK,ERR,kappa1,theta1,sigma1)
                %             xlabel('STK');
                %             ylabel('DIP');
                %             zlabel('RAK');
                %             title('Misfit Error');
                %             shading flat
                %             colorbar
                %             axis equal
                
                % Plot for multiple slices
                %             figure
                %             for j=1:length(dip1)
                %                 subplot(5,4,j)
                %                 slice(STK,DIP,RAK,ERR,stk_min,dip1(j),rak_min)
                %                 view(0,0)
                %                 shading flat
                %             end
                %             figure
                %             for j=1:length(dip1)
                %                 subplot(5,4,j)
                %                 slice(STK,DIP,RAK,ERR,stk1(4*j),dip_min,rak_min)
                %                 view(0,90)
                %                 shading flat
                %             end
                %             figure
                %             for j=1:length(dip1)
                %                 subplot(5,4,j)
                %                 slice(STK,DIP,RAK,ERR,stk_min,dip_min,rak1(2*j))
                %                 view(2)
                %                 shading flat
                %             end
                %--------------------------------------------------------
                figure
                subplot(2,2,1)
                stk_sec=squeeze(ERR(:,find(stk_min == stk1),:));
                stk_omega=squeeze(omega0(:,find(stk_min == stk1),:));
                [x1,y1] = meshgrid(rak1,dip1);
                Z1 = [x1(:) y1(:) stk_sec(:)];
                C1 = [x1(:) y1(:) stk_omega(:)];
                imagesc(stk_sec)
                subplot(2,2,2)
                dip_sec=squeeze(ERR(find(dip_min == dip1),:,:));
                dip_omega=squeeze(omega0(find(dip_min == dip1),:,:));
                [x2,y2] = meshgrid(rak1,stk1);
                Z2 = [x2(:) y2(:) dip_sec(:)];
                C2 = [x2(:) y2(:) dip_omega(:)];
                imagesc(dip_sec)
                subplot(2,2,3)
                rak_sec=squeeze(ERR(:,:,find(rak_min == rak1)));
                rak_omega=squeeze(omega0(:,:,find(rak_min == rak1)));
                [x3,y3] = meshgrid(stk1,dip1);
                Z3 = [x3(:) y3(:) rak_sec(:)];
                C3 = [x3(:) y3(:) rak_omega(:)];
                imagesc(rak_sec)
                
                % figure, surf(x1,y1,stk_sec) % surf(rak1,dip1 , stk_sec)
                
                % sections of opposite solutions
                if iminus==1
                    figure
                    subplot(2,2,1)
                    [a,b]=min(abs(stk1-kappa1)); % plane closest to opposite MT solution
                    stk_sec=squeeze(ERR(:,b,:)); % Error for all dip and rake in that strike plane
                    stk_omega=squeeze(omega0_opp(:,b,:)); % Omega for all dip and rake in that strike plane
                    [x1,y1] = meshgrid(rak1,dip1); % create mesh
                    ZZ1 = [x1(:) y1(:) stk_sec(:)]; % for output as xyz file
                    CC1 = [x1(:) y1(:) stk_omega(:)];
                    imagesc(stk_sec)
                    subplot(2,2,2)
                    [a,b]=min(abs(dip1-theta1));
                    dip_sec=squeeze(ERR(b,:,:));
                    dip_omega=squeeze(omega0_opp(b,:,:));
                    [x2,y2] = meshgrid(rak1,stk1);
                    ZZ2 = [x2(:) y2(:) dip_sec(:)];
                    CC2 = [x2(:) y2(:) dip_omega(:)];
                    imagesc(dip_sec)
                    subplot(2,2,3)
                    [a,b]=min(abs(rak1-sigma1));
                    rak_sec=squeeze(ERR(:,:,b));
                    rak_omega=squeeze(omega0_opp(:,:,b));
                    [x3,y3] = meshgrid(stk1,dip1);
                    ZZ3 = [x3(:) y3(:) rak_sec(:)];
                    CC3 = [x3(:) y3(:) rak_omega(:)];
                    imagesc(rak_sec)
                end
                
                % Will plot kagan angle contour in GMT instead of omega
                %             figure
                %             subplot(2,2,1)
                %             [a,b]=min(abs(stk1-kappa1));
                %             stk_sec=squeeze(ERR(:,b,:));
                %             stk_omega=squeeze(xi0(:,b,:));
                %             [x1,y1] = meshgrid(rak1,dip1);
                %             ZZ1 = [x1(:) y1(:) stk_sec(:)];
                %             CC1 = [x1(:) y1(:) stk_omega(:)];
                %             imagesc(stk_sec)
                %             subplot(2,2,2)
                %             [a,b]=min(abs(dip1-theta1));
                %             dip_sec=squeeze(ERR(b,:,:));
                %             dip_omega=squeeze(xi0(b,:,:));
                %             [x2,y2] = meshgrid(rak1,stk1);
                %             ZZ2 = [x2(:) y2(:) dip_sec(:)];
                %             CC2 = [x2(:) y2(:) dip_omega(:)];
                %             imagesc(dip_sec)
                %             subplot(2,2,3)
                %             [a,b]=min(abs(rak1-sigma1));
                %             rak_sec=squeeze(ERR(:,:,b));
                %             rak_omega=squeeze(xi0(:,:,b));
                %             [x3,y3] = meshgrid(stk1,dip1);
                %             ZZ3 = [x3(:) y3(:) rak_sec(:)];
                %             CC3 = [x3(:) y3(:) rak_omega(:)];
                %             imagesc(rak_sec)
                
            end
            
            % --------------------plotting pdf-----------
            if ifig==1
                figure;
                subplot(3,1,1)
                hold on
                plot(stk1,PROB(idip,:,irak),'-x')
                text(stk1(end),PROB(idip,stk_len,irak),num2str(mw(1)));
                xlabel('STK');
                ylabel('pdf');
                title(sprintf('Dip=%3.0f Rak=%d',dip1(idip),rak1(irak)));
                subplot(3,1,2)
                hold on
                plot(dip1,PROB(:,istk,irak),'-x')
                text(dip1(end),PROB(dip_len,istk,irak),num2str(mw(1)));
                xlabel('DIP');
                ylabel('pdf');
                title(sprintf('Stk=%d Rak=%d',stk1(istk),rak1(irak)));
                subplot(3,1,3)
                hold on
                plot(rak1,squeeze(PROB(idip,istk,:)),'-x')
                text(rak1(end),PROB(idip,istk,rak_len),num2str(mw(1)));
                xlabel('RAK');
                ylabel('pdf');
                title(sprintf('Stk=%d Dip=%3.0f',stk1(istk),dip1(idip)));
            end
            
            % --------------------plotting normalized pdf-----------
            PROBN = PROB/max_sigma;
            if ifig==1
                figure;
                subplot(3,1,1)
                hold on
                plot(stk1,PROBN(idip,:,irak),'-x')
                text(stk1(end),PROBN(idip,stk_len,irak),num2str(mw(1)));
                xlabel('STK');
                ylabel('pdfN');
                title(sprintf('Dip=%3.0f Rak=%d',dip1(idip),rak1(irak)));
                ylim([0 1])
                subplot(3,1,2)
                hold on
                plot(dip1,PROBN(:,istk,irak),'-x')
                text(dip1(end),PROBN(dip_len,istk,irak),num2str(mw(1)));
                xlabel('DIP');
                ylabel('pdfN');
                title(sprintf('Stk=%d Rak=%d',stk1(istk),rak1(irak)));
                ylim([0 1])
                subplot(3,1,3)
                hold on
                plot(rak1,squeeze(PROBN(idip,istk,:)),'-x')
                text(rak1(end),PROBN(idip,istk,rak_len),num2str(mw(1)));
                xlabel('RAK');
                ylabel('pdfN');
                title(sprintf('Stk=%d Dip=%3.0f',stk1(istk),dip1(idip)));
                ylim([0 1])
            end
            % check that sum of probability in complete model space is 1
            pdf_check=sum(sum(sum(PROB)))*dstk*ddip*drak*(pi/180)^2
            
            if multi_run==0
                break
            end
        end
        %% =========================random sample generation=======================
        clear data4 data3 Nplot1 Nplot2 Nplot3 Nplot4
        filename = strcat(evtdir,'log/log_',edep,'_rand');
        if cap_aeic
        switch jj
            case 1
                filename = strcat(evtdir,'cap_aeic_log/M11_grid');
            case 2
                filename = strcat(evtdir,'cap_aeic_log/M12_grid');
            case 3
                filename = strcat(evtdir,'cap_aeic_log/M21_grid');
            case 4
                filename = strcat(evtdir,'cap_aeic_log/M22_grid');
        end
        end
        filename = strcat(evtdir,'norm_log/M11_l1_rand');
        filename = strcat(evtdir,'norm_log/full_l1_rand');
        filename = strcat(evtdir,'/log/log_077_rand');
        filename = strcat(evtdir,'/log/log_023_rand');
        
        fid = fopen(filename);
        
        wt=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f');       % read the log file
        fclose(fid);
        
        rnd_stk=wt{1};
        rnd_dip=wt{2};
        rnd_rak=wt{3};
        rnd_err=wt{4};
        rnd_mw=wt{5};
        N=length(rnd_stk);
        fprintf('Number of random samples = %d\n',N);
        
        [ival,ix] = min(rnd_err);
        stk_min2 = rnd_stk(ix);
        dip_min2 = rnd_dip(ix);
        rak_min2 = rnd_rak(ix);
        M_min0 = dcfaultpar2CMT(stk_min2,dip_min2,rak_min2,0);
        if iminus==1
            [kappa2,theta2,sigma2] =faultpar2minusM(stk_min2,dip_min2,rak_min2); % 180 (opposite fault plane)
        end
        if err_norm==1
            %rnd_err = rnd_err-ival;
            %%%rnd_err = (rnd_err-ival)/min(rnd_err)*misfit_factor;
            rnd_err = log(rnd_err/min(rnd_err))/log(max(rnd_err)/min(rnd_err));
            rnd_err = rnd_err*misfit_factor;
        end
        %rnd_err=rnd_err*facto;
        %==============
        
        %     if iminus == 1
        %         M_min0 = faultpar2CMT(kappa2,theta2,sigma2,0);
        %         M_min = repmat(M_min0,1,N);
        %         Mi1 = faultpar2CMT(rnd_stk, rnd_dip , rnd_rak,0);
        %         [omega1,xi1] = CMT2omega_xi0(Mi1,M_min,0,0);
        %         data2(:,1) = omega1(:);
        %         data2(:,2) = rnd_err(:);
        %         data2 = sortrows(data2,1);
        %         figure,
        %         plot(data2(:,1),data2(:,2),'.')
        %         figure
        %         subplot(2,1,1)
        %         plot_histo(xi1(:),[0:5:120],1)
        %         xlabel('xi (kagan)');
        %         subplot(2,1,2)
        %         [N1,Nplot1]=plot_histo(omega1(:)*pi/180,[0:pi/36:pi],3)     % for plotting
        %         [N3,Nplot3]=plot_histo(omega1(:)*pi/180,[0:pi/180:pi],3)    % for computing uncertainty
        %         [N2,Nplot2]=plot_histo(xi1(:)*pi/180,[0:pi/36:pi*2/3],3)
        %         [N4,Nplot4]=plot_histo(xi1(:)*pi/180,[0:pi/180:pi*2/3],3)
        %         xlabel('omega (tape)');
        %     end
        
        if 1 == 1
            M_min = repmat(M_min0,1,N);
            Mi1 = dcfaultpar2CMT(rnd_stk, rnd_dip , rnd_rak,0);
            [omega1,xi1] = CMT2omega_xi0(Mi1,M_min,0,0);
            data2(:,1) = omega1(:);
            data2(:,2) = rnd_err(:);
            data2 = sortrows(data2,1);
            figure,
            plot(data2(:,1),data2(:,2),'.')
            figure
            subplot(2,1,1)
            plot_histo(xi1(:),[0:5:120],1)
            xlabel('xi (kagan)');
            subplot(2,1,2)
            plot_histo(omega1(:),[0:5:180],1)
            [N1,Nplot1]=plot_histo(omega1(:)*pi/180,[0:pi/36:pi],3)     % for plotting
            [N3,Nplot3]=plot_histo(omega1(:)*pi/180,[0:pi/180:pi],3)    % for computing incertainty
            [N2,Nplot2]=plot_histo(xi1(:)*pi/180,[0:pi/36:pi*2/3],3)
            [N4,Nplot4]=plot_histo(xi1(:)*pi/180,[0:pi/180:pi*2/3],3)
            xlabel('omega (tape)');
        end
        
        [a,b]=(min(180-omega1));
        sol_minus =[rnd_stk(b),rnd_dip(b),rnd_rak(b),rnd_err(b),omega1(b),xi1(b)];
        
        % looking at the trend of xi and omega and misfit (all sorted, no
        % correltation between them)
        %     figure
        %     plot(sort(xi1)/max(xi1),'.')
        %     hold on
        %     plot(sort(omega1)/max(omega1),'.r')
        %     plot(sort(rnd_err)/max(rnd_err),'.g')
        %     legend('kagan','omega','misfit')
        %     title('All normalized and sorted:\n only pattern is useful')
        %     title('All normalized and sorted: only pattern is useful')
        %--------------------------------------------
        
        %     figure
        %     subplot(2,2,1)
        %     plot(rnd_stk,rnd_err,'.')
        %     subplot(2,2,2)
        %     plot(rnd_dip,rnd_err,'.')
        %     subplot(2,2,3)
        %     plot(rnd_rak,rnd_err,'.')
        %     subplot(2,2,4)
        %     plot_histo(rnd_err,linspace(min(err),max(err),20));
        
        
        % Different ways of normalizing
        ptry=(1/K)*exp(-rnd_err);
        ptryN = ptry/max_sigma;
        
        % Old efforts of normalizing (all discarded)
        % if misfit_factor~=1
        %     chance=rand(N,1);
        % else
        %     if true_per==1
        %         true_per=mean(ptryN);
        %         true_per=min(ptryN);
        %         %true_per=.95
        %     end
        %     chance=true_per+(1-true_per)*rand(N,1);
        % end
        
        chance=rand(N,1);
        ikeep = find(chance < ptryN);
        %ikeep=ikeep(1:10000);
        
        % POSTERIOR SAMPLES
        stk_sample=rnd_stk(ikeep);
        dip_sample=rnd_dip(ikeep);
        rak_sample=rnd_rak(ikeep);
        err_sample=rnd_err(ikeep);
        [stk_sample2,dip_sample2,rak_sample2]=dcfaultpar2nodal(stk_sample,dip_sample,rak_sample);
        [stk_min3,dip_min3,rak_min3]=dcfaultpar2nodal(stk_min,dip_min,rak_min);
        
        %scatter3(rnd_stk,rnd_dip,rnd_rak)
        if ifig==1
            figure
            scatter3(stk_sample,dip_sample,rak_sample,'.')
            xlim([0 360]);
            ylim([0 90]);
            zlim([-90 90]);
            xlabel('STK');
            ylabel('DIP');
            zlabel('RAK');
            hold on
            plot3(stk_min,dip_min,rak_min,'m.','Markersize',30)
        end
        % plot cross-section
        stk_ind = find(stk_sample > stk_min-5);
        stk_ind = find(stk_sample(stk_ind) < stk_min+5);
        
        dip_ind = find(dip_sample > dip_min-5);
        dip_ind = find(dip_sample(dip_ind) < dip_min+5);
        
        rak_ind = find(rak_sample > rak_min-5);
        rak_ind = find(rak_sample(rak_ind) < rak_min+5);
        
        figure
        subplot(2,2,1)
        scatter(dip_sample(stk_ind),rak_sample(stk_ind),'.')
        xlabel('DIP');
        ylabel('RAK');
        title(sprintf('stk = %3.2f',stk_min));
        subplot(2,2,2)
        scatter(stk_sample(dip_ind),rak_sample(dip_ind),'.')
        xlabel('STK');
        ylabel('RAK');
        title(sprintf('dip = %3.2f',dip_min));
        subplot(2,2,3)
        scatter(stk_sample(rak_ind),dip_sample(rak_ind),'.')
        xlabel('STK');
        ylabel('DIP');
        title(sprintf('rak = %3.2f',rak_min));
        
        % plot histogram for samples
        figure
        subplot(2,2,1)
        plot_histo(stk_sample,[0:20:360]);
        title('stk\_samples');
        subplot(2,2,2)
        plot_histo(dip_sample,[0:5:90]);
        title('dip\_samples');
        subplot(2,2,3)
        plot_histo(rak_sample,[-90:10:90]);
        title('rak\_samples');
        
        %break
        
        %% ----------------plotting beachballs of posterior samples------------------
        nr = 5;
        nc = 4;
        meca_len = (nr*nc);
        omega_len = length(ikeep);
        try_samp = [stk_sample(1:omega_len) dip_sample(1:omega_len) rak_sample(1:omega_len)];
        M_min_v = repmat(M_min0,1,length(try_samp));
        Mi2 = dcfaultpar2CMT(try_samp(:,1), try_samp(:,2) , try_samp(:,3),0);
        [omega2,xi2] = CMT2omega_xi0(Mi2,M_min_v,0,0);
        
        data3(:,1) = omega2;
        data3(:,2) = [err_sample(1:omega_len)];
        data3 = [data3 try_samp xi2];
        data3 = sortrows(data3,2);
        figure,
        subplot(2,2,1)
        plot(data3(:,1),data3(:,2),'.')
        xlim([0 180]);
        subplot(2,2,2)
        plot(data3(:,6),data3(:,2),'.')
        xlim([0 120]);
        subplot(2,2,3)
        [N7,Nplot7]=plot_histo(omega2(:)*pi/180,[0:pi/180:pi],3);
        [N5,Nplot5]=plot_histo(omega2(:)*pi/180,[0:pi/36:pi],3);
        xlabel('omega (tape)');
        subplot(2,2,4)
        [N8,Nplot8]=plot_histo(xi2(:)*pi/180,[0:pi/180:pi*2/3],3);
        [N6,Nplot6]=plot_histo(xi2(:)*pi/180,[0:pi/36:pi*2/3],3);
        xlabel('xi (kagan)');
        
        % distance in misfit and kagan angle (omega) space
        % normalized
        %     E1= (0.5*((data3(:,1)*pi/180)/max(data3(:,1)*pi/180)).^2 + (data3(:,2)/max(data3(:,2))).^2).^0.5;
        %     E2= (0.5*((data3(:,8)*pi/180)/max(data3(:,8)*pi/180)).^2 + (data3(:,2)/max(data3(:,2))).^2).^0.5;
        %     figure; subplot(1,2,1); plot(E1,'.');
        %     subplot(1,2,2); plot(E2,'.');
        
        % Fixing the number of elements
        Nplot1 = [Nplot1; Nplot1(end)];
        Nplot2 = [Nplot2; Nplot2(end)];
        Nplot3 = [Nplot3; Nplot3(end)];
        Nplot4 = [Nplot4; Nplot4(end)];
        Nplot5 = [Nplot5; Nplot5(end)];
        Nplot6 = [Nplot6; Nplot6(end)];
        Nplot7 = [Nplot7; Nplot7(end)];
        Nplot8 = [Nplot8; Nplot8(end)];
        
        % Finding uncertainity in OMEGA
        bin_plot = pi/36;
        bin_unc = pi/180;
        [X1,Y1] = stairs([0:bin_plot:pi],Nplot1);      % homogeneous for plotting
        [X5,Y5] = stairs([0:bin_plot:pi],Nplot5);      % posterior for plotting
        [X3,Y3] = stairs([0:bin_unc:pi],Nplot3);      % homogeneous for uncertainty
        [X7,Y7] = stairs([0:bin_unc:pi],Nplot7);      % posterior for uncertainty
        % ploting
        dis1 = Nplot5./Nplot1;                          % Likelihood for plotting
        dis1(find(isnan(dis1)==1))=0;
        const1=sum(dis1)*bin_plot;
        dis1 = dis1/const1;
        [XX1,YY1] = stairs([0:bin_plot:pi],dis1);
        % unc
        dis3 = Nplot7./Nplot3;                          % Likelihood for uncertainty
        dis3(find(isnan(dis3)==1))=0;
        const3=sum(dis3)*bin_unc;
        dis3 = dis3/const3;
        [XX3,YY3] = stairs([0:bin_unc:pi],dis3);
        % 95% confidence interval
        cd3=cumsum(dis3*bin_unc);
        iOMEGA=min(find(cd3>.95));
        OMEGA = (iOMEGA)*bin_unc*180/pi
        figure
        subplot(2,2,1)
        plot(X1,Y1)
        hold on; plot(X5,Y5,'g')
        plot(XX1,YY1,'r')
        plot([OMEGA*pi/180 OMEGA*pi/180],[0 4],'k-.')
        xlim([0 pi])
        grid on
        xlabel('omega');
        title(sprintf('uncertainity: %2.0f',OMEGA));
        legend('homogeneous','posterior','likelihood');
        
        % Finding uncertainity in XI
        [X2,Y2] = stairs([0:bin_plot:pi*2/3],Nplot2);     % homogeneous for plotting
        [X6,Y6] = stairs([0:bin_plot:pi*2/3],Nplot6);     % posterior for plotting
        [X4,Y4] = stairs([0:bin_unc:pi*2/3],Nplot4);     % homogeneous for uncertainty
        [X8,Y8] = stairs([0:bin_unc:pi*2/3],Nplot8);     % posterior for uncertainty
        % ploting
        dis2 = Nplot6./Nplot2;                          % Likelihood for plotting
        dis2(find(isnan(dis2)==1))=0;
        const2=sum(dis2)*bin_plot;
        dis2 = dis2/const2;
        [XX2,YY2] = stairs([0:bin_plot:pi*2/3],dis2);
        % unc
        dis4 = Nplot8./Nplot4;                          % Likelihood for uncertainty
        dis4(find(isnan(dis4)==1))=0;
        const4=sum(dis4)*bin_unc;
        dis4 = dis4/const4;
        [XX4,YY4] = stairs([0:bin_unc:pi*2/3],dis4);
        % 95% confidence interval
        cd4=cumsum(dis4*bin_unc);
        iXI=min(find(cd4>.95));
        XI = (iXI)*bin_unc*180/pi
        subplot(2,2,2)
        plot(X2,Y2)
        hold on; plot(X6,Y6,'g')
        plot(XX2,YY2,'r')
        plot([XI*pi/180 XI*pi/180],[0 6],'-.')
        xlim([0 2*pi/3])
        grid on
        xlabel('kagan');
        title(sprintf('uncertainity: %2.0f',XI));
        legend('homogeneous','posterior','likelihood');
                
        %% Plotting P-T axis on beachball
        Msamp=dcfaultpar2CMT(data3(:,3),data3(:,4),data3(:,5));
        [lam,U] = CMTdecom(Msamp);
        Usamp = convertv(1,5,U);  % convert GCMT to south-east-up
        Upasamp = U2pa(Usamp,1,0);
        Upasamp(:,3:4) = [];     % cut the nodal axis
        %Upasamp(end+1,:) = [0 0 0 90];
        Tp = Upasamp(:,1);
        Ta = Upasamp(:,2);
        Pp = Upasamp(:,3);
        Pa = Upasamp(:,4);
        [Txs,Tys] = pa2xy(Tp,Ta);
        [Pxs,Pys] = pa2xy(Pp,Pa);
        
        subplot(2,2,3);
        th=linspace(0,2*pi); plot(cos(th),sin(th),'k'); hold on;
        %plot(Txs,Tys,'.');
        %hold on, plot(Pxs,Pys,'.r');
        scatter(Pxs(1:end),Pys(1:end),50,data3(:,2),'.')  % flip to plot the low misfit ones on top
        scatter(Txs(1:end),Tys(1:end),50,data3(:,2),'.')
        plot(Txs(1),Tys(1),'ko','markersize',7,'linewidth',2,'markerfacecolor','m')
        plot(Pxs(1),Pys(1),'ko','markersize',7,'linewidth',2,'markerfacecolor','m'); axis equal
        
        % samples within 95% CI
        %     samp95=max(find(data3(:,1)<OMEGA));
        %     subplot(2,2,4);
        %     th=linspace(0,2*pi); plot(cos(th),sin(th),'k'); hold on;
        %     plot(Txs(1:samp95),Tys(1:samp95),'.');
        %     hold on, plot(Pxs(1:samp95),Pys(1:samp95),'.r');
        %     plot(Txm(1),Tym(1),'ko','markersize',7,'linewidth',2,'markerfacecolor','c')
        %     plot(Pxm(1),Pym(1),'ko','markersize',7,'linewidth',2,'markerfacecolor','c'); axis equal
        
        %     post_samp=100;
        %     isamples=randperm(length(data3),post_samp);
        %     subplot(2,2,4);
        %     th=linspace(0,2*pi); plot(cos(th),sin(th),'k'); hold on;
        %     scatter(Pxs(isamples),Pys(isamples),50,data3(isamples,2),'.')
        %     scatter(Txs(isamples),Tys(isamples),50,data3(isamples,2),'.')
        %     plot(Txs(1),Tys(1),'ko','markersize',7,'linewidth',2,'markerfacecolor','m')
        %     plot(Pxs(1),Pys(1),'ko','markersize',7,'linewidth',2,'markerfacecolor','m'); axis equal
        
        % PLotting beachball for best solution
        subplot(2,2,4)
        plot_samples([data3(1,3),data3(1,4),data3(1,5)],1,1,1)
        
        pt_rnd = randperm(length(ikeep),50);
        Txs = [Txs(1);  Txs(pt_rnd)];
        Tys = [Tys(1);  Tys(pt_rnd)];
        Pxs = [Pxs(1);  Pxs(pt_rnd)];
        Pys = [Pys(1);  Pys(pt_rnd)];
        
        %% plotting beachballs of posterior samples------------------
        data4 = sortrows(data3(randperm(length(ikeep),20),:),2);
        
        figure
        if 1==1
            plot_samples([data4(:,3),data4(:,4),data4(:,5)],nr,nc,1)
        end
        %=====================================================
        if icov==1
            cov = zeros(3,3);
            m(1,:) = stk_sample - stk_min;
            m(2,:) = dip_sample - dip_min;
            m(3,:) = rak_sample - rak_min;
            for ii=1:length(ikeep)
                cov = cov+m(:,ii)*m(:,ii)';
            end
            cov = cov/length(ikeep);
            figure, imagesc(cov)
            colorbar
        end
        
        %=========================================
        % stk_inc=10.000000 	 cos(dip_inc)=0.111111 	 rak_inc=10.000000
        % K=21.722655 	 max_sigma=0.045354
        % stk_inc=5.000000 	 cos(dip_inc)=0.055556 	 rak_inc=5.000000
        % K=19.848311 	 max_sigma=0.049750
        % stk_inc=3.000000 	 cos(dip_inc)=0.033333 	 rak_inc=3.000000
        % K=19.125331 	 max_sigma=0.051628
        
        %break
        %% ==========OUTPUT for GMT================================
        if igmt==1
            
            plotdir=strcat(gmtdir,eid,'/');
            mkdir(plotdir);
            
            copyfile(strcat(evtdir,eid,'/',strcat('scak','_',edep,'.out')),strcat(plotdir,eid,'.out'));
            copyfile(strcat(evtdir,eid,'/',strcat('scak','_',edep,'.ps')),strcat(plotdir,eid,'.ps'));
            copyfile(strcat(evtdir,eid,'/',strcat('scak','_',edep,'_beach.ps')),strcat(plotdir,eid,'_beach.ps'));
            copyfile(strcat(evtdir,eid,'/',strcat('tactmod','_',edep,'.out')),strcat(plotdir,eid,'.out'));
            copyfile(strcat(evtdir,eid,'/',strcat('tactmod','_',edep,'.ps')),strcat(plotdir,eid,'.ps'));
            copyfile(strcat(evtdir,eid,'/',strcat('tactmod','_',edep,'_beach.ps')),strcat(plotdir,eid,'_beach.ps'));
            copyfile(strcat(evtdir,eid,'/',strcat(eid,'_event.dat')),strcat(plotdir));
            copyfile(strcat(evtdir,eid,'/',strcat(eid,'_station.dat')),strcat(plotdir));
            % file for preparing histogram of stk,dip,rak
            data_gmt = [stk_sample';dip_sample';rak_sample'];
            fid = fopen(strcat(plotdir,eid,'_histo.dat'),'w');
            fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',data_gmt);
            fclose(fid);
            data_gmt2 = [stk_sample2';dip_sample2';rak_sample2'];
            data_gmt2 = [nearest([stk_min3;dip_min3;rak_min3]),[stk_sample2';dip_sample2';rak_sample2']];
            fid = fopen(strcat(plotdir,eid,'_histo2.dat'),'w');
            fprintf(fid,'%3.0f\t%3.0f\t%3.0f\n',data_gmt2);
            fclose(fid);
            % file of posterior samples
            fid = fopen(strcat(plotdir,eid,'_samples.dat'),'w');
            fprintf(fid,'%3.2f\t%3.3f\t%3.3f\t%3.2f\t%3.2f\t%3.2f\n',data4');
            fclose(fid);
            % negative moment tensor solution
            fid = fopen(strcat(plotdir,eid,'_opp_solution.dat'),'w');
            fprintf(fid,'%3.0f\t%3.0f\t%3.0f\t%3.6f\t%3.2f\t%3.2f\n',sol_minus');
            fclose(fid);
            % XYZ files
            fid = fopen(strcat(plotdir,eid,'_stkfile.dat'),'w');
            fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',Z1');
            fclose(fid);
            fid = fopen(strcat(plotdir,eid,'_dipfile.dat'),'w');
            fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',Z2');
            fclose(fid);
            fid = fopen(strcat(plotdir,eid,'_rakfile.dat'),'w');
            fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',Z3');
            fclose(fid);
            if iminus==1
                fid = fopen(strcat(plotdir,eid,'_stkfile_opp.dat'),'w');
                fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',ZZ1');
                fclose(fid);
                fid = fopen(strcat(plotdir,eid,'_dipfile_opp.dat'),'w');
                fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',ZZ2');
                fclose(fid);
                fid = fopen(strcat(plotdir,eid,'_rakfile_opp.dat'),'w');
                fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',ZZ3');
                fclose(fid);
            end
            % XYZ for omega contour
            fid = fopen(strcat(plotdir,eid,'_stkcont.dat'),'w');
            fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',C1');
            fclose(fid);
            fid = fopen(strcat(plotdir,eid,'_dipcont.dat'),'w');
            fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',C2');
            fclose(fid);
            fid = fopen(strcat(plotdir,eid,'_rakcont.dat'),'w');
            fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',C3');
            fclose(fid);
            if iminus==1
                fid = fopen(strcat(plotdir,eid,'_stkcont_opp.dat'),'w');
                fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',CC1');
                fclose(fid);
                fid = fopen(strcat(plotdir,eid,'_dipcont_opp.dat'),'w');
                fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',CC2');
                fclose(fid);
                fid = fopen(strcat(plotdir,eid,'_rakcont_opp.dat'),'w');
                fprintf(fid,'%3.2f\t%3.2f\t%3.2f\n',CC3');
                fclose(fid);
            end
            % (omega and misfit) for all samples
            xi_omega_sort = sortrows([omega1 rnd_err rnd_err],1);
            fid = fopen(strcat(plotdir,eid,'_omega_err.dat'),'w');
            fprintf(fid,'%3.2f\t%3.3f\t%3.3f\n',xi_omega_sort');
            fclose(fid);
            xi_err_sort = sortrows([xi1 rnd_err rnd_err],1);
            fid = fopen(strcat(plotdir,eid,'_xi_err.dat'),'w');
            fprintf(fid,'%3.2f\t%3.3f\t%3.3f\n',xi_err_sort');
            fclose(fid);
            % (omega and misfit) for posterior samples
            fid = fopen(strcat(plotdir,eid,'_omega_err_post.dat'),'w');
            fprintf(fid,'%3.2f\t%3.3f\n',[omega2'; rnd_err(ikeep)']);
            fclose(fid);
            fid = fopen(strcat(plotdir,eid,'_xi_err_post.dat'),'w');
            fprintf(fid,'%3.2f\t%3.3f\n',[xi2'; rnd_err(ikeep)']);
            fclose(fid);
            % pdf (probability density functions) for homogeneous and posterior
            % samples in omega space
            fid = fopen(strcat(plotdir,eid,'_omega_err_dis.dat'),'w');
            fprintf(fid,'%3.3f\t%3.3f\n',[XX1'; YY1']);
            fclose(fid);
            fid = fopen(strcat(plotdir,eid,'_omega_homo_pdf.dat'),'w');
            fprintf(fid,'%3.3f\t%3.3f\n',[X1'; Y1']);
            fclose(fid);
            fid = fopen(strcat(plotdir,eid,'_omega_post_pdf.dat'),'w');
            fprintf(fid,'%3.3f\t%3.3f\n',[X5'; Y5']);
            fclose(fid);
            fid = fopen(strcat(plotdir,eid,'_xi_err_dis.dat'),'w');
            fprintf(fid,'%3.3f\t%3.3f\n',[XX2'; YY2']);
            %fprintf(fid,'%3.3f\t%3.3f\n',[XX4'; YY4']);
            fclose(fid);
            fid = fopen(strcat(plotdir,eid,'_xi_homo_pdf.dat'),'w');
            fprintf(fid,'%3.3f\t%3.3f\n',[X2'; Y2']);
            %fprintf(fid,'%3.3f\t%3.3f\n',[X4'; Y4']);
            fclose(fid);
            fid = fopen(strcat(plotdir,eid,'_xi_post_pdf.dat'),'w');
            fprintf(fid,'%3.3f\t%3.3f\n',[X6'; Y6']);
            %fprintf(fid,'%3.3f\t%3.3f\n',[X8'; Y8']);
            fclose(fid);
            fid = fopen(strcat(plotdir,eid,'_OMEGA.dat'),'w');
            fprintf(fid,'%2.0f',OMEGA);
            %fprintf(fid,'0 0 14 0 1 LB \\@~W\\@~=%2.1f\\@+o\\@+',OMEGA);
            fclose(fid);
            fid = fopen(strcat(plotdir,eid,'_XI.dat'),'w');
            fprintf(fid,'%2.0f',XI);
            fclose(fid);
            % P-T axis
            fid = fopen(strcat(plotdir,eid,'_pxs.dat'),'w');
            fprintf(fid,'%2.4f\t%2.4f\t%2.4f\n',[Pxs'; Pys'; 0 omega2(pt_rnd)']);
            fclose(fid);
            fid = fopen(strcat(plotdir,eid,'_txs.dat'),'w');
            fprintf(fid,'%2.4f\t%2.4f\t%2.4f\n',[Txs'; Tys'; 0 omega2(pt_rnd)']);
            fclose(fid);
        end
end


