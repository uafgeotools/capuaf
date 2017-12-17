eid = '20150912032512711';
Mref = TT2CMT(0,0,1,210,52,-20); % Original solution from nennuc_supp2

ddir = strcat('/home/vipul/CAP/inv/tact/nennuc/',eid,'/mw_dur_save/');
ftag = 'OUTPUT_DIR*';
t = dir(strcat(ddir,ftag));

outfilename = strcat(eid,'_tactmod_020.out');

for ii=1:length(t)
    dname = t(ii).name;
    aa = strsplit(dname,'_');
    bb = strsplit(aa{3},'L');
    dur(ii) = str2num(bb{2});
    bb = strsplit(aa{4},'bp');
    f1(ii) = str2num(bb{1});
    f2(ii) = str2num(bb{2});
    fc(ii) = (f1(ii)+f2(ii))/2;
    filename = strcat(ddir,dname,'/',outfilename);
    [otime,elat,elon,edep,strike(ii),dip(ii),rake(ii),M(:,ii),mw(ii),eid,capmod,capdep,rms,vr(ii),Pwin,Swin,Nstn,Pstn,Sstn,...
    stnm,lampPV,lampPR,lampSV,lampSR,lampSH,corrPV,corrPR,corrSV,corrSR,corrSH,pPol,ePol] = read_capout(filename);
end

omega = CMT2omega(Mref,M);

% CAP adds half extra as rise time (trapezoidal soure-time function)
% see srcfile (inside OUTPUT_DIR)
dur = dur + dur/2;

% Plotting (This is not working at the moment)
%[DUR,FC] = meshgrid(dur,fc);
%MW = griddata(dur,fc,mw,DUR,FC);
%VR = griddata(dur,fc,vr,DUR,FC);
%figure; surf(DUR,FC,MW*1000); %caxis([3.8,4.3])
%set(gca, 'clim', [3.6 -4.4]);
%colormap([0 0 0; jet]);
%colorbar;
%figure; pcolor(DUR,FC,VR); shading interp; colorbar

% Plotting (This works!)
figure; scatter(dur,fc,5000,vr,'s','filled')
xlim([4 48]); colorbar
set(gca,'FontSize',16)
xlabel('source duration');
ylabel('center of the banpass filter (filter width = 60secs)');
title('Variance reduction');
figure; scatter(dur,fc,5000,mw,'s','filled')
xlim([4 48]); colorbar
set(gca,'FontSize',16)
xlabel('source duration');
ylabel('center of the banpass filter (filter width = 60secs)');
title('magnitude (Mw)');
figure; scatter(dur,fc,5000,omega,'s','filled')
xlim([4 48]); colorbar
set(gca,'FontSize',16)
xlabel('source duration');
ylabel('center of the banpass filter (filter width = 60secs)');
title('Omega (angle from the Mref)');
