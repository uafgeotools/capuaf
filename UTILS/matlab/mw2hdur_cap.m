function [durCAPe,durCAP,durGCMT] = mw2hdur_cap(Mw)
% Outputs the source duration for CAP and GCMT
% Input: Mw
% Output:   
%   durCAPe     source duration used by CAP
%   durCAP      source duration calculated by CAP 
%   durGCMT     source duration used by GCMT   
% See cap.pl and m02hdur
% See example below.
%
% Note: If you input source duration using -L(t) flag in CAP then the actual
% source duration is going to be: (1.5 x t)
%
% Might be a good idea to just use the GCMT equation. It works even for
% smaller magnitudes

iplot = 1;
rise = 0.5;

durCAP = floor((10.^((Mw-5)/2)+0.5));  % From cap.pl (uses (int) function of perl)
%durCAP = 10.^((Mw-5)/2+0.5); 

% This reduces to:
% durCAP = 10^((Mw-4)/2);
% durCAP < 1 for Mw < 4

% This section is from cap.pl
%  $dura = 1 if $dura < 1;
%  $dura = 9 if $dura > 9;

for ii=1:length(Mw)
    if durCAP(ii)<1
        durCAPe(ii)=1;
        
    elseif durCAP(ii)>9
        durCAPe(ii)=9;
    else
        durCAPe(ii)=durCAP(ii);
    end
end

riseTime = rise*durCAPe

% Total source function duration after convolution of (riseTime and durCAP)
% This is the full duration
durCAPe = riseTime+durCAPe

%------------------------------------------
% GCMT dura
imag = 2; % 2 for GCMT; 1 for Kanamori 1977 (CAP uses Kanamori)

M0 = mw2m0(imag,Mw);
hdur = m02hdur(M0); % this returns half duration
durGCMT = 2*hdur;

% Plotting
if iplot
    mvec = 0:0.01:10;
    % For CAP
    durvec1 = floor((10.^((mvec-5)/2)+0.5));
    %
    durvec1e = durvec1;
    ivec = find(durvec1e<1);
    durvec1e(ivec)=1;
    ivec = find(durvec1e>9);
    durvec1e(ivec)=9;
    durvec1e = rise*(durvec1e) + durvec1e; % add riseTime length (after convolution)
    % For GCMT
    M0vec = mw2m0(2,mvec);
    durvec2 = 2*m02hdur(M0vec); % full duration
    
    figure
    set(gca,'FontSize',16);
    semilogy(mvec,durvec1e,'--r','Linewidth',2)
    hold on
    semilogy(mvec,durvec2,'--g','Linewidth',2)
    grid on
    xlabel('Mw');
    ylabel('source duraiton')
    plot(Mw,durCAPe,'ok','Markersize',10,'MarkerFaceColor','r')
    plot(Mw,durGCMT,'ok','Markersize',10,'MarkerFaceColor','g')
    title(sprintf('Mw = %2.1f; dura: CAP = %.2fs, GCMT = %.2fs',Mw,durCAPe,durGCMT))
    legend('CAP','GCMT')
end

%==========================================================================
if 0
    for Mw=2:6
        [durCAP,durGCMT,durCAPe] = mw2hdur_cap(Mw);
        disp(sprintf('Mw %5.2f : duration %7.3f%7.3f%7.3f',Mw,durCAP,durGCMT,durCAPe));
    end
end
%==========================================================================

