function [err1,err2,err3]= cap_error(filename)
%
% Plots cross section of misfit error with strike, dip and rake plane
% filename = 5 coloumn log files output by cap (uncomment the line :
%       " fprintf(stdout,"%f\t%f\t%f\t%e\t%f\n", sol.meca.stk, sol.meca.dip,
%       sol.meca.rak, sol.err, m0);" ----(in error function) and save
%       output in a log file
% Log file format:
%       strike    Dip     Rake    misfir_error    M0

 close all

fid = fopen(filename);
% Or remove input(filename) from function declaration and try with these fid's
% fid(1) = fopen('/home/vipul/CAP/inv/tact/log272787_118');
% fid(2) = fopen('/home/vipul/CAP/inv/tact/log272787_198');
% fid(1) = fopen('/home/vipul/CAP/inv/tact/log266672');
% fid(4) = fopen('/home/vipul/CAP/inv/scak/log319605');


wt=textscan(fid,'%f%f%f%f%f');       % read the log file
fclose(fid);

temp=wt{1};
wt{1}=temp(1:(length(temp)-1));

[ival,ix]=min(wt{4});                   % find minimum error location (ix) and its value (ival)

stk=wt{1};
dip=wt{2};
rak=wt{3};
err=wt{4};
m0=wt{5};
k1=1;
k2=1;
k3=1;

% stk(ix);
% dip(ix);
% rak(ix);
% err(ix);
% m0(ix);

st = [0:40:360];
dp = [0:10:90];
rk = [0:20:180]-90;
dipl = linspace(0,90,10);   %   x
stkl = linspace(0,360,37);  %   y
rakl = linspace(-90,90,19); %   z

if rak(ix)~=rk(:)
    rk(:)=rk(:)+10;
end
rk(10)=90;

st(:)=st(:)-min((st(:)-stk(ix)));
st(:)=wrapTo360(st(:));

for jj=1:10
    for ii=1:length(m0)
        if m0(ii)==m0(ix)
            if st(jj)==stk(ii)
                st(jj);
                % dip1(k1)=dip(ii);
                % rak1(k1)=rak(ii);
                err1(k1)=err(ii);
                k1=k1+1;
            end
            if dp(jj)==dip(ii)
                % stk1(k2)=stk(ii);
                % rak2(k2)=rak(ii);
                err2(k2)=err(ii);
                k2=k2+1;
            end
            if rk(jj)==rak(ii)
                % dip2(k3)=dip(ii);
                % stk2(k3)=stk(ii);
                err3(k3)=err(ii);
                k3=k3+1;
            end
        end
    end
    
    [rakg,dipg]=meshgrid(rakl,dipl);
    for ii=1:190
        if  mod(ii,10) ~= 0
            err4(mod(ii,10),floor(ii/10)+1)=err1(ii);    % arrange error according to mesh
        else
            err4(mod(ii,10)+10,floor(ii/10))=err1(ii);   % changes for last row of mesh (multiple of 37)
        end
    end
    subplot(10,3,3*(jj-1)+1)
    surface(rakg,dipg,err4,'LineStyle','none');
    xlim([-90 90])
    ylim([0 90])
    text(-90,70,num2str(stk(jj)),'erasemode','none');
    
    [stkg,rakg]=meshgrid(rakl,stkl);
    for ii=1:703
        if  mod(ii,37) ~= 0
            err5(mod(ii,37),floor(ii/37)+1)=err2(ii);    % arrange error according to mesh
        else
            err5(mod(ii,37)+37,floor(ii/37))=err2(ii);   % changes for last row of mesh (multiple of 37)
        end
    end
    subplot(10,3,3*(jj-1)+2)
    surface(rakg,stkg,err5,'LineStyle','none');
    xlim([0 360])
    ylim([-90 90])
    text(0,70,num2str(dip(jj)),'erasemode','none');
    
    [dipg,stkg]=meshgrid(dipl,stkl);
    for ii=1:370
        if  mod(ii,37) ~= 0
            err6(mod(ii,37),floor(ii/37)+1)=err3(ii);    % arrange error according to mesh
        else
            err6(mod(ii,37)+37,floor(ii/37))=err3(ii);   % changes for last row of mesh (multiple of 37)
        end
    end
    subplot(10,3,3*(jj-1)+3)
    surface(stkg,dipg,err6,'LineStyle','none');
    xlim([0 360])
    ylim([0 90])
    text(0,70,num2str(rak(jj)),'erasemode','none');
    
    k1=1;
    k2=1;
    k3=1;
end

subplot(10,3,28)
xlabel('rake');
ylabel('dip');
subplot(10,3,30)
xlabel('strike');
ylabel('dip');
subplot(10,3,29)
xlabel('strike');
ylabel('rake');
% suptitle([' Minimum at: STK ',num2str(stk(ix)), ' DIP ', int2str(dip(ix)),...
%    ' RAK ',int2str(rak(ix)),' -- Misfit ',num2str(err(ix))]);
colorbar('location','southoutside');
%cap_contour(filename);
%cap_slice(filename);


%==== Example======
if 1==0
    % fid(1) = fopen('/home/vipul/CAP/inv/tact/log272787_118');
    % fid(2) = fopen('/home/vipul/CAP/inv/tact/log272787_198');
    % fid(1) = fopen('/home/vipul/CAP/inv/tact/log266672');
    % fid(4) = fopen('/home/vipul/CAP/inv/scak/log319605');
    % filename= '/home/vipul/CAP/inv/tact/log272787_118';
    filename= '/home/vipul/CAP/inv/cus/log_20080418093700';
    cap_error(filename);
end
