function [t1,t2,err1]=cap_contour(filename)
%
% Plots contour map of error with strike dip and rake at minimum rake.
% Input :  5 coloumn log files output by cap (uncomment the line :
%       " fprintf(stdout,"%f\t%f\t%f\t%e\t%f\n", sol.meca.stk, sol.meca.dip,
%       sol.meca.rak, sol.err, m0);" ----(in error function) and save
%       output in a log file
% Log file format:
%       strike    Dip     Rake    misfir_error    M0


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
j=1;

for ii=1:length(m0)
    if m0(ii)==m0(ix)
        if rak(ii)==rak(ix)
            err1(j)=err(ii);        % extract 'all' (strike dip and error) at minimum rake
            t1(j)=stk(ii);
            t2(j)=dip(ii);
            j=j+1;
        end
    end
end


dip1 = linspace(0,90,10);
stk1 = linspace(0,360,37);

[dip2,stk2]=meshgrid(dip1,stk1);    % create meshgrid


for ii=1:370
    if  mod(ii,37) ~= 0
        err2(mod(ii,37),floor(ii/37)+1)=err1(ii);    % arrange error according to mesh
    else
        err2(mod(ii,37)+37,floor(ii/37))=err1(ii);   % changes for last row of mesh (multiple of 37)
    end
end

figure;
subplot(2,1,2);
surfc(dip2,stk2,err2,'LineStyle','none');
% meshc(dip2,stk2,err2);
% surface(dip2,stk2,err2);
ylabel('strike');
ylim([0 360]);
xlabel('dip');
xlim([0 90]);
zlabel('error');
title(['Minimum at: STK ',num2str(stk(ix)), ' DIP ', int2str(dip(ix)),...
    ' RAK ',int2str(rak(ix)),' -- Misfit ',num2str(err(ix))]);
colorbar('location','eastoutside');
hold on
plot3(stk(ix),dip(ix),err(ix),'k.','MarkerSize', 30)    % Mark minimum error point

%=====================================

xi=linspace(0,90,90);
yi=linspace(0,360,360);
[xii,yii] = meshgrid(xi,yi);
zii = griddata(dip2,stk2,err2,xii,yii);
subplot(2,1,1);
surfc(xii,yii,zii,'LineStyle','none');
ylabel('strike');
ylim([0 360]);
xlabel('dip');
xlim([0 90]);
zlabel('error');

%==========Example=====================
if 1==0
    % fid(1) = fopen('/home/vipul/CAP/inv/tact/log272787_118');
    % fid(2) = fopen('/home/vipul/CAP/inv/tact/log272787_198');
    % fid(1) = fopen('/home/vipul/CAP/inv/tact/log266672');
    % fid(4) = fopen('/home/vipul/CAP/inv/scak/log319605');
    filename= '~/CAP/inv/tact/log272787_118';
    cap_contour(filename);
end