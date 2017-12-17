function [err1,err2,dip2,stk2,rak2] = cap_slice(filename)
%
% Plots cross-section of strike-dip-rake alongwith misfit error sliced at 
% mimimum error point
% Input :  5 coloumn log files output by cap (uncomment this line in cap.c:
%       " fprintf(stdout,"%f\t%f\t%f\t%e\t%f\n", sol.meca.stk, sol.meca.dip,
%       sol.meca.rak, sol.err, m0);" ----(in error function) and save
%       output in a log file
% Log file format:
%       strike    Dip     Rake    misfir_error    M0


% Or remove input(filename) from function declaration and try with these fid's
% fid(1) = fopen('/home/vipul/CAP/inv/tact/log272787_118');
% fid(2) = fopen('/home/vipul/CAP/inv/tact/log272787_198');
% fid(1) = fopen('/home/vipul/CAP/inv/tact/log266672');
% fid(4) = fopen('/home/vipul/CAP/inv/scak/log319605');

fid = fopen(filename);

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
%err=log(ival)./log(err);

stk_min = stk(ix);
dip_min = dip(ix);
rak_min = rak(ix);

for ii=1:length(m0)
    if m0(ii)==m0(ix)
            err1(j)=err(ii);        % extract 'all' (strike dip and error) at minimum rake
            t1(j)=stk(ii);
            t2(j)=dip(ii);
            t3(j)=rak(ii);
            j=j+1;
    end
end

dip1 = linspace(0,90,10);
stk1 = linspace(0,360,37);
rak1 = linspace(-90,90,19);
[stk2,dip2,rak2]=meshgrid(stk1,dip1,rak1);

err2(1:10,1:37,1:19)=0;

for i=1:19
    for j=1:10
        for k=1:37
            err2(j,k,i)=err1(370*(i-1)+37*(j-1)+k);
        end
    end
end


slice(stk2,dip2,rak2,err2,stk_min,dip_min,rak_min)
xlabel('strike');
xlim([0 360]);
ylabel('dip');
ylim([0 90]);
zlabel('rake');
zlim([-90 90]);
title(['Minimum at: STK ',num2str(stk(ix)), ' DIP ', int2str(dip(ix)),...
    ' RAK ',int2str(rak(ix)),' -- Misfit ',num2str(err(ix))]);
colorbar('location','eastoutside');      


%==== Example======
if 1==0
    % fid(1) = fopen('/home/vipul/CAP/inv/tact/log272787_118');
    % fid(2) = fopen('/home/vipul/CAP/inv/tact/log272787_198');
    % fid(1) = fopen('/home/vipul/CAP/inv/tact/log266672');
    % fid(4) = fopen('/home/vipul/CAP/inv/scak/log319605');
    filename= '~/CAP/inv/tact/log272787_118';
    cap_slice(filename);
end