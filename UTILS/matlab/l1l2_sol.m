clear all
close all
iwrite=1;
ddir = '/home/vipul/CAP/inv/scak/MOOS/RESULTS/';
dirs = dir(strcat(ddir,'2*'));

gmtdir = '/home/vipul/gmt/data/cmt/l1l2/';

for ii=1:length(dirs)
    eid{ii} = dirs(ii).name;
    ffname = strcat(ddir,eid{ii},'/summary.dat');
    fid = fopen(ffname);
    w = textscan(fid,'%f %f %f %f %f %f %f %f');
    fclose(fid);
    norm = w{1};
    wt = w{2};
    stk = w{3};
    dip = w{4};
    rak = w{5};
    mag = w{6};
    dep = w{7};
    vr = w{8};
    Mii = dcfaultpar2CMT(stk,dip,rak,0);
    otime(1:4,1) = eid2otime(eid{ii});
    lat_vec=ones(4,1);
    lon_vec=ones(4,1);
    dep_vec=ones(4,1);
    for jj=1:length(norm)
        if norm(jj) == 1 %L1
            if wt(jj)==111
                M1(:,1) = Mii(:,jj);
            elseif wt(jj)==110
                M1(:,2) = Mii(:,jj);
            end
        end
               
        if norm(jj) == 2 %L2
            if wt(jj)==111
                M1(:,3) = Mii(:,jj);
            elseif wt(jj)==110
                M1(:,4) = Mii(:,jj);
            end
        end
    end
    
    [omegadc1(:,ii)] = CMT2omegadc_xi0(repmat(M1(:,1),1,4),M1,0,0);
    
    psmecaname = strcat(gmtdir,eid{ii});
    write_psmeca(psmecaname,otime,lat_vec,lon_vec,dep_vec,M1);
end

if iwrite==1
    fid = fopen(strcat(gmtdir,'Omega_L1'),'w');
    fprintf(fid,'%3.0f\t%3.0f\t%3.0f\t%3.0f\n',omegadc1)
    fclose(fid);
end


% OLD CODE
% filename = '/home/vipul/CAP/inv/scak/MOOS_old/l1l2.dat';
% fid = fopen(filename);
% 
% igmt=1;
% % direcoty for input files for gmt plotting
% gmtdir = '/home/vipul/gmt/data/';
% 
% wt=textscan(fid,'%s %f %f %f %f');       % read the log file
% fclose(fid);
% 
% eids = unique(wt{1});
% stks = wt{2};
% dips = wt{3};
% raks = wt{4};
% Mws = wt{5};
% 
% l1b = [];
% l1f = [];
% l2b = [];
% l2f = [];
% eid = [];
% for ii=1:length(eids)
%     l1b = [l1b;stks((ii-1)*4+1) dips((ii-1)*4+1) raks((ii-1)*4+1) Mws((ii-1)*4+1)];
%     l1f = [l1f;stks((ii-1)*4+2) dips((ii-1)*4+2) raks((ii-1)*4+2) Mws((ii-1)*4+2)];
%     l2b = [l2b;stks((ii-1)*4+3) dips((ii-1)*4+3) raks((ii-1)*4+3) Mws((ii-1)*4+3)];
%     l2f = [l2f;stks((ii-1)*4+4) dips((ii-1)*4+4) raks((ii-1)*4+4) Mws((ii-1)*4+4)];
%     eid = [eid;eids{ii}];
% end
% 
% 
% refdc = l1b;
% o1 = dc2omega(l1b,refdc);
% o2 = dc2omega(l1f,refdc);
% o3 = dc2omega(l2f,refdc);
% o4 = dc2omega(l2b,refdc);
% 
% l1b = [l1b o1];
% l1f = [l1f o2];
% l2f = [l2f o3];
% l2b = [l2b o4];
% 
% if igmt==1
%     fid = fopen(strcat(gmtdir,'l1bfile.dat'),'w');
%     fprintf(fid,'%3.0f\t%3.0f\t%3.0f\t%2.2f\t%3.0f\n',l1b');
%     fclose(fid);
%     fid = fopen(strcat(gmtdir,'l1ffile.dat'),'w');
%     fprintf(fid,'%3.0f\t%3.0f\t%3.0f\t%2.2f\t%3.0f\n',l1f');
%     fclose(fid);
%     fid = fopen(strcat(gmtdir,'l2bfile.dat'),'w');
%     fprintf(fid,'%3.0f\t%3.0f\t%3.0f\t%2.2f\t%3.0f\n',l2b');
%     fclose(fid);
%     fid = fopen(strcat(gmtdir,'l2ffile.dat'),'w');
%     fprintf(fid,'%3.0f\t%3.0f\t%3.0f\t%2.2f\t%3.0f\n',l2f');
%     fclose(fid);
% end