% Be EXTREMELY careful when running this script !!!!!!!!!!!!!!!!!!!!!!
eid_dir = '/home/vipul/CAP/inv/scak/MOOSE/';
%eidlist = {'20090407201255','20090215193500','20090430045457','20090730223910',...
    '20080828231418','20090626164820','20090622192805','20080327230745','20090223000427',...
    '20081228071310','20090317011333','20090124180950','20090414171427','20090524094004'...
    '20070911234634','20080314093821','20071003140612','20071128235703','20070919112226'...
    '20071010180326'};

eidlist = {'20070911234634'};
ifname= 'M11.dat';
ofname= 'M_1.dat';

pfact = 0;
sfact = 1;

for jj=1:length(eidlist)
    eid = eidlist{jj};
    
    filename=strcat(eid_dir,eid,'/',eid,'/',ifname);
    fid = fopen(filename);
    wt = textscan(fid,'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n');
    fclose(fid);
    
    st=wt{1};
    st2=[wt{2}';wt{3}'*pfact;wt{4}'*pfact;wt{5}'*sfact;wt{6}'*sfact;wt{7}'*sfact;wt{8}';wt{9}';wt{10}';wt{11}';wt{12}'];
    
    outputfile = strcat(eid_dir,eid,'/',eid,'/',ofname);
    fid = fopen(outputfile,'w');
    for  ii=1:length(st)
        %fprintf(fid,'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',[wt{1}';wt{2}';wt{3}';wt{4}';wt{5}';wt{6}';wt{7}';wt{8}';wt{9}';wt{10}';wt{11}';wt{12}']);
        %fprintf(fid,'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',st{ii},st2(1,ii),st2(2,ii),st2(3,ii),st2(4,ii),st2(5,ii),st2(6,ii),st2(7,ii),st2(8,ii),st2(9,ii),st2(10,ii),st2(11,ii));
        fprintf(fid,'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',st{ii},st2(:,ii));
    end
    fclose(fid);
end