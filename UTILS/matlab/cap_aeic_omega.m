filename = '/home/vipul/CAP/inv/scak/MOOSE/aeic_cap.dat';
fid = fopen(filename);

wt=textscan(fid,'%s %f %f %f %f');       % read the log file
fclose(fid);

eids = unique(wt{1});
stks = wt{2};
dips = wt{3};
raks = wt{4};
Mws = wt{5};

M11a = [];
M11 = [];
M12 = [];
M21 = [];
M22 = [];
eid = [];
for ii=1:length(eids)
    M11a = [M11a;stks((ii-1)*5+1) dips((ii-1)*5+1) raks((ii-1)*5+1) Mws((ii-1)*5+1)];
    M11 = [M11;stks((ii-1)*5+2) dips((ii-1)*5+2) raks((ii-1)*5+2) Mws((ii-1)*5+2)];
    M12 = [M12;stks((ii-1)*5+3) dips((ii-1)*5+3) raks((ii-1)*5+3) Mws((ii-1)*5+3)];
    M21 = [M21;stks((ii-1)*5+4) dips((ii-1)*5+4) raks((ii-1)*5+4) Mws((ii-1)*5+4)];
    M22 = [M22;stks((ii-1)*5+5) dips((ii-1)*5+5) raks((ii-1)*5+5) Mws((ii-1)*5+5)];
    eid = [eid;eids{ii}];
end

[o1,x1] = dc2omega(M11a,M11);
[o2,x2] = dc2omega(M12,M11);
[o3,x3] = dc2omega(M21,M11);
[o4,x4] = dc2omega(M22,M11);

