close all
clear all
clc
    
dircaprun = '/home/ksmith/REPOSITORIES/capuaf/nenanabasin_MTs/input_files/';
    
a = ls([dircaprun '*caprun']);
files = textscan(a,'%s','delimiter','\n');
files = files{1};
n = length(files);
    
% get proper tags
pre_capruntags = {'K','H','P','p','S','T','D','C','Y','Z','M','m','I','L','R','A','20'};
for cr_i = 1:length(pre_capruntags)-1
    capruntags{cr_i} = ['-' pre_capruntags{cr_i}];
end
capruntags{length(pre_capruntags)} = pre_capruntags{length(pre_capruntags)};
for ii = 1:n
    fid = fopen(files{ii});
    runstr = fgetl(fid)
    C = strsplit(runstr);
    fclose(fid);
    preC = {};
    for c_i = 1:length(C)
	if ~isempty(C{c_i})
    	    preC{c_i} = C{c_i}(1:2);
	end
    end
    [K,~,IC] = intersect(capruntags,preC,'stable');
    K2 = cell(length(capruntags),1);
    K2(IC) = preC(IC);
    
    %[K2,OIC,IC2] = union(capruntags,preC,'stable');
    %K2 = K2(1:length(capruntags));

    if length(preC) == length(capruntags)
	error('stopping')
    end
    orderedC = {C{1},C{IC}};
    fid2 = fopen([dircaprun 'reorder_MT_meta.txt'],'a');
    fprintf(fid2,[strjoin(orderedC) '\n'])
    fclose(fid2);
end

%fopen(files)
%fgetl(fid);
%[~,K,H,P,p,S,T,D,C,Y,Z,M,m,I,L,R,A,eidx] = textread(filename,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s')
%fclose(fid)
