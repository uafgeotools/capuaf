close all
clear all
clc
    
dircaprun = '/home/ksmith/REPOSITORIES/capuaf/nenanabasin_MTs/input_files/';
delete([dircaprun 'ordered_MT_meta.txt'])
a = ls([dircaprun 'unordered_caprun/20*caprun']);
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
    [K,ID,IC] = intersect(capruntags,preC,'stable');
    J2 = cell(length(capruntags),1);
    K2 = cell(length(capruntags),1);
    J2(ID) = capruntags(ID);
    %K2(ID) = preC(IC);
    L_i = find(strcmp(capruntags,'-L'));
    Y_i = find(strcmp(capruntags,'-Y'));
    A_i = find(strcmp(capruntags,'-A'));
    M_i = find(strcmp(capruntags,'-M'));
    EID_i = find(strcmp(capruntags,'20'));
    K2{L_i} = '-L1.0'; 
    K2{Y_i} = '-Y1'; 
    K2{A_i} = '-A0/0/0'; 
    K2(ID) = C(IC);

    if length(preC) == length(capruntags)
	%error('stopping')
    end
    %orderedC = {C{1},C{IC}};
    orderedC = {C{1},K2{:}}

    if 1 == 0
    % use get_meta in path to make this file
    fid2 = fopen([dircaprun 'ordered_MT_meta.txt'],'a');
    fprintf(fid2,[strjoin(orderedC) '\n'])
    fclose(fid2);
    end

    fid3 = fopen([dircaprun K2{EID_i} '_' K2{M_i}(3:end) '_caprun'],'w');
    fprintf(fid3,[strjoin(orderedC) '\n'])
    fclose(fid3);
end

%fopen(files)
%fgetl(fid);
%[~,K,H,P,p,S,T,D,C,Y,Z,M,m,I,L,R,A,eidx] = textread(filename,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s')
%fclose(fid)
