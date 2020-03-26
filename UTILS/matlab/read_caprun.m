function [Kf,Hf,Pf,pf,Sf,Tf,Df,Cf,Yf,Zf,Mfmod,capdep,mf,If,Lf,Rf,Af,eidout] = read_caprun(filenames)
%READ_CAPRUN read the input file command line for CAP
%
% TEMPORARY SCRIPT: This works for only a limited set of cases.
%
% WORK NEEDED:
% 1) flexibility for different types of input lines, many of which will
%    have a different number of input flags
% 2) flexibility for the arbitrary order of input flags that is allowed by CAP
%
% tweaks: remove a -W1 flag, add -Y1 flag (could be removed from all), add a -L0 flag, add a -A0/0/0 flag
%
% Carl Tape, 2020-03-25
%

n = length(filenames);
eidout = repmat(cellstr(''),n,1);

for ii=1:n
    ii
    filename = filenames{ii}
    if ~exist(filename,'file')
        warning('file does not exist');
        continue
    end
    [~,K,H,P,p,S,T,D,C,Y,Z,M,m,I,L,R,A,eidx] = textread(filename,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
    K = char(K); Kf(ii) = str2num(K(3));
    H = char(H); Hf(ii) = str2num(H(3:end));
    Px = strsplit(P{1},'/'); P1 = Px{1}; Pf(ii,:) = [str2num(P1(3)) str2num(Px{2}) str2num(Px{3})];
    p = char(p); pf(ii) = str2num(p(3:end));
    Sx = strsplit(S{1},'/'); S1 = Sx{1}; Sf(ii,:) = [str2num(S1(3)) str2num(Sx{2}) str2num(Sx{3})];
    Tx = strsplit(T{1},'/'); T1 = Tx{1}; Tf(ii,:) = [str2num(T1(3)) str2num(Tx{2})];
    Dx = strsplit(D{1},'/'); D1 = Dx{1}; Df(ii,:) = [str2num(D1(3)) str2num(Dx{2}) str2num(Dx{3})];
    Cx = strsplit(C{1},'/'); C1 = Cx{1}; Cf(ii,:) = [str2num(C1(3)) str2num(Cx{2}) str2num(Cx{3}) str2num(Cx{4})];
    Y = char(Y); Yf(ii) = str2num(Y(3));
    Z = char(Z); Zf(ii) = cellstr(Z(3:end));
    M = strsplit(M{1},'_'); Mmod = M{1}; Mfmod(ii) = cellstr(Mmod(3:end)); capdep(ii) = str2num(M{2});
    mx = strsplit(m{1},'/'); m1 = mx{1}; mf(ii,:) = [str2num(m1(3)) str2num(mx{2}) str2num(mx{3})];
    I = char(I); If(ii) = str2num(I(3:end));
    L = char(L); Lf(ii) = str2num(L(3:end));
    Rx = strsplit(R{1},'/'); R1 = Rx{1}; Rf(ii,:) = [str2num(R1(3)) str2num(Rx{2})];
    Ax = strsplit(A{1},'/'); A1 = Ax{1}; Af(ii,:) = [str2num(A1(3)) str2num(Ax{2}) str2num(Ax{3})];
    eidout(ii) = eidx;
end

%==========================================================================
