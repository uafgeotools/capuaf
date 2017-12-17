function [omega,xi] = dc2omega(dc1,dc2)
% returns the omega angle difference between two fault parameters
% If 2 arrays of fault solutions are eneterd, then the omega returned will
% be between the corresponding enteries (i.e. two arrays must be of equal
% lengths)
% If length(dc1!=1) and length(dc2==1), omega will between all enteries of
% dc1 and dc2

[M1,N1] = size(dc1)
[M2,N2] = size(dc2)

if M2==1
    dc2 = repmat(dc2,M1,1)
end

[MDC1,k1,d1,n1,p1,p2,p3] = dcfaultpar2CMT(dc1(:,1),dc1(:,2),dc1(:,3),0);
[MDC2,k2,d2,n2,p12,p22,p33] = dcfaultpar2CMT(dc2(:,1),dc2(:,2),dc2(:,3),0);

[omega,xi] = CMT2omega_xi0(MDC1,MDC2,0,0);