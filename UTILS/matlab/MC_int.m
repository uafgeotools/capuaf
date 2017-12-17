function [x_mid, int_y, cum_y, N] = MC_int(x,y,dx)
% 2 D monte carlo integration
% INPUT
%   x   x values
%   y   y values
%   dx  step size in x
% OUTPUT - here is the output 
%   x_mid   midpoints in x vectorized array (for plotting)
%   int_y   integration over x
%   cum_y   cumsum(int_y) - cumulative integral
%   N       Number of elements in each bin
%
% P.S. For normalizing purposes check (or make) -> [sum(int_y)*dx = 1]
% 
% Vipul Silwal
% Jul 21, 2015
% 

% sort x and y (Not needed)
[x,indx] =sort(x);
y = y(indx);

% create wedges for integration
x_min = min(x);
x_max = max(x);
x_min = 0;
x_max = pi;
x_vec = x_min:dx:x_max;
x_mid = x_min+dx/2:dx:x_max-dx/2;

% pre-create vectors
N = zeros(length(x_mid),1);
int_y = zeros(length(x_mid),1);

% Monte-carlo integration algorithm
for ii = 1:length(x_vec)-1
    indx = find(and(x >= x_vec(ii),x <= x_vec(ii+1)));
    N(ii) = length(indx); % number of elements in each dx slot
    int_y(ii) = sum(y(indx)); 
    % use this only if N is not distributed homogeneously
    % int_y(ii) = (dx/N(ii))*sum(y(indx));  
end

% cumulative integral 
cum_y = cumsum(int_y);
end
