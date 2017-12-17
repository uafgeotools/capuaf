function [kappa,theta,sigma,F,Mw,delta,gamma,M0,M] = read_CAP_log(filename)
% SCRIPT SHOULD:
% INPUT log file name (grid or random) - saved in mat format (its faster!)
% OUTPUT [kappa,theta,sigma,misfit,Mw,delta,gamma,M0,M]
%
% Vipul Silwal
% 29 July, 2015

% FUTURE : 
% 
filename
load(filename);

% read the parameters
kappa = file(:,1);      % strike
sigma = file(:,3);      % rake
h = cosd(file(:,2));
theta = file(:,2);      % dip
F = file(:,4);
Mw = unique(file(:,5));
delta = file(:,6);
gamma = file(:,7);
M0 = file(:,8);
% get moment tensors
M = file(:,9:14)';
% conevrt M (from Aki&Richard to GCMT)
M = convert_MT(2,1,M);
