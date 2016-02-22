% Plot PageRank solutions paths for netscience dataset
% for different values of the parameter rho.
clear; clc;
load ../../data/netscience-cc; % must change directory to point to location of dataset

% Ensure adjacency matrix is symmetrized and binary
% (can skip if data is pre-processed)
A = A|A';
A = A - diag(diag(A));
n = size(A,1);
addpath ../..; %for ppr_path tools
addpath ../../util; % for set_figure_size.m

%%
seed = 211;

alpha1 = 0.99;
epsmin=1e-5;
rhos = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99];

%%

for j=1:length(rhos),
    rho = rhos(j);
    tic;
    rval = ppr_path_rho(A,seed,'epsmin',epsmin,'rho',rho);
    timet = toc; fprintf( '\nNum eps = %f,  time=%f', length(rval.ep_stats), toc);
end