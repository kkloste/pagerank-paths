%% Set up experiment for many-term senator
% paths, labelled and non.
% call from /code/experiments/senate/

clear; clc;
addpath ../..; % for ppr_path functions
addpath ../../util;

% Get data
file_tag = '3knn';
filename = ['all_senators' file_tag];
load(['../../data/' filename]) ;
n = size(A,1);

%% Run paths

% some seeds with fewer_terms = [50:55]; 
% some seeds with more_terms = [1810:1814];
% Also looked at these: [309, 310, 624, 982];
%taus = [ 3*1e-4, 1e-4, 3*1e-5]; % checking out most useful parameters

try_these = [309,310];
taus = [3*1e-5];
rho = 0.9;
alpha = 0.99;

for j=1:length(try_these),
    for which_tau = 1:length(taus),
        tau = taus(which_tau);
        conds = senate_paths(A, try_these(j), attendance, 'alpha', alpha, 'epsmin', tau, 'rho', rho);
    end
end

fprintf('\nDone with senate trials\n');