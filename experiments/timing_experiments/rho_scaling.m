function rho_scaling(filename, numtrials, grid_size, load_dir)
% rho_scaling(filename, numtrials, rho, load_dir)
%
% Analyze how ppr_paths_rho_mex performs as rho varies
% Value of alpha is 0.99
% value of epsmin is 1e-5
%
% CHANGE default_load_dir IN THIS FILE
% If load_dir is specified, the dataset will be loaded from that directory,
% instead of the default directory.
%
% Call from [project]/experiments

save_directory = '~/ppr-all-eps/results/';
outputname = strcat('rho-scaling-',char(filename));
default_load_dir = '/scratch/dgleich/kyle/';

addpath ../..;       % for mex code
addpath ../../util ;   % for tools

alpha = 0.99;
epsmin = 1e-5;


load_directory = default_load_dir;
if nargin == 4,
    load_directory = load_dir;
end
load(strcat(load_directory,char(filename)));

if nargin < 2, numtrials = 5; end
if nargin < 3, grid_size = 5; end

n = size(A,1);
Gvol = nnz(A);

seeds = randi(n,2*numtrials,1);
seeds = unique(seeds);
seeds = seeds(1:numtrials);

rhos = [0:1:grid_size]./grid_size;
times = zeros(numtrials,grid_size);
conds = zeros(numtrials,grid_size);
num_eps = zeros(numtrials,grid_size);

fprintf(' %s \n', char(filename));

for trial_id=1:numtrials,
    node = seeds(trial_id);
    for grid_id=1:grid_size,
        tic;
        rho = rhos(grid_id);
        [step_stats,ep_stats] = ppr_paths_rho_mex(A,node,alpha,epsmin,rho);
        times(trial_id, grid_id) = toc;
        num_eps(trial_id, grid_id) = size(ep_stats,1);
        conds(trial_id, grid_id) = step_stats(end,8);
        fprintf('%i: %f \t', trial_id, times(trial_id, grid_id));
    end    
    fprintf('\n trial=%d\n', trial_id);
end



save( [save_directory outputname '.mat'] ,'filename','n','Gvol','seeds','conds','num_eps','times','-v7.3');
clear;
exit

% output from ppr_paths_rho_mex:
%
%   ep_stats
%       epsilons
%       conductances
%       cuts
%       volumes
%       setsizes
%       stepnum
%
%   step_stats
%       start index
%       end index
%       node label
%       degree of node pushed
%       size of solution vector
%       size of residual
%       value pushed
%       global best conductance


