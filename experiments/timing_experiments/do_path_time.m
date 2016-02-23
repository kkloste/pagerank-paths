function do_path_time(filename, epsmin, numtrials, rho, load_dir)
% do_path_time(filename, epsmin, numtrials, rho, load_dir)
%
% Analyzing ppr_path for time and number of eps values
% on dataset = input 'filename'.
%
% CHANGE default_load_dir IN THIS FILE
% If load_dir is specified, the dataset will be loaded from that directory,
% instead of the default directory.
%
% Call from [project]/experiments/timing_experiments

save_directory = './results/';
outputname = strcat('path-time-',char(filename));
% default_load_dir = '../../data/';
default_load_dir = '/scratch/dgleich/kyle/';
load_directory = default_load_dir;
if nargin == 5,
    load_directory = load_dir;
end
load(strcat(load_directory,char(filename)));

addpath ../.. ;       % for mex code
addpath ../../util ;   % for tools

if nargin < 4,
    rho = 0;
end


n = size(A,1);
Gvol = nnz(A);

% numtrials = 100;
alpha = 0.99;
% epsmin = 1/3000000;

seeds = randi(n,2*numtrials,1);
seeds = unique(seeds);
seeds = seeds(1:numtrials);

times = zeros(numtrials,2);
conds = zeros(numtrials,2);
num_eps = zeros(numtrials,1);
    fprintf(' %s \n', char(filename));
for trial_id=1:numtrials,
    node = seeds(trial_id);
    alg_id = 1;
    tic;
    [step_stats,ep_stats] = ppr_paths_rho_mex(A,node,alpha,epsmin,rho);
    times(trial_id, alg_id) = toc;
    num_eps(trial_id) = size(ep_stats,1);
    conds(trial_id, alg_id) = step_stats(end,8);
    
    alg_id = 2;
    tic;
    [bestset1,cond1,cut1,vol1] = pprgrow_path_comp(A,node);
    times(trial_id, alg_id) = toc;
    [cond1 cut1 vol1 inds] = cut_cond(A,bestset1);
    conds(trial_id, alg_id) = cond1;
    
    fprintf('%i: %f \t', trial_id, times(trial_id, alg_id));    
end

rhost = num2str(rho);
lastchar = min( 3, length(rhost) );
rhost = rhost(1:lastchar);
save( [save_directory outputname '-rho' rhost '.mat'] ,'filename','rho','n','Gvol','seeds','conds','num_eps','times','-v7.3');
clear;
exit

% output from ppr_paths_mex:
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


