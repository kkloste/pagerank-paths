function grid_v_grow_path(filename, numtrials, load_dir)
% function grid_v_grow_path(filename, numtrials, load_dir)
%
% Compares two ppr-grow, grid, and path versions in runtime
% on dataset = input 'filename'.
%
% CHANGE default_load_dir IN THIS FILE
% If load_dir is specified, the dataset will be loaded from that directory,
% instead of the default directory.
%
% Call from [project]/experiments/timing_experiments

save_directory = './results/';
outputname = strcat('grid_v_grow_path_',char(filename));
% default_load_dir = '../../data/';
default_load_dir = '/scratch/dgleich/kyle/';
load_directory = default_load_dir;
if nargin == 3,
    load_directory = load_dir;
end
load(strcat(load_directory,char(filename)));

addpath ../.. ;       % for mex code
addpath ../../util ;   % for tools

n = size(A,1);
Gvol = nnz(A);

if nargin < 2, numtrials = 10; end

epsmin = 1/3000000;
alpha = 0.99;
rho = 0;

seeds = randi(n,2*numtrials,1);
seeds = unique(seeds);
seeds = seeds(1:numtrials);

thetas = [0.99, 0.817, 0.66];
num_thetas = numel(thetas);
times = zeros(2+num_thetas, numtrials);
conds = zeros(2+num_thetas, numtrials);
vols = zeros(2+num_thetas, numtrials);
cuts = zeros(2+num_thetas, numtrials);


fprintf('     p-cond  .99-cond  .817-cond  .66-cond  path || p-time  .99-time  .817-time  .66-time   path \n');    
for id=1:length(seeds),
fprintf('#%i  ', id);    
    vert = seeds(id);
    
    algorithm_id = 1;
    tic; [bestset1,cond1,cut1,vol1] = pprgrow(A,vert); pprtime = toc;
    [cond1 cut1 vol1 inds] = cut_cond(A,bestset1);
    times(algorithm_id,id) = pprtime;
    conds(algorithm_id,id) = cond1;
    vols(algorithm_id,id) = vol1;
    cuts(algorithm_id,id) = cut1;
    fprintf('%.4f  ', cond1);

    for algorithm_id = 1:3,
        theta = thetas(algorithm_id);
        tic; [bestset1,cond1,cut1,vol1] = ppr_fast_grid_mex(A,vert,alpha, epsmin, theta, 0); fastgridtime = toc;
        fprintf('%.4f  ', cond1);
        % record data for original ppr_grid_mex code
        [cond1 cut1 vol1 inds] = cut_cond(A,bestset1);
        times(1+algorithm_id,id) = fastgridtime;
        conds(1+algorithm_id,id) = cond1;
        vols(1+algorithm_id,id) = vol1;
        cuts(1+algorithm_id,id) = cut1;
    end
    
    algorithm_id = 5;
    tic; [step_stats,ep_stats,bestset1] = ppr_paths_rho_mex(A,vert,alpha,epsmin,rho);  pathtime = toc;
    [cond1 cut1 vol1 inds] = cut_cond(A,bestset1);
    fprintf('%.4f  ', cond1);
    times(algorithm_id,id) = pathtime;
    conds(algorithm_id,id) = cond1;
    vols(algorithm_id,id) = vol1;
    cuts(algorithm_id,id) = cut1;
    num_eps(id) = size( ep_stats, 1 );
    
    fprintf(' || %.3f    %.3f     %.3f    %.3f   %.3f  \n', pprtime, times(2,id), times(3,id), times(4,id), pathtime );
end

save( [save_directory outputname '.mat'] ,'filename','num_eps','n', ...
    'Gvol','seeds','conds','vols','cuts','times','-v7.3');
clear;
exit;

