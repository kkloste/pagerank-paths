function [times conds seeds] = test_ppr_grid(numtrials, filename, eps, alpha, theta, debugflag, newseeds)
% [times conds seeds] = test_ppr_grid(numtrials, filename, eps, alpha, theta, debugflag, newseeds)
% Set debugflag to 1 to turn on messages in the matlab and mex code.
%
% if no inputs are given, then settings are as follows:
% A = netscience-cc
% eps = 1e-5
% alpha = .99
% theta = .8
% debugflag = 0
% numtrials = 100
% newseeds = [1:n]
%
% examples:
% [times conds seeds] = test_ppr_grid(1);
% [times conds seeds] = test_ppr_grid(1000,'pgp-cc',1/3000000,.99,.8);
% [times conds seeds] = test_ppr_grid(1000,'itdk0304-cc',1/3000000,.99,.8);
% [times conds seeds] = test_ppr_grid(1,'netscience-cc',1e-4,.95,.9,2);
%
% To simply call the mex function directly, try
% [bset,cond,cut,vol,eps_stats,sol_track] = ppr_grid_mex(A,set,alpha,eps_min,theta,debugflag);
% [bset,cond,cut,vol,eps_stats,sol_track] = ppr_grid_mex(A,set,alpha,eps_min,theta,debugflag)
% [bset,cond,cut,vol,eps_stats,sol_track] = ppr_grid_mex(A,1,0.99,1e-5,0.8,debugflag)
%

addpath ../util ;   % for loading, and cut-set tools

tic;
% check inputs
if nargin < 1, numtrials = 1; end
if nargin < 2, filename = 'netscience-cc'; end
if nargin < 3, eps  = 1e-5; end
if nargin < 4, alpha = .99; end
if nargin < 5, theta = .8; end
if nargin < 6, debugflag = 0; end

assert(eps > 0 && eps <= 1, 'eps violates 0<eps<=1');
assert(alpha > 0, 'alphat violates alpha>0');

% setup inputs
A = load_graph(filename,'../data'); n = size(A,1);

if nargin < 7,
    seeds = [1:numtrials];
end
if (nargin == 7), seeds = newseeds; end
numtrials = length(seeds);
times = zeros(numtrials,1);
tcuts = zeros(numtrials,1);
tvols = zeros(numtrials,1);
conds = zeros(numtrials,1);
setuptime = toc

tic;
if debugflag >= 1, fprintf('test_ppr_grid_mex: setup time=%f \n', setuptime); end
% Test on random seeds
for trial_num=1:numtrials
    if debugflag>=1, fprintf('test_ppr_grid_mex:  start rand trial=%i \n', trial_num); end
    tic; [bset1 condr cut vol] = ppr_grid_mex(A,seeds(trial_num),alpha, eps, theta, debugflag);
    times(trial_num) = toc;
    if debugflag >= 1, fprintf('test_ppr_grid_mex:  end rand trial=%i \n', trial_num); end

    [cond1 cut1 vol1 inds1] = cut_cond(A,bset1);
    if (nnz(inds1)~=nnz(bset1)),
        fprintf('bset1 has repeated entry on seed=%i , diff=%i\n', seeds(trial_num) , abs(nnz(bset1)-nnz(inds1)));
    end

    if (cond1~=condr),
        fprintf('truecond=%f ~= reportedcond=%f on seed=%i \n', cond1, condr, seeds(trial_num) );
        fprintf('\t truevol=%i ~= reportedvol=%i \n', vol1, vol);
        fprintf('\t truecut=%i ~= reportedcut=%i \n', cut1, cut);
    end
    conds(trial_num) = cond1;
end
averuntime = sum(times)/numtrials, avecond = sum(conds)/numtrials
end