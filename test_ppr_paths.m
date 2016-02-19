function [times conds seeds] = test_ppr_paths(numtrials, filename, eps, alpha, debugflag)
% [times conds seeds] = test_ppr_paths(numtrials, filename, eps, alpha, debugflag)
% Set debugflag to 1 to turn on messages in the matlab and mex code.
%
% If no inputs are given, then settings are as follows:
% A = netscience-cc
% eps = 1e-5
% alpha = .99
% debugflag = 0
%
%   For mex only:
%     [step_stats,eps_stats,bestset] = ppr_paths_mex(A,set,alpha,eps_min,debugflag);
%     [step_stats,eps_stats,bestset] = ppr_paths_mex(A,set);
%

addpath ../util ;   % for loading, and cut-set tools

tic;
% check inputs
if nargin < 1, numtrials = 1; end
if nargin < 2, filename = 'netscience-cc'; end
if nargin < 3, eps  = 1e-4; end
if nargin < 4, alpha = .99; end
if nargin < 5, debugflag = 0; end

assert(eps > 0 && eps <= 1, 'eps violates 0<eps<=1');
assert(alpha > 0, 'alpha violates alpha>0');

A = load_graph(filename,'./data'); n = size(A,1);

seeds = randi(n,numtrials, 1);
times = zeros(numtrials,1);
conds = zeros(numtrials,1);

setuptime = toc
for j=1:length(seeds),
    tic;
    [step_stats,eps_stats,bestset] = ppr_paths_mex(A,seeds(j),alpha, eps, debugflag);
    times(j) = toc;
    
    % Now use info to reconstruct the set of best conductances
    [Y, inds] = min(step_stats(:,8));
    bcond_step = min(inds);
    [cond2 cut2 vol2 inds2] = cut_cond(A,bestset);
    condr = step_stats(bcond_step,8); % get best global cond
    if (cond2~=condr),
        fprintf('truecond=%f ~= reported bcond=%f on seed=%i \n', cond2, condr, seeds(j) );
    end

        
    rank = zeros(n,1);
    for step = 1:bcond_step,
        node = step_stats(step,3);
        Istart = step_stats(step,1);
        Iend = step_stats(step,2);
        if rank(Istart) == step_stats(step,5)+1, % new node if it starts at the end position
            rank(Istart) = node;
        end
        % now swap:
        if (Iend ~= Istart),
            rank(Iend+1:Istart) = rank(Iend:Istart-1);
            rank(Iend) = node;
        end
    end

    % Now solution vector rankings at step bcond_step
    % has been reconstructed.
    ep_ind = find(eps_stats(6,:)==bcond_step);
    bset1 = rank(1:eps_stats(ep_ind,5)); %get setsize of best cond from eps_stats
    condr = eps_stats(ep_ind,2); % get conductance recorded then in eps_stats
    
    % Now test best set -- make sure it matches what ppr_paths has recorded as the actual stats
    [cond1 cut1 vol1 inds1] = cut_cond(A,bset1);
    if (nnz(setdiff(inds1,inds2)) ),
        fprintf('bestset differs from reconstructed cluster seed=%i , diff=%i\n', seeds(j) , nnz(setdiff(inds1,inds2)) );
    end

    conds(j) = cond1;
    if (cond1~=condr),
        fprintf('reconstructed cond=%f ~= reported best cond=%f on seed=%i \n', cond1, condr, seeds(j) );
    end

end

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
%       start indices
%       end indices
%       node labels
%       degree of node pushed
%       size of solution vector
%       size of residual
%       value pushed
%       global best conductance
%
%   Note that you can use start_indices, end_indices, and node_labels
%   to reconstruct the PPR-rankings after every single push of the algorithm.
%   This, in turn, can be used to recosntruct the sets of best conductance
%   after any given step.


