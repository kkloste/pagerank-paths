function rval=ppr_path_rho(A,vert,varargin)
% ppr_path_rho Compute data for the Personalized PageRank solution Path
%
% rval = ppr_path_rho(A,seeds)
%   Computes the PageRank solution path for alpha = 0.99, epsmin=1e-4,
%   and who = 0. Starting from the nodes "seeds"
%
%   If you wish to use weighted seeds, we allow integer weights by listing
%   a seed multiple times. 
%
% rval = ppr_path_rho(A,seeds,'key',value,'key',value...) allows the following
% optional key and value pairs:
% 
%   'alpha', <value> : the value of alpha to use
%   'epsmin', <value> : the smallest value of eps considered
%   'rho', <value> : the value of rho used
%   'degweights', true|{false} : a way to automatically degree weight a set
%      of seeds. Note that this will ignore any repetition among the seeds.
% 

parser = inputParser;
parser.addOptional('epsmin',1e-4);
parser.addOptional('alpha',0.99);
parser.addOptional('rho',0);
parser.addOptional('degweights',false);
parser.addOptional('debug',0);
parser.parse(varargin{:});

%parser.Results

if parser.Results.degweights
    seeds = unique(vert);
    % each seed shows up exactly d(seed) times
    % in the array vert.
    [~,vert] = find(A(:,seeds)); 
    vert = seeds(vert);
end

t0=tic;
[rval.step_stats,rval.ep_stats] = ppr_paths_rho_mex(A,vert,parser.Results.alpha,...
        parser.Results.epsmin, parser.Results.rho, parser.Results.debug);
rval.dt=toc(t0);

epstosteps=rval.ep_stats(:,[1,6]);
n = size(A,1);
rval.solution = @(eps) solution(eps,epstosteps,n,rval.step_stats);


function x=solution(eps,epstosteps,n,steps)
    ind=epstosteps(find(epstosteps(:,1) < eps,1,'first'),2)+1;
    x=accumarray(steps(1:ind,3),steps(1:ind,7),[n,1]);
