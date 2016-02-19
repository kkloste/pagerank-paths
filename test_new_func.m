function [times conds seeds] = test_new_func(numtrials, filename, eps, alpha, theta, debugflag, newseeds)
% [times conds seeds] = test_new_func(numtrials, filename, eps, alpha, theta, debugflag, newseeds)
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
% Test on random seeds
for trial_num=1:numtrials
    tic; [bset1 cond1 cut1 vol1] = ppr_grid_mex(A,seeds(trial_num),alpha, eps, theta, debugflag);

    tic; [bset2 cond2 cut2 vol2] = [NEW FUNCTION](A,seeds(trial_num),alpha, eps, theta, debugflag);
    
    times(trial_num) = toc;

    [condT cutT volT indsT] = cut_cond(A,bset2);
    if (nnz(indsT)~=nnz(bset2)),
        fprintf('new bestset repeats nodes on seed=%i , diff=%i\n', seeds(trial_num) , abs(nnz(bset2)-nnz(indsT)));
    end

    if (cond1~=condr),
        fprintf('new function misreports cond=%f ~= true_cond=%f on seed=%i \n', cond2, condT, seeds(trial_num) );
        fprintf('\t truevol=%i ~= reportedvol=%i \n', vol2, volT);
        fprintf('\t truecut=%i ~= reportedcut=%i \n', cut2, cutT);
    end
    if (cond2~=cond1),
        fprintf('new function conflicts with ppr_grid on seed=%i \n', seeds(trial_num) ); 
        fprintf('\t new cond=%f ~= grid cond=%f \n', cond2, cond1);        
    end
    conds(trial_num) = cond1;
end
averuntime = sum(times)/numtrials, avecond = sum(conds)/numtrials
end