function [x,r,s] = pprfast(A, d, beta, v, tau, varargin)
% PPRFAST Fast Personalized PageRank with 1-norm regularization
% 
% This algorithm is designed for medium scale problems where you can
% factorize the PageRank linear system in order to solve systems of
% equations with it. 
%
% [x,r,s] = pprfast(A, d, beta, v, tau)
%
% Algorithm:
%   if the system is large
%     run 1/(tau*(1-beta)) iterations of the push method
%     if this solve it, quit!
%     if not, see if this was the right non-zero set by solving
%     for the components of the solution given this active set
%     if that still wasn't the solution, then...
%   (in either case)
%     solve the PageRank system with the 1-norm regularization
%     check if the "full" solution is correct
%     find the smallest vector that produces a strictly feasible dual
%     solution
%     reduce the dual solution until we arrive at a primal solution
%
% Optional parameters
%   'tol' the stopping tolerance
%   'force_push' always start with the push procedure
%   'skip_push' always skip the push procedure
%   'quiet' force the method to be quiet
%   'solver' a precomputed PageRank solver (e.g. via factorization)
%   'force_primal' force using the primal ascent method
%   'force_dual' force using the dual descent method
%     if both are forced, then we use the primal

parser = inputParser;
parser.KeepUnmatched = true;
parser.addParamValue('tol', sqrt(eps(1)), @(x) x > 0 && x <= 1);
parser.addParamValue('force_push', false, @islogical);
parser.addParamValue('skip_push', false, @islogical);
parser.addParamValue('force_dual', false, @islogical);
parser.addParamValue('force_primal', false, @islogical);
parser.addParamValue('quiet', false, @islogical);
parser.addParamValue('prsolver', @prsolve);
parser.parse(varargin{:});
tol = parser.Results.tol;
prsolver = parser.Results.prsolver;

flag = 1; % non-convergence flag, 

% two timers, t0 is start, t1 is increment, dt1 is print increment
t0 = tic; t1 = tic; print1 = 1; dt1 = 5;
if parser.Results.quiet 
    dt1 = Inf;
end


d = full(d);

% check for early return
n = size(A,1);
x = zeros(n,1);
r = (1-beta)*v;
s = tau*d - r;

% check for trivial return
if abs(x'*s) < tol && sum(max(r - tau*d,0)) < tol
    flag = 0;
    return;
end

npushes = 1/(tau*(1-beta));

if (npushes < n^2.5 || parser.Results.force_push) && ~parser.Results.skip_push
    % try using push to accelerate
    vs = sparse(v);
    As = sparse(A);
    % -1 is a sentinal here to denote that we actually want to run with
    % tol set to machine eps and we want to run for "1" total iterations.
    [~,~,~,~,x] = pprpush_weighted_mex(As, d, vs, tau, beta, -1-tol);
    r = (1-beta)*v + beta*(A'*(x./d)) - x;
    s = tau*d - r;
    
    if abs(x'*s)/n <= tol && sum(max(r - tau*d,0)) < tol
        % we have the solution!
        return
    end
    
    % otherwise, check if the current set of non-zeros is the right
    % solution
    active = x > 0;
    if sum(active) > 2000
        xactive = prsolveiter(A(active,active),d(active),beta,(1-beta)*v(active) - tau*d(active));
    else
        xactive = prsolve(A(active,active),d(active),beta,(1-beta)*v(active) - tau*d(active));
    end
    x = zeros(n,1);
    x(active) = xactive;
    
    r = (1-beta)*v + beta*(A'*(x./d)) - x;
    s = tau*d - r;
end

% check for trivial return
if abs(x'*s) < tol && sum(max(r - tau*d,0)) < tol
    flag = 0;
    return;
end


if parser.Results.force_primal == true || ...
        (issparse(A) && parser.Results.force_dual == false)
    
    % In this case, we want to use the primal ascent strategy.
    maxit = 10*n; % Kyle has a better bound here.

    rp = max(r - tau*d,0);

    for i=1:maxit
        x = x + rp;
        r = (1-beta)*v + beta*(A'*(x./d)) - x;
        rp = max(r - tau*d,0);
        s = tau*d - r;
        dt = toc(t1);
        if dt > dt1
            if print1 % print the header
                fprintf('  %5s  %11s  %11s  %11s  %8s\n', ...
                    'iter', 'sum rp', 'dual', 'sum x', 'time');
                fprintf('  -----  -----------  -----------  -----------  --------\n');
            end
            fprintf('  %5i  %11.3e  %11.3e  %11.3e  %8.1f\n', ...
                i, sum(rp), x'*s/n, sum(x), toc(t0));
            t1 = tic;
            print1 = 0;
        end
        
        if abs(x'*s)/n < tol && sum(rp) < tol
            flag = 0;
            break;
        end
    end
    
    if flag == 0
        if print1 == 0
            % final print
            fprintf('* %5i  %11.3e  %11.3e  %11.3e  %8.1f\n', ...
                i, sum(rp), x'*s/n, sum(x), toc(t0));
        end
        return
    end
end
    


V = prsolver(A,d,beta,[(1-beta)*v tau*d]);

xs = V(:,1);
z = V(:,2);

% check the full solution
xt = xs - z;
if xt >= 0
    rt = (1-beta)*v + beta*(A'*(xt./d)) - xt;
    st = tau*d - rt;
    if abs(x'*s)/n < tol && sum(max(rt - tau*d,0)) < tol
        x = xt;
        r = rt;
        s = st;
        return;
    end
end

% Run the inverse iteration
maxit = 10*n;

% Do a search to find the smallest term we can add to keep s >= 0
rho = 1; % rho = 1 is always valid
for i=1:10
    % so test rho/2
    xt = min(x + rho/2*xs,xs);
    rt = (1-beta)*v + beta*(A'*(xt./d)) - xt;
    st = tau*d - rt;
    if any(st < 0)
        break; % rho/2 isn't valid, so leave rho as is
    end
    rho = rho/2; % this rho is valid
    r = rt; % save the dual variables for the return
    s = st;
end
x = max(x + rho*xs,xs);
r = (1-beta)*v + beta*(A'*(x./d)) - x; % TODO could be avoided if I was a little more careful
s = tau*d - r;

for i=1:maxit
    % alg: 
    x = x - s;
    x = max(x,0);
    r = (1-beta)*v + beta*(A'*(x./d)) - x;
    s = tau*d - r;
    
    dt = toc(t1);
    if dt > dt1
        if print1 % print the header
            fprintf('  %5s  %11s  %11s  %11s  %8s\n', ...
                'iter', 'sum x', 'dual', 'min s', 'time');
            fprintf('  -----  -----------  -----------  -----------  --------\n');
        end
        fprintf('  %5i  %11.3e  %11.3e  %11.3e  %8.1fs\n', ...
            i, sum(x), x'*s/n, min(s), toc(t0));
        t1 = tic;
        print1 = 0;
    end

    % we maintain strict primal feasibility, so we just need to worry 
    % about slackness 
    if abs(x'*s)/n <= tol
        flag = 0;
        break;
    end
end

if print1 == 0
    if flag == 0
        fprintf('* %5i  %11.3e  %11.3e  %11.3e  %8.1fs\n', ...
            i, sum(x), x'*s/n, min(s), toc(t0));
    else
        fprintf('- %5i  %11.3e  %11.3e  %11.3e  %8.1fs\n', ...
            i, sum(x), x'*s/n, min(s), toc(t0));
    end
end

function X = prsolve(A,d,alpha,B)
% Solve a PageRank system

X = (diag(sparse(d)) - alpha*A) \ B;
for i=1:size(X,2)
    X(:,i) = X(:,i).*d;
end


function X = prsolveiter(A,d,alpha,B)
% Solve a PageRank system iteratively
maxiter = log(10*eps(1))/log(alpha);
X = B;
iD = repmat(1./d,1,size(B,2));
for i=1:maxiter
    Xnew = B + alpha*(A'*(B.*iD));
    errnorm = (1/(1-alpha))*norm(Xnew - X,1);
    X = Xnew;
    if errnorm < 10*eps(1), break; end
end