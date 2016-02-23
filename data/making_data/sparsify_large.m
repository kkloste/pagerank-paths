function [S,thresh] = sparsify_large(W,type,opt)
% SPARSIFY_LARGE Keep only the large entries of a matrix.
% 
% S = sparsify_large(W,type,opt) keeps only the largest entries of a 
% a matrix to form a weighted graph S. There are a few different types 
% of sparification. The parameter option is a type-dependent threshold.
%
% Assumptions on W:
%   W is sqaure.
%   Diagonal entries of W don't matter (and should be excluded from the
%   computations.)
%   Only the non-zero entries of W matter (so we can be sparse if possible)
%
% All operations are applied on the non-zeros of W.
%
%   type='percentile', opt=the percentile level (from 0,100), reasonable
%   choices for the percentile are 99, 95, or 90.
%
%   type='connected', opt=the size of the bisection interval required
%   to terminate.
%
%   type='threshold', opt=the threshold to use to sparisfy
%
%   type='range', opt=the fraction of the range. In this case, we take
%   the entries of the matrix and remap them to [0,1], The optional value
%   is 
%
%   type='below-max', opt=the percentage allowed below the maximum entry
%   so if opt = 0.95, then we'll allow any entry within 5% of the maximum
%   entry.
%
%   type='nearest', opt=the number of nearest neighbors
%   
%   type='neighbor-max', opt=the fraction of the maximum of each neighbor
%   to include.
%
%   type='neighbor-prctile', opt=the percentile of each neighbor list to
%   include. So if opt=95, then we include any entry above the 95th
%   percentile of each neighbor
%
% Note that opt can either be a fraction or a percentage. So opt=50 will be
% the fraction .5 or the percentage 50. There are equivalent commands.
%
% 

% TODO Modify this function to work for rectangular matrices too
% TODO Optimize some of the computations here, but probably 
% not worth it unless we run into scalability issues, and even
% then, more memory is easy-ish

n = size(W,1);

% replace diagonals with NaNs so they drop out of the computation
Wd = W + diag(sparse(NaN*ones(size(W,1),1)));
Wnd = W - diag(diag(W));

switch type
    case 'threshold'
        thresh = opt;
        S = spfun(@(x) x.*(x >= thresh), Wnd);
        
    case 'prctile'
        [~,p] = interpret_opt(opt);
        thresh = prctile(nonzeros(Wd),p);
        S = spfun(@(x) x.*(x >= thresh), Wnd);
        
    case 'below-max'
        f = interpret_opt(opt);
        M = max(nonzeros(Wd));
        thresh = f*M;
        S = spfun(@(x) x.*(x >= thresh), Wnd);
        
    case 'range'
        f = interpret_opt(opt);
        M = max(nonzeros(Wd));
        m = min(nonzeros(Wd));
        thresh = f*(M-m)+m;
        S = spfun(@(x) x.*(x >= thresh), Wnd);
        
    case 'nearest'        
        k = opt;
        [V,I] = sort(Wd,'descend');
        V = V(2:k+1,:);
        thresh = min(V);
        R = I(2:k+1,:);
        C = repmat(1:n,k,1);
        S = sparse(R(:),C(:),1,n,n);
        S = S|S';
        S = S.*Wnd; % get the entries out of S symmetrically.
        
    case 'neighbor-max'
        % find an element within some fraction of the maximum entry in each
        % column
        
        f = interpret_opt(opt);
        M = max(Wd)'; % get the max on each column
        thresh = f*M;
        [wi,wj,wv] = find(Wnd); % find all the entries
        S = sparse(wi,wj,wv>=thresh(wj),n,n);
        S = S|S';
        S = S.*Wnd;
        
    case {'neighbor-percentile','neighbor-prctile'}
        Wd(W==0) = NaN;
        [~,p] = interpret_opt(opt);
        thresh = prctile(Wd,p)';
        [wi,wj,wv] = find(Wnd); % find all the entries
        S = sparse(wi,wj,wv>=thresh(wj),n,n);
        S = S|S';
        S = S.*Wnd;
end

function [f,p] = interpret_opt(opt)
assert(opt > 0);
if opt < 1
    f = opt;
    p = 100*opt;
else
    f = opt/100;
    p = opt;
end


    