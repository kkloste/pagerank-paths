function [cond cut vol inds] = cut_cond(A,set, degrees)
% function [cond cut vol inds] = cut_cond(A,set, degrees)

inds = set(find(set)); % get rid of zeros
inds = unique(inds);
if nnz(inds) == 0,
   cond = 1;
   cut = 0;
   vol = 0;
   return;
end

if nargin < 3, degrees = full(sum(A,2)); end

n = size(A,1);
eS = sparse( inds, 1, 1, n, 1);
interior_vol = full(eS'*(A*eS));
if size(degrees,1) > size(degrees, 2), vol = full(degrees'*eS);
else, vol = full(degrees*eS);
end
totvol = nnz(A);

cut = vol - interior_vol;
temp_vol = min(vol, totvol-vol);
if temp_vol == 0,
    cond = 1;
    return;
end
cond = cut/temp_vol;