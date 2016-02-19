function cutvals = fastcutsweep(A,p)
% FASTCUTSWEEP use a mex code for the cusweep
%
% cutvals = cutsweep(A,p)
%
% cutvals.cut = cut(S)
% cutvals.vol = vol(S)
% cutvals.edges = edges(S) ( = vol(S) - cut(S))
% cutvals.ncut = cut(S)/vol(S) + cut(S)/vol(notS)
% cutvals.conductance = cut(A,B)/min(vol(A),vol(B))
% cutvals.modularity = 
% cutvals.expansion = 
% cutvals.qcut = cut(S)/min(|S|,n-|S|)
%

n = size(A,1);
d = full(sum(spones(A),2));
np = numel(p);

[cuts volG] = graphcutsweep_val(sparse(A),p);

cutvals = struct;
cutvals.cut = cuts(end-1,:)';
cutvals.vol = cuts(end,:)';

cutvals.edges = cutvals.vol - cutvals.cut;

cutvals.conductance = cutvals.cut./min(cutvals.vol,volG-cutvals.vol);
cutvals.ncut = cutvals.cut./cutvals.vol + cutvals.cut./(volG-cutvals.vol);
cutvals.modularity = (1/volG) * (cutvals.edges - (1/volG)*cutvals.vol.^2);
cutvals.qcut = cutvals.cut./min(1:np,n-(1:np))';

