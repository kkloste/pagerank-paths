G = load_graph('netscience-cc');
n = size(G,1);
epsmin=1/(3*1e6);
alpha=0.99;
for i=1:n
    [bestset1,cond1,cut1,vol1] = pprgrow(G,i);
    [bestset2,cond2,cut2,vol2] = ppr_fast_grid_mex(G,i,alpha, epsmin,.99, 0); 
    assert(cond2 <= cond1);
    fprintf('Passed vertex %i\n', i);
end
    