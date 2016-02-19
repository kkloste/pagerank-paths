% test on itdk0304
% test on fbone
% test on twitter?

if exist('A','var') && exist('B','var')
    if nnz(A)/2 ~= 4404989
        clear A
    end
    
    if nnz(B)/2 ~= 607610
        clear B
    end
end

if ~exist('A','var') || ~exist('B','var')
    A=readSMAT('~/data/graph-db/traces/anony-interactions-oneyearA-cc.smat');
    A=A-diag(diag(A));
    
    B=readSMAT('~/data/router/itdk0304-cc.smat');
    B=B-diag(diag(B));
end
    
compile

%%
rng(1); % reset random seed
Graphs = {A,B};
numtrials = 5;
epsmin = 1/3000000;
alpha = 0.99;

for gi=1:numel(Graphs)
    G = Graphs{gi};
    n = size(G,1);
    Gvol = nnz(G);
    
    seeds = randi(n,2*numtrials,1);
    seeds = unique(seeds);
    seeds = seeds(1:numtrials);
    
    curtimes = [];
    newtimes = [];
    

    for id=1:length(seeds)
        vert = seeds(id);
    
        tic; [bestset1,cond1,cut1,vol1] = pprgrow(A,vert); reftime = toc;
        tic; [bset2,cond2,cut2,vol2] = ppr_grid_mex(A,vert,alpha, epsmin, .99, 0); curtimes(end+1) = toc/reftime;
        tic; [bset2,cond2,cut2,vol2] = ppr_grid_mex(A,vert,alpha, epsmin, .817, 0); curtimes(end+1) = toc/reftime;
        tic; [bset2,cond2,cut2,vol2] = ppr_grid_mex(A,vert,alpha, epsmin, .66, 0); curtimes(end+1) = toc/reftime;

        tic; [bset2,cond2,cut2,vol2] = ppr_fast_grid_mex(A,vert,alpha, epsmin, .99, 0); newtimes(end+1) = toc/reftime;
        tic; [bset2,cond2,cut2,vol2] = ppr_fast_grid_mex(A,vert,alpha, epsmin, .817, 0); newtimes(end+1) = toc/reftime;
        tic; [bset2,cond2,cut2,vol2] = ppr_fast_grid_mex(A,vert,alpha, epsmin, .66, 0); newtimes(end+1) = toc/reftime;
        
        fprintf('%6s    %5.2f  %5.2f  %5.2f    %5.2f  %5.2f  %5.2f\n', ' ', ...
            curtimes(end-2:end), newtimes(end-2:end));
    end
    curtimes = reshape(curtimes,3,numtrials)';
    newtimes = reshape(newtimes,3,numtrials)';
    
    fprintf('%6s    %5.2f  %5.2f  %5.2f    %5.2f  %5.2f  %5.2f\n', '80\%', prctile(curtimes, 80), prctile(newtimes, 80));
    fprintf('%6s    %5.2f  %5.2f  %5.2f    %5.2f  %5.2f  %5.2f\n', '50\%', prctile(curtimes, 50), prctile(newtimes, 50));
    fprintf('%6s    %5.2f  %5.2f  %5.2f    %5.2f  %5.2f  %5.2f\n', 'mean', mean(curtimes), mean(newtimes));
    fprintf('%6s    %5.2f  %5.2f  %5.2f    %5.2f  %5.2f  %5.2f\n', '20\%', prctile(curtimes, 20), prctile(newtimes, 20));
    
end
    