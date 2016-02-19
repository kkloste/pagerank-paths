%% load the livejournal graph
A = readSMAT('~/data/dsi/ljournal-2008-wcc.smat');
A = A|A';
A = A - diag(diag(A));

%%
rng(1);
ntrials = 100;
tfrac = zeros(ntrials,1);
for i=1:ntrials
    [~,~,~,~,dts] = pprgrow(A,randi(size(A,1),1));
    tfrac(i) = dts(end)/sum(dts);
end
[median(tfrac) mean(tfrac)]

%%
% On a Core i7-990X, with livejournal, I got 0.3618 as the mean fraction
% of time spent in the final run.