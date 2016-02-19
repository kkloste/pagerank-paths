% call from /code/experiments/senate/
%
% This script used to identify node's terms and location
% in the network, e.g. core v periphery

clc; clear;
load ../../data/all_senators3knn_layouts.mat
load ../../data/all_senators3knn.mat
addpath ../..;

%%
num_terms = size(attendance,2);

terms309 = find(attendance(309,:));
other_senators = find( sum( attendance(:,terms309), 2 ) );

spattn = spones(attendance);
terms_per_sen = sum(spattn, 2);

single_termers = find( terms_per_sen == 1 );

single_termers_concurrent309 = intersect( single_termers, other_senators );
seeds = single_termers_concurrent309;
length(seeds)

seeds = [309,310];
seed_terms = attendance(seeds,:)';


%%
%taus = [3*1e-4, 1e-4, 3*1e-5];
taus = [3*1e-4];
layout_num = 4;
xy = layouts{layout_num};

% Colors
cmappart = hot(9);
cmap = @() colormap(flipud(cmappart(2:6,:)));
crange = [-3,0];
cmap();
colors=colormap;
ncolors = size(colors,1);
%
for seed_num=1:length(seeds),
    clf
    
    seed = seeds(seed_num);

    [px,py] = gplot(A,xy);
    gdraw = @() plot(px,py,'k-','LineWidth',0.4,'Color',0.55*[1,1,1]);

    sdraw = @() plot(xy(S,1),xy(S,2),'ko','MarkerSize',7);
    bigmarkers = 15;
    smallmarkers = 2;
    figsize = [2.75,2.75];

    for which_tau=1:length(taus),
        tau = taus(which_tau);
        rval = ppr_path_rho(A,seed,'epsmin',tau,'alpha',0.99,'rho',0.7);
        n = size(A,1);
        xfinal = accumarray(rval.step_stats(:,3),rval.step_stats(:,7),[n,1]);
        xinds = find(xfinal);

        gdraw(); hold on;

        % draw senator's co-term-members
        sen_inds = find( attendance*seed_terms(:,seed_num) );
        scatter(xy(sen_inds,1), xy(sen_inds,2), 5, [0,0,0], 'filled'); cmap();
        plot(xy(seed,1),xy(seed,2),'ko','MarkerSize',10);
        plot(xy(309,1),xy(309,2),'ko','MarkerSize',10);
        axis off;
print(gcf, ['./senate_figures/senate-layout-eps', num2str(log10(1/tau)) ,'-s', num2str(seed), '.png'], '-dpng','-r600');
    end
end
%%
        