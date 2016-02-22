% call from /experiments/senate/
clc; clear;
load ../../data/all_senators3knn_layouts.mat
load ../../data/all_senators3knn.mat
addpath ../..; % for ppr_path scripts

%%
seeds = [309,310];
num_terms = size(attendance,2);
seed_terms = attendance(seeds,:)';

taus = [3*1e-4, 1e-4, 3*1e-5];
layout_num = 4; % This layout gives the most insightful display
xy = layouts{layout_num};

% Colors
cmappart = hot(9);
cmap = @() colormap(flipud(cmappart(2:6,:)));
crange = [-3,0];
cmap();
colors=colormap;
ncolors = size(colors,1);
%%
for seed_num=1:length(seeds),
    clf
    
    seed = seeds(seed_num);

    [px,py] = gplot(A,xy);
    gdraw = @() plot(px,py,'k-','LineWidth',0.2,'Color',0.55*[1,1,1]);

    sdraw = @() plot(xy(S,1),xy(S,2),'ko','MarkerSize',3);
    bigmarkers = 1;
    smallmarkers = 0.5;
    figsize = [2.75,2.75];

    for which_tau=1:length(taus),
        tau = taus(which_tau);
        rval = ppr_path_rho(A,seed,'epsmin',tau,'alpha',0.99,'rho',0.9);
        n = size(A,1);
        xfinal = accumarray(rval.step_stats(:,3),rval.step_stats(:,7),[n,1]);
        xinds = find(xfinal);

        gdraw(); hold on;
        
        % set color
        cvals = zeros( length(xinds), 3);
        for cj=1:length(cvals),
            cval = log10(xfinal(xinds(cj)));
            cind = round(((cval-(crange(1)))/(crange(2)-crange(1)))*(ncolors-1)+1);
            cind = max(cind,1); % cap to range
            cind = min(cind,ncolors);
            cvals(cj,:) = colors(cind,:);
        end

        scatter(xy(xinds,1),xy(xinds,2),bigmarkers,cvals,'filled'); cmap();
        
%         % If you want to draw senator's co-term-members:
%         sen_inds = find( attendance*seed_terms(:,seed_num) );
%         scatter(xy(sen_inds,1), xy(sen_inds,2), 5, [0,0,0], 'filled'); cmap();
        plot(xy(seed,1),xy(seed,2),'ko','MarkerSize',1.5);
        axis off;
        eps_str = num2str(log10(1/tau));
        eps_str = eps_str( [1, 3:end] );
        print(gcf, ['./senate_figures/senate-layout-eps', eps_str ,'-s', num2str(seed), '.png'], '-dpng','-r600');
    end
end
fprintf('\n Done with senate_layout.\n');