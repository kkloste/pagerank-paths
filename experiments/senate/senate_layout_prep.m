% Call from /experiments/senate
% Used to identify which graph layout was the best
%
%
clc; clear; clf;
load ../../data/all_senators3knn_layouts.mat

%seeds = [308,309,310];
seeds = [309,310];
layout_nums = [1,2,3,4];

for j=1:length(layout_nums),
    for seed_num=1:length(seeds),
        clf
        
        layout_num = layout_nums(j);
        xy = layouts{layout_num};
        seed = seeds(seed_num);

        cmappart = hot(9);
        cmap = @() colormap(flipud(cmappart(2:6,:)));
        [px,py] = gplot(A,xy);
        gdraw = @() plot(px,py,'k-','LineWidth',0.4,'Color',0.55*[1,1,1]);

        sdraw = @() plot(xy(S,1),xy(S,2),'ko','MarkerSize',7);
        bigmarkers = 15;
        smallmarkers = 2;
        figsize = [2.75,2.75];
        %figarea = [0,220,-600,-120];

        d = sum(A,2);
        v1 = zeros(size(A,1),1); v1(seed) = 1; v1 = v1 / sum(v1); v1 = sparse(v1);
        x = (diag(d) - 0.99*A) \ (v1);
        gdraw(); hold on; 
        scatter(xy(:,1),xy(:,2),bigmarkers,log10(x),'filled'); cmap();
        plot(xy(find(v1),1),xy(find(v1),2),'ko','MarkerSize',7);
        axis off;
        print(gcf, ['./senate_figures/senate-layout-', num2str(layout_num) ,'-s', num2str(seed), '.png'], '-dpng','-r600');
    end
end