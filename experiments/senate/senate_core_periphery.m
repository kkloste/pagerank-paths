% call from /code/exper_new/senate/
clc; clear;
load ../../data/all_senators3knn_layouts.mat
load ../../data/all_senators3knn.mat

num_terms = size(attendance,2);
%seed_terms = attendance(seeds,:)';
n = size(A,1);

layout_num = 4;
xy = layouts{layout_num};

cmappart = hot(9);
cmap = @() colormap(flipud(cmappart(2:6,:)));
colors=colormap();
%colors = cmappart(flipud(cmappart(1:30,:)));
ncolors = size(colors,1);
%%
    clf

    [px,py] = gplot(A,xy);
    gdraw = @() plot(px,py,'k-','LineWidth',0.4,'Color',0.55*[1,1,1]);

    %sdraw = @() plot(xy(S,1),xy(S,2),'ko','MarkerSize',7);
    bigmarkers = 15;
    smallmarkers = 2;
    figsize = [2.75,2.75];

    gdraw(); hold on;

    % set color
    cvals = zeros(n,3);
    W = spones(attendance);
    number_terms = full( sum( W, 2 ) );
    maxnum = log(max(number_terms));
    
    for cj=1:length(cvals),
        cind = ceil(ncolors*log(number_terms(cj))/maxnum)
        cind = max(cind,1); % cap to range
        cind = min(cind,ncolors);
        cvals(cj,:) = colors(cind,:);
    end
    
    scatter(xy(:,1),xy(:,2),bigmarkers,cvals,'filled'); cmap();

    % draw senator's co-term-members
%     sen_inds = find( attendance*seed_terms(:,seed_num) );
%     scatter(xy(sen_inds,1), xy(sen_inds,2), 5, [0,0,0], 'filled'); cmap();
%     plot(xy(seed,1),xy(seed,2),'ko','MarkerSize',10);
    axis off;
    print(gcf, ['./senate_figures/senate-core-periphery.png'], '-dpng','-r600');
