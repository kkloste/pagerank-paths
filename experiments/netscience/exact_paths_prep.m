% Call from /experiments/netscience/;
%
% This script was used in preparation for the netscience exact_paths experiment,
% to determine which seed and which value of rho to use.
clear; clc;
load ../../data/netscience-cc
n = size(A,1);
addpath ../..;
addpath ../../util;

A = A|A';
A = A - diag(diag(A));

alpha = 0.99;
epsmin=1e-5;

rhos = [0,0.1,0.25,0.5,0.6,0.7,0.8,0.9];
rho = 0.1; % 0.25, 0.5, 0.6, 0.9
seed = 211; % 211 is used in paper hand-picked for a good image
% If you want more than 1 seed:
% neighb = find(A(:,seed)); neighb = union( neighb, seed );
% seeds = neighb( randi(length(neighb), 3,1) );
seeds = seed;

num_rhos = length(rhos);
for j=1:num_rhos
    rho = rhos(j);
%
    rval = ppr_path_rho(A,seeds,'epsmin',epsmin,'alpha',alpha,'rho',rho);

    % find the set of non-zeros in final vector and build a local index
    xfinal = accumarray(rval.step_stats(:,3),rval.step_stats(:,7),[n,1]);
    xnnz = find(xfinal);
    xinds = zeros(n,1);
    xinds(xnnz) = 1:numel(xnnz);
    numel(xnnz);

    %
    X = zeros(numel(xnnz),size(rval.ep_stats,1));

    newconds = [];
    mincond = Inf;
    bsetthresh = zeros(size(rval.ep_stats,1),1);
    for i=1:size(rval.ep_stats,1)
        ep = rval.ep_stats(i,1);
        step = rval.ep_stats(i,6)+1;
        xvec=accumarray(rval.step_stats(1:step,3),rval.step_stats(1:step,7),[n,1]);
        X(:,i) = xvec(xnnz);
        [~,xperm] = sort(X(:,i),'descend');
        % find the value of x such that we are in/out of the best set
        bsetthresh(i) = X(xperm(rval.ep_stats(i,5)),i);
        curcond = rval.ep_stats(i,2);
        if curcond < mincond
            newconds(end+1,:) = [ep,curcond];
            mincond = curcond;
        end
    end

    x_paths = zeros(n,1);
    x_paths(xnnz) = X(:,end);

    %

    %
    % Colormap for paths
    clf;
    cmappart = hot(9);
    cmap = @() colormap(flipud(cmappart(2:6,:)));

    hold on;
    xvalind = find(rval.ep_stats(:,1) < 1e-4,1,'first');
    hs=loglog((1./rval.ep_stats(:,1)), X','color','k');
    plot(1./rval.ep_stats(:,1),bsetthresh,'k','LineWidth',2);
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    crange = [-3,0];
    cmap();
    colors=colormap;
    ncolors = size(colors,1);
    for i=1:numel(hs)
        cval = log10(X(i,xvalind));
        cind = round(((cval-(crange(1)))/(crange(2)-crange(1)))*(ncolors-1)+1);
        cind = max(cind,1); % cap to range
        cind = min(cind,ncolors);
        set(hs(i),'Color',colors(cind,:),'LineWidth',0.3);
    end

    xl = [1e1,1e5];
    yl = [1e-5,1e0];
    line([xl(1),1./epsmin],[(1-rho)./xl(1),(1-rho)*epsmin],'Color','k','LineWidth',1);
    xlabel('1/\epsilon');
    ylabel('Degree normalized PageRank');
    box off;

    xlim(xl);
    set(gca,'XTick',[10,100,1000,10000,100000]);
    set_figure_size([3.5,3]);
    dummy = num2str(rho);
    if length(dummy)>1, dummy = dummy(3:end);
    else dummy = '0';
    end
    title( sprintf('rho = %.2f', rho) );
    print(gcf, ['./netsci-path-rho', dummy,'-s', num2str(seed), '.png'],'-dpng','-r600');
end
