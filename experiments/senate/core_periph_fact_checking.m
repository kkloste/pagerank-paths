% call from /code/exper_new/senate/
clc; clear;
load ../../data/all_senators3knn_layouts.mat
load ../../data/all_senators3knn.mat

num_terms = size(attendance,2);
n = size(A,1);

layout_num = 4;
xy = layouts{layout_num};

cmappart = hot(9);
cmap = @() colormap(flipud(cmappart(2:6,:)));
colors=colormap();
ncolors = size(colors,1);
%%
seeds = [310];
seed_terms = attendance(seeds,:)';


for j=1:10,
    clf

    [px,py] = gplot(A,xy);
    gdraw = @() plot(px,py,'k-','LineWidth',0.4,'Color',0.55*[1,1,1]);

    bigmarkers = 15;
    smallmarkers = 2;
    figsize = [2.75,2.75];

    gdraw(); hold on;

    % set color

    sen_inds = find(attendance(:,j));
    scatter(xy(sen_inds,1), xy(sen_inds,2), bigmarkers, [0,0,0], 'filled'); cmap();
    %plot(xy(seed,1),xy(seed,2),'ko','MarkerSize',10);

    axis off;
   print(gcf, ['./senate_figures/senate-layout-term' num2str(j) '.png'], '-dpng','-r600');
end

%% NOW PATH PLOTS

seed = 310;
taus = [1e-2, 3*1e-3, 1e-3]; 

    for which_tau = 1:length(taus),
        tau = taus(which_tau);
        [bestset, cond] = ppr_fast_grid_mex(A, seed, 0.99, tau, 0.66);

                epsmin = 1e-5;
                alpha = 0.99;
                rho = 0.9;
                seed = 743479;

                n = size(A,1);
                rval = ppr_path_rho(A,seed,'epsmin',epsmin,'rho',rho,'alpha',alpha);

                % find the set of non-zeros and build a local index
                xfinal = accumarray(rval.step_stats(:,3),rval.step_stats(:,7),[n,1]);
                xnnz = find(xfinal);
                xinds = zeros(n,1);
                xinds(xnnz) = 1:numel(xnnz);
                %%

                X = zeros(numel(xnnz),size(rval.ep_stats,1));
                newconds = [];
                mincond = Inf;
                bsetthresh = zeros(size(rval.ep_stats,1),1);
                mincond_which_ep = 1;
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
                        mincond_which_ep = i;
                    end
                end
                cond = mincond;
                xsolvec = zeros(n,1);
                xsolvec(xnnz) = X(:,mincond_which_ep);
                [~,xperm] = sort(xsolvec,'descend');
                bestcutvals = cutsweep(A,xperm); % find the best conductance set
                [ind, ~] = min( bestcutvals.conductance );
                bset_inds = xperm(1:ind);

                lasteps = Inf;
                for i=size(newconds,1):-1:1
                    curcond=newconds(i,2);
                    curep = 1./newconds(i,1);
                    if curep<lasteps/2;

                    end
                end


cond
    clf
    [px,py] = gplot(A,xy);
    gdraw = @() plot(px,py,'k-','LineWidth',0.4,'Color',0.55*[1,1,1]);

    bigmarkers = 15;
    smallmarkers = 2;
    figsize = [2.75,2.75];

    gdraw(); hold on;

    % set color

    sen_inds = bestset;% find(attendance(:,j));
    scatter(xy(sen_inds,1), xy(sen_inds,2), bigmarkers, [0,0,1], 'filled'); cmap();
    %plot(xy(seed,1),xy(seed,2),'ko','MarkerSize',10);

    axis off;

   print(gcf, ['./senate_figures/senate-layout-bset-tau' num2str(j) '.png'], '-dpng','-r600');
end

