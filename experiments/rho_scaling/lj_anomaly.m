% Explaining the livejournal anomaly
clear; clc;

% CHANGE DIRECTORY
load ../../../ljournal-2008;
assert( max(diag(A)) == 0, 'not hollow');
assert( nnz(A - A') == 0, 'not sym');

addpath ../..;
addpath ../../util; % for cutsweep


%% rho = 0 paths


epsmin = 1e-5;
alpha = 0.99;
rho = 0;
seed = 743479;


fprintf('\nFirst approximate soln paths with rho=%.2f\n', rho);

n = size(A,1);
rval = ppr_path_rho(A,seed,'epsmin',epsmin,'rho',rho,'alpha',alpha);

% find the set of non-zeros and build a local index
xfinal = accumarray(rval.step_stats(:,3),rval.step_stats(:,7),[n,1]);
xnnz = find(xfinal);
xinds = zeros(n,1);
xinds(xnnz) = 1:numel(xnnz);

fprintf('Done. Number of epsilon values explored in our approximate soln paths, %d\n', length(rval.ep_stats) );
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
bestcutvals = cutsweep(A,xperm(1:length(xnnz)); % find the best conductance set


%%
bset_line_width = 1.5;
cmappart = hot(9);
cmap = @() colormap(flipud(cmappart(2:6,:)));

clf; hold on;
xvalind = find(rval.ep_stats(:,1) < epsmin,1,'first');
hs=loglog((1./rval.ep_stats(:,1)), X','color','k');
plot(1./rval.ep_stats(:,1),bsetthresh,'k','LineWidth',bset_line_width);
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

xl = xlim;
xl = [xl(1), 1/epsmin];
yl = ylim;
line([xl(1),1./epsmin],[(1-rho)./xl(1),(1-rho)*epsmin],'Color','k','LineWidth',1);
title('LiveJournal, \rho = 0');
ylabel('Degree-normalized PageRank');
box off;

lasteps = Inf;
for i=size(newconds,1):-1:1
    curcond=newconds(i,2);
    curep = 1./newconds(i,1);
    if curep<lasteps/2;
        line([1./newconds(i,1) 1./newconds(i,1)],[yl(1) (1-rho)*newconds(i,1)],'LineWidth',0.5,'Color','b');
        if curep > 10^4
            ht=text(curep,newconds(i,1), sprintf('\\phi = %.3f',curcond),...
                'Rotation',90,'VerticalAlignment','top','HorizontalAlignment','left',...
                'FontSize',8);
        else
            ht=text(curep,(yl(1)*newconds(i,1))^(1/2), sprintf('\\phi = %.3f',curcond),...
                'Rotation',90,'VerticalAlignment','top','HorizontalAlignment','center',...
                'FontSize',8);
        end
        lasteps = curep;
    end
end

xlim(xl);
set(gca,'XTick',[10,100,1000,10000,100000]);
set(gca,'XTickLabel','');
set_figure_size([3.5,2.5]);
print(gcf,['./images/lj-seed' num2str(seed) '-rho0.png'],'-dpng','-r600');

%% Plot conductance for rho = 0
fprintf('Now, plot conductances. \n');

clf;
plot( 1./rval.ep_stats(:,1), rval.ep_stats(:,2) )

set(gca,'XScale','log');
set(gca,'YScale','log');

xlabel('1/\epsilon');
ylabel('Best \phi');
box off;

xlim(xl);
ylim([ min( rval.ep_stats(:,2)*0.9 ), 1] );
yl_cond = ylim;
set(gca,'XTick',[10,100,1000,10000,100000]);
set_figure_size([3.5,1.5]);

yl = ylim;
line([1e-5,1e-5],[yl(1), yl(2)],'Color','k','LineWidth',1);

print(gcf,['./images/lj-seed' num2str(seed) '-rho0-cond.png'],'-dpng','-r600');

fprintf('\n lj rho=0 done.\n');


%% rho = 0.9 paths

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


%%

cmappart = hot(9);
cmap = @() colormap(flipud(cmappart(2:6,:)));

clf; hold on;
xvalind = find(rval.ep_stats(:,1) < epsmin,1,'first');
hs=loglog((1./rval.ep_stats(:,1)), X','color','k');
plot(1./rval.ep_stats(:,1),bsetthresh,'k','LineWidth',bset_line_width);
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

xl = xlim;
xl = [xl(1), 1/epsmin];
yl = ylim;
line([xl(1),1./epsmin],[(1-rho)./xl(1),(1-rho)*epsmin],'Color','k','LineWidth',1);
title('LiveJournal, \rho = 0.9');
ylabel('Degree-normalized PageRank');
box off;

lasteps = Inf;
for i=size(newconds,1):-1:1
    curcond=newconds(i,2);
    curep = 1./newconds(i,1);
    if curep<lasteps/2;
        line([1./newconds(i,1) 1./newconds(i,1)],[yl(1) (1-rho)*newconds(i,1)],'LineWidth',0.5,'Color','b');
        if curep > 10^4
            ht=text(curep,newconds(i,1), sprintf('\\phi = %.3f',curcond),...
                'Rotation',90,'VerticalAlignment','top','HorizontalAlignment','left',...
                'FontSize',8);
        else
            ht=text(curep,(yl(1)*newconds(i,1))^(1/2), sprintf('\\phi = %.3f',curcond),...
                'Rotation',90,'VerticalAlignment','top','HorizontalAlignment','center',...
                'FontSize',8);
        end
        lasteps = curep;
    end
end

xlim(xl);
set(gca,'XTick',[10,100,1000,10000,100000]);
set(gca,'XTickLabel','');
set_figure_size([3.5,2.5]);
print(gcf,['./images/lj-seed' num2str(seed) '-rho9.png'],'-dpng','-r600');

%% Plot conductance for rho = 0
fprintf('Now, plot conductances. \n');

clf;
plot( 1./rval.ep_stats(:,1), rval.ep_stats(:,2) )

set(gca,'XScale','log');
set(gca,'YScale','log');

xlabel('1/\epsilon');
ylabel('Best \phi');
box off;

xlim(xl);
% ylim([ min( rval.ep_stats(:,2)*0.9 ), 1] );
ylim(yl_cond);
set(gca,'XTick',[10,100,1000,10000,100000]);
set_figure_size([3.5,1.5]);
yl = ylim;
line([1e-5,1e-5],[yl(1), yl(2)],'Color','k','LineWidth',1);

print(gcf,['./images/lj-seed' num2str(seed) '-rho9-cond.png'],'-dpng','-r600');

fprintf('\n lj rho=0.9 done.\n');

%% rho = 0.9 paths

epsmin = 5e-6;
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


%%

cmappart = hot(9);
cmap = @() colormap(flipud(cmappart(2:6,:)));

clf; hold on;
xvalind = find(rval.ep_stats(:,1) < epsmin,1,'first');
hs=loglog((1./rval.ep_stats(:,1)), X','color','k');
plot(1./rval.ep_stats(:,1),bsetthresh,'k','LineWidth',bset_line_width);
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

xl = xlim;
xl = [xl(1), 1/epsmin];
yl = ylim;
line([xl(1),1./epsmin],[(1-rho)./xl(1),(1-rho)*epsmin],'Color','k','LineWidth',1);
title('LiveJournal, \rho=0.9, small \epsilon');
ylabel('Degree-normalized PageRank');
box off;

lasteps = Inf;
for i=size(newconds,1):-1:1
    curcond=newconds(i,2);
    curep = 1./newconds(i,1);
    if curep<lasteps/2;
        line([1./newconds(i,1) 1./newconds(i,1)],[yl(1) (1-rho)*newconds(i,1)],'LineWidth',0.5,'Color','b');
        if curep > 10^4
            ht=text(curep,newconds(i,1), sprintf('\\phi = %.3f',curcond),...
                'Rotation',90,'VerticalAlignment','top','HorizontalAlignment','left',...
                'FontSize',8);
        else
            ht=text(curep,(yl(1)*newconds(i,1))^(1/2), sprintf('\\phi = %.3f',curcond),...
                'Rotation',90,'VerticalAlignment','top','HorizontalAlignment','center',...
                'FontSize',8);
        end
        lasteps = curep;
    end
end

xlim(xl);
set(gca,'XTick',[10,100,1000,10000,100000]);
set(gca,'XTickLabel','');
set_figure_size([3.5,2.5]);
print(gcf,['./images/lj-seed' num2str(seed) '-rho9fulleps.png'],'-dpng','-r600');

%% Plot conductance for rho = 0.9, full eps
fprintf('Now, plot conductances. \n');

clf;
plot( 1./rval.ep_stats(:,1), rval.ep_stats(:,2) )

set(gca,'XScale','log');
set(gca,'YScale','log');

xlabel('1/\epsilon');
ylabel('Best \phi');
box off;

xlim(xl);
% ylim([ min( rval.ep_stats(:,2)*0.9 ), 1] );
ylim(yl_cond);
set(gca,'XTick',[10,100,1000,10000,100000]);
set_figure_size([3.5,1.5]);
yl = ylim;
line([1e-5,1e-5],[yl(1), yl(2)],'Color','k','LineWidth',1);

print(gcf,['./images/lj-seed' num2str(seed) '-rho9fulleps-cond.png'],'-dpng','-r600');

fprintf('\n lj rho=0.9, full eps done.\n');