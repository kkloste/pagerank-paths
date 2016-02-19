function [cond, bestcutvals] = ppr_pathplot_rho(A, seed, varargin)
% [cond, bestcutvals] = ppr_pathplot_rho(A, seeds, varargin)
%
% Uses our ppr_paths_rho to plot solution paths for seed node with
% given input parameters in varargin.
%
% rval = ppr_path_rho(A,seeds)
%   Computes the PageRank solution path for alpha = 0.99, epsmin=1e-4,
%   and who = 0. Starting from the nodes "seeds"
%
%   If you wish to use weighted seeds, we allow integer weights by listing
%   a seed multiple times. 
%
% rval = ppr_path_rho(A,seeds,'key',value,'key',value...) allows the following
% optional key and value pairs:
% 
%   'alpha', <value> : the value of alpha to use
%   'epsmin', <value> : the smallest value of eps considered
%   'rho', <value> : the value of rho used
%
% Default values:
% alpha = 0.99;
% epsmin=1e-4;
% rho = 0.7;

parser = inputParser;
parser.addOptional('epsmin',1e-4);
parser.addOptional('alpha',0.99);
parser.addOptional('rho',0.7);
parser.parse(varargin{:});

epsmin = parser.Results.epsmin;
alpha = parser.Results.alpha;
rho = parser.Results.rho;

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
bestcutvals = cutsweep(A,xperm);


%%

cmappart = hot(9);
cmap = @() colormap(flipud(cmappart(2:6,:)));

clf; hold on;
xvalind = find(rval.ep_stats(:,1) < epsmin,1,'first');
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

xl = xlim;
xl = [xl(1), 1/epsmin];
yl = ylim;
line([xl(1),1./epsmin],[(1-rho)./xl(1),(1-rho)*epsmin],'Color','k','LineWidth',1);
xlabel('1/\epsilon');
ylabel('Degree normalized PageRank');
box off;

for i=1:size(newconds,1)
    line([1./newconds(i,1) 1./newconds(i,1)],[yl(1) (1-rho)*newconds(i,1)],'LineWidth',0.5,'Color','b');
end
lasteps = Inf;
for i=size(newconds,1):-1:1
    curcond=newconds(i,2);
    curep = 1./newconds(i,1);
    if curep<lasteps/2;
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
set_figure_size([3.5,3]);
% print(gcf,'./figures/fbA-path-rho.png','-dpng','-r600');

