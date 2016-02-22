% Plot PageRank solutions paths for facebook dataset
clear; clc; 
load ../../data/fbA; % must change directory to point to location of dataset

% Ensure adjacency matrix is symmetrized and binary
% (can skip if data is pre-processed)
A = A|A';
A = A - diag(diag(A));
n = size(A,1);
addpath ../..; %for ppr_path tools
addpath ../../util; % for set_figure_size

seed = 62;
epsmin=1e-5;
rho = 0.9;

tic;
rval = ppr_path_rho(A,seed,'epsmin',epsmin,'rho',rho);
timet = toc;
fprintf('\nSoln paths (with rho = %.2f) computed in %f \n', rho, timet );
fprintf('\nMin conductance found: %f,  num eps=%d \n', min(rval.ep_stats(:,2)), size(rval.ep_stats,1) );


%%
% find the set of non-zeros and build a local index

fprintf('\nPaths computed; number of eps values: %d', size( rval.ep_stats, 1 ) );

xfinal = accumarray(rval.step_stats(:,3),rval.step_stats(:,7),[n,1]);
xnnz = find(xfinal);
xinds = zeros(n,1);
xinds(xnnz) = 1:numel(xnnz);

%% determine subset of eps values, to speed things up
% comment out this section of code to compute exact paths.
% For faster termination, uncomment this section
% old_eps = 1;
% diff = zeros(size(rval.ep_stats,1), 1);
% for j=1:size(rval.ep_stats,1),
%     new_eps = rval.ep_stats(j,1);
%     diff(j) = abs(old_eps-new_eps)/new_eps;
%     old_eps = new_eps;
% end
% smallest_eps_gap = 50*epsmin ; % for faster testing
% % smallest_eps_gap = epsmin ; % for good approximate path plot
% eps_inds = find( diff >= smallest_eps_gap ); 
% 
% rval.ep_stats = rval.ep_stats(eps_inds,:);

% fprintf('\nBegin plotting subset; number of eps values: %d, max used: %d', size( rval.ep_stats, 1 ), max(eps_inds)  );

%% Find best conductance info
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
%%

cmappart = hot(9);
cmap = @() colormap(flipud(cmappart(2:6,:)));

clf; hold on;
hs=loglog((1./rval.ep_stats(:,1)), X','color','k');
plot(1./rval.ep_stats(:,1),bsetthresh,'k','LineWidth',2);
set(gca,'XScale','log');
set(gca,'YScale','log');
crange = [-3,0];
cmap();
colors=colormap;
ncolors = size(colors,1);
for i=1:numel(hs)
    cval = log10(X(i,end));
    cind = round(((cval-(crange(1)))/(crange(2)-crange(1)))*(ncolors-1)+1);
    cind = max(cind,1); % cap to range
    cind = min(cind,ncolors);
    set(hs(i),'Color',colors(cind,:),'LineWidth',0.3);
end

% Construct diagonal line
xl = xlim;
yl = ylim;
line([xl(1),1./epsmin],[(1-rho)./xl(1),(1-rho)*epsmin],'Color','k','LineWidth',1);
title('Facebook');
ylabel('Degree-normalized PageRank');
box off;

% Insert vertical blue lines marking best conductance sets
lasteps = Inf;
for i=size(newconds,1):-1:1
    curcond=newconds(i,2);
    curep = 1./newconds(i,1);
    if curep<lasteps/2;
        line([1./newconds(i,1) 1./newconds(i,1)],[yl(1) (1-rho)*newconds(i,1)],'LineWidth',0.5,'Color','b');
        if curep > 10^4
            ht=text(curep,yl(1), sprintf('\\phi = %.3f',curcond),...
                'Rotation',90,'VerticalAlignment','top','HorizontalAlignment','left',...
                'FontSize',9);
        else
            ht=text(curep,yl(1), sprintf('\\phi = %.3f',curcond),...
                'Rotation',90,'VerticalAlignment','top','HorizontalAlignment','left',...
                'FontSize',9);
        end
        lasteps = curep;
    end
end

xlim(xl);
set(gca,'XTick',[10,100,1000,10000,100000]);
set(gca,'XTickLabel','');
set_figure_size([3.5,2.5]);
% set_figure_size([3.5,3]);

% Label image according to value of rho used
dummy = num2str(rho);
if length(dummy)>1, dummy = dummy(3:end);
else dummy = '0'; end
print(gcf,['./figures/fbA-path-rho', dummy ,'.png'],'-dpng','-r600');
fprintf('\n paths done.\n');

%% Plot conductance for our paths
fprintf('Now, plot conductances. \n');

clf;
plot( 1./rval.ep_stats(:,1), rval.ep_stats(:,2) )

set(gca,'XScale','log');
set(gca,'YScale','log');

%title( );
xlabel('1/\epsilon');
ylabel('Best \phi');
box off;

xlim(xl);
ylim([ min( rval.ep_stats(:,2)*0.9 ), 1] );
set(gca,'XTick',[10,100,1000,10000,100000]);
set_figure_size([3.5,1.5]);

print(gcf,['./figures/fbA-cond-rho', dummy ,'.png'],'-dpng','-r600');

fprintf('\n prpath_fbA done.\n');
