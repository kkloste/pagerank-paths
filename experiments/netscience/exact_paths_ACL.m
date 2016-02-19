% Call from [project]/experiments/netscience/
clear; clc;
addpath ../.. % for ppr_path_rho.m
addpath ../../util % for cutsweep info
% Get dataset
load ../../data/netscience-cc;

A = A |A';
A = A - diag(diag(A));
n = size(A,1);

alpha = 0.99;
epsmin=1e-5;
rho = 0.9;
seed = 211; % hand-picked for a good image
seeds = seed;

fprintf('\nFirst approximate soln paths with rho=%.2f\n', rho);

rval = ppr_path_rho(A,seeds,'epsmin',epsmin,'alpha',alpha,'rho',rho);

% find the set of non-zeros in final vector and build a local index
xfinal = accumarray(rval.step_stats(:,3),rval.step_stats(:,7),[n,1]);
xnnz = find(xfinal);
xinds = zeros(n,1);
xinds(xnnz) = 1:numel(xnnz);

fprintf('Done. Number of epsilon values explored in our approximate soln paths, %d\n', length(rval.ep_stats) );

fprintf('Next, extract subset of those epsilon values, such that consecutive values are at least 1e-5 apart.\n');
% find subset of eps values, taking values that are at least 1e-5 apart.
eps_inds = [];
old_eps = 1;
diff = zeros(size(rval.ep_stats,1), 1);
for j=1:size(rval.ep_stats,1),
    new_eps = rval.ep_stats(j,1);
    diff(j) = abs(old_eps-new_eps);
    if abs(old_eps-new_eps)/new_eps>=3*1e-4,
        eps_inds = [eps_inds; j];
    end
    old_eps = new_eps;
end

%%

eps_subset = rval.ep_stats(eps_inds);
X = zeros(numel(xnnz),size(eps_inds,1));

for i=1:size(eps_subset,1)
    ep = eps_subset(i);% rval.ep_stats(i,1);
    step = rval.ep_stats( eps_inds(i) ,6)+1;
    xvec=accumarray(rval.step_stats(1:step,3),rval.step_stats(1:step,7),[n,1]);
    X(:,i) = xvec(xnnz);
    [~,xperm] = sort(X(:,i),'descend');
end

x_paths = zeros(n,1);
x_paths(xnnz) = X(:,end);

fprintf('Done. Number of epsilons in subset, %d\n', length(eps_inds) );

fprintf('Next, plot results for our approximate paths. \n');

%% Colormap for paths
clf;
cmappart = hot(9);
cmap = @() colormap(flipud(cmappart(2:6,:)));

hold on;
hs=loglog((1./eps_subset), X','color','k');

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

xl = xlim;
yl = ylim;
xl = [xl(1), 1e5/3];
yl = [1e-4, yl(2)];

%title(sprintf('Our \rho=%.2f approximate paths', rho) );
title('Our \rho = 0.9 approximate paths' );
% xlabel('1/\epsilon');
ylabel('Degree-normalized PageRank');
box off;

ylim(yl)
xlim(xl);
% set(gca,'XTick',[10,100,1000,10000,100000]);
% set_figure_size([3.5,3]);
set(gca,'XTickLabel','');
set_figure_size([3.5,2.5]);

print(gcf, ['../images/netscience-our-path', num2str(seed), '.png'],'-dpng','-r600');  

%% Plot conductance for our paths
fprintf('Now, plot conductances for our approximate paths. \n');

clf;
plot( 1./eps_subset, rval.ep_stats(eps_inds,2) )

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

print(gcf, ['../images/netscience-our-cond', num2str(seed), '.png'],'-dpng','-r600');


%%
% Now do ACL-paths
fprintf('\nNext, compute exact solution paths (rho=1.0) on the subset of epsilon values. \n');

% Set up ACL inputs
d = full(sum(A,2));
vS = zeros(n,1);
vS(seeds) = 1;
 
taus = eps_subset*(1-alpha);

num_taus = numel(taus);
v1 = d.*vS; 
volS = full(sum(v1));
v1 = v1 / volS; v1 = sparse(v1);

X2 = zeros( length(d), num_taus) ;
infdists = zeros(num_taus,1);
minconds = zeros(num_taus,1);
thresh = sqrt(eps(1));
for taui = 1:num_taus
    tau = taus(taui);

    z = pprfast(A,d,alpha,v1,tau,'tol',sqrt(eps(1)),'skip_push',true);
    x = (z./d)./(1-alpha);
    xnnzs( taui ) = nnz(x);

    X2(:,taui) = x;
    infdists(taui) = norm( X(:,taui) - x,'inf' );
    % check how far away exact paths are from our paths
    
    % now get conductance information
    [~,permut] = sort( x, 'descend' );
    cutvals = fastcutsweep(A, permut);
    minconds(taui) = min(cutvals.conductance);
end

X = X2;
x_acl = X(:,end);

% fprintf('\nCheck path sums: %f  %f\n',  sum(d.*x_acl)*(1-alpha), sum(d.*x_paths)*(1-alpha) );
fprintf('Check inf-norm distance between paths at end: %f\n', norm( x_acl - x_paths, 'inf' )*(1-alpha) );
fprintf('Check inf-norm distance between paths, global: %f\n', max(infdists)*(1-alpha) );

%% ACL PATHS PLOT

fprintf('Finally, plot results from exact paths.\n');

xfinal = sum(X,2);
xnnz = find(xfinal);
xinds = zeros(n,1);
xinds(xnnz) = 1:numel(xnnz);
X = X( xnnz, :);

clf; hold on;
hs=loglog(((1-alpha)./taus), X','color','k');
set(gca,'XScale','log');
set(gca,'YScale','log');
ncolors = size(colors,1);
for i=1:numel(hs)
    cval = log10(X(i,end));
    cind = round(((cval-(crange(1)))/(crange(2)-crange(1)))*(ncolors-1)+1);
    cind = max(cind,1); % cap to range
    cind = min(cind,ncolors);
    set(hs(i),'Color',colors(cind,:),'LineWidth',0.3);
end

title('Exact paths');
% xlabel('1/\epsilon');
ylabel('Degree-normalized PageRank');
box off;

xlim(xl);
ylim(yl);
% set(gca,'XTick',[10,100,1000,10000,100000]);
% set_figure_size([3.5,3]);
set(gca,'XTickLabel','');
set_figure_size([3.5,2.5]);

print(gcf,['../images/netscience-ACL-path', num2str(seed), '.png'],'-dpng','-r600');

%% Plot conductance for our paths
fprintf('... and plot conductances for exact paths. \n');

clf;
plot( ((1-alpha)./taus) , minconds );

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

print(gcf, ['../images/netscience-exact-cond', num2str(seed), '.png'],'-dpng','-r600');  


