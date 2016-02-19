% Plot PageRank solutions paths and ground-truth labels
% for USPS digits dataset
clear; clc;
fname = 'usps_3nn';
load(['../../data/' fname '.mat']); % must change directory to point to location of dataset
A = G; clear G;

% Ensure adjacency matrix is symmetrized and binary
% (can skip if data is pre-processed)
A = A|A';
A = A - diag(diag(A));
n = size(A,1);

addpath ../..; %for ppr_path tools
addpath ../../util; % for set_figure_size.m


%% Choose parameters and inputs

num_seeds = 4;
digit = 5;
node_labels = labels; clear labels;
digit_label = find( node_labels == digit );

alpha1 = 0.99;
epsmin=1e-4;
rho = 0.9;

% So we can label image according to value of rho used
dummyrho = num2str(rho);
if length(dummyrho)>1, dummyrho = dummyrho(3:end);
else dummyrho = '0'; end

%% GET SEEDS
% get num_seeds random nodes of the same label (i.e. digit)

% % Used to find interesting looking seed sets:
% % temp_ind = randi( length(digit_label), 2*num_seeds, 1);
% % temp_ind = unique(temp_ind);
% % temp_ind = temp_ind( 1:num_seeds );
% % temp_ind = [146   335   559   579];
% % temp_ind = 28    40    83   236

temp_ind = [30   160   326   374];
seeds = digit_label( temp_ind );
tic;
rval = ppr_path_rho(A,seeds,'epsmin',epsmin, 'degweights', true, 'rho', rho);
min(rval.ep_stats(:,2));
toc

%%
% find the set of non-zeros and build a local index
xfinal = accumarray(rval.step_stats(:,3),rval.step_stats(:,7),[n,1]);
xnnz = find(xfinal);
xinds = zeros(n,1);
xinds(xnnz) = 1:numel(xnnz);

% Find best conductance info
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

%% Draw path plots
%%% individually first, then as subplots
clf; hold on;
hs=loglog((1./rval.ep_stats(:,1)), X','color','k');
plot(1./rval.ep_stats(:,1),bsetthresh,'k','LineWidth',1);
set(gca,'XScale','log');
set(gca,'YScale','log');
crange = [-3,0];
cmappart = hot(9);
cmap = @() colormap(flipud(cmappart(2:6,:)));
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

% Draw in diagonal line
xl = xlim;
yl = ylim;
line([xl(1),1./epsmin],[(1-rho)./xl(1),(1-rho)*epsmin],'Color','k','LineWidth',1);
%xlabel('1/\epsilon');
ylabel('Degree-normalized PageRank');
title('USPS-digits');
box off;

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
ylim([1e-5, 3*1e-1]);
%set(gca,'XTick',[10,100,1000,10000,100000]);
set(gca,'XTickLabel','');
set_figure_size([3.5,2.5]);
print(gcf,strcat( './figures/', fname, '-phi-rho', dummyrho,  '-', num2str(digit),  '.png'),'-dpng','-r600');

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

print(gcf,strcat( './figures/', fname, '-cond-rho', dummyrho,  '-', num2str(digit),  '.png'),'-dpng','-r600');


%% Draw path plots with labels
clf; hold on;
hs=loglog((1./rval.ep_stats(:,1)), X','color','k');
plot(1./rval.ep_stats(:,1),bsetthresh,'k','LineWidth',1);
set(gca,'XScale','log');
set(gca,'YScale','log');
crange = [-3,0];
cmap();

% create labelling color map
colors=colormap(winter(10));
colors = colors( [1,6,2,7,3,8,4,9,5,10], :);
ncolors = size(colors,1);
colors(digit,:) = [ 1,0,0 ]; % assign color of chosen label to be hot red

for i=1:numel(hs)
    i_digit_label = node_labels( xnnz(i) );
    if i_digit_label ~= digit, i_digit_label = 8; end % forces all outliers to be same color
    set(hs(i),'Color',colors(i_digit_label,:),'LineWidth',0.3);
end

% Draw in diagonal line
xl = xlim;
yl = ylim;
line([xl(1),1./epsmin],[(1-rho)./xl(1),(1-rho)*epsmin],'Color','k','LineWidth',1);
%xlabel('1/\epsilon');
ylabel('Degree-normalized PageRank');
title('USPS-digits, labeled');
box off;

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
ylim([1e-5, 3*1e-1]);
%set(gca,'XTick',[10,100,1000,10000,100000]);
set(gca,'XTickLabel','');
set_figure_size([3.5,2.5]);

fprintf('\n paths done.\n');
print(gcf,strcat( './figures/', fname, '-phi-rho', dummyrho,  '-', num2str(digit), '-labels.png'),'-dpng','-r600');

fprintf('\n prpath_%s done.\n', fname);