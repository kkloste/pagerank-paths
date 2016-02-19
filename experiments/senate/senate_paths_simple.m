function [newconds] = senate_paths_simple(A, seed, attendance, varargin)
% [newconds] = senate_paths_simple(A, seed, attendance, varargin)
% Optional paramters: 'epsmin', 'alpha', 'rho'
%
% "simple" means that this omits the vertical blue lines that indicate sets
% of best conductance.

parser = inputParser;
parser.addOptional('epsmin',1e-4);
parser.addOptional('alpha',0.99);
parser.addOptional('rho',0);
parser.parse(varargin{:});
epsmin = parser.Results.epsmin;
alpha = parser.Results.alpha;
rho = parser.Results.rho;

%% Compute pagerank
tic;
rval = ppr_path_rho(A,seed,'epsmin',epsmin,'alpha',alpha,'rho',rho);
fprintf('\nPageRank done in time=%f',toc);
%% Make path plot
% find the set of non-zeros in final vector and build a local index
n = size(A,1);
xfinal = accumarray(rval.step_stats(:,3),rval.step_stats(:,7),[n,1]);
xnnz = find(xfinal);
xinds = zeros(n,1);
xinds(xnnz) = 1:numel(xnnz);
numel(xnnz);

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

% Colormap for paths
clf;
cmappart = hot(9);
cmap = @() colormap(flipud(cmappart(2:6,:)));

hold on;
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
yl = ylim;
line([xl(1),1./epsmin],[(1-rho)./xl(1),(1-rho)*epsmin],'Color','k','LineWidth',1);

taus = [ 3*1e-4, 1e-4, 3*1e-5];
for j =1:length(taus),
    tau = taus(j);
    line([1/tau,1/tau],[yl(1),(1-rho)*tau],'Color','b','LineWidth',0.7);
end

xlabel('1/\epsilon');
ylabel('Degree normalized PageRank');
box off;

xlim([10,1/epsmin]);
set(gca,'XTick',[10,100,1000,10000,100000]);
set_figure_size([3.5,3]);
dummy = num2str(rho);
if length(dummy)>1, dummy = dummy(3:end);
else dummy = '0';
end
%title( sprintf('US-Senate paths, seed %d', seed ));
print(gcf, ['./senate_figures/senate-s', num2str(seed), '-rho', dummy,'.png'], '-dpng','-r600');

%% Now label

% cmap();
% % create labelling color map
% colors =[ 1 0 0; ...    %red
%         1 0.5 0; ...%
%         0.25 1 0.75; ...% 
%           0 0 1]; % blue
% ncolors = size(colors,1);
% [ sen_indices, ~, node_ids ] = find(attendance);
% 
% for i=1:numel(hs)
%     % find terms of the senator corresponding to this node
%     node_id = xnnz(i);
%     temp_ind = find( node_ids == node_id );
%     sen_id = sen_indices(temp_ind);
%     which_terms = find( attendance(sen_id, :) );
% 
%     num_overlap = length( intersect([seed_terms-1:seed_terms+1], which_terms) );
%     seed_overlap = length( intersect(seed_terms, which_terms) );
%     cind = 4; % no overlap -> blue
%     if num_overlap >= 1, cind = 3; end
%     if seed_overlap >= 1, cind = 1; end
%     
%     set(hs(i),'Color',colors(cind,:),'LineWidth',0.3);
% end    

