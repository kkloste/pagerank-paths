function [newconds] = senate_paths(A, seed, attendance, varargin)
% [newconds] = senate_paths(A, seed, attendance, varargin)
% Optional paramters: 'epsmin', 'alpha', 'rho'
%

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
fprintf('\nPageRank paths done in time=%f\n',toc);
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
    cval = log10(X(i,end));
    cind = round(((cval-(crange(1)))/(crange(2)-crange(1)))*(ncolors-1)+1);
    cind = max(cind,1); % cap to range
    cind = min(cind,ncolors);
    set(hs(i),'Color',colors(cind,:),'LineWidth',0.3);
end

xl = xlim;
yl = ylim;
line([xl(1),1./epsmin],[(1-rho)./xl(1),(1-rho)*epsmin],'Color','k','LineWidth',1);

% taus = [ 3*1e-4, 1e-4, 3*1e-5];
% for j =1:length(taus),
%     tau = taus(j);
%     line([1/tau,1/tau],[yl(1),(1-rho)*tau],'Color','b','LineWidth',0.7);
% end

xlabel('1/\epsilon');
ylabel('Degree-normalized PageRank');
box off;

% for i=1:size(newconds,1)
%     line([1./newconds(i,1) 1./newconds(i,1)],[yl(1) (1-rho)*newconds(i,1)],'LineWidth',0.5,'Color','b');
% end
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

xlim([xl(1),1/epsmin]);
set(gca,'XTick',[10,100,1000,10000,100000]);
set_figure_size([3.5,3]);
dummy = num2str(rho);
if length(dummy)>1, dummy = dummy(3:end);
else dummy = '0';
end
%title( sprintf('US-Senate paths, seed %d', seed ));
print(gcf, ['./senate_figures/senate-s', num2str(seed), '-rho', dummy,'.png'], '-dpng','-r600');

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
set(gca,'XTick',[10,100,1000,10000,100000]);
set_figure_size([3.5,2]);

print(gcf,['./senate_figures/senate-s', num2str(seed), '-cond-rho', dummy,'.png'],'-dpng','-r600');

end % END FUNCTION