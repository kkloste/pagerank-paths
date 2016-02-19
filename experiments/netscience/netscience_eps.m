% Call from /experiments/netscience/;
clear; clc; clf;
load ../../data/xy_1
n = size(A,1);
addpath ../..; % for ppr_path.m function

%%
cmappart = hot(9);
cmap = @() colormap(flipud(cmappart(2:6,:)));
[px,py] = gplot(A,xy);
gdraw = @() plot(px,py,'k-','LineWidth',0.4,'Color',0.55*[1,1,1]); 
sdraw = @() plot(xy(S,1),xy(S,2),'ko','MarkerSize',7);
bigmarkers = 15;
smallmarkers = 2;
figsize = [2.75,2.75];
figarea = [0,220,-600,-120];

%%
rho = 0; % also try rho=9.0
epsmin = 1e-5;
alpha = 0.99;
vert = 211; % chosen as an illustrative example

%%
tic;
rval = ppr_path_rho(A,vert,'epsmin',epsmin,'alpha',alpha,'rho',rho);
toc;

tols=[1e-2,1e-3,1e-4];
clf;

for i=1:numel(tols)
    subplot(1,numel(tols),i);
    %x = rval.solution(tols(i));
    % incorporate the solution function here to get step info
    
    ind=find(rval.ep_stats(:,1) < tols(i),1,'first');
    step=rval.ep_stats(ind,6)+1;
    x=accumarray(rval.step_stats(1:step,3),rval.step_stats(1:step,7),[n,1]);

    
    f = x > 0;
    cla;
    gdraw(); hold on; 
    plot(xy((vert),1),xy((vert),2),'ko','MarkerSize',7);
    scatter(xy(f,1),xy(f,2),bigmarkers,log10(x(f)),'filled'); cmap();
    caxis([-3,0]);
    set_figure_size(figsize);
    axis(figarea);
    axis off;
    hold off;
end

%%
% Now save these figures individually

for i=1:numel(tols)
    clf;
    %x = rval.solution(tols(i));
    % incorporate the solution function here to get step info
    
    ind=find(rval.ep_stats(:,1) < tols(i),1,'first');
    step=rval.ep_stats(ind,6)+1;
    x=accumarray(rval.step_stats(1:step,3),rval.step_stats(1:step,7),[n,1]);

    
    f = x > 0;
    cla;
    gdraw(); hold on; 
    plot(xy((vert),1),xy((vert),2),'ko','MarkerSize',7);
    scatter(xy(f,1),xy(f,2),bigmarkers,log10(x(f)),'filled'); cmap();
    caxis([-3,0]);
    set_figure_size(figsize);
    axis(figarea);
    axis off;
    hold off;
    
    print(sprintf('../images/netscience-rho%.2f-%g',rho,log10(tols(i))),'-dpng','-r600','-painters');
end


%% Prep work that helps in picking a vertex
% We considered these as seeds: 151, 211
%
% clf;
% d = sum(A,2);
% v1 = zeros(size(A,1),1); v1(211) = 1; v1 = v1 / sum(v1); v1 = sparse(v1);
% x = (diag(d) - 0.99*A) \ (v1);
% gdraw(); hold on; 
% plot(xy(find(v1),1),xy(find(v1),2),'ko','MarkerSize',7);
% scatter(xy(:,1),xy(:,2),bigmarkers,log10(x),'filled'); cmap();
% set_figure_size(figsize);
% axis(figarea);
% axis off;
% hold off;

