% Generate plot displaying how the runtime of ppr_paths_rho_mex
% scales as rho goes from 0 to 1.
% Other paramter values were constant in these experiments:
% alpha = 0.99
% epsmin = 1e-5

% Call from [project]/plotting/

clear; clc; clf;

load_directory = '../experiments/rho_scaling/results/';
image_directory = '../images/';

files = {
    'rho-scaling-youtube.mat',...    
'rho-scaling-twitter-2010.mat',...
'rho-scaling-ljournal-2008.mat',...
'rho-scaling-hollywood-2009.mat', ... 
'rho-scaling-friendster.mat', ...
'rho-scaling-fbA.mat', ...
    'rho-scaling-dblp.mat'};
filenames = {
'youtube', 'twitter', 'ljournal', 'hollywood', 'friendster', ...
'fbA', 'dblp'
};


% Y = prctile(X,p) returns percentiles of the values in X for p in [0,100].

% RUNTIMES
linehandles = [];
hold all
for j=1:length(files),
    fname = char(files(j));
    load([load_directory, fname]);
    flname = char(filenames(j));
    
    grid_size = size(times,2);
    rhos = [0:1:grid_size-1]./grid_size;

    dat = prctile( times, 50 );
    linehandles(end+1) = plot( rhos, dat, 'LineWidth',1 );
%     prctlA = 25;
%     prctlB = 75;

end
%legend(linehandles, filenames, 'Location', 'Northwest');
xlim([0,1]);
%legend(filenames, 'Location', 'Northwest');
%legend boxoff;
legend boxoff;
xlabel('rho');
ylabel('Runtime (sec)');
set(gca,'yscale','log');

addpath ../util;
set_figure_size([3.5,3]);
print(gcf,strcat(image_directory, 'rho-scaling-times.eps'),'-depsc2');

%%
% CONDS
clf;
hold all

which_prctA = 100;
which_prctB = 90;
linehandles = [];

for j=1:length(files),
    fname = char(files(j));
    load([load_directory, fname]);
    flname = char(filenames(j));
    
    grid_size = size(times,2);
    rhos = [0:1:grid_size-1]./grid_size;

    minconds = min(conds')';
   
    cond_diff = (conds - minconds*ones(1, size(conds,2)));
    dat50 = prctile( cond_diff, 50 );
    datA = prctile( cond_diff, which_prctA );
    datB = prctile( cond_diff, which_prctB );
    %dat = min( conds );
    %dat = max( conds );
    linehandles(end+1) = plot( rhos, datA, 'LineWidth',1 );
    %linehandles(end+1) = plot( rhos, datA, '.-', 'Color', get(linehandles(end),'Color') );
    
end

xlim([0,1]);
legend(linehandles, filenames, 'Location', 'Northwest');
legend boxoff;
xlabel('rho');
ylabel(sprintf('Conductance gap'));
%set(gca,'yscale','log');

addpath ../util;
set_figure_size([3.5,3]);
print(gcf,strcat(image_directory, 'rho-scaling-conductance.eps'),'-depsc2');

%% Number of eps

clf;
hold all
for j=1:length(files),
    fname = char(files(j));
    load([load_directory, fname]);
    flname = char(filenames(j));
    
    grid_size = size(times,2);
    rhos = [0:1:grid_size-1]./grid_size;

    dat = prctile( num_eps, 50 );
    plot( rhos, dat, 'LineWidth',1 );
%     prctlA = 25;
%     prctlB = 75;

end

xlim([0,1]);
%legend(filenames, 'Location', 'Northwest');
%legend boxoff;
xlabel('rho');
ylabel('Number of different \epsilon explored');
set(gca,'yscale','log');

addpath ../util;
set_figure_size([3.5,3]);
%print(gcf,strcat('rho-scaling-numeps.png'),'-dpng','-r600');
print(gcf,strcat(image_directory, 'rho-scaling-numeps.eps'),'-depsc2');
