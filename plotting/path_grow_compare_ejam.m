clear; clc;

load_directory = '../experiments/timing_experiments/results/';
%suffix_rho = '-rho0.9'; % other option is '-rho0'
suffix_rho = '-rho0'; % other option is '-rho0'

files = { 'path-time-itdk0304', ...
'path-time-dblp',...
'path-time-youtube', ...
'path-time-fb-one', ...
'path-time-fbA', ...
'path-time-ljournal-2008', ...
'path-time-hollywood-2009',...
'path-time-twitter-2010',...
'path-time-friendster'
};

files2 = { 'pathgrow-many-itdk0304', ...
'pathgrow-many-dblp',...
'pathgrow-many-youtube', ...
'pathgrow-many-fb-one', ...
'pathgrow-many-fbA', ...
'pathgrow-many-ljournal-2008', ...
'pathgrow-many-hollywood-2009',...
'pathgrow-many-twitter-2010',...
'pathgrow-many-friendster'
};

for j=1:length(files),
     files{j} = [char(files{j}) char(suffix_rho) '.mat'];
     files2{j} = [files2{j} suffix_rho '.mat'];
     disp([ files{j} '  '  files2{j} ])
end
fprintf('\n\n');
%


filenames = {
'itdk0304', 'dblp', 'youtube', 'fb-one', 'fbA', 'ljournal', ...
'hollywood', 'twitter', 'friendster'
};

% Y = prctile(X,p) returns percentiles of the values in X for p in [0,100].

conductances_pprpaths1 = [];
conductances_singlediff = [];
conductances_pprpaths2 = [];
conductances_manydiff = [];

fprintf('\n TABLE 2:  comparing ppr_paths to (1) a single diffusion, and (2) 10,000 pprgrows \n');
fprintf('\n filename  & $\\#$ \\epscu  & single & path & multi &  best ratio \\\\ \n');
for j=1:length(files),
    fname = char(files(j));
    load([load_directory fname]);
    flname = char(filenames(j));
    
    % GET MEDIAN  TIMES, NUM_EPS
    pprpath_times = times(:,1);
    pprpath_numeps = num_eps;
    pprpath_conds = conds(:,1);
    prctl = 50;
    prctlB = 1;

    conductances_pprpaths1 = [conductances_pprpaths1; pprpath_conds];
    
    singlediff_times = times(:,2);
    singlediff_conds = conds(:,2);
    
    conductances_singlediff = [conductances_singlediff; singlediff_conds];

    % get many-diffusion information
    fname = char(files2(j));
    load(['../results/' fname]);

    manydiff_times = times(:,2);
    manydiff_conds = conds(:,2);
    
    conductances_pprpaths2 = [ conductances_pprpaths2; conds(:,1) ];
    conductances_manydiff = [ conductances_manydiff; manydiff_conds ];
    
fprintf('%10.10s  & %2.0f  & %6.2f & %6.2f & %6.2f  & %6.2f & %6.2f & %6.2f & %.1f & %.1f & %.1f & %6.2f \\\\ \n', ... 
        flname, prctile(pprpath_numeps,prctl), ...
        prctile(singlediff_times,25), prctile(singlediff_times,prctl), prctile(singlediff_times,75), ...
        prctile(pprpath_times,25), prctile(pprpath_times,prctl), prctile(pprpath_times,75), ...
        prctile(manydiff_times,25), prctile(manydiff_times,prctl), prctile(manydiff_times,75), ...
         1./(prctile( pprpath_conds./singlediff_conds,prctlB)) );    
end




diff = - (conductances_pprpaths1-conductances_singlediff);
avediff = sum(diff)/length(diff);
top95 = prctile(diff,95);
top90 = prctile(diff,90);
top80 = prctile(diff,80);
top70 = prctile(diff,70);
top60 = prctile(diff,60);

maxdiff = max(diff);
mindiff = min(diff);
%prctile(condtot2,prctl)

fprintf('\n\n  given percentiles of conductance difference, e.g. abs(path-grow) \n');
fprintf('  top95  & top90  & top80 & top70 & top60 & avediff  &  maxdiff & mindiff \\\\ \n');
fprintf(' %6.3f  & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f \\\\ \n', top95  , top90  , top80 , top70 , top60 , avediff  ,  maxdiff, mindiff);

diff = -(conductances_pprpaths2-conductances_manydiff);
avediff = sum(diff)/length(diff);
top95 = prctile(diff,95);
top90 = prctile(diff,90);
top80 = prctile(diff,80);
top70 = prctile(diff,70);
top60 = prctile(diff,60);

maxdiff = max(diff);
mindiff = min(diff);

fprintf('\n\n  Now compare path vs many-diff \n');
fprintf('  top95  & top90  & top80 & top70 & top60 & avediff  &  maxdiff & mindiff \\\\ \n');
fprintf(' %6.3f  & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f \\\\ \n', top95  , top90  , top80 , top70 , top60 , avediff  ,  maxdiff, mindiff);
