clear
clc


files = { 
'path-time-fb-one-cc', ...
'path-time-ljournal-prepped', ...
'path-time-hollywood-2009-cc',...
'path-time-friendster'};


files2 = { 
'pathgrow-many-fb-one-cc', ...
'pathgrow-many-ljournal-prepped', ...
'pathgrow-many-hollywood-2009-cc',...
'pathgrow-many-friendster'};

filenames = {
'fb-one', 'ljournal', ...
'hollywood', 'friendster'
};

% Y = prctile(X,p) returns percentiles of the values in X for p in [0,100].

condtot1 = [];
condtot2 = [];

    fprintf('filename  & $\\#$ \\epscu  & Single & Path & Appx. Path &  Single cond & 95 best Path/single cond \\\\ \n');
for j=1:length(files),
    fname = char(files(j));
    load(['../results/' fname]);
    flname = char(filenames(j));

    
    
    % GET MEDIAN  TIMES, NUM_EPS
    dat_pt = times(:,1);
    dat_eps = num_eps;
    dat_pc = conds(:,1);
    prctl = 50;
    prctlB = 1;

    condtot1 = [condtot1; dat_pc];
    
    dat_gt = times(:,2);
    dat_gc = conds(:,2);
    
    condtot2 = [condtot2; dat_gc];

    % get 
    fname = char(files2(j));
    load(['../results/' fname]);

    dat_gtm = times(:,2);
    
    %         name  times   num_eps
%    fprintf('%10.10s  & %2.0f  & %6.2f & %6.2f & %4f & %6.3f  \\\\ \n', ... 
%        flname, prctile(dat2,prctl), prctile(dat1,prctl), prctile(dat21,prctl), prctile(dat3,prctl), prctile( abs(dat23) ,prctl) );
%    fprintf('      ave  & %2.0f  & %6.3f & %6.3f & %6.3f & %6.3f  \\\\ \n', ... 
%       sum(dat2)/100, sum(dat1)/100, sum(dat21)/100, sum(dat3)/100, sum(abs(dat23))/100 );
    fprintf('%10.10s  & %2.0f  & %6.2f & %6.2f & %.f & %6.2f & %6.2f \\\\ \n', ... 
        flname, prctile(dat_eps,prctl), prctile(dat_gt,prctl), prctile(dat_pt./dat_gt,prctl), prctile(dat_gtm./dat_gt,prctl), ...
        prctile(dat_gc,prctl), prctile( dat_pc./dat_gc ,prctlB) );
    
end



diff = (condtot1-condtot2);
avediff = sum(abs(diff))/length(diff);
top95 = prctile(abs(diff),95);
top90 = prctile(abs(diff),90);
top80 = prctile(abs(diff),80);
top70 = prctile(abs(diff),70);
top60 = prctile(abs(diff),60);

maxdiff = max(abs(diff));
%prctile(condtot2,prctl)

fprintf('\n\n  given percentiles of conductance difference, e.g. abs(path-grow) \n');
fprintf('top95  & top90  & top80 & top70 & top60 & avediff  &  maxdiff \\\\ \n');
fprintf(' %6.3f  & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f \\\\ \n', top95  , top90  , top80 , top70 , top60 , avediff  ,  maxdiff);
