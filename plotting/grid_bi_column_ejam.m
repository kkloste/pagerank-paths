% This script creates LaTeX formatted tables
% displaying the times and conductances for ppr-grid and pprgrow.
%
% Last altered 11/29/2015

clear; clc; 
cd ~/Dropbox/ppr-all/code/plotting

files = {
    'grid_v_grow_on_itdk0304.mat', ...
    'grid_v_grow_on_dblp.mat', ...
    'grid_v_grow_on_youtube.mat',... 
    'grid_v_grow_on_fb-one.mat', ...
    'grid_v_grow_on_fbA.mat', ...
    'grid_v_grow_on_ljournal-2008.mat',...
    'grid_v_grow_on_hollywood-2009.mat', ... 
    'grid_v_grow_on_twitter-2010.mat',...
    'grid_v_grow_on_friendster.mat', ...
};

filenames = {
'itdk0304', 'dblp', 'youtube', 'fb-one', 'fbA', 'ljournal', 'hollywood', ...
    'twitter', 'friendster'
};

% Y = prctile(X,p) returns percentiles of the values in X for p in [0,100].

fprintf('\n First do ratios of conductances (This is Table 4)  \n \n');
for j=1:length(files),
    fname = char(files(j));
    load(['../results/' fname]);
    flname = char(filenames(j));
    %name & pgrow & pgrid & 64 & 32 \\

    % do conds, 25%
    dat = conds';
    dat1 = dat(:,1);
    dat2 = dat1./dat(:,2);
    dat3 = dat1./dat(:,3);
    dat4 = dat1./dat(:,4);
    
    prctlA = 25;
    prctlB = 75;
    %fprintf('%11.11s  &  %.3f  &  %.3f  &  %.3f  &  %.3f  &  %.3f  &  %.3f  &  %.3f  \\\\ \n', flname,...
    fprintf('%11.11s  &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f  \\\\ \n', flname,...
        prctile(dat1,50) , prctile(dat4,prctlA) , prctile(dat4,prctlB), prctile(dat3,prctlA), ...
        prctile(dat3,prctlB) , prctile(dat2,prctlA) , prctile(dat2,prctlB)  );
end

% for j=1:length(files),
%     fname = char(files(j));
%     load(['../results/' fname]);
%     flname = char(filenames(j));
%     %name & pgrow & pgrid & 64 & 32 \\
% 
%     % do conds, 25%
%     dat = conds';
%     dat1 = dat(:,1);
%     dat2 = dat1./dat(:,2);
%     dat3 = dat1./dat(:,3);
%     dat4 = dat1./dat(:,4);
%     
%     prctlA = 25;
%     prctlB = 75;
%     %fprintf('%11.11s  &  %.3f  &  %.3f  &  %.3f  &  %.3f  &  %.3f  &  %.3f  &  %.3f  \\\\ \n', flname,...
%     fprintf('%11.11s  &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f  \\\\ \n', flname,...
%          prctile(dat4,prctlA) , prctile(dat4,prctlB), prctile(dat3,prctlA), ...
%         prctile(dat3,prctlB) , prctile(dat2,prctlA) , prctile(dat2,prctlB)  );
% end


% fprintf('\n Now do times \n \n');
% for j=1:length(files),
%     fname = char(files(j));
%     load(['../results/' fname]);% load(fname);
%     flname = char(filenames(j));
%     %name & pgrow & pgrid & 64 & 32 \\
% 
%     % do conds, 25%
%     dat = times';
%     dat1 = dat(:,1);
%     dat2 = dat(:,2);%./dat1;
%     dat3 = dat(:,3);%./dat1;
%     dat4 = dat(:,4);%./dat1;
%     
%     prctlA = 25;
%     prctlB = 75;
%     fprintf('%11.11s  &  %.2f &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f  \\\\ \n', flname,...
%         prctile(dat1,prctlA), prctile(dat1,prctlB) , prctile(dat4,prctlA) , prctile(dat4,prctlB),...
%         prctile(dat3,prctlA), prctile(dat3,prctlB) , prctile(dat2,prctlA) , prctile(dat2,prctlB)  );
% end

fprintf('\n Now do ratios of times (This is Table 3) \n \n');
for j=1:length(files),
    fname = char(files(j));
    load(['../results/' fname]);% load(fname);
    flname = char(filenames(j));
    %name & pgrow & pgrid & 64 & 32 \\

    % do conds, 25%
    dat = times';
    dat1 = dat(:,1);
    dat2 = dat(:,2)./dat1;
    dat3 = dat(:,3)./dat1;
    dat4 = dat(:,4)./dat1;
    
    prctlA = 25;
    prctlB = 75;
    fprintf('%11.11s  &  %.2f &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f  \\\\ \n', flname,...
        prctile(dat1,prctlA), prctile(dat1,prctlB) , prctile(dat4,prctlA) , prctile(dat4,prctlB),...
        prctile(dat3,prctlA), prctile(dat3,prctlB) , prctile(dat2,prctlA) , prctile(dat2,prctlB)  );
end

% RUNTIME comparison of our \pgrid with \pgrow.
% % & & &time &ratio &time ratio& &time ratio\\
%  & \multicolumn{2}{>{}l}{time (sec.)} & \multicolumn{2}{>{\hspace{-6pt}}l}{time ratio} & \multicolumn{2}{>{\hspace{-6pt}}l}{time ratio} & \multicolumn{2}{>{\hspace{-6pt}}l}{time ratio} \\
% %
% Data & \multicolumn{2}{>{}l}{\pgrow } & \multicolumn{2}{>{\hspace{-6pt}}l}{\pgrid$N=32$} & \multicolumn{2}{>{\hspace{-6pt}}l}{\pgrid $N = 64$} & \multicolumn{2}{>{\hspace{-6pt}}l}{\pgrid $N = 1256$} \\
% & 25 & 75 &
% 25 & 75 & 
%  25 & 75 &
%  25 & 75 \\
% \midrule
%     itdk0304  &  6.42 &  9.06  &  0.56  &  0.60  &  0.61  &  0.64  &  1.12  &  1.18  \\ 
%        dblp  &  4.80 &  7.83  &  0.53  &  0.63  &  0.59  &  0.69  &  1.24  &  1.44  \\ 
%   fb-one  &  1.46 &  1.96  &  0.33  &  0.38  &  0.44  &  0.51  &  3.65  &  4.33  \\ 
%         fbA  &  0.55 &  0.74  &  0.46  &  0.54  &  0.62  &  0.70  &  5.93  &  6.71  \\ 
%    ljournal  &  0.83 &  1.29  &  0.43  &  0.54  &  0.57  &  0.72  &  4.57  &  5.95  \\ 
%   hollywood  &  0.33 &  1.07  &  0.36  &  0.52  &  0.46  &  0.66  &  3.33  &  5.15  \\ 
%     twitter  &  0.17 &  0.44  &  0.41  &  0.45  &  0.55  &  0.62  &  4.41  &  5.47  \\ 
%  friendster  &  0.31 &  0.44  &  0.39  &  0.45  &  0.52  &  0.60  &  4.09  &  4.58  \\ 
%  \bottomrule
%  \vspace{1pt}
% \end{tabularx}
% Columns 2 and 3 display the 25th and 75th
% percentile runtimes for ppr-grow (in seconds). The other columns display the median over
% the 100 trials of the ratios of the runtimes of ppr-grid (using the indicated parameter
% setting) with the runtime of ppr-grow on the same node.


% Data & \texttt{grow} & \multicolumn{2}{>{\hspace{-6pt}}l}{$N=32$} & \multicolumn{2}{>{\hspace{-6pt}}l}{$N = 64$} & \multicolumn{2}{>{\hspace{-6pt}}l}{$N = 1256$} \\
% &  & 25 & 75 & 
%  {25} & {75} &
%  {25} & {75} \\
% \midrule
%    itdk0304  &  0.06  &  1.00  &  1.00  &  1.00  &  1.01  &  1.00  &  1.02  \\ 
%        dblp  &  0.07  &  1.00  &  1.00  &  1.00  &  1.00  &  1.00  &  1.01  \\ 
%   fb-one-cc  &  0.37  &  1.06  &  1.16  &  1.10  &  1.26  &  1.18  &  1.37  \\ 
%         fbA  &  0.56  &  1.00  &  1.05  &  1.00  &  1.06  &  1.00  &  1.09  \\ 
%    ljournal  &  0.32  &  1.00  &  1.01  &  1.00  &  1.01  &  1.00  &  1.01  \\ 
%   hollywood  &  0.25  &  1.00  &  1.03  &  1.00  &  1.04  &  1.00  &  1.05  \\ 
%     twitter  &  0.81  &  1.00  &  1.00  &  1.00  &  1.00  &  1.00  &  1.00  \\ 
%  friendster  &  0.85  &  1.00  &  1.00  &  1.00  &  1.01  &  1.00  &  1.01  \\  

% \caption{Conductance comparison of our \pgrid with \pgrow.
% Column 2 displays the median of the conductances found by \pgrow
