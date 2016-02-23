digits=1:10;

nn = 3; %number of nearest neighbors
%nn = 10; %number of nearest neighbors


%% Load the data
mainpath = fileparts(which('isorank')); % addpath containing /data/
load(fullfile(mainpath, 'data', 'usps.mat'));

%% digit subset
patterns = [];
labels = [];
for d = digits,
    d
    filt = train_labels(d,:)==1;
    labels = [labels d*ones(1,sum(filt))];
    patterns = [patterns train_patterns(:,filt)];
    filt = test_labels(d,:)==1;
    labels = [labels d*ones(1,sum(filt))];
    patterns = [patterns test_patterns(:,filt)];
end
reallabels = labels;
[~,~,labels] = unique(reallabels);
rad = 2.5;
Dist = pdist2(patterns',patterns','euclidean');
K = exp(-Dist.^2/(2*rad^2));
G = K - diag(diag(K)); % remove diagonal
% the next steps are something I'd rather exclude
G(G<sqrt(eps(1))) = 0;
df = sum(G,2) >= sqrt(eps(1));
G = G(df,df);
labels = labels(df);
nmissing = size(K,1) - sum(df);
%% sparsify the graph
    G = (sparsify_large(G,'nearest',nn));
    G = spones(G);
    
    save(['usps_' num2str(nn) 'nn.mat'],'G','labels');