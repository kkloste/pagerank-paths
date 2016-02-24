

%% Load the data
load('usps.mat');

digits=1:10;
% nn = 3; 
% nn = 10; %number of nearest neighbors
nn = 30;

%% digit subset
patterns = [];
labels = [];
for d = digits,
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
    A = (sparsify_large(G,'nearest',nn));
    A = spones(A);
    
    save(['usps_' num2str(nn) 'nn.mat'],'A','labels');