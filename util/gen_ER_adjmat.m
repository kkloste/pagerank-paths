function A = gen_ER_adjmat(n,p)
%% Generating an Erdos-Renyi graph  
% The Erdos-Renyi model is connect _n_ vertices with probability _p_. 
% Let's begin by setting these parameters for an example.

% n = 150;
% p = 0.01;

%%
% With these two parameters, we can instantiate the graph. The 
% variable |G| is the adjacency matrix for the graph. However,
% the first step doesn't treat edges symmetrically. The last
% two operations fix this and yield a symmetric adjacency matrix.

rand('seed',100); % reseed so you get a similar picture
G = rand(n,n) < p;
G = sparse(triu(G,1));
A = G + G'; 