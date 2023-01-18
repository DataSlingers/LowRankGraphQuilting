function [S_obs,Omega,Obs]=observation_pattern(p,K,o,seed)
% Generate random observation patterns (indices V_1,...,V_K).
% Inputs:
%%% p: total number of features.
%%% K: number of observation blocks. 
%%% o: size of each observed block.
%%% seed: random seed.
% Outputs: 
%%% S_obs: cells listing features in each observation subset.
%%% Omega: list of unique covariance matrix entries that are part of the
%%%% observed set.
%%% Obs: binary matrix, denotes if corresponding entry in empirical
%%%% covariance is in observed set.

% Total amount of overlap.
overlap=o*K-p;
rng(seed);
% Create the K blocks
%%% Make sure overlap sizes of each block are as even as possible.
ind=randperm(p);S_obs=cell(K,1);
overlap_size=floor(overlap/(K-1));
for i=1:K
    if i<K
        S_obs{i}=ind(((o-overlap_size)*(i-1)+1):((o-overlap_size)*(i-1)+o));
    else
        S_obs{i}=ind(((o-overlap_size)*(K-1)+1):p);
    end
end
% Create list of unique entries in covariance matrix that are observed
Omega=zeros(0);
for i=1:length(S_obs)
    Omega=cat(1,Omega,reshape((S_obs{i}-1)*p+S_obs{i}',length(S_obs{i})^2,1));
end
Omega=unique(Omega);Obs=zeros(p,p);Obs(Omega)=1;
end