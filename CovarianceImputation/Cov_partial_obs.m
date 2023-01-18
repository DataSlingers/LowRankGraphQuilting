function [obs_covs]=Cov_partial_obs(p,n,K,S_obs,Omega,X)
obs_covs = zeros(p, p);
for i = 1:K
    X_sub = X(S_obs{i}, S_obs{i});
    obs_covs(S_obs{i}, S_obs{i}) = X_sub'*X_sub/n;
end
