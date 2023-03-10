function[Sigma_obs, Sigma_n, S_obs, Omega, X] = patch_simulator(Sigma, n, ...
    K, o, seed, output_folder)

% For empirical studies. 
% Takes true underlying covariance matrix, creates fully observed and 
% partially observed empirical covariance matrices.
% Inputs:
%%% Sigma: True covariance matrix.
%%% n: number of data points for empirical covariance matrix.
%%% K: number of subsets that features are observed in.
%%% o: number of features per block.
%%%% K * o must be greater than number of features.
%%% seed: starting seed.
%%% output_folder: output folder for saving results
% Outputs: 
%%% Sigma_obs: partially observed empirical covariance matrix.
%%% Sigma_n: empirical covariance matrix, sampled from true covariance
%%%% matrix.
%%% S_obs: cell of lists containing features in each observation subset.
%%%% Must have some overlap.
%%% Omega: list of unique covariance matrix entries that are part of the
%%%% observed set.
%%% X: random data generated from true underlying covariance matrix, used
%%%% to create empirical covariance matrix.

% How many features are there?
[p, ~] = size(Sigma);
% Create blocks of observed features.
[S_obs,Omega,Obs]=observation_pattern(p,K,o,seed);
% Simulate data for empirical covariance matrix.
[U,Lambda,~]=svd(Sigma);
Z=normrnd(0,1,n,p);
X=Z*U*(sqrt(Lambda))*U';
% Calculate empirical covariance matrix.
Sigma_n=X'*X/n;
% Create partially observed covariance matrix.
Sigma_obs=zeros(p,p);
Sigma_obs(Omega)=Sigma_n(Omega);
mkdir(output_folder);
% Save results.
writematrix(Sigma_obs, strcat(output_folder, 'ObsCov.csv'));
writematrix(Omega, strcat(output_folder, 'Omega.csv'));
writecell(S_obs, strcat(output_folder, 'Obs.csv'));