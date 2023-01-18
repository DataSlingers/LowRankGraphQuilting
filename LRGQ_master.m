
addpath(genpath('./'))

%%% Output folder for results
output_folder = '/D/ExampleOutput/';
%%% Maximum rank to estimate for in LRF and BSVD methods.
max_rank = 20;

%% Simulate graph structure

% Simulation parameters: 
%%% Number of data points for empirical covariance matrix.
n = 2000;
%%% Number of features.
p = 100;
%%% Number of blocks.
K = 2;
%%% Number of features per block; K * o must be greater than p.
o = 60;
%%% Low-rankness of simulated graph.
%%%% 1/eigg is the leading eigenvalue of the covariance matrix.
eigv = 0.1;
%%% Graph type.
graph_type = 'block';
%%% Seed for random generation.
seed = 1000;

% Create true cov mat and inverse cov mat with desired graph structure.
[Theta,Sigma]=graph_generator(p,seed,eigv,output_folder,graph_type);

%% Generate partially observed empirical covariance matrix

% Create empirical partially observed cov mat, full empirical cov mat,
%%% known observation pattern, raw data.
%%% S_obs denotes features in each block. Must have some overlap.
%%% Omega denotes to unique covariance matrix entries that are observed. 
[Sigma_obs, Sigma_n, S_obs, Omega, X] = patch_simulator(Sigma, n, K, o, seed, output_folder);

%% Run covariance imputation.
% Create cell of cells for each LRCCmethod containing imputed covariance matrices 
%%% for different hyperparameter values
[output] = LRGQ_imp(Sigma_obs, S_obs, max_rank, output_folder);


%% Run graph estimation.
% Saves estimated graphs from imputed covariance matrices in
% /GraphEstimation subfolder in output_folder.
% Uses saved imputed covariance matrices in /Imputed subfolder in
% output_folder as produced by LRGQ_imp.m.

R_call = ['R CMD BATCH ', '"--args ', output_folder ,'" ' ,'./GraphEstimation/GraphEstimationR.R ', ...
    output_folder, 'ROutput.Rout '];
system(R_call);

%% Evaluate imputation performance
% Saves best Frobenius norm in output folder as ImputationError.csv.
% Uses outputs from patch_simulator.m, LRGQ_imp.m.
[imp_err] = ImpAnalysis(output, Sigma_obs, Sigma, Sigma_n, S_obs, output_folder);

%% Evaluate graph performance
% Saves best F-1 score in output folder as F1Results.csv.

R_call = ['R CMD BATCH ', '"--args ', output_folder ,'" ' ,'./SimAnalysis/GraphAnalysis.R ', ...
    output_folder, 'RAnalysis.Rout '];
system(R_call);