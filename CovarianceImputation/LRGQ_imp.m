function[output] = LRGQ_imp(Sigma_obs, S_obs, max_rank, output_folder)

% Runs all of the low-rank graph quilting covariance imputation methods.
% Inputs:
%%% Sigma_obs: matrix, observed parts of covariance matrix; should be 0 in all
%%%% non-observed pairwise entries.
%%% S_obs: cell, list of features in each block.
%%% max_rank: integer, maximum rank to estimate for in LRF and BSVD methods.
%%% output_folder: string, desired location to save imputed covariance matrices.
% Outputs:
%%% output: cell, cells for each LRCCmethod containing imputed covariance matrices 
%%%% for different hyperparameter values.

[p, ~] = size(Sigma_obs);
max_rank = min(p, max_rank);
r_list=1:max_rank;
foldname = [output_folder, 'Imputed/'];
mkdir(foldname);
Omega=zeros(0);
for i=1:length(S_obs)
    Omega=cat(1,Omega,reshape((S_obs{i}-1)*p+S_obs{i}',length(S_obs{i})^2,1));
end
Omega=unique(Omega);

output=cell(6,1);output{1,1}=cell(length(r_list),1);
output{2,1}=output{1,1};output{3,1}=output{1,1};output{4,1}=output{1,1};
maxiter=100;tol=1e-5;eta_1=1/max(eig(Sigma_obs));eta_2=1/max(eig(Sigma_obs))/p;
for jj=1:length(r_list)
    %low-rank, svd+rotation
    [Sigma_hat_BSVDLR,U_init1]=svd_rot(Sigma_obs,S_obs,r_list(jj));
    output{1,1}{jj}=Sigma_hat_BSVDLR;
    %low-rank, gradient descent
    [U,c,Sigma_hat_LRFLR] = gd_planted(p,Omega,Sigma_obs,U_init1,0,0,maxiter,tol,eta_1,eta_2);
    output{2,1}{jj}=Sigma_hat_LRFLR;
    %planted model, svd+rotation with c estimated by the median
    med_c=median(diag(Sigma_obs)); 
    [Sigma_hat_BSVDALR,U_init2]=svd_rot(Sigma_obs-med_c*eye(p),S_obs,r_list(jj));
    output{3,1}{jj}=Sigma_hat_BSVDALR;
    %planted model, gradient descent
    [U,c,Sigma_hat_LRFALR] = gd_planted(p,Omega,Sigma_obs,U_init2,med_c,1,maxiter,tol,eta_1,eta_2);
    output{4,1}{jj}=Sigma_hat_LRFALR;
end

%nuclear norm minimization
lambda_list=[0 0.0003 0.001 0.003 0.01 0.03 0.1];
output{5,1}=cell(length(lambda_list),1);output{6,1}=cell(length(lambda_list),1);
maxiter=1000;
for jj=1:length(lambda_list)
    [L1,c1,Sigma_hat_NNLR,output_nucmin1] = nucmin(p,Omega,Sigma_obs(Omega),lambda_list(jj),0,maxiter,eta_1,eta_2);
    output{5,1}{jj}=Sigma_hat_NNLR;
    [L2,c2,Sigma_hat_NNALR,output_nucmin2] = nucmin(p,Omega,Sigma_obs(Omega),lambda_list(jj),1,maxiter,eta_1,eta_2);
    output{6,1}{jj}=Sigma_hat_NNALR;
end

%output imputed Sigma for fitting graphical lasso
method_list = ["BSVDgq_LR","LRFgq_LR","BSVDgq_ALR","LRFgq_ALR","NNgq_LR","NNgq_ALR"];
nsetting = [repmat(length(r_list),4,1);repmat(length(lambda_list),2,1)];
for i=1:6
    dirname = strcat(foldname, sprintf('%s/', method_list(i)));
    mkdir(dirname);
    for jj=1:nsetting(i)
        filename = strcat(dirname, sprintf('%s_%d.csv', method_list(i), jj));
        writematrix(output{i,1}{jj},filename);
    end
end