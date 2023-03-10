function[imp_err_all] = ImpAnalysis(output, Sigma_obs, Sigma, Sigma_n, ...
    S_obs, output_folder)

% Saves best imputation error for each method in output folder as
% ImputationError.csv.
% Inputs:
%%% output: cell, cells for each LRCCmethod containing imputed covariance matrices 
%%%% for different hyperparameter values.
%%% Sigma_obs: partially observed covariance matrix; should be 0 in all
%%%% non-observed pairwise entries.
% Outputs: 
%%% imp_err_all: Best imputation error, with index, for each method.
%%%% Compares to true underlying cov mat and empirical cov mat from
%%%% regenerated data.

imp_err_all=cell(6,2);
r_list=1:length(output{1});
for jj=1:length(r_list)
    %low-rank, svd+rotation
    imp_err_all{1,1}=cat(1,imp_err_all{1,1},norm(output{1}{jj}-Sigma_n,'Fro'));
    imp_err_all{1,2}=cat(1,imp_err_all{1,2},norm(output{1}{jj}-Sigma,'Fro'));
    %low-rank, gradient descent
    imp_err_all{2,1}=cat(1,imp_err_all{2,1},norm(output{2}{jj}-Sigma_n,'Fro'));
    imp_err_all{2,2}=cat(1,imp_err_all{2,2},norm(output{2}{jj}-Sigma,'Fro'));
    %planted model, svd+rotation 
    imp_err_all{3,1}=cat(1,imp_err_all{3,1},norm(output{3}{jj}-Sigma_n,'Fro'));
    imp_err_all{3,2}=cat(1,imp_err_all{3,2},norm(output{3}{jj}-Sigma,'Fro'));
    %planted model, gradient descent
    imp_err_all{4,1}=cat(1,imp_err_all{4,1},norm(output{4}{jj}-Sigma_n,'Fro'));
    imp_err_all{4,2}=cat(1,imp_err_all{4,2},norm(output{4}{jj}-Sigma,'Fro'));
end

%nuclear norm minimization
lambda_list=1:length(output{5});
for jj=1:length(lambda_list)
    imp_err_all{5,1}=cat(1,imp_err_all{5,1},norm(output{5}{jj}-Sigma_n,'Fro'));
    imp_err_all{5,2}=cat(1,imp_err_all{5,2},norm(output{5}{jj}-Sigma,'Fro'));
    imp_err_all{6,1}=cat(1,imp_err_all{6,1},norm(output{6}{jj}-Sigma_n,'Fro'));
    imp_err_all{6,2}=cat(1,imp_err_all{6,2},norm(output{6}{jj}-Sigma,'Fro'));
end

%Compare imputation error
imp_err_0 = zeros(0);
dd1 = [];
dd2 = [];
for i=1:6
    [a1, d1] = min(imp_err_all{i,1});
    [a2, d2] = min(imp_err_all{i,2});
    dd1 = [dd1 d1];
    dd2 = [dd2 d2];
    imp_err_0 = cat(1,imp_err_0,[min(imp_err_all{i,1}) min(imp_err_all{i,2})]);
end
dd1 = [dd1 0]; dd2 = [dd2 0];
method_list2 = ["BSVDgq, LR","LRFgq, LR","BSVDgq, ALR","LRFgq, ALR","NNgq, LR","NNgq, ALR", "Zero"];
imp_err = cat(1,imp_err_0,[norm(Sigma_obs-Sigma_n,'fro') norm(Sigma_obs-Sigma,'fro')]);
imp_err = cat(2, method_list2', dd1', dd2', imp_err);
filename = strcat(output_folder, 'ImputationError.csv');
imp_err = array2table(imp_err);
imp_err.Properties.VariableNames = ["Model","BestIndexOriginal","BestIndexData","BestFrobOriginal","BestFrobData"];
writetable(imp_err,filename);
