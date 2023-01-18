%% main codes for simulations (data generation+low-rank imputation)
%% generate data from GGM with certain graph structures
%% how would SVT+rotation+glasso work?
function[imp_err] = GQ_lr_imp_sim(K, o, seed, n, p, eigg)
graph_type = 'star';
%generate graph
[Theta,Sigma]=graph_generator(p,seed,eigg,graph_type);
[U,Lambda,~]=svd(Sigma);
plot(sort(eig(Sigma),'descend'))
%generate data
Z=normrnd(0,1,n,p);
% X=Z*U*sqrt(Lambda)*U';
X=Z*U*(sqrt(Lambda))*U';
%observation pattern
[S_obs,Omega,Obs]=observation_pattern(p,K,o,seed);
dirname = sprintf('SimStudyStar/p_%d_n_%d_K_%d_o_%d_eig_%d_seed_%d_gtype_%s/',...
            p,n,K,o,eigg,seed,graph_type);     
mkdir(dirname);

writematrix(Omega,sprintf('SimStudyStar/p_%d_n_%d_K_%d_o_%d_eig_%d_seed_%d_gtype_%s/Omega.csv',...
            p,n,K,o,eigg,seed,graph_type));
        
for i=1:K
    filename = sprintf('SimStudyStar/p_%d_n_%d_K_%d_o_%d_eig_%d_seed_%d_gtype_%s/Obs_%d.csv',...
            p,n,K,o,eigg,seed,graph_type,i);
    writematrix(S_obs{i},filename);
end

Sigma_n=X'*X/n;plot(sort(eig(Sigma_n),'descend'));
filename = sprintf('SimStudyStar/p_%d_n_%d_K_%d_o_%d_eig_%d_seed_%d_gtype_%s/TruePrecision.csv',...
            p,n,K,o,eigg,seed,graph_type);
writematrix(Theta,filename)

filename = sprintf('SimStudyStar/p_%d_n_%d_K_%d_o_%d_eig_%d_seed_%d_gtype_%s/FullCov.csv',...
            p,n,K,o,eigg,seed,graph_type);
writematrix(Sigma_n,filename)
norm(Sigma_n-Sigma,'fro')


%Covariance calculated from the observed data
Sigma_obs=zeros(p,p);Sigma_obs(Omega)=Sigma(Omega);
%Sigma_obs=Cov_partial_obs(p,n,K,S_obs,Omega,X);
norm(Sigma_obs-Sigma_n,'fro')
norm((Sigma_obs-Sigma_n).*Obs,'fro') 
norm((Sigma_obs-Sigma).*Obs,'fro') 
norm(Sigma_n.*Obs,'fro')
norm((Sigma_obs-Sigma_n).*(1-Obs),'fro')
filename = sprintf('SimStudyStar/p_%d_n_%d_K_%d_o_%d_eig_%d_seed_%d_gtype_%s/ObsCov.csv',...
            p,n,K,o,eigg,seed,graph_type);
writematrix(Sigma_obs,filename)


%impute Sigma_n using several methods and models
r_list=1:20;
output=cell(6,3);output{1,1}=cell(length(r_list),1);
output{2,1}=output{1,1};output{3,1}=output{1,1};output{4,1}=output{1,1};
maxiter=100;tol=1e-5;eta_1=1/max(eig(Sigma_obs));eta_2=1/max(eig(Sigma_obs))/p;
for jj=1:length(r_list)
    %low-rank, svd+rotation
    [Sigma_hat_svdrot1,U_init1]=svd_rot(Sigma_obs,S_obs,r_list(jj));
    output{1,1}{jj}=Sigma_hat_svdrot1;
    output{1,2}=cat(1,output{1,2},norm(Sigma_hat_svdrot1-Sigma_n,'Fro'));
    output{1,3}=cat(1,output{1,3},norm(Sigma_hat_svdrot1-Sigma,'Fro'));
    %low-rank, gradient descent
    [U,c,Sigma_hat_gd1] = gd_planted(p,Omega,Sigma_obs,U_init1,0,0,maxiter,tol,eta_1,eta_2);
    output{2,1}{jj}=Sigma_hat_gd1;
    output{2,2}=cat(1,output{2,2},norm(Sigma_hat_gd1-Sigma_n,'Fro'));
    output{2,3}=cat(1,output{2,3},norm(Sigma_hat_gd1-Sigma,'Fro'));
    %planted model, svd+rotation with c estimated by the median
    med_c=median(diag(Sigma_obs)); 
    [Sigma_hat_svdrot2,U_init2]=svd_rot(Sigma_obs-med_c*eye(p),S_obs,r_list(jj));
    output{3,1}{jj}=Sigma_hat_svdrot2;
    output{3,2}=cat(1,output{3,2},norm(Sigma_hat_svdrot2-Sigma_n,'Fro'));
    output{3,3}=cat(1,output{3,3},norm(Sigma_hat_svdrot2-Sigma,'Fro'));
    %planted model, gradient descent
    [U,c,Sigma_hat_gd2] = gd_planted(p,Omega,Sigma_obs,U_init2,med_c,1,maxiter,tol,eta_1,eta_2);
    output{4,1}{jj}=Sigma_hat_gd2;
    output{4,2}=cat(1,output{4,2},norm(Sigma_hat_gd2-Sigma_n,'Fro'));
    output{4,3}=cat(1,output{4,3},norm(Sigma_hat_gd2-Sigma,'Fro')); 
end

%nuclear norm minimization
addpath('/supporting_codes_for_nuc_min');
lambda_list=[0 0.0003 0.001 0.003 0.01 0.03 0.1];
output{5,1}=cell(length(lambda_list),1);output{6,1}=cell(length(lambda_list),1);
maxiter=1000;
for jj=1:length(lambda_list)
    [L1,c1,Sigma_hat_nucmin1,output_nucmin1]  = nucmin(p,Omega,Sigma_obs(Omega),lambda_list(jj),0,maxiter,eta_1,eta_2);
    output{5,1}{jj}=Sigma_hat_nucmin1;
    output{5,2}=cat(1,output{5,2},norm(Sigma_hat_nucmin1-Sigma_n,'Fro'));
    output{5,3}=cat(1,output{5,3},norm(Sigma_hat_nucmin1-Sigma,'Fro')); 
    [L2,c2,Sigma_hat_nucmin2,output_nucmin2]  = nucmin(p,Omega,Sigma_obs(Omega),lambda_list(jj),1,maxiter,eta_1,eta_2);
    output{6,1}{jj}=Sigma_hat_nucmin2;
    output{6,2}=cat(1,output{6,2},norm(Sigma_hat_nucmin2-Sigma_n,'Fro'));
    output{6,3}=cat(1,output{6,3},norm(Sigma_hat_nucmin2-Sigma,'Fro')); 
end

%output imputed Sigma for fitting graphical lasso
method_list = ["lr_SVD+rot","lr_gd","plt_SVD+rot","plt_gd","lr_nucmin","plt_nucmin"];
nsetting = [repmat(length(r_list),4,1);repmat(length(lambda_list),2,1)];
for i=1:6
    for jj=1:nsetting(i)
        filename = sprintf('SimStudyStar/p_%d_n_%d_K_%d_o_%d_eig_%d_seed_%d_gtype_%s/%s_%d.csv',...
            p,n,K,o,eigg,seed,graph_type,method_list(i),jj);
        writematrix(output{i,1}{jj},filename);
    end
end

%Compare imputation error
imp_err = zeros(0);
dd1 = [];
dd2 = [];
for i=1:6
    [a1, d1] = min(output{i,2});
    [a2, d2] = min(output{i,3});
    dd1 = [dd1 d1];
    dd2 = [dd2 d2];
    imp_err = cat(1,imp_err,[min(output{i,2}) min(output{i,3})]);
end
dd1 = [dd1 0]; dd2 = [dd2 0];
method_list2 = ["lr_SVD+rot","lr_gd","plt_SVD+rot","plt_gd","lr_nucmin","plt_nucmin", "zero"];
imp_err2 = cat(1,imp_err,[norm(Sigma_obs-Sigma_n,'fro') norm(Sigma_obs-Sigma,'fro')]);
imp_err2 = cat(2, method_list2', dd1', dd2', imp_err2);
%2 values for each of the 7 methods
filename = sprintf('SimStudyStar/p_%d_n_%d_K_%d_o_%d_eig_%d_seed_%d_gtype_%s/imp_err.csv',...
            p,n,K,o,eigg,seed,graph_type);
writematrix(imp_err2,filename)
end

