clear
clc
close all
rng(1)

format long

n = 1000;

r = 5;

U_truth = orth(randn(n,r));

V_truth = orth(randn(n,r));

kappa = 1000;
S_truth = diag(linspace(kappa, 1, r));

M = U_truth * S_truth * V_truth';

%%  construct the observed matrix
p= 0.2;        %probability of success
A=rand(n,n);
A=(A<p);


sigma = 1 * 10^-3;
E = normrnd(0,sigma, [n,n]);
M_obs = (M+E) .* A;
lambda = 2.5 *  sqrt(n*p)  * sigma;
% sqrt(n * p) * sigma


disp(n^2 * p / (2 * n * r * log(n)))
%%  proximal

Z_prox = M;
T = 1000;
eta = 1;

error_prox = zeros(T,1);
for t = 1:T
    
    Z_new = SVT(Z_prox - eta * (Z_prox - M - E) .* A, lambda * eta);
    grad_norm = norm((Z_new - Z_prox) / eta, 'fro') / norm(M,'fro');
    Z_prox = Z_new;
    error_prox(t) = 1/2 * norm(Z_prox .* A -M_obs, 'fro')^2 + lambda * sum(svd(Z_prox));
%     [t grad_norm/(1e-8) rank(Z_prox)]
    
    if mod(t,1) == 0
        disp([t, grad_norm, rank(Z_prox), norm(Z_prox - M) / norm(M)])
    end
    if grad_norm < 1e-8
        break
    end
end
error_cvx = norm(Z_prox - M, 'fro') / norm(M,'fro')





