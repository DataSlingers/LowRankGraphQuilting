function [L,c,Sigma,output]  = nucmin(p,Omega,b,lambda,update_c,T,eta_1,eta_2)
%% Proximal gradient descent for solving nuclear norm minimization
%% p is dimension of Sigma; Omega is the non-missing indices of Sigma;
%% b=Sigma(Omega); update_c=1 iff the spiked model is considered (o.w., low-rank assumption). T is the maximum number of iterations;
%% eta_1 and eta_2 are initial step sizes for the low-rank component L (corresponding to the S in the writeup) and the constant term c

L = zeros(p,p); c = 0;
output = zeros(0,0); I = eye(p);
Sigma = L + c * I;
error = 1 / 2 * norm(Sigma(Omega) - b,2)^2 + lambda * sum(svd(L));
for t = 1:T
    grad_L = zeros(p,p);
    grad_L(Omega) = Sigma(Omega) - b;
    if update_c
        grad_c = sum((Sigma(Omega)-b).*I(Omega));
    end
    if t > 1
        Delta_L = L - L_prev;
        inprod = Delta_L.*(grad_L - grad_L_prev);
        denom1 = sum(inprod(:));
        num1 = (norm(Delta_L,'fro'))^2; 
        eta1_temp = num1 / denom1;
        if ~isnan(eta1_temp)&&(eta1_temp>=1e-4)
            eta_1=eta1_temp;
        end
        if update_c
            denom2 = (c - c_prev)*(grad_c - grad_c_prev);
            num2 = (c - c_prev)^2;
            eta2_temp = num2 / denom2;
            if ~isnan(eta2_temp)&&(eta2_temp>=1e-7)
                eta_2=eta2_temp;
            end
        end
    end
    L_prev = L; error_prev = error; 
    Sigma_prev = Sigma;
    if update_c
        c_prev = c;
    end
    accept = false;
    while ~accept
        L_temp = L_prev - eta_1 * grad_L;
        L = SVT(L_temp, lambda * eta_1);
        if update_c
            c = c_prev - eta_2 * grad_c;
            Sigma = L + c * I;
        else
            Sigma = L;
        end
        error = 1 / 2 * norm(Sigma(Omega) - b,2)^2 + lambda * sum(svd(L));
        if update_c
            accept = (error <= error_prev)||((eta_1<1e-4)&&(eta_2<1e-7));
            eta_2 = eta_2/2;
        else
            accept = (error <= error_prev)||(eta_1<1e-4);
        end
        eta_1 = eta_1/2;
    end
    grad_L_prev = grad_L;
    if update_c
        grad_c_prev = grad_c;
    end
    grad_norm = norm(Sigma - Sigma_prev, 'fro') / eta_1 / norm(b,2);
    output = cat(1,output,[error, grad_norm, rank(L)]);
    if grad_norm < 1e-8
        break
    end
end
end