function[L,c,Sigma,output] = gd_planted(p,Omega,b,U_init,med,update_c,T,tol,eta_1,eta_2)
    c = 0;
    output = zeros(0,0); I = eye(p);
    L = U_init;
    Sigma = L*L';
    spiked_mat = update_c * c.*I;
    error = norm(Sigma(Omega) + spiked_mat(Omega) - b(Omega),2)^2;
    for t = 1:T
        grad_L = zeros(p,p);
        upL = (Sigma + spiked_mat - b) * L;
        grad_L = upL;
        if update_c
            upC = (Sigma + c.*I - b)*I;
            grad_c = sum(upC, 'all');
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
            L = L_temp;
            if update_c
                c = c_prev - eta_2 * grad_c;
            end
            Sigma = L*L';
            spiked_mat = update_c * c.*I;
            error = norm(Sigma(Omega) + spiked_mat(Omega) - b(Omega),2)^2;
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