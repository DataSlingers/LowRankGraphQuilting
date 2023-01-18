function Z = SVT(X, tau)
    [U,S,V] = mexsvd(X);
    S = max(0, S-tau);
    Z = U * S * V';
end