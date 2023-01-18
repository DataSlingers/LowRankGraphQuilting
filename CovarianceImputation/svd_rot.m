function [SigmaEst,C] = svd_rot(Sigma,S_obs,r)
%The SVD+rotation algorithm. 
%The missing entries in Sigma are 0; 
%S_obs is a cell array with size K, S_obs{k} is the observed indicies in
%kth batch. r is the specified rank.
K = length(S_obs);[p,~]=size(Sigma);
C=zeros(p,r);
[V,D,~]=svd(Sigma(S_obs{1},S_obs{1}));
[~,I]=sort(diag(D),'descend');
C(S_obs{1},:)=V(:,I(1:r))*sqrt(D(I(1:r),I(1:r)));
past_obs=S_obs{1};
for k=2:K
    overlap1=intersect(past_obs,S_obs{k});
    [V,D,~]=svd(Sigma(S_obs{k},S_obs{k}));
    [~,I]=sort(diag(D),'descend');
    C_temp=V(:,I(1:r))*sqrt(D(I(1:r),I(1:r)));
    overlap2=zeros(length(overlap1),1);
    for i=1:length(overlap2)
        overlap2(i)=find(S_obs{k}==overlap1(i));
    end
    [U,~,V]=svd(C_temp(overlap2,:)'*C(overlap1,:));
    W=U*V';
    C_overlap_temp=C(overlap1,:);
    C(S_obs{k},:)=C_temp*W;
    C(overlap1,:)=C_overlap_temp;
    past_obs=union(past_obs,S_obs{k});
end
SigmaEst=C*C';
end