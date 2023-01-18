# Runs the MADgq algorithm of Vinci et. al., 2019.

library(huge)

R <- vector(mode = "list", length = K)
L <- vector(mode = "list", length = K)
thrs1=0
val <- abs(Theta_full[-Omega])
thrs2=0.0005
for(k in 1:K){
  R[[k]]=Theta_est_temp[obs_groups[[k]],obs_groups[[k]]]-Theta_est_temp[obs_groups[[k]],-obs_groups[[k]]] %*%
    solve(Theta_est_temp[-obs_groups[[k]],-obs_groups[[k]]]) %*% Theta_est_temp[-obs_groups[[k]],obs_groups[[k]]]
  L[[k]] <- abs(R[[k]])>thrs1 &abs(R[[k]])<thrs2
}
sum(L[[1]])
H <- NULL
for(i in 1:p){
  ind <- lapply(obs_groups,function(v){which(v==i)})
  check <- TRUE
  for(k in 1:K){
    if(length(ind[[k]])>0){
      check <- check&(sum(L[[k]][ind[[k]],])>0)
    }
  }
  if(check){
    H <- c(H,i)
  }
}
g_reco <-matrix(rep(0,p^2),p,p);g_reco[H,]=1;g_reco[,H]=1;
g_reco[Omega] <- Theta_est_temp[Omega]!=0
TP_RECO <- sum(g_reco[-Omega]==1&Theta_full[-Omega]!=0)#TP
FP_RECO <- sum(g_reco[-Omega]==1&Theta_full[-Omega]==0)#FP
FN_RECO <- sum(g_reco[-Omega]==0&Theta_full[-Omega]!=0)#FN

# gvf1 <- (2*TP_RECO/(2*TP_RECO+FP_RECO+FN_RECO))#F1:0.1105217

      