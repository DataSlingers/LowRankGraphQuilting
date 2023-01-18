#----------------------------------------------------#
library(MTS)
#----------------------------------------------------#
proj_cov <- function(S, N, epsilon, mu, init_B, init_Lambda,tol,maxiter){
  ## projection of sample covariance matrix upon positive semi-definite cone w.r.t. weighted l_infty norm
  # S: sample covariance matrix (entry-wise estimation)
  # N: matrix of pair-wise sample sizes
  # init_B and init_Lambda need to be symmetric to ensure the symmetry of Sigma
  B <- init_B;
  Lambda <- init_Lambda;
  loss <- NULL
  for(i in 1:maxiter){
    if(i>1){
      Sigma_prev <- Sigma;
    }
    Sigma <- B + S + mu * Lambda;
    eig_out <- eigen(Sigma);
    Sigma <- eig_out[[2]]%*%diag(pmax(eig_out[[1]],epsilon))%*%t(eig_out[[2]])
    A <- Sigma - mu * Lambda - S;
    B_prev <- B;
    B <- A - proj_wt_l1(A,sqrt(N),mu/2);
    Lambda_prev <- Lambda;
    Lambda <- Lambda - (Sigma - B - S)/mu;
    loss <- c(loss,max(abs((Sigma - S)*sqrt(N))))
    if(i>1){
      chg_per <- norm(Sigma - Sigma_prev,"F")/norm(Sigma,"F")
      if(chg_per <= tol){
        break;
      }
    }
  }
  return(list(projected_S=Sigma,resid=B,loss=loss,iteration_num=i,chg_per=chg_per))
}
#----------------------------------------------------#
proj_wt_l1 <- function(A,omega,r){
  ## projection on weighted l1 ball
  # omega: weight matrix
  # r: radius
  z <- as.vector(A/omega);
  sort_ind <- order(z,decreasing = TRUE)
  a <- cumsum((as.vector(abs(omega*A)))[sort_ind])
  b <- cumsum((omega^2)[sort_ind])
  j <- length(z);
  for(i in 1:length(z)){
    j_new <- sum(z[sort_ind]>(a[j]-r)/b[j])
    if(j_new == j){
      break;
    }else{
      j <- j_new;
    }
  }
  c <- (a[j]-r)/b[j];
  proj_A <- pmax(abs(A) - c*omega, 0 )*sign(A)
  return(proj_A)
}



# epsilon <- 0.01;mu <- 1;init_B <- matrix(rep(0,p^2),p,p);init_Lambda <- matrix(rep(1,p^2),p,p);
#tol <- 0.001; maxiter <- 500; N <- matrix(rep(1,p^2),p,p);