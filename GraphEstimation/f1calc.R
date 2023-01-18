# Calculates F1 scores, TPR, FPR, FDR, FNR for edge selection.
# Inputs: pxp matrix of estimated graph, pxp matrix of true graph.
# Outputs: list of scores.
f1calc <- function(test, true){
  tp <- 0
  fp <- 0
  fn <- 0
  tn <- 0
  for(ii in 1:(nrow(test) - 1)){
    for(jj in (ii + 1):ncol(test)){
      if(true[ii, jj] & test[ii, jj]){
        tp <- tp + 1
      }
      if(!true[ii, jj] & !test[ii, jj]){
        tn <- tn + 1
      }
      if(!true[ii, jj] & test[ii, jj]){
        fp <- fp + 1
      }
      if(true[ii, jj] & !test[ii, jj]){
        fn <- fn + 1
      }
    }
  }
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  fdr <- fp / (tp + fp)
  fnr <- fn / (tp + fn)
  prec <- tp / (tp + fp)
  rec <- tp / (tp + fn)
  f1 <- 2 * prec * tpr / (prec + tpr)
  return(list(tpr = tpr, fpr = fpr, fdr = fdr, prec = prec, f1 = f1))
}


f1calc_vec <- function(test, true){
  tp <- 0
  fp <- 0
  fn <- 0
  tn <- 0
  for(ii in 1:length(test)){
    if(true[ii] & test[ii]){
      tp <- tp + 1
    }
    if(!true[ii] & !test[ii]){
      tn <- tn + 1
    }
    if(!true[ii] & test[ii]){
      fp <- fp + 1
    }
    if(true[ii] & !test[ii]){
      fn <- fn + 1
    }
  }
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  fdr <- fp / (tp + fp)
  fnr <- fn / (tp + fn)
  prec <- tp / (tp + fp)
  rec <- tp / (tp + fn)
  f1 <- 2 * prec * tpr / (prec + tpr)
  return(list(tpr = tpr, fpr = fpr, fdr = fdr, prec = prec, f1 = f1))
}
