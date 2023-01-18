args=(commandArgs(TRUE))
fold_wd = args[[1]]
print(fold_wd)

source("./GraphEstimation/f1calc.R", local = TRUE)

# Saves best F-1 score for each method in output folder as F1Results.csv.
# Inputs:
## fold_wd: directory location of all the estimated graphs matrices.
### Assumes that the imputed covariance matrices are in \GraphEstimates subfolder of fold_wd.
### Assumes that the imputed covariance matrices are the only files in \GraphEstimates subfolder of fold_wd.
graph_analysis <- function(fold_wd){
  # Find files
  all_files <- list.files(paste0(fold_wd, "/GraphEstimates"), full.names = TRUE)
  true_precision <- read.csv(paste0(fold_wd, "/TruePrecision.csv"), 
                             header=FALSE, stringsAsFactors=FALSE)
  f1_list <- list()
  # For each method, load estimate.
  for(ii in 1:length(all_files)){
    fn <- list.files(all_files[ii], full.names = TRUE)
    f1_list[[ii]] <- numeric()
    # Calculate F-1 score for each estimated graph for each hyperparameter value.
    for(kk in fn){
      est_graph <- read.csv(kk, header = TRUE, row.names = 1)
      f1_list[[ii]] <- c(f1_list[[ii]], f1calc(as.matrix(est_graph), as.matrix(true_precision))$f1)
    }
  }
  
  # Find best estimate by F-1 score for each method.
  ff1 <- data.frame("BestIndex" = c(unlist(lapply(f1_list, which.max))),
                    "BestF1Score" = c(unlist(lapply(f1_list, max, na.rm = TRUE))))
  rownames(ff1) <- data.frame(strsplit(all_files, "/"))[length(strsplit(all_files, "/")[[1]]), ]
  write.csv(ff1, paste0(fold_wd, "/F1Results.csv"))
}

# Run function
graph_analysis(fold_wd)
