#
# DPClust core algorithm
#

subclone.dirichlet.gibbs <- function(mutCount, WTCount, totalCopyNumber = array(1, dim(mutCount)), normalCopyNumber = array(2, dim(mutCount)), copyNumberAdjustment = array(1, dim(mutCount)), C = 30, cellularity = rep(1, ncol(mutCount)), iter = 1000, conc_param = 1, cluster_conc = 10) {
  if (is.null(copyNumberAdjustment)) {
    copyNumberAdjustment <- array(1, dim(mutCount))
  }

  print("Converting data to matrices for Rcpp...")
  # Ensure inputs are in the correct format for Rcpp
  mutCount <- as.matrix(mutCount)
  WTCount <- as.matrix(WTCount)
  totalCopyNumber <- as.matrix(totalCopyNumber)
  normalCopyNumber <- as.matrix(normalCopyNumber)
  copyNumberAdjustment <- as.matrix(copyNumberAdjustment)
  cellularity <- as.vector(cellularity)

  print(paste("Calling Rcpp Gibbs sampler for", nrow(mutCount), "mutations and", iter, "iterations..."))
  # Call C++ implementation
  res <- subclone_dirichlet_gibbs_cpp(mutCount, WTCount, totalCopyNumber, normalCopyNumber, copyNumberAdjustment, C, cellularity, iter, conc_param, cluster_conc)
  print("Gibbs sampler completed.")

  return(res)
}
