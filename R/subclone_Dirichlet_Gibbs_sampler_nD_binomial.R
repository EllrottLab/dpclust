#
# DPClust core algorithm
#

subclone.dirichlet.gibbs <- function(mutCount, WTCount, totalCopyNumber = array(1, dim(mutCount)), normalCopyNumber = array(2, dim(mutCount)), copyNumberAdjustment = array(1, dim(mutCount)), C = 30, cellularity = rep(1, ncol(mutCount)), iter = 1000, conc_param = 1, cluster_conc = 10, keep_aux_fields = FALSE, num_threads = NA_integer_, stored_iters = integer(0)) {
  if (is.null(copyNumberAdjustment)) {
    copyNumberAdjustment <- array(1, dim(mutCount))
  }

  log_info("Converting data to matrices for Rcpp...")
  ensure_numeric_matrix <- function(x) {
    if (!is.matrix(x)) {
      x <- as.matrix(x)
    }
    if (!is.double(x)) {
      storage.mode(x) <- "double"
    }
    x
  }
  # Ensure inputs are in the correct format for Rcpp while avoiding unnecessary copies.
  mutCount <- ensure_numeric_matrix(mutCount)
  WTCount <- ensure_numeric_matrix(WTCount)
  totalCopyNumber <- ensure_numeric_matrix(totalCopyNumber)
  normalCopyNumber <- ensure_numeric_matrix(normalCopyNumber)
  copyNumberAdjustment <- ensure_numeric_matrix(copyNumberAdjustment)
  if (!is.double(cellularity)) {
    storage.mode(cellularity) <- "double"
  }

  log_info(paste("Calling Rcpp Gibbs sampler for", nrow(mutCount), "mutations and", iter, "iterations..."))
  # Call C++ implementation
  cpp_threads <- if (is.na(num_threads)) as.integer(-1) else as.integer(num_threads)
  res <- subclone_dirichlet_gibbs_cpp(mutCount, WTCount, totalCopyNumber, normalCopyNumber, copyNumberAdjustment, C, cellularity, iter, conc_param, cluster_conc, keep_aux_fields, cpp_threads, as.integer(stored_iters), log_info)
  if (!keep_aux_fields) {
    # mutBurdens is large and not used by the current downstream pipeline.
    res$mutBurdens <- NULL
  }
  log_info("Gibbs sampler completed.")

  return(res)
}
