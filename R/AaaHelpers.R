# Small internal helpers to detect placeholder NA fields in dataset/clustering objects.
.is_na_sentinel <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(TRUE)
  }
  if (is.atomic(x) && length(x) == 1 && is.na(x)) {
    return(TRUE)
  }
  FALSE
}

.has_value <- function(x) {
  !.is_na_sentinel(x)
}

.has_assignment_likelihoods <- function(clustering) {
  if (!("all.assignment.likelihoods" %in% names(clustering))) {
    return(FALSE)
  }
  all.assignment.likelihoods <- clustering$all.assignment.likelihoods
  if (!.has_value(all.assignment.likelihoods)) {
    return(FALSE)
  }
  if (is.null(dim(all.assignment.likelihoods))) {
    return(FALSE)
  }
  nrow(all.assignment.likelihoods) > 0 && ncol(all.assignment.likelihoods) > 0
}

.dpclust_num_stored_iters <- function(no.iters, no.iters.burn.in, thin_s_i) {
  if (!thin_s_i) {
    return(as.integer(no.iters))
  }
  stored_iters <- (no.iters.burn.in + 1):no.iters
  stored_iters <- stored_iters[stored_iters != 1]
  if (length(stored_iters) > 1000) {
    stored_iters <- floor(no.iters.burn.in + (1:1000) * (no.iters - no.iters.burn.in) / 1000)
  }
  as.integer(length(unique(stored_iters)))
}

.estimate_dpclust_memory_gb <- function(no.muts, no.samples, no.iters, no.iters.burn.in, max.considered.clusters, thin_s_i, keep_aux_fields) {
  bytes_per_double <- 8
  bytes_per_int <- 4
  stored_iters <- .dpclust_num_stored_iters(no.iters, no.iters.burn.in, thin_s_i)

  bytes_pi_h <- as.numeric(no.iters) * as.numeric(max.considered.clusters) * as.numeric(no.samples) * bytes_per_double
  bytes_v_h <- as.numeric(no.iters) * as.numeric(max.considered.clusters) * bytes_per_double
  bytes_alpha <- as.numeric(no.iters) * bytes_per_double
  bytes_s_i <- as.numeric(stored_iters) * as.numeric(no.muts) * bytes_per_int

  # Input matrices are typically pre-existing in R, but can be temporarily duplicated across the R<->C++ bridge.
  bytes_inputs <- 5 * as.numeric(no.muts) * as.numeric(no.samples) * bytes_per_double
  bytes_working <- as.numeric(no.muts) * as.numeric(no.samples) * bytes_per_double # burden_denom_inv

  bytes_aux <- 0
  if (keep_aux_fields) {
    bytes_aux <- 2 * as.numeric(no.muts) * as.numeric(no.samples) * bytes_per_double
  }

  bytes_total <- bytes_pi_h + bytes_v_h + bytes_alpha + bytes_s_i + bytes_inputs + bytes_working + bytes_aux
  gb_raw <- bytes_total / (1024^3)

  list(
    total_gb = gb_raw,
    raw_gb = gb_raw,
    stored_iters = stored_iters,
    components_gb = c(
      pi_h = bytes_pi_h / (1024^3),
      S_i = bytes_s_i / (1024^3),
      V_h_alpha = (bytes_v_h + bytes_alpha) / (1024^3),
      input_bridge = bytes_inputs / (1024^3),
      working = bytes_working / (1024^3),
      aux = bytes_aux / (1024^3)
    )
  )
}

.memory_guard_plan <- function(no.muts, no.samples, no.iters, no.iters.burn.in, max.considered.clusters, thin_s_i, keep_aux_fields, memory_limit_gb, verbose = TRUE) {
  estimate <- .estimate_dpclust_memory_gb(
    no.muts = no.muts,
    no.samples = no.samples,
    no.iters = no.iters,
    no.iters.burn.in = no.iters.burn.in,
    max.considered.clusters = max.considered.clusters,
    thin_s_i = thin_s_i,
    keep_aux_fields = keep_aux_fields
  )

  if (is.na(memory_limit_gb)) {
    return(list(thin_s_i = thin_s_i, keep_aux_fields = keep_aux_fields, estimate = estimate, limited = FALSE))
  }
  if (!is.finite(memory_limit_gb) || memory_limit_gb <= 0) {
    stop("memory_limit_gb must be a positive finite number when provided.")
  }

  if (verbose) {
    log_info(sprintf("Memory guard estimate: %.2f GB (limit %.2f GB)", estimate$total_gb, memory_limit_gb))
  }

  if (estimate$total_gb <= memory_limit_gb) {
    return(list(thin_s_i = thin_s_i, keep_aux_fields = keep_aux_fields, estimate = estimate, limited = FALSE))
  }

  detail_msg <- sprintf(
    paste0(
      "Estimated peak memory %.2f GB exceeds memory_limit_gb %.2f GB with default memory-saving settings. ",
      "Estimated components (GB): pi.h=%.2f, S.i=%.2f, V.h+alpha=%.2f, input_bridge=%.2f, working=%.2f, aux=%.2f. ",
      "Try reducing no.iters, max.considered.clusters, or num_muts_sample."
    ),
    estimate$total_gb, memory_limit_gb,
    estimate$components_gb["pi_h"], estimate$components_gb["S_i"], estimate$components_gb["V_h_alpha"],
    estimate$components_gb["input_bridge"], estimate$components_gb["working"], estimate$components_gb["aux"]
  )
  stop(detail_msg)
}
