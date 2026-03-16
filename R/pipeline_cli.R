#' DPClust Command Line Interface
#' @description Parses command line arguments and executes the DPClust pipeline.
#' @export
dpclust_cli <- function() {
  options(error = function() {
    # Get the raw calls
    calls <- sys.calls()

    msg <- sprintf("Fatal Error: %s\n\n--- Call Stack ---", geterrmessage())
    for (i in seq_along(calls)) {
      msg <- paste(msg, sprintf("[%2d] %s", i, deparse(calls[[i]], width.cutoff = 500)[1]), sep = "\n")
    }
    log_failure(msg)
    quit(save = "no", status = 1)
  })

  options(show.error.messages = TRUE)
  options(keep.source = TRUE)
  options(width = 10000)
  options(warn = 1) # Print warnings immediately
  options(bitmapType = "cairo")
  options(rgl.useNULL = TRUE)

  option_list <- list(
    # Core Analysis & Sample Info
    optparse::make_option(c("-r", "--run_sample"), type = "integer", default = NULL, help = "Sample index to run (from input file)", metavar = "integer"),
    optparse::make_option(c("-d", "--data_path"), type = "character", default = NULL, help = "Path to where dpinput data files are stored", metavar = "path"),
    optparse::make_option(c("-o", "--outputdir"), type = "character", default = getwd(), help = "Directory where the output is saved [default: %default]", metavar = "path"),
    optparse::make_option(c("-i", "--input"), type = "character", default = NULL, help = "Design file (TSV) with sample information", metavar = "path"),
    optparse::make_option(c("-p", "--prefix"), type = "character", default = NULL, help = "Optional prefix for output files (e.g. snv or indel)", metavar = "string"),

    # Clustering Parameters
    optparse::make_option(c("-a", "--analysis_type"), type = "character", default = "nd_dp", help = "Analysis type to run [default: %default]", metavar = "string"),
    optparse::make_option(c("--iterations"), type = "integer", default = 2000, help = "Number of MCMC iterations [default: %default]", metavar = "integer"),
    optparse::make_option(c("--burnin"), type = "integer", default = 1000, help = "Number of burnin iterations [default: %default]", metavar = "integer"),
    optparse::make_option(c("--mut_assignment_type"), type = "integer", default = 1, help = "Mutation assignment method [default: %default]", metavar = "integer"),
    optparse::make_option(c("--num_muts_sample"), type = "integer", default = 50000, help = "Downsampling threshold [default: %default]", metavar = "integer"),
    optparse::make_option(c("--min_muts_cluster"), type = "integer", default = -1, help = "Min mutations per cluster [default: %default]", metavar = "integer"),
    optparse::make_option(c("--min_frac_muts_cluster"), type = "numeric", default = 0.01, help = "Min fraction mutations per cluster [default: %default]", metavar = "numeric"),
    optparse::make_option(c("--bin_size"), type = "double", default = NULL, help = "Binsize for multi-dimensional density", metavar = "double"),
    optparse::make_option(c("--seed"), type = "integer", default = 123, help = "Random seed [default: %default]", metavar = "integer"),

    # Advanced / Behavior
    optparse::make_option(c("--species"), type = "character", default = "human", help = "Species (e.g. human, mouse) [default: %default]", metavar = "string"),
    optparse::make_option(c("--is_male"), type = "logical", default = NULL, help = "Explicitly set sex (autodetected if NULL)", metavar = "boolean"),
    optparse::make_option(c("--sample_snvs_only"), type = "logical", default = TRUE, help = "Only sample SNVs [default: %default]", metavar = "boolean"),
    optparse::make_option(c("--generate_cluster_ordering"), type = "logical", default = FALSE, help = "Generate phylogenetic ordering [default: %default]", metavar = "boolean"),
    optparse::make_option(c("--assign_sampled_muts"), type = "logical", default = TRUE, help = "Assign mutations removed during sampling [default: %default]", metavar = "boolean"),
    optparse::make_option(c("-k", "--keep_temp_files"), action = "store_true", default = FALSE, help = "Keep intermediate files"),

    # Hardware/Resources
    optparse::make_option(c("--num_threads"), type = "integer", default = NA, help = "Number of CPU threads [default: auto-detected]", metavar = "integer"),
    optparse::make_option(c("--memory_limit_gb"), type = "numeric", default = NA, help = "RAM limit in GB [default: auto-detected]", metavar = "numeric"),

    # Logging & Debug (Battenberg style)
    optparse::make_option(c("--verbose_logging"), type = "logical", default = FALSE, action = "store_true"),
    optparse::make_option(c("--logging_path"), type = "character", default = ".")
  )

  # Parse arguments
  parser <- optparse::OptionParser(option_list = option_list, description = "DPClust Pipeline CLI")
  opt <- optparse::parse_args(parser)

  log_setup(opt$logging_path, opt$verbose_logging)

  # Validation
  if (is.null(opt$run_sample) || is.null(opt$data_path) || is.null(opt$input)) {
    optparse::print_help(parser)
    stop("Missing required arguments: -r, -d, and -i are mandatory.", call. = FALSE)
  }

  log_info(strrep("=", 120))
  log_info("DPCLUST CLI: EXECUTION PARAMETERS")
  log_info(strrep("=", 120))

  opt_names <- sort(names(opt))
  for (name in opt_names) {
    if (name == "help") next
    val <- opt[[name]]
    log_info(sprintf("%-40s : %s", name, paste(val, collapse = ", ")))
  }
  log_info(strrep("=", 120))

  # Resource Auto-Detection
  if (is.na(opt$num_threads)) {
    opt$num_threads <- .get_available_cores()
    log_info(paste("Auto-detected", opt$num_threads, "CPU cores available."))
  }

  if (is.na(opt$memory_limit_gb)) {
    opt$memory_limit_gb <- .get_system_memory_gb()
    if (!is.na(opt$memory_limit_gb)) {
      log_info(paste("Auto-detected", round(opt$memory_limit_gb, 2), "GB system memory."))
    }
  }

  run_dpclust_pipeline(
    run_sample = opt$run_sample,
    data_path = opt$data_path,
    outputdir = opt$outputdir,
    input = opt$input,
    prefix = opt$prefix,
    analysis_type = opt$analysis_type,
    iterations = opt$iterations,
    burnin = opt$burnin,
    mut_assignment_type = opt$mut_assignment_type,
    num_muts_sample = opt$num_muts_sample,
    min_muts_cluster = opt$min_muts_cluster,
    min_frac_muts_cluster = opt$min_frac_muts_cluster,
    seed = opt$seed,
    species = opt$species,
    is_male = opt$is_male,
    sample_snvs_only = opt$sample_snvs_only,
    generate_cluster_ordering = opt$generate_cluster_ordering,
    assign_sampled_muts = opt$assign_sampled_muts,
    keep_temp_files = opt$keep_temp_files,
    num_threads = opt$num_threads,
    memory_limit_gb = opt$memory_limit_gb
  )
}

#' Run the DPClust pipeline for a specific sample
#' @export
run_dpclust_pipeline <- function(run_sample, data_path, outputdir = getwd(), input = NULL,
                                 prefix = NULL, analysis_type = "nd_dp",
                                 iterations = 2000, burnin = 1000,
                                 mut_assignment_type = 1, num_muts_sample = 50000,
                                 min_muts_cluster = -1, min_frac_muts_cluster = 0.01,
                                 seed = 123, species = "human", is_male = NULL,
                                 sample_snvs_only = TRUE, generate_cluster_ordering = FALSE,
                                 assign_sampled_muts = TRUE,
                                 keep_temp_files = FALSE, num_threads = NA,
                                 memory_limit_gb = NA) {
  options(bitmapType = "cairo")
  options(rgl.useNULL = TRUE)

  # Resolve all user-supplied paths to absolute paths immediately.
  # This is essential for container environments (e.g. Singularity with --pwd)
  # where the working directory is set externally and setwd() must never be relied upon.
  if (!is.null(input)) input <- normalizePath(input, mustWork = FALSE)
  if (!is.null(data_path)) data_path <- normalizePath(data_path, mustWork = FALSE)
  outputdir <- normalizePath(outputdir, mustWork = FALSE)

  # Check input
  if (is.null(input)) stop("Input design file must be provided.")
  if (!file.exists(input)) stop(paste("Input file not found:", input))

  # Ensure output directory exists and is writable
  if (!dir.exists(outputdir)) {
    log_info(paste("Output directory not found. Creating path:", outputdir))
    dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(outputdir)) {
      stop(paste("Failed to create output directory:", outputdir))
    }
    log_info("Successfully created output directory.")
  } else {
    log_info(paste("Valid output directory found at:", outputdir))
  }

  if (file.access(outputdir, 2) != 0) {
    stop(paste("ERROR: Output directory exists but is NOT writable:", outputdir))
  }
  log_info("Write-access verified for output directory.")

  # Parse the input file
  sample2purity <- read.table(input, header = T, stringsAsFactors = F)
  if (run_sample > length(unique(sample2purity$sample))) {
    stop(paste("run_sample index", run_sample, "exceeds number of unique samples", length(unique(sample2purity$sample))))
  }

  samplename <- unique(sample2purity$sample)[run_sample]
  datafiles <- sample2purity[sample2purity$sample == samplename, ]$datafile
  subsamples <- sample2purity[sample2purity$sample == samplename, ]$subsample
  cellularity <- sample2purity[sample2purity$sample == samplename, ]$cellularity

  # Sex detection
  if (is.null(is_male)) {
    if ("sex" %in% colnames(sample2purity)) {
      is_male <- (sample2purity[sample2purity$sample == samplename, ]$sex == "male")[1]
    } else {
      is_male <- TRUE
    }
  }

  if ("mutphasing" %in% colnames(sample2purity)) {
    mutphasingfiles <- sample2purity[sample2purity$sample == samplename, ]$mutphasing
  } else {
    mutphasingfiles <- NA
  }

  if ("cndatafile" %in% colnames(sample2purity)) {
    cndatafiles <- sample2purity[sample2purity$sample == samplename, ]$cndatafile
  } else {
    cndatafiles <- NA
  }

  log_info(paste("Running:", samplename))
  log_info(paste("Working dir:", outputdir))

  # Setup parameters
  run_params <- make_run_params(
    no.iters = iterations,
    no.iters.burn.in = burnin,
    mut.assignment.type = mut_assignment_type,
    num_muts_sample = num_muts_sample,
    is.male = is_male,
    min_muts_cluster = min_muts_cluster,
    min_frac_muts_cluster = min_frac_muts_cluster,
    species = species,
    assign_sampled_muts = assign_sampled_muts,
    keep_temp_files = keep_temp_files,
    generate_cluster_ordering = generate_cluster_ordering,
    num_threads = num_threads,
    memory_limit_gb = memory_limit_gb,
    sample.snvs.only = sample_snvs_only,
    remove.snvs = FALSE
  )

  datpath <- if (is.null(data_path)) "" else data_path
  sample_params <- make_sample_params(datafiles, cellularity, is_male, samplename, subsamples, mutphasingfiles, datpath = datpath, cndatafiles = cndatafiles)
  advanced_params <- make_advanced_params(seed)

  # Run clustering
  RunDP(
    analysis_type = analysis_type,
    run_params = run_params,
    sample_params = sample_params,
    advanced_params = advanced_params,
    outdir = outputdir
  )

  # Save run parameters
  .save_pipeline_parameters(
    samplename, prefix, outputdir, run_sample, analysis_type,
    iterations, burnin, mut_assignment_type, num_muts_sample,
    seed, assign_sampled_muts, keep_temp_files,
    min_muts_cluster, min_frac_muts_cluster
  )

  # Prefix renaming logic
  if (!is.null(prefix)) {
    .rename_output_with_prefix(samplename, outputdir, iterations, prefix)
  }

  log_info("DPClust pipeline completed.")
}

.get_available_cores <- function() {
  if (Sys.info()[["sysname"]] == "Linux") {
    # Check for cgroup CPU quota
    if (file.exists("/sys/fs/cgroup/cpu.max")) {
      cpu_max <- readLines("/sys/fs/cgroup/cpu.max", n = 1)
      parts <- strsplit(cpu_max, " ")[[1]]
      if (parts[1] != "max") {
        quota <- as.numeric(parts[1])
        period <- as.numeric(parts[2])
        return(max(1, floor(quota / period)))
      }
    }
    # Legacy cgroup v1
    if (file.exists("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")) {
      quota <- as.numeric(readLines("/sys/fs/cgroup/cpu/cpu.cfs_quota_us", n = 1))
      period <- as.numeric(readLines("/sys/fs/cgroup/cpu/cpu.cfs_period_us", n = 1))
      if (quota > 0) {
        return(max(1, floor(quota / period)))
      }
    }
    # Fallback to nproc if available
    nproc <- suppressWarnings(system("nproc", intern = TRUE))
    if (length(nproc) > 0) {
      return(as.integer(nproc))
    }
  }
  return(parallel::detectCores())
}

.get_system_memory_gb <- function() {
  if (Sys.info()[["sysname"]] == "Linux") {
    if (file.exists("/sys/fs/cgroup/memory.max")) {
      mem <- readLines("/sys/fs/cgroup/memory.max", n = 1)
      if (mem != "max") {
        return(as.numeric(mem) / (1024^3))
      }
    }
    if (file.exists("/sys/fs/cgroup/memory/memory.limit_in_bytes")) {
      mem <- readLines("/sys/fs/cgroup/memory/memory.limit_in_bytes", n = 1)
      if (as.numeric(mem) < 1e15) {
        return(as.numeric(mem) / (1024^3))
      }
    }
  }

  if (Sys.info()[["sysname"]] == "Darwin") {
    mem <- system("sysctl -n hw.memsize", intern = TRUE)
    return(as.numeric(mem) / (1024^3))
  } else if (Sys.info()[["sysname"]] == "Linux") {
    mem <- system("grep MemTotal /proc/meminfo | awk '{print $2}'", intern = TRUE)
    return(as.numeric(mem) / (1024^2))
  }
  return(NA)
}

.save_pipeline_parameters <- function(samplename, prefix, outdir, run_sample, analysis_type,
                                      iterations, burnin, mut_assignment_type, num_muts_sample,
                                      seed, assign_sampled_muts, keep_temp_files,
                                      min_muts_cluster, min_frac_muts_cluster) {
  param_file <- file.path(outdir, paste0(samplename, "__", if (is.null(prefix)) "" else paste0(prefix, "__"), "dpclust_run_parameters.tsv"))

  info_files <- Sys.glob(file.path(outdir, "*__allDirichletProcessInfo.txt"))
  input_loci_count <- if (length(info_files) == 1) length(readLines(info_files[1], warn = FALSE)) - 1 else NA

  params_df <- data.frame(
    parameter = c(
      "timestamp", "dpclust_version", "samplename", "run_sample", "analysis_type",
      "iterations", "burnin", "mut_assignment_type", "num_muts_sample", "seed",
      "assign_sampled_muts", "keep_temp_files", "min_muts_cluster",
      "min_frac_muts_cluster", "prefix", "input_loci_count"
    ),
    value = c(
      as.character(Sys.time()), as.character(packageVersion("DPClust")), samplename,
      run_sample, analysis_type, iterations, burnin, mut_assignment_type,
      num_muts_sample, seed, assign_sampled_muts, keep_temp_files,
      min_muts_cluster, min_frac_muts_cluster, if (is.null(prefix)) "NA" else prefix,
      input_loci_count
    ),
    stringsAsFactors = FALSE
  )

  write.table(params_df, file = param_file, sep = "\t", quote = FALSE, row.names = FALSE)

  best_cluster_files <- list.files(path = outdir, pattern = "bestClusterInfo\\.txt$", full.names = TRUE)
  if (length(best_cluster_files) > 0) {
    bestClusterInfo <- read.table(best_cluster_files[1], header = TRUE)
    log_info(paste0(
      "\n", strrep("=", 30), "\nBest cluster information\n", strrep("=", 30), "\n",
      paste(utils::capture.output(print(bestClusterInfo)), collapse = "\n")
    ))
    write.table(bestClusterInfo, file = param_file, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

.rename_output_with_prefix <- function(samplename, outdir, iterations, prefix) {
  files_with_iters <- list.files(path = outdir, pattern = paste0(iterations, "iters"), full.names = TRUE)
  for (f in files_with_iters) {
    new_f <- sub("_([0-9]+iters)", paste0("_", prefix, "_\\1"), f)
    if (f != new_f) file.rename(f, new_f)
  }

  png1 <- file.path(outdir, paste0(samplename, "_DirichletProcessplot_with_cluster_locations.png"))
  png2 <- file.path(outdir, paste0(samplename, "_DirichletProcessplot_with_cluster_locations_2.png"))

  if (file.exists(png1)) file.rename(png1, file.path(outdir, paste0(samplename, "_", prefix, "_DirichletProcessplot_with_cluster_locations.png")))
  if (file.exists(png2)) file.rename(png2, file.path(outdir, paste0(samplename, "_", prefix, "_DirichletProcessplot_with_cluster_locations_2.png")))
}
