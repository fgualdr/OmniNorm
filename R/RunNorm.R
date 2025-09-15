#' @title RunNorm
#' @description Executes the complete normalization pipeline using mixture models fitted to log-ratio distributions between samples.
#'
#' @param mat_path Input data matrix. Can be:
#'   - A file path to a .csv, .tsv, or .txt file with rownames in the first column
#'   - A dense matrix or data.frame
#'   - A sparse dgCMatrix object (from the Matrix package)
#'
#' @param row_name_index Integer indicating the column index of rownames in the input file (default = 1). Ignored if input is a matrix or dgCMatrix.
#' @param design_path Path to the design matrix. Must contain columns: Sample_ID, Sample_Condition, and Sample_Replicate. Layers in conditions are defined by underscores (_).
#' @param fix_reference Character vector of Sample_IDs to use as normalization reference. Can also be:
#'   - "random" to randomly select ~30% of samples
#'   - A numeric fraction (e.g. 0.3) to select top-N by intensity
#' @param n_pop Number of populations expected in sample-vs-reference fits (default = 1). Set to NULL to auto-detect.
#' @param n_pop_reference Number of populations for reference-vs-reference fits (default = 1). NULL for auto.
#' @param saving_path Output directory to save results (required if Save_results = TRUE).
#' @param sigma_times Multiplier of standard deviation to define invariant interval around the dominant mode (default = 1).
#' @param dist_family Distribution family name to use for fitting. Passed to mixsmsn (default = "Skew.normal").
#' @param Norm_plot Logical. Whether to generate and save QC plots (default = TRUE). Large datasets may be subsampled for plotting.
#' @param Save_results Logical. Whether to save outputs to disk (default = TRUE).
#' @param BiocParam Optional BiocParallel backend. If NULL, defaults to MulticoreParam(workers = cpus).
#' @param return_full Optional whether to return full models per
#' @param cpus Integer. Used if BiocParam is NULL to set number of parallel workers.
#'
#' @return A named list with:
#'   - scaling_factors: a data.frame of computed scaling coefficients
#'   - norm_mat: a normalized sparse matrix (dgCMatrix)
#'   - If return_full = TRUE, also includes model_list: compact model fits
#'
#' @export
#'
#' @import BiocParallel
#' @importFrom data.table as.data.table data.table := setnames rbindlist CJ
#' @importFrom tools file_path_sans_ext
#' @importFrom Matrix colSums rowMeans Diagonal


RunNorm <- function(mat_path,
                    design_path,
                    fix_reference = NULL,
                    row_name_index = 1,
                    n_pop = 1,
                    n_pop_reference = 1,
                    saving_path = NULL,
                    sigma_times = 1,
                    dist_family = "Skew.normal",
                    Norm_plot = TRUE,
                    Save_results = TRUE,
                    return_full = FALSE,
                    BiocParam = NULL,
                    cpus = 2) {

  ## ---- guards ----
  if (isTRUE(Save_results) && is.null(saving_path)) {
    stop("To save plots/results, please provide 'saving_path'.")
  }
  if (!is.null(n_pop_reference) && n_pop_reference > 3) stop("n_pop_reference must be between 1 and 3.")
  if (!is.null(n_pop) && n_pop > 3) stop("n_pop must be between 1 and 3.")

  ## ---- ensure BiocParam ----
  # before creating BPPARAM
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }

  ## ---- import inputs ----
  cat("|| Import matrix\n")
  mat <- parse_input_matrix(mat_path, row_name_index = row_name_index)
  cat("|| DONE\n")

  cat("|| Import Design table\n")
  design <- parse_design(design_path)
  w <- which(design$Sample_ID %in% colnames(mat))
  if (length(w) == 0) stop("No Sample_IDs from design table match matrix column names.")
  design <- design[w, ]
  mat    <- mat[, design$Sample_ID, drop = FALSE]
  cat("|| DONE\n")

  ## ---- resolve fix_reference ----
  if (is.null(fix_reference)) {
    fix_reference <- design$Sample_ID
    message("No fix_reference provided. Using all samples.")
  } else if (is.character(fix_reference) && length(fix_reference) == 1 && fix_reference == "random") {
    set.seed(123456)
    n <- min(50L, max(1L, floor(nrow(design) * 0.3)))
    fix_reference <- sample(design$Sample_ID, n)
    message("Randomly selected ", n, " samples as reference.")
  } else if (is.numeric(fix_reference) && length(fix_reference) == 1 && fix_reference < 1) {
    message("Using top ", fix_reference * 100, "% columns by total counts as reference.")
    counts  <- Matrix::colSums(mat)
    top_idx <- order(counts, decreasing = TRUE)[1:round(ncol(mat) * fix_reference)]
    fix_reference <- colnames(mat)[top_idx]
    if (length(fix_reference) > 50) {
      warning("More than 50 samples selected as reference; limiting to 50 for performance.")
      fix_reference <- fix_reference[1:50]
    }
  } else if (is.numeric(fix_reference) && all(fix_reference %in% seq_len(ncol(mat)))) {
    fix_reference <- colnames(mat)[fix_reference]
  } else if (!all(fix_reference %in% design$Sample_ID)) {
    stop("Invalid fix_reference: must be 'random', a fraction, numeric indices, or Sample_IDs present in design.")
  }

  ## ---- outputs to fill ----
  pairs_scaling <- stats::setNames(numeric(0), character(0))
  pairs_mean_g  <- stats::setNames(numeric(0), character(0))
  model_list    <- NULL   # optional, only if return_full

  if (length(fix_reference) > 1) {

    cat("|| Multi-reference mode: ", length(fix_reference), " reference samples ||\n")

    design_ref <- design[design$Sample_ID %in% fix_reference, , drop = FALSE]
    mat_ref    <- mat[, design_ref$Sample_ID, drop = FALSE]

    ## ----- A) REFERENCE-REFERENCE (bplapply) -----
    pairs_df <- generate_pairwise_comparisons(design_ref$Sample_ID)
    if (nrow(pairs_df) == 0) stop("No pairs generated among reference samples.")

    cat("|| Reference-reference: ", nrow(pairs_df), " pairs (bplapply) ||\n")

    pair_results <- BiocParallel::bplapply(
      split(pairs_df, seq_len(nrow(pairs_df))),
      function(task)
        fun_refref(task, mat_ref, n_pop_reference, sigma_times, dist_family, Norm_plot),
      BPPARAM = BiocParam
    )

    pair_results <- Filter(Negate(is.null), pair_results)
    if (!length(pair_results)) stop("All reference pair fits failed.")

    dt <- data.table::rbindlist(lapply(pair_results, \(x)
      data.table::data.table(sample = x$sample, reference = x$reference, scaling = x$scaling)
    ), use.names = TRUE)
    dt <- dt[is.finite(scaling)]
    pairscombo <- dt[, .(av_scaling = mean(scaling, na.rm = TRUE)), by = sample]
    data.table::setnames(pairscombo, "sample", "samples")

    ## build average_reference as weighted mean of scaled refs
    mm <- mat[, pairscombo$samples, drop = FALSE]
    w  <- pairscombo$av_scaling; names(w) <- pairscombo$samples
    k  <- length(w)
    avg_ref <- as.numeric(mm %*% w) / k
    names(avg_ref) <- rownames(mm)
    
    # remove unused objects:
    rm(mat_ref, design_ref, pairs_df, mm, w, pairscombo, pair_results)
    gc(FALSE)

    ## ----- B) SAMPLE vs AVERAGE_REFERENCE (bplapply) -----
    cat("|| Sample-average_reference (bplapply) ||\n")
    avg_ref_vec <- avg_ref[rownames(mat)]
    sample_ids  <- design$Sample_ID
    cat("Number of samples: ", length(sample_ids), "\n")

    avg_results <- BiocParallel::bplapply(
      as.list(sample_ids),
      function(sname)
        fun_avgref(sname, mat, avg_ref_vec, n_pop, sigma_times, dist_family, Norm_plot),
      BPPARAM = BiocParam
    )

    avg_results <- Filter(Negate(is.null), avg_results)
    if (!length(avg_results)) stop("All sample vs average_reference fits failed.")

    sc_dt <- data.table::rbindlist(lapply(avg_results, \(x)
      data.table::data.table(sample = x$sample, reference = x$reference, scaling = x$scaling)
    ), use.names = TRUE)
    sc_dt <- sc_dt[data.table::CJ(sample = design$Sample_ID), on = "sample"]  # align to design
    pairs_scaling <- stats::setNames(sc_dt$scaling, sc_dt$sample)

    ## Free memory
    rm( sc_dt)  # free memory

    ## Optional: plotting (throttle to avoid explosion)
    if (isTRUE(Norm_plot)) {
      cat("|| Plotting pair models ||\n")

      # map: sample -> model (may be NULL)
      models_by_sample <- stats::setNames(
        lapply(avg_results, `[[`, "model"),
        vapply(avg_results, `[[`, character(1), "sample")
      )

      plot_samples <- names(models_by_sample)
      plot_samples <- plot_samples[!vapply(models_by_sample, is.null, logical(1))]
      plot_samples <- intersect(plot_samples, colnames(mat))

      avg_ref_vec_plot <- avg_ref_vec  # local copy for workers

      BiocParallel::bplapply(
        as.list(plot_samples),
        function(sname) {
          mdl <- models_by_sample[[sname]]
          if (is.null(mdl)) return(NULL)

          s <- as.numeric(mat[, sname, drop = FALSE])
          r <- avg_ref_vec_plot

          keep <- is.finite(r) & is.finite(s) & (r + s) > 0
          if (!any(keep)) return(NULL)

          rn <- rownames(mat)[keep]                 # <- add rownames
          X2 <- cbind(r[keep], s[keep])  # small dense 2-col matrix
          colnames(X2) <- c("average_reference", sname)
          rownames(X2) <- rn                        # <- set rownames explicitly

          plot_pair_model(list(
            model          = mdl,
            mat            = X2,                    # small dense 2-col matrix w/ rownames
            sample_name    = sname,
            reference_name = "average_reference",
            saving_path    = saving_path,
            n_pop          = mdl$n_pop
          ))
          NULL
        },
        BPPARAM = BiocParam
      )
    }
    ## Free memory
    rm(avg_results)  # free memory

  } else {

    ## ----- C) SINGLE REFERENCE (bplapply) -----
    ref <- fix_reference[1]
    cat("|| Single reference mode (bplapply): ", ref, " ||\n")

    sample_ids <- setdiff(design$Sample_ID, ref)
    r_full <- as.numeric(mat[, ref, drop = FALSE])  # slice once on master
    keep_r <- is.finite(r_full)

    cat("Number of samples: ", length(sample_ids), "\n")

    fix_results <- BiocParallel::bplapply(
      as.list(sample_ids),
      function(sname)
        fun_single(sname, mat, ref, r_full, keep_r, n_pop, sigma_times, dist_family, Norm_plot),
      BPPARAM = BiocParam
    )
    fix_results <- Filter(Negate(is.null), fix_results)
    if (!length(fix_results)) stop("All sample vs fixed-reference fits failed.")

    sc_dt <- data.table::rbindlist(lapply(fix_results, \(x)
      data.table::data.table(sample = x$sample, reference = x$reference, scaling = x$scaling)
    ), use.names = TRUE)
    sc_dt <- rbind(sc_dt, data.table::data.table(sample = ref, reference = ref, scaling = 1.0))
    sc_dt <- sc_dt[data.table::CJ(sample = design$Sample_ID), on = "sample"]
    pairs_scaling <- stats::setNames(sc_dt$scaling, sc_dt$sample)

    ## Optional: plotting (throttle to avoid explosion)
    if (isTRUE(Norm_plot)) {
      cat("|| Plotting pair models ||\n")

      # map: sample -> model (may be NULL)
      models_by_sample <- stats::setNames(
        lapply(fix_results, `[[`, "model"),
        vapply(fix_results, `[[`, character(1), "sample")
      )

      plot_samples <- names(models_by_sample)
      plot_samples <- plot_samples[!vapply(models_by_sample, is.null, logical(1))]
      plot_samples <- intersect(plot_samples, colnames(mat))

      ref_vec_plot <- r_full  # local copy for workers

      BiocParallel::bplapply(
        as.list(plot_samples),
        function(sname) {
          mdl <- models_by_sample[[sname]]
          if (is.null(mdl)) return(NULL)

          s <- as.numeric(mat[, sname, drop = FALSE])
          r <- ref_vec_plot

          keep <- is.finite(r) & is.finite(s) & (r + s) > 0
          if (!any(keep)) return(NULL)

          rn <- rownames(mat)[keep]                 # <- add rownames
          X2 <- cbind(r[keep], s[keep])  # small dense 2-col matrix
          colnames(X2) <- c("average_reference", sname)
          rownames(X2) <- rn                        # <- set rownames explicitly

          plot_pair_model(list(
            model          = mdl,
            mat            = X2,                    # small dense 2-col matrix w/ rownames
            sample_name    = sname,
            reference_name = "average_reference",
            saving_path    = saving_path,
            n_pop          = mdl$n_pop
          ))
          NULL
        },
        BPPARAM = BiocParam
      )
    }

  }

  ## ---- apply scaling (backend-aware, no unnecessary materialization) ----
  cat("|| Apply scaling to the matrix ||\n")
  sample_ids <- design$Sample_ID
  scaling    <- pairs_scaling[sample_ids]

  ## In-memory sparse: slice + scale in place (fast)
  mat_subset <- mat[, sample_ids, drop = FALSE]
  norm_mat   <- scale_sparse_columns(mat_subset, scaling)
  backend    <- "sparse"

  ## ---- prepare design_scaling ----
  design_df <- as.data.frame(design)
  rownames(design_df) <- design_df$Sample_ID
  design_df$scaling <- as.numeric(NA)
  
  design_scaling <- design_df[sample_ids, , drop = FALSE]
  design_scaling$scaling <- scaling

  ## ---- save (size- and backend-aware) ----
  if (isTRUE(Save_results)) {
    file_name <- if (is.character(mat_path)) basename(tools::file_path_sans_ext(mat_path)) else "matrix"

    ## 1) Always write the normalization parameters
    design_path_out <- file.path(saving_path, "Normalisation_Parameters.txt")
    write.table(design_scaling, file = design_path_out, sep = "\t", quote = FALSE, row.names = FALSE)
    message("Saved normalization parameters: ", design_path_out)

    if (!inherits(norm_mat, "dgCMatrix")) {
      ## Realize and coerce to sparse to keep size down
      norm_mat <- as(Matrix::Matrix(as.matrix(norm_mat), sparse = TRUE), "dgCMatrix")
    }
    rds_path <- file.path(saving_path, paste0(file_name, "_normalized.rds"))
    saveRDS(norm_mat, file = rds_path, compress = "xz")
    message("Saved normalized sparse matrix as RDS: ", rds_path)

    ## Optional dense TXT for small outputs (kept conservative)
    if (ncol(norm_mat) <= 2000) {
      txt_path <- file.path(saving_path, paste0(file_name, "_normalized.txt"))
      write.table(as.matrix(norm_mat), file = txt_path, sep = "\t", quote = FALSE, col.names = NA)
      message("Saved normalized matrix as TXT (dense): ", txt_path)
    } else {
      message("Matrix has >2000 columns; skipping dense TXT export.")
    }
  }
  
  ## ---- return (lightweight) ----
  Result <- list(
    scaling_factors = design_scaling,
    norm_mat        = norm_mat
  )
  if (isTRUE(return_full)) Result$model_list <- model_list
  class(Result) <- "GMSN"
  return(Result)

}


