#' @title noceros
#' @description Count or retrieve non-zero elements from a numeric vector.
#'
#' @param x Numeric vector.
#' @param num Logical. If `TRUE`, returns the count; if `FALSE`, returns indices (default = TRUE).
#' @param k Threshold. Elements greater than `k` are considered non-zero (default = 0).
#'
#' @return Either a count or vector of indices of non-zero elements.
#' @export

noceros <- function(x, num = TRUE, k = 0) {
  nn <- length(which(x > k))
  if (num) {
    nn
  } else {
    if (nn > 0) {
      which(x > k)
    } else {
      NULL
    }
  }
}

#' @title quantile_normalisation
#' @description Apply quantile normalization to a data frame.
#'
#' @param df A numeric data frame or matrix to normalize.
#'
#' @return A quantile-normalized data frame.
#' @export

quantile_normalisation <- function(df) {
  df_rank <- apply(df, 2, rank, ties.method = "min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  index_to_mean <- function(my_index, my_mean) {
    return(my_mean[my_index])
  }
  df_final <- apply(df_rank, 2, index_to_mean, my_mean = df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

#' @title median_of_ratios_normalization
#' @description Median-of-ratios normalization method (as in DESeq2).
#'
#' @param data A data frame or matrix of counts.
#'
#' @return A normalized data frame.
#' @export

median_of_ratios_normalization <- function(data) {
  # take the log
  log_data <- log(data)
  # find the psuedo-references per sample by taking the geometric mean
  log_data <- dplyr::mutate(log_data, pseudo_reference = rowMeans(log_data))
  log_data <- tibble::rownames_to_column(log_data, "gene")
  log_data <- dplyr::filter(log_data, .data$pseudo_reference != "-Inf")
  # the last columns is the pseudo-reference column
  pseudo_column <- ncol(log_data)
  # where to stop before the pseudo column
  before_pseduo <- pseudo_column - 1
  # find the ratio of the log data to the pseudo-reference e.g. the difference of log of the ratio
  ratios <- sweep(log_data[, 2:before_pseduo], 1, log_data[, pseudo_column], "-")
  # find the median of the ratios
  sample_medians <- apply(ratios, 2, stats::median)
  # convert the median to a scaling factor
  scaling_factors <- exp(sample_medians)
  # use scaling factors to scale the original data
  manually_normalized <- sweep(data, 2, scaling_factors, "/")
  return(manually_normalized)
}

#' @title uqua
#' @description Upper-quartile normalization with optional length correction.
#'
#' @param datos Numeric matrix of read counts.
#' @param long Numeric vector of feature lengths. If 1000, no correction is applied.
#' @param lc Length correction exponent. Default = 1.
#' @param k Replacement for zeros to avoid log issues. Default = 0.
#'
#' @return A normalized numeric data frame.
#' @export

uqua <- function(datos, long = 1000, lc = 1, k = 0) {
  sinceros <- function(datos, k) {
    datos0 <- datos
    if (is.null(k)) {
      mini0 <- min(datos[noceros(datos, num = FALSE, k = 0)])
      kc <- mini0 / 2
      datos0[datos0 == 0] <- kc
    } else {
      datos0[datos0 == 0] <- k
    }
    datos0
  }
  # lc: Length correction. Expression is divided by long^lc. lc can be any real number.
  L <- long^lc
  datos0 <- sinceros(datos, k)
  if (ncol(as.matrix(datos)) > 1) {
    sumatot <- rowSums(datos)
    supertot <- sum(sumatot)
    counts0 <- which(sumatot == 0)
    if (length(counts0) > 0) {
      datitos <- datos[-counts0, ]
    } else {
      datitos <- datos
    }
    q3 <- apply(datitos, 2, stats::quantile, probs = 0.75)
    d <- q3 * supertot / sum(q3)
    datos.norm <- (t(t(datos0) / d) * 10^9) / L
  } else {
    datos.norm <- datos0 / L
  }
  stats::na.omit(datos.norm)
}

#' @title tmm
#' @description Trimmed Mean of M-values (TMM) normalization with optional length correction.
#'
#' @param datos Numeric matrix of read counts.
#' @param long Numeric vector of feature lengths.
#' @param lc Length correction exponent.
#' @param k Replacement for zeros.
#' @param refColumn Column used as reference.
#' @param logratioTrim Trim amount for log-ratios.
#' @param sumTrim Trim amount for absolute levels.
#' @param doWeighting Whether to apply binomial weighting.
#' @param Acutoff Cutoff for "A" values before trimming.
#'
#' @return A normalized numeric data frame.
#' @export

tmm <- function(datos, long = 1000, lc = 1, k = 0, refColumn = NULL,
                logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE,
                Acutoff = -1e10) {
  sinceros <- function(datos, k) {
    datos0 <- datos
    if (is.null(k)) {
      mini0 <- min(datos[noceros(datos, num = FALSE, k = 0)])
      kc <- mini0 / 2
      datos0[datos0 == 0] <- kc
    } else {
      datos0[datos0 == 0] <- k
    }
    datos0
  }
  # lc: Length correction. Expression is divided by long^lc. lc can be any real number.
  L <- long^lc
  datos0 <- sinceros(datos, k)
  if (ncol(as.matrix(datos)) > 1) {
    fk <- edgeR::calcNormFactors(as.matrix(datos),
      method = "TMM", refColumn = refColumn,
      logratioTrim = logratioTrim, sumTrim = sumTrim,
      doWeighting = doWeighting, Acutoff = Acutoff
    )
    datos.norm <- (t(t(datos0) * fk) * 10^3) / L
  } else {
    datos.norm <- datos0 / L
  }
  stats::na.omit(datos.norm)
}

#' @title VarianceStabilized
#' @description Applies a variance stabilization transformation to a normalized matrix.
#'
#' @param norm_mat A numeric matrix to be variance-stabilized.
#'
#' @return A variance-stabilized matrix.
#' @export

VarianceStabilized <- function(norm_mat) {
  ncounts <- norm_mat
  xg <- sinh(seq(asinh(0), asinh(max(ncounts)), length.out = 1000))[-1]
  xim <- 1
  baseVarsAtGrid <- xg * xg^2 + xim * xg
  integrand <- 1 / sqrt(baseVarsAtGrid)
  splf <- stats::splinefun(
    asinh((xg[-1] + xg[-length(xg)]) / 2),
    cumsum(
      (xg[-1] - xg[-length(xg)]) *
        (integrand[-1] + integrand[-length(integrand)]) / 2
    )
  )
  h1 <- stats::quantile(rowMeans(ncounts), .95)
  h2 <- stats::quantile(rowMeans(ncounts), .999)
  eta <- (log2(h2) - log2(h1)) / (splf(asinh(h2)) - splf(asinh(h1)))
  xi <- log2(h1) - eta * splf(asinh(h1))
  tc <- sapply(colnames(norm_mat), function(clm) {
    eta * splf(asinh(ncounts[, clm])) + xi
  })
  rownames(tc) <- rownames(norm_mat)
  return(tc)
}

#' @title convert_to_numeric
#' @description Converts input to numeric format, preserving sparse matrix structure.
#'
#' @param x A matrix or data frame.
#'
#' @return A numeric matrix or sparse dgCMatrix.
#' @export

convert_to_numeric <- function(x) {
  is_dgCMatrix <- function(x) inherits(x, "dgCMatrix")
  if (is_dgCMatrix(x)) {
    return(Matrix::drop0(x))  # Remove explicit zeros, keep sparse
  } else {
    return(apply(x, 2, function(col) as.numeric(as.character(col))))
  }
}

#' @title load_sparse_matrix_optimized
#' @description Fast loader that returns a sparse dgCMatrix from a delimited text file.
#'
#' @param file_path Path to the input matrix file.
#' @param row_name_index Column index for rownames. Default = 1.
#' @param sep Field separator (default tab).
#'
#' @return A dgCMatrix with rownames set.
#' @export

load_sparse_matrix_optimized <- function(file_path, row_name_index = 1, sep = "\t") {
  cat("|| Reading matrix using data.table::fread...\n")
  df <- data.table::fread(file_path, sep = sep, data.table = FALSE)

  # Extract rownames and remove that column
  rownames_vec <- df[[row_name_index]]
  df <- df[, -row_name_index, drop = FALSE]

  # Convert all columns to numeric (safely)
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))

  # Convert to dense matrix (temporary)
  dense_mat <- as.matrix(df)
  rownames(dense_mat) <- rownames_vec

  # Convert to sparse dgCMatrix
  sparse_mat <- as(Matrix::Matrix(dense_mat, sparse = TRUE), "dgCMatrix")

  cat("|| Sparse matrix created (dgCMatrix), dim: ", paste(dim(sparse_mat), collapse = " x "), "\n")
  return(sparse_mat)
}

#' @title parse_input_matrix
#' @description Parses and converts a file path or R object into a sparse dgCMatrix.
#'
#' @param mat_path File path or matrix-like object.
#' @param row_name_index Column index for rownames (only used if `mat_path` is a file).
#'
#' @return A dgCMatrix.
#' @export

parse_input_matrix <- function(mat_path, row_name_index = 1) {
  is_dgCMatrix <- function(x) inherits(x, "dgCMatrix")
  if (is_dgCMatrix(mat_path)) {
    cat("|| Input is already a dgCMatrix\n")
    return(mat_path)
  }

  if (is.matrix(mat_path) || is.data.frame(mat_path)) {
    cat("|| Converting dense object to sparse dgCMatrix\n")
    mat_dense <- as.matrix(mat_path)
    return(as(Matrix::Matrix(mat_dense, sparse = TRUE), "dgCMatrix"))
  }

  if (is.character(mat_path) && file.exists(mat_path)) {
    # Handle RDS input
    if (grepl("\\.rds$", mat_path, ignore.case = TRUE)) {
      cat("|| Loading matrix from .rds file\n")
      mat <- readRDS(mat_path)
      return(parse_input_matrix(mat, row_name_index = row_name_index))  # re-parse
    }

    # Determine separator
    ext <- tools::file_ext(mat_path)
    sep <- switch(tolower(ext),
      csv = ",",
      tsv = "\t",
      txt = "\t",
      "\t"  # default fallback
    )
    cat("|| Loading matrix from file:", mat_path, "with sep =", shQuote(sep), "\n")
    return(load_sparse_matrix_optimized(mat_path, row_name_index = row_name_index, sep = sep))
  }

  stop("Unsupported input type for mat_path. Must be a dgCMatrix, dense matrix/data.frame, or valid file path (.tsv, .csv, .rds, .qs).")
}

#' @title parse_design
#' @description Parses a design matrix from a file or data.table.
#'
#' @param design_path File path or in-memory data.table.
#'
#' @return A validated design data.table with Sample_ID as rownames.
#' @export

parse_design <- function(design_path) {
  if (is.character(design_path) && file.exists(design_path)) {
    design <- data.table::fread(design_path)
  } else {
    design <- data.table::as.data.table(design_path)
  }
  if (!all(c("Sample_ID", "Sample_Condition", "Sample_Replicate") %in% colnames(design))) {
    stop("Design table must include Sample_ID, Sample_Condition, and Sample_Replicate columns.")
  }
  rownames(design) <- design$Sample_ID
  return(design)
}

#' @title Compute Log Ratios
#' @description Compute log(reference/sample) for a given pair of columns, efficiently for dgCMatrix.
#' @param reference Column name or integer index of the reference.
#' @param sample    Column name or integer index of the sample.
#' @param mat       Numeric matrix or sparse Matrix::dgCMatrix.
#' @return Numeric vector of finite log-ratios.
#' @export

compute_log_ratios <- function(reference, sample, mat) {
  # Fast column extraction for dense or dgCMatrix
  ref  <- as.numeric(mat[, reference, drop = FALSE])
  samp <- as.numeric(mat[, sample,    drop = FALSE])
  ratio <- log(ref / samp)
  ratio <- ratio[is.finite(ratio)]
  return(ratio)
}

#' @title Compute Scaling Factor
#' @description Computes the pseudo-coverage ratio between a reference and sample using model cutoffs to define valid log-ratio ranges.
#' @param mat Matrix or sparse matrix of values
#' @param reference_name Column name of the reference
#' @param sample_name Column name of the sample
#' @param model Model object with interval
#' @return Numeric scaling factor (sample to reference)
#' @export
compute_scaling_factor <- function(mat, reference_name, sample_name, model) {

  log_ratios <- log((mat[, reference_name]) / (mat[, sample_name]))
  in_bounds <- which(is.finite(log_ratios) & 
                      log_ratios >= model$interval["lb"] & 
                      log_ratios <= model$interval["ub"]
  )

  mm <- mat[in_bounds, c(sample_name, reference_name), drop = FALSE]

  pseudo_cov <- colSums(mm)
  rat <- pseudo_cov[reference_name] / pseudo_cov[sample_name]
  return(rat)
}

#' @title Scale Sparse Columns
#' @description Efficiently scale columns of a sparse dgCMatrix by a numeric vector.
#' @param mat A dgCMatrix to scale.
#' @param scalars A numeric vector of scaling factors, one per column.
#' @return A scaled dgCMatrix with each column multiplied by the corresponding scalar.
#' @export
scale_sparse_columns <- function(mat, scalars) {
  cat("|| Scaling sparse matrix columns with scalars", "\n")
  stopifnot(length(scalars) == ncol(mat))
  is_dgCMatrix <- function(x) inherits(x, "dgCMatrix")
  if (is_dgCMatrix(mat)) {
    mat@x <- mat@x * rep(scalars, diff(mat@p))
    return(mat)
  } else {
    return(t(t(mat) * scalars))
  }
}

#' @title Generate Pairwise Comparisons
#' @description Generate all pairwise combinations (both directions) with sample indices.
#' @param sample_ids Character vector of sample IDs
#' @return A data.frame with columns: samples, reference, sample_index, reference_index
generate_pairwise_comparisons <- function(sample_ids) {
  n <- length(sample_ids)
  
  # Get upper triangle indices (i < j)
  idx <- which(upper.tri(matrix(1, n, n)), arr.ind = TRUE)
  
  # Original direction (A/B)
  all_pairs <- data.frame(
    samples = sample_ids[idx[, 2]],
    reference = sample_ids[idx[, 1]],
    sample_index = idx[, 2],
    reference_index = idx[, 1],
    stringsAsFactors = FALSE
  )

  # Optional: sort by indices if needed
  all_pairs <- all_pairs[order(all_pairs$reference_index, all_pairs$sample_index), ]
  rownames(all_pairs) <- NULL
  
  return(all_pairs)
}

#' @title Reference-Reference Normalization
#' @description Computes scaling factor for a sample against a reference using pairwise log ratios.
#' @param task A task object containing 'reference' and 'samples'.
#' @param mat_ref A matrix containing reference data.
#' @param n_pop_reference Number of populations in the reference.
#' @param sigma_times Scaling factor for sigma in the model.
#' @param dist_family Distribution family for the model.
#' @param Norm_plot Logical indicating whether to return the model for plotting.
#' @return A list with 'sample', 'reference', 'scaling' factor, and optionally 'model'.
#' @export
fun_refref <- function(task, mat_ref, n_pop_reference, sigma_times, dist_family, Norm_plot = FALSE) {
  ref  <- task$reference
  samp <- task$samples

  X2 <- mat_ref[, c(ref, samp), drop = FALSE]
  keep <- Matrix::rowSums(X2) > 0L
  if (!any(keep)) return(NULL)
  X2 <- X2[keep, , drop = FALSE]

  r <- as.numeric(X2[, ref,  drop = FALSE])
  s <- as.numeric(X2[, samp, drop = FALSE])

  ratio <- log(r) - log(s)
  ratio <- ratio[is.finite(ratio)]
  if (!length(ratio)) return(NULL)

  fit <- try(pair_fit(list(
    ratio       = ratio,
    n_pop       = n_pop_reference,
    sigma_times = sigma_times,
    dist_family = dist_family,
    index       = 0L
  )), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit)) return(NULL)

  lr <- log(r) - log(s)
  ok <- is.finite(lr) & lr >= fit$interval["lb"] & lr <= fit$interval["ub"]
  if (!any(ok)) return(NULL)

  sc <- sum(r[ok]) / sum(s[ok])  # reference / sample
  list(sample = samp, reference = ref, scaling = sc,mean_g    = unname(fit$interval["mean_g"]),
       model  = if (isTRUE(Norm_plot)) fit else NULL)
}

#' @title Average Reference Normalization
#' @description Computes scaling factor for a sample against an average reference from multiple samples.
#' @param sname Sample name to normalize.
#' @param mat Matrix containing the data.
#' @param avg_ref_vec Numeric vector of average reference values.
#' @param n_pop Number of populations in the reference.
#' @param sigma_times Scaling factor for sigma in the model.
#' @param dist_family Distribution family for the model.
#' @param Norm_plot Logical indicating whether to return the model for plotting.
#' @return A list with 'sample', 'reference', 'scaling' factor, and optionally 'model'.
#' @export
fun_avgref <- function(sname, mat, avg_ref_vec, n_pop, sigma_times, dist_family, Norm_plot = FALSE) {
  s <- as.numeric(mat[, sname, drop = FALSE])
  r <- avg_ref_vec
  keep <- is.finite(r) & is.finite(s) & (r + s) > 0
  if (!any(keep)) return(NULL)

  ratio <- log(r[keep]) - log(s[keep])
  ratio <- ratio[is.finite(ratio)]
  if (!length(ratio)) return(NULL)

  fit <- try(pair_fit(list(
    ratio       = ratio,
    n_pop       = n_pop,
    sigma_times = sigma_times,
    dist_family = dist_family,
    index       = 0L
  )), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit)) return(NULL)

  lr <- log(r[keep]) - log(s[keep])
  ok <- is.finite(lr) & lr >= fit$interval["lb"] & lr <= fit$interval["ub"]
  if (!any(ok)) return(NULL)

  sc <- sum(r[keep][ok]) / sum(s[keep][ok])
  list(sample = sname, reference = "average_reference", scaling = sc,mean_g    = unname(fit$interval["mean_g"]),
       model  = if (isTRUE(Norm_plot)) fit else NULL)
}

#' @title Single Reference Normalization
#' @description Computes scaling factor for a sample against a single reference sample.
#' @param sname Sample name to normalize.
#' @param mat Matrix containing the data.
#' @param ref Reference sample name.
#' @param r_full Numeric vector of reference values.
#' @param keep_r Logical vector indicating which reference values to keep.
#' @param n_pop Number of populations in the reference.
#' @param sigma_times Scaling factor for sigma in the model.
#' @param dist_family Distribution family for the model.
#' @param Norm_plot Logical indicating whether to return the model for plotting.
#' @return A list with 'sample', 'reference', 'scaling' factor, and optionally 'model'.
#' @export
fun_single <- function(sname, mat, ref, r_full, keep_r, n_pop, sigma_times, dist_family, Norm_plot = FALSE) {
  s <- as.numeric(mat[, sname, drop = FALSE])
  keep <- keep_r & is.finite(s) & (r_full + s) > 0
  if (!any(keep)) return(NULL)

  ratio <- log(r_full[keep]) - log(s[keep])
  ratio <- ratio[is.finite(ratio)]
  if (!length(ratio)) return(NULL)

  fit <- try(pair_fit(list(
    ratio       = ratio,
    n_pop       = n_pop,
    sigma_times = sigma_times,
    dist_family = dist_family,
    index       = 0L
  )), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit)) return(NULL)

  lr <- log(r_full[keep]) - log(s[keep])
  ok <- is.finite(lr) & lr >= fit$interval["lb"] & lr <= fit$interval["ub"]
  if (!any(ok)) return(NULL)
  
  sc <- sum(r_full[keep][ok]) / sum(s[keep][ok])
  list(sample = sname, reference = ref, scaling = sc,mean_g    = unname(fit$interval["mean_g"]),
       model  = if (isTRUE(Norm_plot)) fit else NULL)
}

utils::globalVariables(c(
  ".",           # the dot symbol in data.table
  ":=",          # needed sometimes, though not always flagged
  "scaling",     # any column names used inside data.table expressions
  "sample",
  "reference",
  "samples",
  "av_scaling"
))
