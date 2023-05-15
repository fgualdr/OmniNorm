#' RunNorm function to execute the whole normalisation flow.
#'
#' @param mat_path path to the matrix (n * m having features i.e. genes, proteins, chipseq peaks in the raws and individual samples in the columns)
#' @param design_path design matrix describing the samples. design matrix must be composed by:  Sample_ID, Sample_Condition and Sample_Replicate columns (layers in Sample_conditions will be defined by "_" separator)
#' @param fix_reference vector of Sample_ID names to be used as initial reference for normalisation. If random specified it select randomly 30% of the samples (To be improoved considering global coverage and intensity in order to avoid failed samples) (default NULL)
#' @param row_name_index column index where to find the rownames of the matrix (default 1)
#' @param n_pop n° of populations expected (default NULL)
#' @param n_pop_reference n° of expected populations in the initial within references normalisation (default NULL)
#' @param saving_path full path where results are saved default to the working directory
#' @param sigma_times number of sigmas defining lower and upper bounds to include features (i.e. genes, peaks, proteins ,..) considered to be unchanged relative to the mean distribution
#' @param dist_family distribution family see package mixsmsn (default "Skew.normal")
#' @param Norm_plot Whether or not to save individual plots during normalisation (defulat TRUE)
#' @param Save_results whether to save individual plots and results (default TRUE)
#' @param BiocParam BiocParallel param object to be specified as : "param <- BiocParallel::MulticoreParam(workers=2,progressbar = TRUE)" (default NULL)
#'
#' @returns Return a list which includes the fields: scaling_factors; ori_mat; norm_mat and model_list
#' @export
#'
#' @import BiocParallel
#' @importFrom dplyr mutate_all
#' @importFrom data.table as.data.table
#' @import tools
#' @importFrom rlang .data
#' @import methods
#'
#' @examples
#' \dontrun{
#' # Typical usage
#' param <- BiocParallel::MulticoreParam(workers = 2, progressbar = TRUE)
#'
#' Result <- RunNorm(mat,
#'   deg_design,
#'   fix_reference = "random",
#'   row_name_index = 1,
#'   saving_path = out_dir,
#'   n_pop = 1, BiocParam = param
#' )
#' }
#'
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
                    BiocParam = NULL) {
  if (is.null(saving_path) & Save_results == TRUE) {
    stop("ERROR: to plot the Normalisation pair plots, a path to a general output folder is required")
  }

  if (n_pop_reference > 3 | is.null(n_pop_reference)) {
    stop("n_pop_reference must be between 1 and 3")
  }

  if (n_pop > 3 | is.null(n_pop)) {
    stop("n_pop must be between 1 and 3")
  }

  # Import the tables
  cat("|| Import matrix\n")
  # If character is a path, if data.frame just assign
  if (is.character(mat_path)) {
    mat <- read.delim(mat_path, sep = "\t", stringsAsFactors = FALSE, row.names = row_name_index)
  } else {
    mat <- mat_path
  }
  cat("|| DONE\n")
  # Import the Design Table
  cat("|| Import Design table\n")
  if (is.character(design_path)) {
    design <- read.delim(design_path, sep = "\t", stringsAsFactors = FALSE)
  } else {
    design <- design_path
  }
  rownames(design) <- design$Sample_ID
  w <- which(design$Sample_ID %in% colnames(mat))
  mat <- mat[, design$Sample_ID]
  mat <- dplyr::mutate_all(mat, function(x) as.numeric(as.character(x)))
  design <- design[w, ]
  cat("|| DONE\n")
  if (dim(design)[1] == 0) {
    stop(
      "No samples in design table is included in the data matrix\n
    Check design and matrix sample format\n
    Design table should include: Sample_ID, Sample_Condition, Sample_Replicate\n
    Sample_ID in the design table and columns names in data matrix must match\n
    see example tables and vignette."
    )
  }
  # Add filter in case user decide on a fixed reference or set of reference -
  if (!is.null(fix_reference)) {
    # Select the indexes matching the reference (which should be a vector of sample IDs)
    cat("Selecting reference samples for normalisation:\n")
    wrev <- which(fix_reference %in% rownames(design))
    if (length(wrev) != length(fix_reference)) {
      if (fix_reference == "random") {
        set.seed(123456)
        tot_samp <- dim(design)[1]
        # By random we always select 30% with a maximum of 50 samples - is samples are less than 50 we will by default select all samples:
        n <- round(tot_samp / 100 * 30, 0)
        if (n > 50) {
          n <- 50
        } else {
          if (tot_samp <= 50) {
            n <- tot_samp
          }
        }
        r_index <- sample(1:tot_samp, n, replace = FALSE)
        fix_reference <- design$Sample_ID[r_index]
        if (tot_samp <= 50) {
          fix_reference <- design$Sample_ID
        }
        cat("Selected ", length(fix_reference), " random reference samples\n")
      } else {
        stop("ERROR: fix_reference must be a vector of Samples_ID all part of the Sample_ID column in the design table")
      }
    }else{
      fix_reference = fix_reference[wrev]
    }
  } else {
    cat("Warning : no reference set was specified\n
      For large number of sample pair-wise computation of Skew-Normal distribution is computational intensive.\n
      Make sure to set up BiocParam accordingly and control memory availability.
      For 100 samples ~5000 fittings are goin to be computed.\n")
    fix_reference <- rownames(design)
  }
  # Consider that if the reference is only one we skip the reference normalization and we go straigth to the by_sample norm:
  if(length(fix_reference) > 1){
    # if multiple reference samples are selected we need to compute the pair-wise skew-normal distribution across them
    design_ref <- as.data.frame(as.matrix(design[fix_reference, ]), stringsAsFactors = FALSE)
    mat_sample <- as.data.frame(as.matrix(mat[, rownames(design_ref)]), stringsAsFactors = FALSE)
    # pair_wise data.frame set up
    sn <- 1:length(design_ref$Sample_ID)
    names(sn) <- design_ref$Sample_ID
    pairs <- t(combn(design_ref$Sample_ID, 2, FUN = NULL, simplify = TRUE))
    pairs <- as.data.frame(rbind(pairs, pairs[, 2:1]))
    colnames(pairs) <- c("samples", "reference")
    pairs$sample_index <- sn[pairs$samples]
    pairs$reference_index <- sn[pairs$reference]
    pairs <- pairs[order(pairs$reference_index), ]
    pairs <- pairs[order(pairs$sample_index), ]
    rownames(pairs) <- NULL
    # Select upper_triangle
    pairs <- pairs[which(pairs$reference_index > pairs$sample_index), ]
    cat(" || Compute Mixture of Skew-Normal Distributions ||\n")
    if (!is.null(BiocParam)) {
      model_list <- list()
      rat_l <- lapply(1:nrow(pairs), function(x) {
        ratio <- log(mat_sample[, as.character(pairs[x, "reference"])] / mat_sample[, as.character(pairs[x, "samples"])])
        # Only consider finite values
        ratio <- ratio[is.finite(ratio)]
        if (length(ratio) == 0) {
          cat("ERROR!!!")
        } # Check!
        return(list("ratio" = ratio, "n_pop" = n_pop_reference, "sigma_times" = sigma_times, "dist_family" = dist_family, "index" = x))
      })
      model_list <- tryCatch(
        {
          BiocParallel::bplapply(rat_l, pair_fit, BPPARAM = BiocParam)
        },
        error = identity
      )
    } else {
      model_list <- list()
      model_list <- lapply(1:nrow(pairs), function(x) {
        # get log2-ratio of the selected pair
        ratio <- log(mat_sample[, as.character(pairs[x, "reference"])] / mat_sample[, as.character(pairs[x, "samples"])])
        # Only consider finite values
        ratio <- ratio[is.finite(ratio)]
        ll <- list("ratio" = ratio, "n_pop" = n_pop_reference, "sigma_times" = sigma_times, "dist_family" = dist_family, "index" = x)
        model <- pair_fit(ll)
        return(model)
      })
    }
    cat("|| Computing models Done\n")
    # Scaling and create Av. Sample, Append to mat[,rownames(design)]
    pairs$scaling <- NA
    scaling_ll <- lapply(1:length(model_list), function(x) {
      reference_name <- as.character(pairs[x, "reference"])
      sample_name <- as.character(pairs[x, "samples"])
      model <- model_list[[x]]
      ww <- which(log(as.numeric(as.character(mat[, reference_name])) / as.numeric(as.character(mat[, sample_name]))) >= model$interval["lb"] &
        log(as.numeric(as.character(mat[, reference_name])) / as.numeric(as.character(mat[, sample_name]))) <= model$interval["ub"])
      mm <- mat[ww, c(sample_name, reference_name)]
      pseudo_cov <- colSums(mm)
      rat <- pseudo_cov[2] / pseudo_cov[1]
      names(rat) <- x
      return(rat)
    })
    pairs$scaling <- unlist(scaling_ll)
    pairs_add <- pairs[, c(2, 1, 4, 3, 5)]
    colnames(pairs_add) <- colnames(pairs)
    pairs_add$scaling <- 1 / pairs$scaling
    pairs_add <- pairs_add[order(pairs_add$reference_index), ]
    pairs_add <- pairs_add[order(pairs_add$sample_index), ]
    self <- cbind(names(sn), names(sn), 1:length(sn), 1:length(sn), rep(1, length(sn)))
    colnames(self) <- colnames(pairs_add)
    pairs_add <- rbind(pairs_add, self)
    pairscombo <- rbind(pairs, pairs_add)
    pairscombo <- pairscombo[order(pairscombo$reference_index), ]
    pairscombo <- pairscombo[order(pairscombo$sample_index), ]
    pairscombo <- data.table::as.data.table(pairscombo)
    pairscombo$scaling <- as.numeric(as.character(pairscombo$scaling))
    pairscombo <- split(pairscombo, pairscombo$samples)
    pairscombo <- lapply(pairscombo, function(x) {
      mean(x$scaling)
    })
    pairscombo <- data.table::as.data.table(cbind(names(pairscombo), unlist(pairscombo)))
    colnames(pairscombo) <- c("samples", "av_scaling")
    mm <- mat[, as.character(pairscombo$samples)]
    mat_append <- mat[, rownames(design)]
    # Add the average reference to the mat_append
    mat_append$average_reference <- rowMeans(as.matrix(mm) %*% diag(pairscombo$av_scaling))
    # Compute Parallel model fit pairwise every sample to the mean_reference
    # pair_wise data.frame set up
    sn <- 1:length(design$Sample_ID)
    names(sn) <- design$Sample_ID
    pairs <- as.data.frame(cbind(design$Sample_ID, rep("average_reference", length(sn))))
    colnames(pairs) <- c("samples", "reference")
    rownames(pairs) <- NULL
    cat("|| Compute Mixture of Skew-Normal Distributions between samples and Average reference ||\n")

  } else {

    # if fix_reference == 1 we want to scale it to one reference sample only therefore we set average_reference to the fix_reference sample
    mat_append <- mat[, rownames(design)]
    # Compute Parallel model fit pairwise every sample to the mean_reference
    # pair_wise data.frame set up
    samp_left = design$Sample_ID[!grepl(fix_reference,design$Sample_ID)]
    sn <- 1:length(samp_left)
    names(sn) <- samp_left
    pairs <- as.data.frame(cbind(samp_left, rep(fix_reference, length(sn))))
    colnames(pairs) <- c("samples", "reference")
    rownames(pairs) <- NULL
    cat("|| Compute Mixture of Skew-Normal Distributions between samples and Average reference ||\n")

  }

  if (!is.null(BiocParam)) {
    model_list <- list()
    rat_l <- lapply(1:nrow(pairs), function(x) {
      ratio <- log(mat_append[, as.character(pairs[x, "reference"])] / mat_append[, as.character(pairs[x, "samples"])])
      # Only consider finite values
      ratio <- ratio[is.finite(ratio)]
      return(list("ratio" = ratio, "n_pop" = n_pop, "sigma_times" = sigma_times, "dist_family" = dist_family, "index" = x))
    })
    model_list <- BiocParallel::bplapply(rat_l, pair_fit, BPPARAM = BiocParam)
  } else {
    model_list <- list()
    model_list <- lapply(1:nrow(pairs), function(x) {
      cat(x, "\n")
      # get log2-ratio of the selected pair
      ratio <- log(mat_append[, as.character(pairs[x, "reference"])] / mat_append[, as.character(pairs[x, "samples"])])
      # Only consider finite values
      ratio <- ratio[is.finite(ratio)]
      ll <- list("ratio" = ratio, "n_pop" = n_pop, "sigma_times" = sigma_times, "dist_family" = dist_family, "index" = x)
      model <- pair_fit(ll)
      return(model)
    })
  }

  cat("|| Done.\n")

  pairs$scaling <- NA
  pairs$mean_g <- NA

  scaling_ll <- lapply(1:length(model_list), function(x) {
    reference_name <- as.character(pairs[x, "reference"])
    sample_name <- as.character(pairs[x, "samples"])
    model <- model_list[[x]]
    ww <- which(log(mat_append[, reference_name] / mat_append[, sample_name]) >= model$interval["lb"] &
      log(mat_append[, reference_name] / mat_append[, sample_name]) <= model$interval["ub"])
    pseudo_cov <- colSums(mat_append[ww, c(sample_name, reference_name)])
    rat <- pseudo_cov[2] / pseudo_cov[1]
    names(rat) <- x
    return(rat)
  })

  mean_g <- lapply(1:length(model_list), function(x) {
    reference_name <- as.character(pairs[x, "reference"])
    sample_name <- as.character(pairs[x, "samples"])
    model <- model_list[[x]]
    return(exp(model$interval["mean_g"]))
  })

  pairs$scaling <- unlist(scaling_ll)
  pairs$mean_g <- unlist(mean_g)

  # Plot the normalisation graphs relative to the average_reference
  cat("|| Plot pairs with average reference\n")
  if (Norm_plot == TRUE) {
    plot_l <- lapply(1:length(model_list), function(x) {

      reference_name <- as.character(pairs[x, "reference"])
      sample_name <- as.character(pairs[x, "samples"])
      model <- model_list[[x]]
      ll <- list("model" = model, "mat_append" = mat_append, "sample_name" = sample_name, "reference_name" = reference_name, "saving_path" = saving_path, "n_pop" = model$n_pop)
      return(ll)
    })
    if (!is.null(BiocParam)) {
      BiocParallel::bplapply(plot_l, plot_pair_model, BPPARAM = BiocParam)
    } else {
      lapply(plot_l, plot_pair_model)
    }
  }

  # if fix_reference == 1 we add "1" to the fix_reference sample
  if(length(fix_reference) == 1){
    pairs = rbind(pairs, c(fix_reference, fix_reference, 1, 1))
  }
  
  # Save the scaling factor to file "filling in the design table":
  norm_mat <- mat[, rownames(design)]
  norm_mat <- sweep(norm_mat[, pairs$samples], 2, as.numeric(as.character(pairs[, "scaling"])), "*")
  # Save the normalised table
  design_scaling <- design
  design_scaling$scaling <- pairs$scaling
  design_scaling$mean_g <- pairs$mean_g

  if (Save_results == TRUE) {
    if (is.character(mat_path)) {
      file_name <- tools::file_path_sans_ext(mat_path)
      file_name <- basename(file_name)
    } else {
      file_name <- "matrix"
    }
    write.table(design_scaling, file = paste0(saving_path, "/Normalisation_Parameters.txt"), sep = "\t", quote = FALSE, col.names = NA)
    write.table(norm_mat, file = paste0(saving_path, "/", file_name, "_normalized.txt"), sep = "\t", quote = FALSE, col.names = NA)
  }

  Result <- list(scaling_factors = design_scaling, ori_mat = mat[, rownames(design)], norm_mat = norm_mat, model_list = model_list)
  # Save to oath the scaled tables as "Norm_table"
  class(Result) <- "GMSN"
  return(Result)
}
