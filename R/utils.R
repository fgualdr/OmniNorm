#' noceros To count number of non-zero elements
#'
#' @param x vector to check zeros
#' @param num if vector is numeric TRUE FALSE
#' @param k default is 0
#'
#' @returns return non zero elements
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

#' quantile_normalisation
#'
#' @param df data frame to apply quantile norm
#'
#' @returns a data frame  normalized
#' @export
#'

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

#' median_of_ratios_normalization
#'
#' @param data data frame to apply quantile norm
#'
#' @importFrom dplyr mutate filter
#' @importFrom tibble rownames_to_column
#' @importFrom stats median
#'
#' @returns a data frame  normalized
#' @export
#'

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

#' Upper quartile normalisation
#'
#' @param datos Matrix containing the read counts for each sample.
#' @param long Numeric vector containing the length of the features. If long == 1000, no length correction is applied (no matter the value of parameter lc).
#' @param lc Correction factor for length normalization. This correction is done by dividing the counts vector by (length/1000)^lc.
#'            If lc = 0, no length correction is applied. By default, lc = 1 for RPKM and lc = 0 for the other methods.
#' @param k Counts equal to 0 are changed to k in order to avoid indeterminations when applying logarithms, for instance. By default, k = 0.
#'
#' @importFrom stats quantile na.omit
#'
#' @returns a data frame  normalized
#' @export
#'

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

## TMM: Trimmed Mean of M values normalization (Robinson & Oshlack, 2010)
#' Upper quartile normalisation
#'
#' @param datos Matrix containing the read counts for each sample.
#' @param long Numeric vector containing the length of the features. If long == 1000, no length correction is applied (no matter the value of parameter lc).
#' @param lc Correction factor for length normalization. This correction is done by dividing the counts vector by (length/1000)^lc.
#'            If lc = 0, no length correction is applied. By default, lc = 1 for RPKM and lc = 0 for the other methods.
#' @param k Counts equal to 0 are changed to k in order to avoid indeterminations when applying logarithms, for instance. By default, k = 0.
#' @param refColumn	Column to use as reference (only needed for tmm function).
#' @param logratioTrim	Amount of trim to use on log-ratios ("M" values) (only needed for tmm function).
#' @param sumTrim	Amount of trim to use on the combined absolute levels ("A" values) (only needed for tmm function).
#' @param doWeighting	Logical, whether to compute (asymptotic binomial precision) weights (only needed for tmm function).
#' @param Acutoff	Cutoff on "A" values to use before trimming (only needed for tmm function).
#'
#' @importFrom edgeR calcNormFactors
#' @importFrom stats na.omit
#'
#' @returns a data frame  normalized
#' @export
#'

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

#' Variance stabilisation
#'
#' @param norm_mat data frame to apply Variance stabilisation
#'
#' @importFrom stats quantile splinefun
#'
#' @returns Variance stabilized data frame
#' @export
#'
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
