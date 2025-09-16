#' @title plot_pair_model
#' @description Plots the fitted distributions and saves diagnostic PNGs of each normalization pair.
#'
#' This function expects a list with:
#' - model: a compact GMSNfit object from pair_fit()
#' - mat_append: a matrix or dgCMatrix with at least sample_name and reference_name columns
#' - sample_name: character string
#' - reference_name: character string
#' - saving_path: character string
#'
#' The function will save a PNG under saving_path/Norm_pair_plots/.
#'
#' @param ll A list containing ratio information and parameters, including: "ratio", "n_pop", "sigma_times", "dist_family", and "index".
#'
#' @return Invisibly returns NULL. Plots are saved to disk under saving_path/Norm_pair_plots/.
#'
#' @export
#'
#' @import mixsmsn
#' @import grDevices
#' @import graphics
#' @import utils
#'

plot_pair_model <- function(ll) {
  # Sparse/dense-safe column extractor

  dSN <- utils::getFromNamespace("dSN", "mixsmsn")

  model          <- ll[["model"]]        # GMSNfit (compact list)
  mat            <- ll[["mat"]]   # 2-column slice is ideal
  sample_name    <- ll[["sample_name"]]
  reference_name <- ll[["reference_name"]]
  path           <- ll[["saving_path"]]
  n_pop          <- ll[["n_pop"]]

  if (is.null(path)) stop("saving_path is required")
  outdir <- file.path(path, "Norm_pair_plots")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # Pull columns safely
  ref  <- as.numeric(mat[rownames(mat), reference_name, drop = FALSE])
  samp <- as.numeric(mat[rownames(mat), sample_name, drop = FALSE])

  ratio <- log(ref / samp)
  ratio <- ratio[is.finite(ratio)]

  if (!length(ratio)) {
    warning(sprintf("plot_pair_model: no finite ratios for %s vs %s", sample_name, reference_name))
    return(invisible(NULL))
  }

  # X-grid (use model hint if present)
  if (!is.null(model$ratio_range)) {
    x_min <- model$ratio_range["min"]; x_max <- model$ratio_range["max"]
    if (!is.finite(x_min) || !is.finite(x_max) || x_min >= x_max) {
      x_min <- min(ratio); x_max <- max(ratio)
    }
  } else { x_min <- min(ratio); x_max <- max(ratio) }
  xx <- seq(x_min, x_max, length.out = 1000L)

  # Highlight best component
  best <- which.max(model$score)
  cols <- rep("black", n_pop)
  if (length(best)) cols[best] <- "red"

  # scatter bounds
  both <- c(ref, samp)
  mmax <- max(both[is.finite(both)], na.rm = TRUE)
  nz   <- both[is.finite(both) & both != 0]
  mmin <- if (length(nz)) min(nz, na.rm = TRUE) else 1

  lb <- unname(model$interval["lb"]); ub <- unname(model$interval["ub"])
  ww <- which(is.finite(log(ref / samp)) & log(ref / samp) >= lb & log(ref / samp) <= ub)

  fname <- sprintf("Sample_%s_Reference_%s.png", sample_name, reference_name)
  fpath <- file.path(outdir, fname)

  png(fpath, width = 2000, height = 2000, res = 600, type = "cairo-png")
    on.exit(grDevices::dev.off(), add = TRUE)

    # smaller margins + smaller text
    par(mar = c(4.2, 4.2, 3, 2), 
        cex.axis = 0.7,   # shrink axis numbers
        cex.lab  = 0.8,   # shrink axis labels
        cex.main = 0.9)   # shrink titles

    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 
                    2, 2, 3, 3, 2, 2, 3, 3), 4, 4, byrow = TRUE))

    ## Panel 1: histogram + mixture/components (density scale)
    h <- hist(ratio, breaks = 100, plot = FALSE)
    plot(h, freq = FALSE, col = "grey85", border = FALSE,
        xlab = "log(reference / sample)", 
        ylab = "Density",
        main = paste0("Mixture fit: ", n_pop, " populations"),
        xlim = c(x_min, x_max), 
        ylim = c(0, max(h$density, na.rm = TRUE) * 1.2),
        axes = TRUE)

    # components with thinner lines
    for (j in seq_len(n_pop)) {
      comp <- model$pii[j] * dSN(xx, model$mu[j], model$sigma2[j], model$shape[j])
      lines(xx, comp, col = cols[j], lwd = if (j == best) 1.2 else 0.8)
    }
    mix <- rowSums(sapply(seq_len(n_pop), function(j)
      model$pii[j] * dSN(xx, model$mu[j], model$sigma2[j], model$shape[j])))
    lines(xx, mix, col = "blue", lwd = 1)
    abline(v = c(lb, ub), col = "red", lty = 2, lwd = 0.8)

    ## Panel 2: original scatter
    plot(ref, samp, log = "xy", pch = 20, cex = 0.2,
        xlab = "Reference", ylab = "Sample",
        xlim = c(mmin, mmax), ylim = c(mmin, mmax), main = "Original")
    if (length(ww)) points(ref[ww], samp[ww], pch = 20, cex = 0.2, col = "red")
    abline(a = 0, b = 1, col = "red", lwd = 0.8)

    ## Panel 3: scaled scatter using exp(mean_g)
    scale_g <- exp(unname(model$interval["mean_g"]))
    plot(ref, samp * scale_g, log = "xy", pch = 20, cex = 0.2,
        xlab = "Reference", ylab = "Sample * exp(mean_g)",
        xlim = c(mmin, mmax), ylim = c(mmin, mmax), main = "Scaled")
    abline(a = 0, b = 1, col = "red", lwd = 0.8)

  invisible(fpath)
}
