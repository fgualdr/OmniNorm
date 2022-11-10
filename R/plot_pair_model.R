#' plot_pair_model Plots fitted distributions to the ratios.
#'
#' @param ll List composed by: "ratio" for the selected pair, "n_pop" i.e. n° of expected populations,
#'          "sigma_times" n° of sigmas to select invariant observations, "dist_family" distribution famly see mixsmsn, ...
#'
#' @returns Returns automatically pdf saved within "saving_path"/Norm_pair_plots/ of the fitting and QC
#' @export
#'
#' @import mixsmsn
#' @import grDevices
#' @import graphics 
#' @import utils
#'
#' @examples
#' \dontrun{
#' # Typical usage
#' ll <- list(
#'   "ratio" = ratio, "n_pop" = n_pop_reference,
#'   "sigma_times" = sigma_times, "dist_family" = dist_family, "index" = x
#' )
#' plot_pair_model(ll)
#' }
#'
plot_pair_model <- function(ll) {
  # Import generate functions:
  smsnmean <- function(mu, sigma2, shape) {
    xi <- mu
    omega <- sqrt(sigma2)
    alpha <- shape
    C <- sqrt(2 / pi)
    delta <- alpha / sqrt(1 + alpha^2)
    mean <- xi + omega * delta * C
    mean
  }
  dSN <- utils::getFromNamespace("dSN", "mixsmsn")

  model <- ll[["model"]]
  mat <- ll[["mat_append"]]
  sample_name <- ll[["sample_name"]]
  reference_name <- ll[["reference_name"]]
  path <- ll[["saving_path"]]
  n_pop <- ll[["n_pop"]]
  # Save_pair plots:
  if (!is.null(path)) {
    save <- paste0(path, "/Norm_pair_plots/")
    dir.create(save)
  } else {
    stop("ERROR: to plot the Normalisation pair plots a path to a general output folder is required")
  }

  ratio <- log(mat[, reference_name] / mat[, sample_name])
  # Only consider finite values
  ratio <- ratio[is.finite(ratio)]
  xx <- seq(min(ratio), max(ratio), (max(ratio) - min(ratio)) / (length(ratio) * 100))
  # We produce grey lines for all but the selected population:
  # To select the reference population we can two main criteria plus an additional:
  # 1) Small error
  # 2) comprise the largest population
  score <- model$score
  # Define colours e.g. reference population in "red"
  myColors <- rep("black", n_pop)
  myColors[which.max(score)] <- "red"
  samp_apirs <- paste0("Sample_", sample_name, "_Reference_", reference_name, ".pdf")
  nm <- paste0(save, "/", samp_apirs)
  pdf(nm)
  par(mar = c(2, 2, 2, 2))
  layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 2, 2, 3, 3), 4, 4, byrow = TRUE))
  # par(mfrow = c(1, 3), pty="s")
  # Model fitting plot
  mixsmsn::mix.hist(ratio, model, col.hist = "grey", border = "white", cex = 0.5, cex.lab = .5, cex.axis = .5, cex.main = .5, cex.sub = .5)
  mixsmsn::mix.lines(ratio, model, cex = 0.5, cex.lab = .5, cex.axis = .5, cex.main = .5, cex.sub = .5, col = "blue")
  dens <- matrix(NA_real_, ncol = n_pop, nrow = length(xx))
  dens.list <- vector("list", n_pop)
  means <- vector("numeric", n_pop)
  modes <- vector("numeric", n_pop)
  dens.max <- vector("numeric", n_pop)
  n_calls <- vector("numeric", n_pop)
  for (j in 1:n_pop) {
    dens <- dSN(xx, model$mu[j], model$sigma2[j], model$shape[j])
    dens.list[[j]] <- data.frame(x = xx, y = dens)
    means[j] <- smsnmean(model$mu[j], model$sigma2[j], model$shape[j])
    dens.max[j] <- which.max(dens.list[[j]]$y)
    modes[j] <- dens.list[[j]]$x[dens.max[j]]
    lines(xx, model$pii[j] * dens, type = "l", col = myColors[j], lwd = 1)
    lines(segments(means[j], 0, means[j], model$pii[j] * dSN(model$means[j], model$mu[j], model$sigma2[j], model$shape[j]), lty = 2, lwd = 1, col = myColors[j]))
    lines(segments(modes[j], 0, modes[j], model$pii[j] * dens[dens.max], col = myColors[j], lwd = 1))
  }
  # Scatter with selected population
  mmax <- as.numeric(c(mat[, reference_name], mat[, sample_name]))
  mmax <- max(mmax[is.finite(mmax)])
  mmin <- as.numeric(c(mat[, reference_name], mat[, sample_name]))
  mmin <- mmin[is.finite(mmin)]
  mmin <- min(mmin[mmin != 0])
  ##
  ww <- which(log(mat[, reference_name] / mat[, sample_name]) >= model$interval["lb"] &
    log(mat[, reference_name] / mat[, sample_name]) <= model$interval["ub"])
  # Ori mat
  plot(mat[, reference_name], mat[, sample_name], xlab = "Average reference", ylab = "Selected Sample", cex = 0.5, cex.lab = .5, cex.axis = .5, cex.main = .5, cex.sub = .5, log = "xy", xlim = c(mmin, mmax), ylim = c(mmin, mmax), pch = 20, asp = 1)
  points(mat[ww, reference_name], mat[ww, sample_name], log = "xy", xlim = c(mmin, mmax), ylim = c(mmin, mmax), cex = 0.5, col = "red")
  abline(coef = c(0, 1), col = "red")
  abline(v = 100)
  abline(h = 100)

  # Scaled mat
  plot(mat[, reference_name], mat[, sample_name] * (exp(model$interval["mean_g"])), xlab = "Average reference", ylab = "Selected Sample", cex = 0.5, cex.lab = .5, cex.axis = .5, cex.main = .5, cex.sub = .5, log = "xy", xlim = c(mmin, mmax), ylim = c(mmin, mmax), pch = 20, asp = 1)
  abline(coef = c(0, 1), col = "red")
  abline(v = 100)
  abline(h = 100)

  dev.off()
}
