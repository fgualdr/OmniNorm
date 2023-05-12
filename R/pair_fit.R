#' pair_fit execute pair wise fitting of Skewed normal on the freq. distribution of the ratios.
#'
#' @param ll List composed by: "ratio" for the selected pair, "n_pop" i.e. n° of expected populations,
#'          "sigma_times" n° of sigmas to select invariant observations, "dist_family" distribution famly see mixsmsn, ...
#'
#' @returns Return EM algorithm output for FM-SMSN distributions with additional fields: n_pop,means,modes,fre  q,freq_dist,sigma_norm
#' @export
#'
#' @import mixsmsn
#' @import fitdistrplus
#'
#' @examples
#' \dontrun{
#' # Typical usage
#' ll <- list(
#'   "ratio" = ratio, "n_pop" = n_pop_reference,
#'   "sigma_times" = sigma_times, "dist_family" = dist_family, "index" = x
#' )
#' model <- pair_fit(ll)
#' }
#'
pair_fit <- function(ll) {
  # Import/create dedicated functions
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
  d.mixedSN <- utils::getFromNamespace("d.mixedSN", "mixsmsn")

  ratio <- ll[["ratio"]]
  sigma_times <- ll[["sigma_times"]]
  family <- ll[["dist_family"]]
  n_pop <- ll[["n_pop"]]

  if (is.null(sigma_times)) {
    sigma_times <- 1
  }

  if (is.null(n_pop)) {
    # Discover best n° of peaks:
    model_search <- mixsmsn::smsn.search(ratio, 3,
      g.min = 1, g.max = 3,
      family = family, criteria = "bic",
      error = 0.00001, iter.max = 100,
      calc.im = FALSE, uni.Gama = FALSE, kmeans.param = NULL
    )
    n_pop <- length(unique(model_search$best.model$group))
  }

  # Compute the mixed Skewed model
  model <- mixsmsn::smsn.mix(ratio,
    nu = 3,
    g = n_pop,
    get.init = TRUE,
    criteria = TRUE,
    iter.max = 1000,
    error = 0.00001,
    group = TRUE,
    family = family,
    calc.im = FALSE
  )
  model$n_pop <- n_pop
  xx <- seq(min(ratio), max(ratio), (max(ratio) - min(ratio)) / (length(ratio) * 100))
  dens <- matrix(NA_real_, ncol = n_pop, nrow = length(xx))
  means <- vector("numeric", n_pop)
  modes <- vector("numeric", n_pop)
  dens.max <- vector("numeric", n_pop)
  n_calls <- vector("numeric", n_pop)
  sigmas <- vector("numeric", n_pop)
  for (j in 1:n_pop) {
    dens[, j] <- dSN(xx, model$mu[j], model$sigma2[j], model$shape[j])
    means[j] <- smsnmean(model$mu[j], model$sigma2[j], model$shape[j])
    dens.max[j] <- which.max(dens[, j])
    modes[j] <- xx[dens.max[j]]
    sk <- ratio[model$group == j]
    f <- fitdistrplus::fitdist(data = sk, distr = "norm", keepdata = FALSE, discrete = FALSE, fix.arg = list("mean" = modes[j]))
    sigmas[j] <- f$estimate["sd"]
  }
  global_dens <- d.mixedSN(xx, model$pii, model$mu, model$sigma2, model$shape)
  # Define the freq_sum as the mode of the global distribution
  freq_sum <- xx[which.max(global_dens)]
  # Decorate the model
  model$means <- means
  model$modes <- modes
  model$freq <- apply(dens, 2, max)
  model$freq_dist <- abs(modes - freq_sum)
  model$sigma_norm <- sigmas
  # Select the reference population based on score: i.e. 
  # highest frequency and smaller error:
  # consider that up to 3 populations are possible
  # We sum 3 scores:
  # 1 - we give higher score to the population with higher frequency
  # 2 - we give higher score to the population with smaller distance from the freq_sum which is the mode of the global distribution
  # if we have one population with freq_dist = 0, then we have a perfect fit and we give a score of 1
  # 3 - we give higher score to the population with smaller sigma_norm
  score1 <- (model$freq) / (max(model$freq)) # higher the better # min max(model$freq)
  score2 <- max(model$freq_dist) / model$freq_dist # smaller the better
  score2 <- score2 / max(score2)
  # if all freq_dist are 0, then we have a perfect fit
  if(length(model$freq_dis)==1){
    if (model$freq_dist == 0) {
      score2 <- 1
    }
  }
  score3 <- max(model$sigma_norm) / model$sigma_norm # smaller the better
  score3 <- score3 / max(score3)
  # We sum them all up 
  model$score <- score1 + score2 + score3
  # we fix to the mode to get the sd of the best fit normal
  interval <- c(
    model$means[which.max(model$score)],
    model$modes[which.max(model$score)],
    model$modes[which.max(model$score)],
    model$sigma_norm[which.max(model$score)],
    model$modes[which.max(model$score)] - (sigma_times * model$sigma_norm[which.max(model$score)]),
    model$modes[which.max(model$score)] + (sigma_times * model$sigma_norm[which.max(model$score)])
  )
  names(interval) <- c("mean_sn", "mode_sn", "mean_g", "sd", "lb", "ub")
  model$interval <- interval
  return(model)
}
