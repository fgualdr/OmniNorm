#' @title pair_fit
#' @description Executes pairwise fitting of skewed normal mixtures on log-ratio distributions.
#'
#' @param ll A list with the following elements:
#'   - `ratio`: numeric vector of log-ratios
#'   - `n_pop`: (optional) integer, expected number of populations (default = auto-detect)
#'   - `sigma_times`: numeric, multiplier for standard deviation to define invariant region
#'   - `dist_family`: character string, distribution family for `mixsmsn`
#'   - `index`: character or numeric identifier (used for tracking/logging)
#' @param do_gc whether to clean memory via gc
#'
#' @return A list containing the fitted model parameters, including `n_pop`, `means`, `modes`, `freq`, `freq_dist`, `sigma_norm`, `score`, and `interval`.
#'
#' @export
#'
#' @import mixsmsn
#' @import fitdistrplus
#'

pair_fit <- function(ll,do_gc = FALSE) {
  # helpers
  smsnmean <- function(mu, sigma2, shape) {
    xi <- mu; omega <- sqrt(sigma2); alpha <- shape
    xi + omega * (alpha / sqrt(1 + alpha^2)) * sqrt(2/pi)
  }
  dSN       <- utils::getFromNamespace("dSN", "mixsmsn")
  d.mixedSN <- utils::getFromNamespace("d.mixedSN", "mixsmsn")

  # unpack
  ratio       <- ll[["ratio"]]
  sigma_times <- ll[["sigma_times"]]; if (is.null(sigma_times)) sigma_times <- 1
  fam         <- ll[["dist_family"]]
  n_pop       <- ll[["n_pop"]]
  index       <- ll[["index"]]

  # clean
  ratio <- ratio[is.finite(ratio)]
  if (length(ratio) < 100 || stats::sd(ratio, na.rm = TRUE) == 0) {
    warning(sprintf("pair_fit[%s]: invalid ratio; forcing n_pop = 1", index))
    n_pop <- 1L
  }

  # infer n_pop if needed
  fallback_reason <- NULL
  if (is.null(n_pop)) {
    ms <- tryCatch({
      mixsmsn::smsn.search(
        ratio, 3, 
        g.min = 1, g.max = 3,
        family = fam, criteria = "bic",
        error = 0.00001, iter.max = 500,
        calc.im = FALSE, uni.Gama = FALSE, kmeans.param = NULL
      )
    }, error = function(e) {
      sprintf("pair_fit[%s]: smsn.search failed - %s", index, conditionMessage(e))
    })
    if (is.character(ms) || is.null(ms) || is.null(ms$best.model$group)) {
      fallback_reason <- if (is.character(ms)) ms else
        sprintf("pair_fit[%s]: smsn.search returned no groups", index)
      warning(paste0(fallback_reason, " - fallback to n_pop = 1"))
      n_pop <- 1L
    } else {
      n_pop <- max(1L, length(unique(ms$best.model$group)))
    }
  }

  # fit mixture
  fit <- tryCatch({
    mixsmsn::smsn.mix(
      ratio,
      nu = 3,
      g = n_pop,
      get.init = TRUE,
      criteria = FALSE,
      iter.max = 500,
      error = 0.00001,
      group = TRUE,
      family = fam,
      calc.im = FALSE
    )
  }, error = function(e) {
    warning(sprintf("pair_fit[%s]: smsn.mix failed - %s", index, conditionMessage(e)))
    NULL
  })
  if (is.null(fit)) return(NULL)
  if (any(!is.finite(fit$mu)) || any(!is.finite(fit$sigma2)) || any(!is.finite(fit$shape))) {
    warning(sprintf("pair_fit[%s]: non-finite parameters", index))
    return(NULL)
  }

  # build a small grid (no need to store it)
  grid_len <- min(2000L, max(1000L, length(ratio)))
  xx_min <- min(ratio, na.rm = TRUE)
  xx_max <- max(ratio, na.rm = TRUE)
  xx <- seq(xx_min, xx_max, length.out = grid_len)

  # per-component summaries
  means  <- rep(NA_real_, n_pop)
  modes  <- rep(NA_real_, n_pop)
  sigmas <- rep(NA_real_, n_pop)
  freq   <- rep(NA_real_, n_pop)  # max density per comp

  for (j in seq_len(n_pop)) {
    sk <- ratio[fit$group == j]
    if (length(sk) < 3 || stats::sd(sk, na.rm = TRUE) == 0 || any(!is.finite(sk))) {
      warning(sprintf("pair_fit[%s]: component %d skipped (flat/invalid)", index, j))
      next
    }
    dens <- tryCatch(
      dSN(xx, fit$mu[j], fit$sigma2[j], fit$shape[j]),
      error = function(e) rep(NA_real_, grid_len)
    )
    freq[j]  <- suppressWarnings(max(dens, na.rm = TRUE))
    means[j] <- smsnmean(fit$mu[j], fit$sigma2[j], fit$shape[j])
    modes[j] <- xx[which.max(dens)]

    fa <- tryCatch({
      fitdistrplus::fitdist(sk, "norm", keepdata = FALSE, fix.arg = list(mean = modes[j]))
    }, error = function(e) NULL)
    if (!is.null(fa)) sigmas[j] <- unname(fa$estimate["sd"])
  }

  valid <- which(is.finite(means) & is.finite(modes) & is.finite(sigmas) & is.finite(freq))
  if (!length(valid)) {
    warning(sprintf("pair_fit[%s]: all components invalid", index))
    return(NULL)
  }

  # compute global density for freq_dist (no need to store dens)
  global_dens <- d.mixedSN(xx, fit$pii[valid], fit$mu[valid], fit$sigma2[valid], fit$shape[valid])
  freq_sum <- xx[which.max(global_dens)]
  freq_dist <- abs(modes[valid] - freq_sum)

  # scores
  score1 <- (freq[valid]) / max(freq[valid])
  score2 <- ifelse(freq_dist == 0, 1, max(freq_dist) / freq_dist); score2 <- score2 / max(score2)
  score3 <- max(sigmas[valid]) / sigmas[valid];                  score3 <- score3 / max(score3)
  score  <- score1 + score2 + score3
  best   <- which.max(score)

  interval <- c(
    mean_sn = means[valid][best],
    mode_sn = modes[valid][best],
    mean_g  = modes[valid][best],
    sd      = sigmas[valid][best],
    lb      = modes[valid][best] - (sigma_times * sigmas[valid][best]),
    ub      = modes[valid][best] + (sigma_times * sigmas[valid][best])
  )

  # COMPACT return
  out <- list(
    # params sufficient to redraw curves
    mu     = unname(fit$mu[valid]),
    sigma2 = unname(fit$sigma2[valid]),
    shape  = unname(fit$shape[valid]),
    pii    = unname(fit$pii[valid]),
    n_pop  = length(valid),

    # diagnostics
    means      = unname(means[valid]),
    modes      = unname(modes[valid]),
    sigma_norm = unname(sigmas[valid]),
    freq       = unname(freq[valid]),
    freq_dist  = unname(freq_dist),
    score      = unname(score),

    # scaling interval
    interval = interval,

    # logging
    search_fallback = !is.null(fallback_reason),
    fallback_reason = fallback_reason,

    # tiny hint for plotting grid
    ratio_range = c(min = xx_min, max = xx_max)
  )

  # aggressively drop big locals before returning (helps peak mem in long chains)
  rm(ratio, fit, xx, means, modes, sigmas, freq, freq_dist, score, score1, score2, score3, valid)
  if (do_gc) gc(FALSE)

  return(out)
}
