#' Bootstrap confidence interval for Spearman's rho
#'
#' Nonparametric BCa bootstrap: resamples paired (x, y) observations and
#' recomputes Spearman's rho on each resample. No normality or asymptotic
#' variance-formula assumptions, unlike a Fisher z approach.
#'
#' @param data data.frame containing the x and y columns
#' @param x,y quoted names of the columns in `data` to correlate, e.g. "MS1"
#' @param conf.level confidence level, default 0.95
#' @param R number of bootstrap replicates, default 2000
#' @return data.frame with r, n, se (bootstrap SD), ci.lower, ci.upper, p.value
#'   (p.value is from the standard asymptotic cor.test, not the bootstrap itself —
#'   a bootstrap resample preserves the observed association, so it can't be used
#'   to test rho = 0; the CI and the p-value are deliberately from different methods)
spearman_boot_ci <- function(data, x, y, conf.level = 0.95, R = 2000) {
  if (!all(c(x, y) %in% names(data))) {
    stop("x and y must be column names present in data")
  }
  
  dat <- data.frame(x = data[[x]], y = data[[y]])
  if (any(!complete.cases(dat))) {
    stop("x and y contain NA values; remove or impute them before bootstrapping")
  }
  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("package 'boot' is required: install.packages('boot')")
  }
  
  stat <- function(d, i) cor(d$x[i], d$y[i], method = "spearman")
  
  boot_out <- boot::boot(dat, stat, R = R)
  ci <- boot::boot.ci(boot_out, conf = conf.level, type = "bca")
  bca <- ci$bca[4:5]  # lower, upper
  
  p_value <- suppressWarnings(
    cor.test(dat$x, dat$y, method = "spearman")$p.value
  )
  
  data.frame(r = boot_out$t0, n = nrow(dat),
             se = sd(boot_out$t),
             ci.lower = bca[1], ci.upper = bca[2],
             p.value = p_value)
}
