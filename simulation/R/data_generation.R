###############################################################################
# Manuscript-reference data generation
###############################################################################


#' Generate clustered data from the manuscript reference design.
#'
#' Generates one trial sample and one auxiliary sample using the same
#' data-generating mechanism used to validate the PopART estimator. Trial and
#' auxiliary sample sizes are rounded to multiples of the number of clusters.
#'
#' @param m Integer number of clusters.
#' @param n_trial Integer nominal trial sample size.
#' @param n_auxiliary Integer nominal auxiliary sample size.
#' @param p_resp Numeric response probability setting. The reference design uses the
#'   original fixed value `0.5`; changing it would require deriving a new
#'   response-model intercept.
#' @param p_cens Numeric censoring probability setting. The reference design uses the
#'   original fixed value `0.3`; changing it would require deriving a new
#'   censoring-model intercept.
#' @param seed Integer random seed for the full Monte Carlo replicate.
#'
#' @return A list with:
#' \itemize{
#'   \item `dat`: combined trial/auxiliary data frame.
#'   \item `eta0`, `eta1`: true potential outcome means in the trial sample.
#'   \item `n_trial`, `n_auxiliary`: rounded sample sizes actually used.
#' }
#' @export
generate_reference_data <- function(m, n_trial, n_auxiliary, p_resp, p_cens,
                                    seed) {
  load_reference_study_packages()
  n_aux <- n_auxiliary
  m <- suppressWarnings(as.integer(m))
  n_trial <- suppressWarnings(as.integer(n_trial))
  n_aux <- suppressWarnings(as.integer(n_aux))
  seed <- suppressWarnings(as.integer(seed))
  if (is.na(m) || m < 2L || m > 20L || m %% 2L != 0L) {
    stop("m must be an even integer between 2 and 20.", call. = FALSE)
  }
  if (is.na(n_trial) || is.na(n_aux) ||
      n_trial < m || n_aux < m) {
    stop("n_trial and n_aux must each be at least m.", call. = FALSE)
  }
  if (is.na(seed) || seed < 0L) {
    stop("seed must be a nonnegative integer.", call. = FALSE)
  }
  p_resp <- suppressWarnings(as.numeric(p_resp))
  p_cens <- suppressWarnings(as.numeric(p_cens))
  if (length(p_resp) != 1L || !is.finite(p_resp) ||
      abs(p_resp - 0.5) > 1e-12) {
    stop(
      "The manuscript reference design requires p_resp = 0.5.",
      call. = FALSE
    )
  }
  if (length(p_cens) != 1L || !is.finite(p_cens) ||
      abs(p_cens - 0.3) > 1e-12) {
    stop(
      "The manuscript reference design requires p_cens = 0.3.",
      call. = FALSE
    )
  }

  R_int <- 0.4449358
  C_int <- -0.8711302

  n_trial <- m * round(n_trial / m)
  n_aux <- m * round(n_aux / m)

  set.seed(seed)

  XX1 <- c(
    0.92, 0.78, 0.07, -1.99, 0.62, -0.06, -0.16, -1.47, -0.48, 0.42,
    1.36, -0.10, 0.39, -0.05, -1.38, -0.41, -0.39, -0.06, 1.10, 0.76
  )
  AA <- sample(rep(c(0, 1), times = m / 2), replace = FALSE)

  trial_dat <- data.frame(id = seq_len(n_trial)) %>%
    mutate(
      clust = rep(seq_len(m), length = n_trial),
      X1 = XX1[clust],
      A = AA[clust],
      W1 = rbinom(n_trial, 1, 0.5),
      W2 = rnorm(n_trial, 0, 1),
      W3 = rnorm(n_trial, 0, 1),
      pR = plogis(R_int - 2 * A * W1 + 0.2 * W2 * W3^2),
      R = rbinom(n_trial, 1, pR),
      pC = plogis(C_int - 0.25 * A + 0.25 * W1 + 0.2 * W2^2 * W3),
      C = rbinom(n_trial, 1, pC),
      C = ifelse(R == 1, C, 0),
      mu0 = plogis(-1 + 2 * W1),
      mu1 = plogis(0 - 1 * W1 + 0.2 * sin(W2) - 0.2 * W3^2 + 0.25 * X1),
      Y0 = rbinom(n_trial, 1, mu0),
      Y1 = rbinom(n_trial, 1, mu1),
      Y = (1 - A) * Y0 + A * Y1,
      Y = ifelse(R == 1 & C == 0, Y, 0),
      wt = 1
    )

  aux_dat <- data.frame(id = seq_len(n_aux) + n_trial) %>%
    mutate(
      clust = rep(seq_len(m), length = n_aux),
      X1 = XX1[clust],
      A = AA[clust],
      W1 = rbinom(n_aux, 1, 0.75),
      W2 = rnorm(n_aux, 0, 1),
      W3 = rnorm(n_aux, 0, 1),
      wt = ifelse(W1 == 0, 0.75, 0.25),
      R = 0,
      C = 0,
      Y = 0
    )
  aux_dat$wt <- aux_dat$wt / mean(aux_dat$wt)

  dat <- bind_rows(
    trial_dat %>% select(-c(pR, pC, mu0, mu1, Y0, Y1)) %>% mutate(S = 1),
    aux_dat %>% mutate(S = 0)
  )

  result <- list(
    dat = dat,
    eta0 = mean(trial_dat$mu0),
    eta1 = mean(trial_dat$mu1),
    n_trial = n_trial,
    n_aux = n_aux
  )
  result$n_auxiliary <- result$n_aux
  result
}
