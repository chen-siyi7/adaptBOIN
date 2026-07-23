#' CRM skeleton sensitivity
#'
#' @param scenarios_idx Scenario indices.
#' @param n_sim Trials per cell.
#' @param params Parameter list.
#' @return A list of results, each tagged with \code{skel}.
#' @export
crm_sensitivity <- function(scenarios_idx = c(1, 3, 4, 5, 6, 8),
                            n_sim = 2000L, params = adapt_params()) {
  skels <- list(
    calibrated   = c(0.05, 0.12, 0.25, 0.40, 0.55, 0.70),
    shifted_up   = c(0.10, 0.17, 0.30, 0.45, 0.60, 0.75),
    shifted_down = c(0.01, 0.07, 0.20, 0.35, 0.50, 0.65))
  atbl <- adapt_table(params); res <- list()
  for (sk in names(skels)) {
    pp <- params; pp$crm_skel <- skels[[sk]]
    for (sc in scenarios_idx) {
      r <- run_scenario(sc, "crm", pp, n_sim, atbl); r$skel <- sk
      res[[length(res) + 1]] <- r
    }
  }
  res
}

#' Bernstein fallback threshold sensitivity
#'
#' @param thresholds Effective sample size thresholds.
#' @param n_sim Trials per cell.
#' @param params Parameter list.
#' @return A list of results, each tagged with \code{thr}.
#' @export
ess_sensitivity <- function(thresholds = c(50, 75, 100, 150, 200),
                            n_sim = 2000L, params = adapt_params()) {
  atbl <- adapt_table(params); res <- list()
  for (thr in thresholds) {
    pp <- params; pp$ESS_thr <- thr
    for (sc in monotone_idx) {
      r <- run_scenario(sc, "boin_bern", pp, n_sim, atbl); r$thr <- thr
      res[[length(res) + 1]] <- r
    }
  }
  res
}

#' Sample-size sensitivity
#'
#' @param Nvals Maximum sample sizes.
#' @param designs Design names.
#' @param n_sim Trials per cell.
#' @param params Parameter list.
#' @return A list of results, each tagged with \code{N_val}.
#' @export
n_sensitivity <- function(Nvals = c(15L, 20L, 25L, 30L),
                          designs = c("adaptive_iso", "boin"),
                          n_sim = 2000L, params = adapt_params()) {
  res <- list()
  for (N in Nvals) {
    pp <- params; pp$N <- as.integer(N); atbl <- adapt_table(pp)
    for (sc in monotone_idx) for (dn in designs) {
      r <- run_scenario(sc, dn, pp, n_sim, atbl); r$N_val <- N
      res[[length(res) + 1]] <- r
    }
  }
  res
}

#' Bernstein degree and concentration sensitivity
#'
#' @param Ks Bernstein degrees.
#' @param alphas Dirichlet concentrations. Values below one concentrate mass at
#'   the simplex vertices; values above one concentrate toward the centroid.
#' @param n_sim Trials per cell.
#' @param params Parameter list.
#' @return A list of results, each tagged with \code{K} and \code{alpha}.
#' @export
bern_KA_sensitivity <- function(Ks = c(3L, 4L, 5L), alphas = c(0.5, 1.0, 2.0),
                                n_sim = 2000L, params = adapt_params()) {
  res <- list()
  for (K in Ks) for (a in alphas) {
    pp <- params; pp$K <- as.integer(K); pp$alpha_d <- a
    atbl <- adapt_table(pp)
    for (sc in monotone_idx) {
      r <- run_scenario(sc, "boin_bern", pp, n_sim, atbl)
      r$K <- K; r$alpha <- a
      res[[length(res) + 1]] <- r
    }
  }
  res
}

#' Overdose-aversion sweep
#'
#' Rebuilds the decision table for each multiplier and measures the resulting
#' accuracy and safety tradeoff.
#'
#' @param k_overs Overdose-aversion multipliers. Note that 3 is inadmissible and
#'   is included only for completeness.
#' @param n_sim Trials per cell.
#' @param params Parameter list.
#' @return A list of results, each tagged with \code{k_over}.
#' @export
loss_sensitivity <- function(k_overs = c(1.0, 1.5, 2.0, 3.0),
                             n_sim = 2000L, params = adapt_params()) {
  res <- list()
  for (k in k_overs) {
    pp <- params; pp$loss <- make_loss(k_over = k)
    atbl <- adapt_table(pp)
    for (sc in monotone_idx) {
      r <- run_scenario(sc, "adaptive_iso", pp, n_sim, atbl); r$k_over <- k
      res[[length(res) + 1]] <- r
    }
  }
  res
}
