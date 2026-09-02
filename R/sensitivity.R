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

#' gBOINS calibration sensitivity at a target rate of 0.25
#'
#' The gBOINS paper reports binary-endpoint calibration constants for target
#' rates 0.20 and 0.30, but not 0.25. This study brackets the primary
#' interpolation by using the two published endpoint calibrations and their
#' midpoint on the multiplier scale. The same deterministic trial seeds are
#' used for all three calibrations within each scenario.
#'
#' @param scenarios_idx Scenario indices.
#' @param n_sim Trials per cell.
#' @param params Parameter list.
#' @return A data frame of accuracy and safety operating characteristics.
#' @export
gboins_calibration_sensitivity <- function(
    scenarios_idx = seq_along(scenarios), n_sim = 2000L,
    params = adapt_params()) {
  multipliers <- c(published_020 = 1.05, midpoint_025 = 1.075,
                   published_030 = 1.10)
  out <- list()
  for (cal in names(multipliers)) {
    pp <- params
    pp$gb_c1 <- log(multipliers[[cal]])
    pp$gb_c2 <- pp$gb_c1 / 3
    atbl <- adapt_table(pp)
    for (sc in scenarios_idx) {
      r <- run_scenario(sc, "gboins", pp, n_sim, atbl)
      out[[length(out) + 1L]] <- data.frame(
        calibration = cal, multiplier = multipliers[[cal]],
        c1 = pp$gb_c1, c2 = pp$gb_c2, sc = sc,
        pcs = round(r$pcs, 1), pct_over = round(r$pct_over_sel, 1),
        pod = round(r$pod, 1), mean_dlts = round(r$mean_dlts, 1),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
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
  kmax <- kappa_max(params)
  res <- list()
  for (k in k_overs) {
    pp <- params; pp$loss <- make_loss(k_over = k)
    # Values at or above kappa_max() violate the no-DLT condition. They are
    # included deliberately to show the failure, so the warning that
    # adapt_table() would raise is replaced by a single informative message.
    inadmissible <- k >= kmax
    if (inadmissible)
      message(sprintf(
        "  k_over = %.2f exceeds the admissible bound %.2f; included to show the failure.",
        k, kmax))
    atbl <- adapt_table(pp, check = !inadmissible)
    for (sc in monotone_idx) {
      r <- run_scenario(sc, "adaptive_iso", pp, n_sim, atbl)
      r$k_over <- k; r$admissible <- !inadmissible
      res[[length(res) + 1]] <- r
    }
  }
  res
}

#' Tie-convention sensitivity
#'
#' Reruns the primary grid under both tie conventions. Because the convention is
#' applied consistently to the definition of the true MTD and to every
#' end-of-trial estimator, switching it changes both what the design does and
#' what counts as correct, so the two columns are not directly comparable as
#' estimates of the same quantity. The comparison shows how sensitive reported
#' accuracy is to a convention that is often left unstated.
#'
#' @param n_sim Trials per cell.
#' @param params Parameter list; \code{tie_high} is overridden.
#' @param designs Design names.
#' @param scenarios_idx Scenario indices.
#' @return A data frame with accuracy under each convention and the difference.
#' @examples
#' \dontrun{
#' tie_sensitivity(n_sim = 2000L)
#' }
#' @export
tie_sensitivity <- function(n_sim = 2000L, params = adapt_params(),
                            designs = c("adaptive_iso", "boin", "aboin", "gboins"),
                            scenarios_idx = monotone_idx) {
  out <- list()
  for (th in c(TRUE, FALSE)) {
    pp <- params; pp$tie_high <- th
    atbl <- adapt_table(pp)
    for (sc in scenarios_idx) for (dn in designs) {
      r <- run_scenario(sc, dn, pp, n_sim, atbl)
      out[[length(out) + 1]] <- data.frame(
        sc = sc, sc_name = scenarios[[sc]]$name,
        mtd = mtd_of(scenarios[[sc]]$pi, pp$phi_tgt, th),
        tie = if (th) "higher" else "lower", design = dn,
        pcs = round(r$pcs, 1), pct_over = round(r$pct_over_sel, 1),
        pod = round(r$pod, 1), stringsAsFactors = FALSE)
    }
  }
  d <- do.call(rbind, out)
  hi <- d[d$tie == "higher", c("sc", "design", "pcs")]
  lo <- d[d$tie == "lower",  c("sc", "design", "pcs")]
  names(hi)[3] <- "pcs_higher"; names(lo)[3] <- "pcs_lower"
  m <- merge(hi, lo, by = c("sc", "design"))
  m$difference <- round(m$pcs_lower - m$pcs_higher, 1)
  m[order(m$sc, m$design), ]
}
