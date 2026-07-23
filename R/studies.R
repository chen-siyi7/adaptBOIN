#' Steepness-gradient scenarios
#'
#' The MTD is fixed at the third dose with toxicity exactly at the target, and
#' the gap to the next dose increases. Isolating steepness as a single factor
#' makes a trend easier to interpret than a set of unrelated scenarios.
#'
#' @format A list of lists with \code{name}, \code{gap}, \code{pi} and \code{mtd}.
#' @export
steep_scenarios <- local({
  gaps <- c(0.05, 0.10, 0.15, 0.25, 0.40)
  lapply(gaps, function(g) {
    d4 <- 0.25 + g
    pi <- c(0.03, 0.08, 0.25, d4, min(d4 + 0.20, 0.92), min(d4 + 0.33, 0.97))
    list(name = sprintf("gap %.2f", g), gap = g, pi = pi, mtd = 3L)
  })
})

#' Offset-MTD scenarios
#'
#' Configurations in which no dose sits at exactly the target rate, which is the
#' generic case on a discrete dose grid. The third is an exact tie and is scored
#' with the acceptable-selection criterion as well as the strict one.
#'
#' @format A list of lists with \code{name} and \code{pi}.
#' @export
offset_scenarios <- list(
  list(name = "MTD below target (0.22)",
       pi = c(0.05, 0.10, 0.22, 0.42, 0.58, 0.72)),
  list(name = "MTD above target (0.28)",
       pi = c(0.04, 0.09, 0.28, 0.45, 0.60, 0.75)),
  list(name = "Symmetric straddle (0.20/0.30)",
       pi = c(0.05, 0.10, 0.20, 0.30, 0.48, 0.62))
)

#' Steepness-gradient study
#'
#' @param n_sim Trials per cell.
#' @param params Parameter list.
#' @param designs Design names to compare.
#' @return A data frame.
#' @export
steep_curve_study <- function(n_sim = 2000L, params = adapt_params(),
                              designs = c("adaptive_iso", "boin", "gboins")) {
  atbl <- adapt_table(params); out <- list()
  for (i in seq_along(steep_scenarios)) {
    s <- steep_scenarios[[i]]
    base <- 700000L + 1000L * i          # shared across designs
    for (dn in designs) {
      r <- run_scenario_pi(s$pi, dn, params, n_sim, atbl, base)
      out[[length(out) + 1]] <- data.frame(
        scen = s$name, gap = s$gap, mtd = s$mtd, design = dn,
        pcs = round(r$pcs, 1), pct_over = round(r$pct_over_sel, 1),
        pod = round(r$pod, 1), dlts = round(r$mean_dlts, 2),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

#' Offset-MTD study across sample sizes
#'
#' @param Nvals Maximum sample sizes.
#' @param n_sim Trials per cell.
#' @param params Parameter list.
#' @param designs Design names to compare.
#' @return A data frame including strict and acceptable-selection PCS.
#' @export
straddle_N_study <- function(Nvals = c(30L, 45L, 60L, 90L), n_sim = 2000L,
                             params = adapt_params(),
                             designs = c("adaptive_iso", "boin", "gboins")) {
  out <- list()
  for (N in Nvals) {
    pp <- params; pp$N <- as.integer(N)
    atbl <- adapt_table(pp)
    for (i in seq_along(offset_scenarios)) {
      s <- offset_scenarios[[i]]
      acc <- acceptable_set(s$pi)
      base <- 800000L + 1000L * i + N    # shared across designs
      for (dn in designs) {
        r <- run_scenario_pi(s$pi, dn, pp, n_sim, atbl, base)
        sd <- selection_dist(r)
        out[[length(out) + 1]] <- data.frame(
          scen = s$name, N = N, mtd = mtd_of(s$pi), design = dn,
          pcs = round(r$pcs, 1), pcs_acc = round(sum(sd$dose[acc]), 1),
          pct_over = round(r$pct_over_sel, 1), pod = round(r$pod, 1),
          stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, out)
}

#' Elimination-rule comparison
#'
#' Compares the standard target-based rule against an equivalence-bound rule for
#' every scenario whose true MTD lies above the target, where the two differ.
#' Seeds are shared across designs and across rules, so the contrast is paired.
#'
#' @param n_sim Trials per cell.
#' @param params Parameter list.
#' @param designs Design names to compare.
#' @param Nvals Maximum sample sizes; defaults to the value in \code{params}.
#' @return A data frame including \code{pct_mtd_elim}, the percentage of trials
#'   in which the true MTD was eliminated.
#' @export
elim_rule_study <- function(n_sim = 2000L, params = adapt_params(),
                            designs = c("adaptive_iso", "boin", "gboins"),
                            Nvals = NULL) {
  scs <- lapply(affected_scenarios(), function(i)
    list(name = sprintf("Sc%d %s", i, scenarios[[i]]$name),
         pi = scenarios[[i]]$pi))
  scs <- c(scs, list(list(name = "Offset MTD above target (0.28)",
                          pi = offset_scenarios[[2]]$pi)))
  rules <- list(
    standard   = list(phi_elim = params$phi_tgt, elim_a0 = 1.0, elim_b0 = 1.0),
    equiv_phi2 = list(phi_elim = params$phi2, elim_a0 = params$a0,
                      elim_b0 = params$b0))
  if (is.null(Nvals)) Nvals <- params$N
  out <- list()
  for (N in Nvals) for (i in seq_along(scs)) {
    s <- scs[[i]]
    base <- 950000L + 1000L * i + N
    for (rn in names(rules)) {
      pp <- params; pp$N <- as.integer(N)
      pp[names(rules[[rn]])] <- rules[[rn]]
      atbl <- adapt_table(pp)
      for (dn in designs) {
        r <- run_scenario_pi(s$pi, dn, pp, n_sim, atbl, base)
        out[[length(out) + 1]] <- data.frame(
          scen = s$name, N = N, mtd = mtd_of(s$pi),
          pi_mtd = s$pi[mtd_of(s$pi)], elim_rule = rn, design = dn,
          pcs = round(r$pcs, 1),
          pct_mtd_elim = if (!is.null(r$pct_mtd_elim))
            round(r$pct_mtd_elim, 1) else NA_real_,
          pct_over = round(r$pct_over_sel, 1),
          pod = round(r$pod, 1), dlts = round(r$mean_dlts, 2),
          stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, out)
}
