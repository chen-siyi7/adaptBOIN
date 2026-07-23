#' Per-dose allocation diagnostic
#'
#' Loops \code{\link{one_trial}} to recover where patients actually went and how
#' selection related to allocation, which \code{\link{run_scenario_cpp}} does not
#' expose. This separates two explanations for poor selection: escalation
#' drifting off the correct dose, versus correct allocation followed by a failure
#' at the end-of-trial step.
#'
#' @param pi True DLT probabilities.
#' @param design_name Design name.
#' @param N Maximum sample size.
#' @param n_sim Number of simulated trials.
#' @param params Parameter list.
#' @param base_seed Base seed.
#' @return A data frame with one row per dose, including mean allocation,
#'   observed rate, selection percentage, and accuracy conditional on the correct
#'   dose having received the most patients.
#' @export
diagnose_alloc <- function(pi, design_name, N, n_sim = 2000L,
                           params = adapt_params(), base_seed = 900000L) {
  pp <- params; pp$N <- as.integer(N)
  atbl <- adapt_table(pp)
  J <- pp$J
  mtd <- mtd_of(pi)
  nmat <- matrix(0L, n_sim, J); smat <- matrix(0L, n_sim, J)
  sel  <- integer(n_sim)
  for (i in seq_len(n_sim)) {
    tr <- one_trial(pi, design_codes[[design_name]], atbl, pp,
                    as.integer(base_seed + i))
    nmat[i, ] <- as.integer(tr$n_vec)
    smat[i, ] <- as.integer(tr$s_vec)
    sel[i]    <- as.integer(tr$mtd_sel)   # 0 = no recommendation
  }
  argmax_n <- apply(nmat, 1, which.max)
  correct_alloc <- argmax_n == mtd
  data.frame(
    design = design_name, N = N, mtd = mtd, dose = 1:J,
    mean_n = round(colMeans(nmat), 2),
    mean_rate = round(colSums(smat) / pmax(colSums(nmat), 1), 3),
    pct_visited = round(100 * colMeans(nmat > 0), 1),
    pct_selected = round(100 * vapply(1:J, function(j) mean(sel == j), 0), 1),
    pct_norec = round(100 * mean(sel == 0), 1),
    pct_alloc_max = round(100 * vapply(1:J, function(j) mean(argmax_n == j), 0), 1),
    pcs_given_correct_alloc = round(100 * mean(sel[correct_alloc] == mtd), 1),
    pct_correct_alloc = round(100 * mean(correct_alloc), 1),
    stringsAsFactors = FALSE)
}

#' Allocation diagnostic for the offset-MTD scenario
#'
#' Runs \code{\link{diagnose_alloc}} across sample sizes and designs for the
#' scenario whose MTD lies above the target, where accuracy is non-monotone in
#' sample size under the standard elimination rule.
#'
#' @param n_sim Number of simulated trials.
#' @param params Parameter list.
#' @param Nvals Maximum sample sizes.
#' @param designs Design names.
#' @return A data frame.
#' @export
diagnose_offset_anomaly <- function(n_sim = 2000L, params = adapt_params(),
                                    Nvals = c(30L, 60L, 90L),
                                    designs = c("adaptive_iso", "boin",
                                                "gboins")) {
  pi <- offset_scenarios[[2]]$pi
  out <- list()
  for (N in Nvals) for (dn in designs) {
    message(sprintf("  diagnosing %s at N=%d", dn, N))
    out[[length(out) + 1]] <- diagnose_alloc(pi, dn, N, n_sim, params,
                                             base_seed = 900000L + N * 10L)
  }
  do.call(rbind, out)
}
