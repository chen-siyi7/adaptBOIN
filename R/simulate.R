#' Simulate one built-in scenario
#'
#' @param sc_idx Scenario index, 1 to 8.
#' @param design_name Design name; see \code{\link{design_codes}}.
#' @param params Parameter list.
#' @param n_sim Number of simulated trials.
#' @param atbl Optional pre-built decision table.
#' @return A list of operating characteristics with scenario metadata attached.
#' @examples
#' r <- run_scenario(1, "boin", n_sim = 200L)
#' r$pcs
#' @export
run_scenario <- function(sc_idx, design_name, params = adapt_params(),
                         n_sim = 2000L, atbl = NULL) {
  if (is.null(atbl)) atbl <- adapt_table(params)
  base_seed <- 1000L * sc_idx + design_codes[[design_name]] * 100L + 1L
  res <- run_scenario_cpp(scenarios[[sc_idx]]$pi, design_codes[[design_name]],
                          atbl, params, as.integer(n_sim), as.integer(base_seed))
  res$sc_idx   <- sc_idx
  res$sc_name  <- scenarios[[sc_idx]]$name
  res$true_mtd <- scenarios[[sc_idx]]$mtd
  res$design   <- design_name
  res
}

#' Simulate an arbitrary dose-toxicity curve
#'
#' Unlike \code{\link{run_scenario}} this takes the curve directly, which the
#' targeted studies use. Passing the same \code{base_seed} to several designs
#' gives a partially paired comparison, since trials then start from a common
#' random stream and diverge only once decisions differ.
#'
#' @param pi True DLT probabilities.
#' @param design_name Design name.
#' @param params Parameter list.
#' @param n_sim Number of simulated trials.
#' @param atbl Pre-built decision table.
#' @param base_seed Base seed.
#' @return A list of operating characteristics.
#' @export
run_scenario_pi <- function(pi, design_name, params, n_sim, atbl, base_seed) {
  r <- run_scenario_cpp(pi, design_codes[[design_name]], atbl, params,
                        as.integer(n_sim), as.integer(base_seed))
  r$design <- design_name
  r
}

#' Simulate a grid of scenarios and designs
#'
#' @param params Parameter list.
#' @param n_sim Number of simulated trials per cell.
#' @param designs Design names to run.
#' @param scenarios_idx Scenario indices to run.
#' @param n_cores Number of cores. Uses fork-based parallelism, so values above
#'   one have no effect on Windows.
#' @param verbose Print progress.
#' @return A list of results, one per scenario and design.
#' @examples
#' res <- run_all(n_sim = 100L, designs = c("boin", "adaptive_iso"),
#'                scenarios_idx = 1:2)
#' summarise_results(res)
#' @export
run_all <- function(params = adapt_params(), n_sim = 2000L,
                    designs = names(design_codes), scenarios_idx = 1:8,
                    n_cores = 1L, verbose = FALSE) {
  atbl   <- adapt_table(params)
  combos <- expand.grid(sc = scenarios_idx, design = designs,
                        stringsAsFactors = FALSE)
  worker <- function(i) run_scenario(combos$sc[i], combos$design[i],
                                     params, n_sim, atbl)
  if (n_cores > 1L && .Platform$OS.type == "unix" &&
      requireNamespace("parallel", quietly = TRUE)) {
    parallel::mclapply(seq_len(nrow(combos)), worker, mc.cores = n_cores)
  } else {
    if (n_cores > 1L) message("Multi-core requires a unix fork; running serially.")
    out <- vector("list", nrow(combos))
    for (i in seq_len(nrow(combos))) {
      if (verbose) message(sprintf("  [%d/%d] sc=%d %s", i, nrow(combos),
                                   combos$sc[i], combos$design[i]))
      out[[i]] <- worker(i)
    }
    out
  }
}

#' Collapse simulation results into a data frame
#'
#' @param results List returned by \code{\link{run_all}}.
#' @return A data frame, one row per scenario and design.
#' @export
summarise_results <- function(results) {
  do.call(rbind, lapply(results, function(r) data.frame(
    sc = r$sc_idx, sc_name = r$sc_name, true_mtd = r$true_mtd, design = r$design,
    pcs = round(r$pcs, 1), pod = round(r$pod, 1), pcs_se = round(r$pcs_se, 2),
    ess_med = if (!is.null(r$ess_median)) round(r$ess_median) else NA_real_,
    ess_p10 = if (!is.null(r$ess_p10)) round(r$ess_p10) else NA_real_,
    pct_fb = round(r$pct_fallback, 1), mean_dlts = round(r$mean_dlts, 1),
    pct_over = round(r$pct_over_sel, 1), pct_under = round(r$pct_under_sel, 1),
    pct_mtd_elim = if (!is.null(r$pct_mtd_elim)) round(r$pct_mtd_elim, 1)
                   else NA_real_,
    stringsAsFactors = FALSE)))
}

#' Final dose-selection distribution
#'
#' @param res A single result from \code{\link{run_scenario}}.
#' @return A list with \code{no_rec}, the percentage of trials ending without a
#'   recommendation, and \code{dose}, the percentage selecting each dose.
#' @export
selection_dist <- function(res) {
  md <- as.numeric(res$mtd_dist)
  n  <- res$n_sim
  list(no_rec = 100 * md[1] / n, dose = 100 * md[-1] / n)
}
