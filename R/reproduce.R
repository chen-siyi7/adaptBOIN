#' Regenerate every table and figure in the manuscript
#'
#' Runs the primary grid, the sensitivity analyses, the targeted comparator
#' studies and the elimination-rule comparison, writing CSV files and, if
#' \pkg{ggplot2} is installed, PDF figures.
#'
#' The full run at 2000 trials is slow, chiefly because the two Bernstein
#' configurations draw 2000 importance samples per trial. Start with a small
#' \code{n_sim} to confirm the pipeline works.
#'
#' @param out_dir Output directory; created if absent.
#' @param n_sim Number of simulated trials per cell.
#' @param n_cores Number of cores. Fork-based, so values above one have no
#'   effect on Windows.
#' @param do_figures Write figures.
#' @param do_studies Run the steepness, offset-MTD and elimination studies.
#' @param params Parameter list.
#' @return Invisibly, a list holding every result object.
#' @examples
#' \dontrun{
#' reproduce_paper(n_sim = 200L)                     # quick check
#' reproduce_paper(n_sim = 2000L, n_cores = 6L)      # full run
#' }
#' @export
reproduce_paper <- function(out_dir = "adaptboin_output", n_sim = 2000L,
                            n_cores = 1L, do_figures = TRUE,
                            do_studies = TRUE, params = adapt_params()) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  wr <- function(x, f) utils::write.csv(x, file.path(out_dir, f),
                                        row.names = FALSE)

  message("Primary grid (", n_sim, " trials per cell) ...")
  results <- run_all(params, n_sim = n_sim, n_cores = n_cores, verbose = TRUE)
  df <- summarise_results(results)
  wr(df, "sim_results.csv")
  saveRDS(results, file.path(out_dir, "results.rds"))

  message("Tables ...")
  wr(boundary_rates(params),   "table_boundaries.csv")
  wr(boundaries_by_loss(params = params), "table_loss_boundaries.csv")
  wr(table_scenarios(),        "table_scenarios.csv")
  wr(table_selection(results, 1), "table_selection_sc1.csv")
  for (sc in which(vapply(scenarios, function(s) is_tie(s$pi), logical(1))))
    wr(table_tie_pcs(results, sc), sprintf("table_tiePCS_sc%d.csv", sc))
  t6 <- table_safety(df)
  wr(t6$over, "table_overdose.csv")
  wr(t6$dlts, "table_dlts.csv")
  wr(t6$pod,  "table_pod.csv")

  message("Sensitivity analyses ...")
  crm_res <- crm_sensitivity(n_sim = n_sim, params = params)
  wr(table_crm(crm_res), "table_crm.csv")
  ess_res <- ess_sensitivity(n_sim = n_sim, params = params)
  wr(do.call(rbind, lapply(ess_res, function(r) data.frame(
    sc = r$sc_idx, thr = r$thr, pcs = round(r$pcs, 1),
    fallback = round(r$pct_fallback, 1)))), "ess_sensitivity.csv")
  ka_res <- bern_KA_sensitivity(n_sim = n_sim, params = params)
  wr(do.call(rbind, lapply(ka_res, function(r) data.frame(
    K = r$K, alpha = r$alpha, sc = r$sc_idx, pcs = round(r$pcs, 1)))),
    "bern_KA_sensitivity.csv")
  n_res <- n_sensitivity(n_sim = n_sim, params = params)
  wr(do.call(rbind, lapply(n_res, function(r) data.frame(
    sc = r$sc_idx, design = r$design, N = r$N_val, pcs = round(r$pcs, 1)))),
    "n_sensitivity.csv")
  tie_res <- tie_sensitivity(n_sim = n_sim, params = params)
  wr(tie_res, "tie_sensitivity.csv")
  loss_res <- loss_sensitivity(n_sim = n_sim, params = params)
  wr(do.call(rbind, lapply(loss_res, function(r) data.frame(
    k_over = r$k_over, admissible = r$admissible, sc = r$sc_idx,
    pcs = round(r$pcs, 1), pct_over = round(r$pct_over_sel, 1),
    pod = round(r$pod, 1)))),
    "loss_sensitivity.csv")

  steep_df <- NULL; str_df <- NULL; elim_df <- NULL
  if (do_studies) {
    message("Steepness gradient ...")
    steep_df <- steep_curve_study(n_sim = n_sim, params = params)
    wr(steep_df, "study_steep_curve.csv")
    message("Offset MTD across sample sizes ...")
    str_df <- straddle_N_study(n_sim = n_sim, params = params)
    wr(str_df, "study_offset_mtd_byN.csv")
    message("Elimination-rule comparison ...")
    elim_df <- elim_rule_study(n_sim = n_sim, params = params,
                               Nvals = c(30L, 60L, 90L))
    wr(elim_df, "study_elim_rule.csv")
  }

  if (do_figures && requireNamespace("ggplot2", quietly = TRUE)) {
    message("Figures ...")
    sv <- function(p, f, w = 7, h = 5)
      ggplot2::ggsave(file.path(out_dir, f), p, width = w, height = h)
    sv(fig_ess(params, n_sim),      "fig_ess.pdf")
    sv(fig_pcs(df),                 "fig_pcs.pdf")
    sv(fig_decomposition(df),       "fig_decomposition.pdf")
    sv(fig_nsens(n_res),            "fig_nsensitivity.pdf", w = 8, h = 6)
    sv(fig_loss(loss_res),          "fig_loss_tradeoff.pdf")
    if (!is.null(steep_df)) sv(fig_steep(steep_df), "fig_steepness.pdf", w = 8)
    if (!is.null(str_df))   sv(fig_straddle(str_df), "fig_offset_mtd.pdf",
                               w = 10, h = 6)
    if (!is.null(elim_df))  sv(fig_elim(elim_df[elim_df$N == 90, ]),
                               "fig_elim_rule.pdf", w = 10, h = 6)
  } else if (do_figures) {
    message("ggplot2 not installed; skipping figures.")
  }

  message("Done. Outputs in ", normalizePath(out_dir))
  invisible(list(df = df, results = results, crm = crm_res, ess = ess_res,
                 ka = ka_res, nsens = n_res, loss = loss_res, tie = tie_res,
                 steep = steep_df, straddle = str_df, elim = elim_df))
}
