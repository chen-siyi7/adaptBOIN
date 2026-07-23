.need_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for the figure functions. ",
         "Install it with install.packages(\"ggplot2\").", call. = FALSE)
}

#' Importance-sampling effective sample size and fallback rate
#'
#' @param params Parameter list.
#' @param n_sim Trials per scenario.
#' @return A ggplot object.
#' @export
fig_ess <- function(params = adapt_params(), n_sim = 2000L) {
  .need_ggplot()
  atbl <- adapt_table(params)
  d <- do.call(rbind, lapply(seq_along(scenarios), function(sc) {
    r <- run_scenario(sc, "boin_bern", params, n_sim, atbl)
    data.frame(sc = sc, ess_med = r$ess_median, ess_p10 = r$ess_p10,
               fallback = r$pct_fallback)
  }))
  d$sc_lab <- factor(paste0("Sc ", d$sc), levels = paste0("Sc ", d$sc))
  scale <- max(d$ess_med, na.rm = TRUE) / 100
  ggplot2::ggplot(d, ggplot2::aes(.data$sc_lab)) +
    ggplot2::geom_col(ggplot2::aes(y = .data$ess_med), fill = "steelblue4",
                      width = 0.6) +
    ggplot2::geom_col(ggplot2::aes(y = .data$ess_p10), fill = "lightblue",
                      width = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = .data$fallback * scale, group = 1),
                       color = "firebrick", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 100, linetype = "dotted") +
    ggplot2::scale_y_continuous(
      name = "IS effective sample size",
      sec.axis = ggplot2::sec_axis(~ . / scale, name = "% fallback")) +
    ggplot2::labs(x = NULL) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Accuracy by design and scenario
#'
#' @param df Data frame from \code{\link{summarise_results}}.
#' @return A ggplot object.
#' @export
fig_pcs <- function(df) {
  .need_ggplot()
  d <- df[df$design %in% primary_designs, ]
  d$sc_lab <- factor(paste0("Sc ", d$sc, " ", d$sc_name),
                     levels = rev(unique(paste0("Sc ", d$sc, " ", d$sc_name))))
  ggplot2::ggplot(d, ggplot2::aes(.data$pcs, .data$sc_lab,
                                  color = .data$design, shape = .data$design)) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::labs(x = "Probability of correct selection (%)", y = NULL,
                  color = "Design", shape = "Design") +
    ggplot2::theme_minimal(base_size = 11)
}

#' Decomposition of accuracy differences relative to BOIN
#'
#' @param df Data frame from \code{\link{summarise_results}}.
#' @return A ggplot object.
#' @export
fig_decomposition <- function(df) {
  .need_ggplot()
  base <- df[df$design == "boin" & df$sc %in% monotone_idx,
             c("sc", "sc_name", "pcs")]
  names(base)[3] <- "boin"
  ai <- df[df$design == "adaptive_iso" & df$sc %in% monotone_idx, c("sc", "pcs")]
  bb <- df[df$design == "boin_bern"   & df$sc %in% monotone_idx, c("sc", "pcs")]
  m <- merge(merge(base, stats::setNames(ai, c("sc", "ai")), by = "sc"),
             stats::setNames(bb, c("sc", "bb")), by = "sc")
  long <- rbind(
    data.frame(sc = m$sc, sc_name = m$sc_name,
               component = "Boundaries (Adpt+Iso - BOIN)", value = m$ai - m$boin),
    data.frame(sc = m$sc, sc_name = m$sc_name,
               component = "Bernstein (BOIN+Bern - BOIN)", value = m$bb - m$boin))
  long$sc_lab <- factor(paste0("Sc ", long$sc, " ", long$sc_name),
                        levels = rev(unique(paste0("Sc ", long$sc, " ",
                                                   long$sc_name))))
  ggplot2::ggplot(long, ggplot2::aes(.data$value, .data$sc_lab,
                                     fill = .data$component)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
                      width = 0.6) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::labs(x = "PCS difference vs BOIN (percentage points)", y = NULL,
                  fill = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

#' Accuracy against sample size
#'
#' @param n_res List from \code{\link{n_sensitivity}}.
#' @param scs Scenario indices to show.
#' @return A ggplot object.
#' @export
fig_nsens <- function(n_res, scs = c(1, 6, 3, 5)) {
  .need_ggplot()
  d <- do.call(rbind, lapply(n_res, function(r)
    data.frame(sc = r$sc_idx, sc_name = r$sc_name, design = r$design,
               N = r$N_val, pcs = r$pcs)))
  d <- d[d$sc %in% scs, ]
  d$panel <- factor(paste0("Sc ", d$sc, " ", d$sc_name),
                    levels = unique(paste0("Sc ", d$sc, " ", d$sc_name)))
  ggplot2::ggplot(d, ggplot2::aes(.data$N, .data$pcs, color = .data$design,
                                  linetype = .data$design)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(~ panel, scales = "free_y") +
    ggplot2::labs(x = "Maximum sample size N",
                  y = "Probability of correct selection (%)") +
    ggplot2::theme_minimal(base_size = 11)
}

#' Overdose-aversion frontier
#'
#' @param loss_res List from \code{\link{loss_sensitivity}}.
#' @return A ggplot object.
#' @export
fig_loss <- function(loss_res) {
  .need_ggplot()
  d <- do.call(rbind, lapply(loss_res, function(r) data.frame(
    sc = r$sc_idx, sc_name = r$sc_name, k_over = r$k_over,
    pcs = r$pcs, over = r$pct_over_sel)))
  d$sc_lab <- factor(paste0("Sc ", d$sc, " ", d$sc_name))
  ggplot2::ggplot(d, ggplot2::aes(.data$over, .data$pcs, group = .data$sc_lab,
                                  color = .data$sc_lab)) +
    ggplot2::geom_path(alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(size = .data$k_over)) +
    ggplot2::labs(x = "Overdose selection (%)", y = "PCS (%)",
                  color = "Scenario", size = "Aversion") +
    ggplot2::theme_minimal(base_size = 11)
}

#' Steepness gradient
#'
#' @param steep_df Data frame from \code{\link{steep_curve_study}}.
#' @return A ggplot object.
#' @export
fig_steep <- function(steep_df) {
  .need_ggplot()
  long <- rbind(
    data.frame(steep_df[c("gap", "design")], metric = "PCS (%)",
               value = steep_df$pcs),
    data.frame(steep_df[c("gap", "design")], metric = "POD (%)",
               value = steep_df$pod))
  ggplot2::ggplot(long, ggplot2::aes(.data$gap, .data$value,
                                     color = .data$design,
                                     shape = .data$design)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::labs(x = "Toxicity gap above the MTD", y = NULL,
                  color = "Design", shape = "Design") +
    ggplot2::theme_minimal(base_size = 11)
}

#' Offset-MTD behaviour across sample sizes
#'
#' @param str_df Data frame from \code{\link{straddle_N_study}}.
#' @return A ggplot object.
#' @export
fig_straddle <- function(str_df) {
  .need_ggplot()
  long <- rbind(
    data.frame(str_df[c("scen", "N", "design")], metric = "PCS (%)",
               value = str_df$pcs),
    data.frame(str_df[c("scen", "N", "design")], metric = "POD (%)",
               value = str_df$pod))
  ggplot2::ggplot(long, ggplot2::aes(.data$N, .data$value,
                                     color = .data$design,
                                     shape = .data$design)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 2) +
    ggplot2::facet_grid(metric ~ scen, scales = "free_y") +
    ggplot2::labs(x = "Maximum sample size N", y = NULL,
                  color = "Design", shape = "Design") +
    ggplot2::theme_minimal(base_size = 10)
}

#' Elimination rule comparison
#'
#' @param elim_df Data frame from \code{\link{elim_rule_study}}.
#' @return A ggplot object.
#' @export
fig_elim <- function(elim_df) {
  .need_ggplot()
  long <- rbind(
    data.frame(elim_df[c("scen", "design", "elim_rule")], metric = "PCS (%)",
               value = elim_df$pcs),
    data.frame(elim_df[c("scen", "design", "elim_rule")],
               metric = "True MTD eliminated (%)", value = elim_df$pct_mtd_elim))
  ggplot2::ggplot(long, ggplot2::aes(.data$design, .data$value,
                                     fill = .data$elim_rule)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75),
                      width = 0.65) +
    ggplot2::facet_grid(metric ~ scen, scales = "free_y") +
    ggplot2::labs(x = NULL, y = NULL, fill = "Elimination rule") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                   legend.position = "bottom")
}
