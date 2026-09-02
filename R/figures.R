.need_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for the figure functions. ",
         "Install it with install.packages(\"ggplot2\").", call. = FALSE)
}

# Restrained Nature-style system used for every manuscript figure. The palette
# follows the muted NPG family; categorical distinctions are also encoded by
# shape or line type so that figures remain interpretable in grayscale.
.nature_palette <- c(
  blue = "#3C5488", vermillion = "#E64B35", green = "#00A087",
  purple = "#8491B4", orange = "#F39B7F", sky = "#4DBBD5",
  black = "#202124"
)

.design_labels <- c(
  adaptive_iso = "Adaptive+Iso", boin = "BOIN", aboin = "aBOIN",
  crm = "CRM", mtpi2 = "mTPI-2", gboins = "gBOINS"
)

.design_colors <- c(
  "Adaptive+Iso" = .nature_palette[["blue"]],
  "BOIN" = .nature_palette[["black"]],
  "aBOIN" = .nature_palette[["vermillion"]],
  "CRM" = .nature_palette[["sky"]],
  "mTPI-2" = .nature_palette[["purple"]],
  "gBOINS" = .nature_palette[["green"]]
)

.design_shapes <- c(
  "Adaptive+Iso" = 16, "BOIN" = 17, "aBOIN" = 15,
  "CRM" = 3, "mTPI-2" = 4, "gBOINS" = 18
)

.nature_theme <- function(base_size = 13) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.45, colour = "#202124"),
      axis.ticks = ggplot2::element_line(linewidth = 0.45, colour = "#202124"),
      axis.ticks.length = grid::unit(1.8, "mm"),
      axis.text = ggplot2::element_text(colour = "#202124"),
      axis.title = ggplot2::element_text(colour = "#202124"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", colour = "#202124"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      legend.key.width = grid::unit(5, "mm"),
      legend.key.height = grid::unit(5, "mm"),
      legend.spacing.x = grid::unit(2, "mm"),
      legend.spacing.y = grid::unit(1.5, "mm"),
      plot.margin = ggplot2::margin(5, 8, 5, 5)
    )
}

#' Importance-sampling effective sample size and fallback rate
#'
#' @param params Parameter list.
#' @param n_sim Trials per scenario.
#' @param df Optional summary data frame from \code{\link{summarise_results}};
#'   supplying it avoids rerunning the BOIN+Bern simulations.
#' @return A ggplot object.
#' @export
fig_ess <- function(params = adapt_params(), n_sim = 2000L, df = NULL) {
  .need_ggplot()
  if (is.null(df)) {
    atbl <- adapt_table(params)
    d <- do.call(rbind, lapply(seq_along(scenarios), function(sc) {
      r <- run_scenario(sc, "boin_bern", params, n_sim, atbl)
      data.frame(sc = sc, ess_med = r$ess_median, ess_p10 = r$ess_p10,
                 fallback = r$pct_fallback)
    }))
  } else {
    d <- df[df$design == "boin_bern", c("sc", "ess_med", "ess_p10", "pct_fb")]
    names(d)[4] <- "fallback"
  }
  d$sc_lab <- factor(paste0("Sc ", d$sc), levels = paste0("Sc ", d$sc))
  scale <- max(d$ess_med, na.rm = TRUE) / 100
  ggplot2::ggplot(d, ggplot2::aes(.data$sc_lab)) +
    ggplot2::geom_col(ggplot2::aes(y = .data$ess_med),
                      fill = .nature_palette[["blue"]],
                      width = 0.6) +
    ggplot2::geom_col(ggplot2::aes(y = .data$ess_p10),
                      fill = .nature_palette[["sky"]],
                      width = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = .data$fallback * scale, group = 1),
                       color = .nature_palette[["vermillion"]],
                       linewidth = 0.9) +
    ggplot2::geom_point(ggplot2::aes(y = .data$fallback * scale),
                        color = .nature_palette[["vermillion"]], size = 3.0,
                        stroke = 0.8) +
    ggplot2::geom_hline(yintercept = 100, linetype = "dashed", linewidth = 0.55) +
    ggplot2::scale_y_continuous(
      name = "IS effective sample size",
      sec.axis = ggplot2::sec_axis(~ . / scale, name = "% fallback")) +
    ggplot2::labs(x = NULL) +
    .nature_theme(base_size = 14)
}

#' Accuracy by design and scenario
#'
#' @param df Data frame from \code{\link{summarise_results}}.
#' @return A ggplot object.
#' @export
fig_pcs <- function(df) {
  .need_ggplot()
  d <- df[df$design %in% primary_designs, ]
  d$design_lab <- factor(unname(.design_labels[d$design]),
                         levels = unname(.design_labels[primary_designs]))
  d$sc_lab <- factor(paste0("Sc ", d$sc, " ", d$sc_name),
                     levels = rev(unique(paste0("Sc ", d$sc, " ", d$sc_name))))
  # Use fixed display lanes within each scenario. Horizontal position remains
  # the observed PCS; the small vertical offsets only prevent overplotting.
  design_offset <- stats::setNames(
    seq(0.42, -0.42, length.out = length(levels(d$design_lab))),
    levels(d$design_lab)
  )
  d$display_y <- as.numeric(d$sc_lab) +
    unname(design_offset[as.character(d$design_lab)])
  scenario_y <- seq_along(levels(d$sc_lab))

  ggplot2::ggplot(d, ggplot2::aes(.data$pcs, .data$display_y,
                                  color = .data$design_lab,
                                  shape = .data$design_lab)) +
    ggplot2::geom_hline(
      yintercept = scenario_y, colour = "#D9D9D9", linewidth = 0.35
    ) +
    ggplot2::geom_point(size = 3.8, stroke = 1.05) +
    ggplot2::scale_color_manual(values = .design_colors, drop = FALSE) +
    ggplot2::scale_shape_manual(values = .design_shapes, drop = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = seq(20, 80, by = 10), limits = c(20, 80),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      breaks = scenario_y, labels = levels(d$sc_lab),
      limits = c(0.35, length(scenario_y) + 0.65),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(x = "Probability of correct selection (%)", y = NULL,
                  color = "Design", shape = "Design") +
    .nature_theme(base_size = 15) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 13.5),
      axis.title.x = ggplot2::element_text(size = 15),
      panel.grid.major.x = ggplot2::element_line(
        colour = "#ECECEC", linewidth = 0.35
      ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = ggplot2::element_text(size = 13.5),
      legend.title = ggplot2::element_text(size = 13.5, face = "bold")
    )
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
               component = "Boundaries (Adaptive+Iso - BOIN)", value = m$ai - m$boin),
    data.frame(sc = m$sc, sc_name = m$sc_name,
               component = "Bernstein (BOIN+Bern - BOIN)", value = m$bb - m$boin))
  long$sc_lab <- factor(paste0("Sc ", long$sc, " ", long$sc_name),
                        levels = rev(unique(paste0("Sc ", long$sc, " ",
                                                   long$sc_name))))
  ggplot2::ggplot(long, ggplot2::aes(.data$value, .data$sc_lab,
                                     fill = .data$component)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
                      width = 0.6) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.55) +
    ggplot2::scale_fill_manual(values = c(
      "Boundaries (Adaptive+Iso - BOIN)" = .nature_palette[["blue"]],
      "Bernstein (BOIN+Bern - BOIN)" = .nature_palette[["vermillion"]]
    )) +
    ggplot2::labs(x = "PCS difference vs BOIN (percentage points)", y = NULL,
                  fill = NULL) +
    .nature_theme() +
    ggplot2::theme(legend.position = "bottom")
}

#' Accuracy against sample size
#'
#' @param n_res List from \code{\link{n_sensitivity}}.
#' @param scs Scenario indices to show.
#' @return A ggplot object.
#' @export
fig_nsens <- function(n_res, scs = c(6, 9, 1, 5)) {
  .need_ggplot()
  d <- do.call(rbind, lapply(n_res, function(r)
    data.frame(sc = r$sc_idx, sc_name = r$sc_name, design = r$design,
               N = r$N_val, pcs = r$pcs)))
  d <- d[d$sc %in% scs, ]
  d$design_lab <- factor(unname(.design_labels[d$design]),
                         levels = c("Adaptive+Iso", "BOIN"))
  d$panel <- factor(paste0("Sc ", d$sc, " ", d$sc_name),
                    levels = unique(paste0("Sc ", d$sc, " ", d$sc_name)))
  ggplot2::ggplot(d, ggplot2::aes(.data$N, .data$pcs,
                                  color = .data$design_lab,
                                  linetype = .data$design_lab,
                                  shape = .data$design_lab)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 3.0, stroke = 0.8) +
    ggplot2::scale_color_manual(values = .design_colors) +
    ggplot2::scale_shape_manual(values = .design_shapes) +
    ggplot2::scale_linetype_manual(values = c("Adaptive+Iso" = "solid",
                                              "BOIN" = "dashed")) +
    ggplot2::facet_wrap(~ panel, scales = "free_y") +
    ggplot2::labs(x = "Maximum sample size N",
                  y = "Probability of correct selection (%)") +
    ggplot2::labs(color = "Design", linetype = "Design", shape = "Design") +
    .nature_theme(base_size = 14)
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
    ggplot2::geom_path(linewidth = 0.9) +
    ggplot2::geom_point(ggplot2::aes(size = .data$k_over), stroke = 0.9) +
    ggplot2::scale_color_manual(values = stats::setNames(
      unname(.nature_palette[c("blue", "vermillion", "green", "purple",
                               "orange", "sky", "black")]), levels(d$sc_lab))) +
    ggplot2::scale_size_continuous(range = c(2.8, 4.8), breaks = c(1, 1.5, 2)) +
    ggplot2::labs(x = "Overdose selection (%)", y = "PCS (%)",
                  color = "Scenario", size = expression(kappa)) +
    .nature_theme(base_size = 14)
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
  long$design_lab <- factor(unname(.design_labels[long$design]),
                            levels = c("Adaptive+Iso", "BOIN", "gBOINS"))
  ggplot2::ggplot(long, ggplot2::aes(.data$gap, .data$value,
                                     color = .data$design_lab,
                                     shape = .data$design_lab)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 3.0, stroke = 0.8) +
    ggplot2::scale_color_manual(values = .design_colors) +
    ggplot2::scale_shape_manual(values = .design_shapes) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::labs(x = "Toxicity gap above the MTD", y = NULL,
                  color = "Design", shape = "Design") +
    .nature_theme(base_size = 14)
}

#' Offset-MTD behaviour across sample sizes
#'
#' @param str_df Data frame from \code{\link{straddle_N_study}}.
#' @return A ggplot object.
#' @export
fig_straddle <- function(str_df) {
  .need_ggplot()
  panel_labels <- c(
    "MTD above target (0.28)" = "Above target (0.28)",
    "MTD below target (0.22)" = "Below target (0.22)",
    "Symmetric straddle (0.20/0.30)" = "Straddle (0.20/0.30)"
  )
  long <- rbind(
    data.frame(str_df[c("scen", "N", "design")], metric = "PCS (%)",
               value = str_df$pcs),
    data.frame(str_df[c("scen", "N", "design")], metric = "POD (%)",
               value = str_df$pod))
  long$design_lab <- factor(unname(.design_labels[long$design]),
                            levels = c("Adaptive+Iso", "BOIN", "gBOINS"))
  ggplot2::ggplot(long, ggplot2::aes(.data$N, .data$value,
                                     color = .data$design_lab,
                                     shape = .data$design_lab)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 3.0, stroke = 0.8) +
    ggplot2::scale_color_manual(values = .design_colors) +
    ggplot2::scale_shape_manual(values = .design_shapes) +
    ggplot2::facet_grid(
      metric ~ scen, scales = "free_y",
      labeller = ggplot2::labeller(scen = panel_labels)
    ) +
    ggplot2::labs(x = "Maximum sample size N", y = NULL,
                  color = "Design", shape = "Design") +
    .nature_theme(base_size = 18)
}

#' Elimination rule comparison
#'
#' @param elim_df Data frame from \code{\link{elim_rule_study}}.
#' @return A ggplot object.
#' @export
fig_elim <- function(elim_df) {
  .need_ggplot()
  panel_labels <- c(
    "Offset MTD above target (0.28)" = "Above target (0.28)",
    "Sc1 Standard sigmoidal" = "Sc 1 Standard sigmoid",
    "Sc4 Early toxicity" = "Sc 4 Early toxicity"
  )
  long <- rbind(
    data.frame(elim_df[c("scen", "design", "elim_rule")], metric = "PCS (%)",
               value = elim_df$pcs),
    data.frame(elim_df[c("scen", "design", "elim_rule")],
               metric = "MTD eliminated (%)", value = elim_df$pct_mtd_elim))
  long$design_lab <- factor(unname(.design_labels[long$design]),
                            levels = c("Adaptive+Iso", "BOIN", "gBOINS"))
  long$metric <- factor(long$metric,
                        levels = c("PCS (%)", "MTD eliminated (%)"))
  long$elim_lab <- factor(long$elim_rule,
                          levels = c("standard", "equiv_phi2"),
                          labels = c("Standard", "Equivalence"))
  ggplot2::ggplot(long, ggplot2::aes(.data$design_lab, .data$value,
                                     fill = .data$elim_lab)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75),
                      width = 0.65) +
    ggplot2::facet_grid(
      metric ~ scen, scales = "free_y",
      labeller = ggplot2::labeller(scen = panel_labels)
    ) +
    ggplot2::scale_fill_manual(values = c(
      "Standard" = .nature_palette[["vermillion"]],
      "Equivalence" = .nature_palette[["blue"]])) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Elimination rule") +
    .nature_theme(base_size = 18) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                   legend.position = "bottom")
}
