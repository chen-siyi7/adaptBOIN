#!/usr/bin/env Rscript

# Regenerate the manuscript figures from the archived 2,000-trial outputs.
# This changes only presentation; it does not rerun any simulation.

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
pkg <- if (file.exists(file.path(root, "DESCRIPTION"))) {
  root
} else {
  file.path(root, "adaptBOIN")
}
out <- file.path(pkg, "revision_output")

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("ggplot2 is required")
primary_designs <- c("adaptive_iso", "boin", "aboin", "crm", "mtpi2", "gboins")
monotone_idx <- c(1:6, 9L)
source(file.path(pkg, "R", "figures.R"), local = FALSE)

pdf_device <- function(filename, width, height, ...) {
  if (capabilities("aqua")) {
    grDevices::quartz(file = filename, type = "pdf", width = width,
                      height = height, family = "Arial")
  } else {
    grDevices::cairo_pdf(filename = filename, width = width, height = height,
                         family = "sans", bg = "white")
  }
}

rd <- function(name) utils::read.csv(file.path(out, name), check.names = FALSE)
sim <- rd("sim_results.csv")
ns <- rd("n_sensitivity.csv")
loss <- rd("loss_sensitivity.csv")
steep <- rd("study_steep_curve.csv")
offset <- rd("study_offset_mtd_byN.csv")
elim <- rd("study_elim_rule.csv")

scenario_name <- tapply(sim$sc_name, sim$sc, function(x) unique(x)[1])
name_for <- function(sc) unname(scenario_name[as.character(sc)])

ns_list <- lapply(seq_len(nrow(ns)), function(i) list(
  sc_idx = ns$sc[i], sc_name = name_for(ns$sc[i]), design = ns$design[i],
  N_val = ns$N[i], pcs = ns$pcs[i]
))
loss_list <- lapply(seq_len(nrow(loss)), function(i) list(
  sc_idx = loss$sc[i], sc_name = name_for(loss$sc[i]),
  k_over = loss$k_over[i], pcs = loss$pcs[i],
  pct_over_sel = loss$pct_over[i]
))

plots <- list(
  fig_ess.pdf = fig_ess(params = list(), df = sim),
  fig_pcs.pdf = fig_pcs(sim),
  fig_decomposition.pdf = fig_decomposition(sim),
  fig_nsensitivity.pdf = fig_nsens(ns_list),
  fig_loss_tradeoff.pdf = fig_loss(loss_list),
  fig_steepness.pdf = fig_steep(steep),
  fig_offset_mtd.pdf = fig_straddle(offset),
  fig_elim_rule.pdf = fig_elim(elim[elim$N == 90, ])
)

dims <- list(
  fig_ess.pdf = c(7, 5), fig_pcs.pdf = c(7.4, 7.5),
  fig_decomposition.pdf = c(7, 5), fig_nsensitivity.pdf = c(8, 6),
  fig_loss_tradeoff.pdf = c(7, 5), fig_steepness.pdf = c(8, 5),
  fig_offset_mtd.pdf = c(10, 6), fig_elim_rule.pdf = c(10, 6)
)

for (name in names(plots)) {
  ggplot2::ggsave(
    filename = file.path(out, name), plot = plots[[name]],
    width = dims[[name]][1], height = dims[[name]][2], units = "in",
    bg = "white", device = pdf_device
  )
  message("WROTE ", name)
}

message("All figures regenerated from archived outputs in ", out)
