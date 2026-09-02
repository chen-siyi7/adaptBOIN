#!/usr/bin/env Rscript

# Verify every numerical coordinate in the manuscript plots against the
# archived CSV outputs used to build the figures.

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

rd <- function(name) utils::read.csv(file.path(out, name), check.names = FALSE)
check_values <- function(label, got, expected, tol = 1e-8) {
  got <- sort(as.numeric(got))
  expected <- sort(as.numeric(expected))
  if (length(got) != length(expected) ||
      any(abs(got - expected) > tol))
    stop("Mismatch in ", label)
  message("PASS ", label, " (", length(expected), " coordinates)")
}
check_xy <- function(label, got, x, y, digits = 8) {
  key <- function(a, b) sprintf(paste0("%.", digits, "f|%.", digits, "f"), a, b)
  observed <- sort(key(got$x, got$y))
  expected <- sort(key(x, y))
  if (!identical(observed, expected)) stop("Mismatch in ", label)
  message("PASS ", label, " (", length(expected), " points)")
}

sim <- rd("sim_results.csv")
ns <- rd("n_sensitivity.csv")
loss <- rd("loss_sensitivity.csv")
steep <- rd("study_steep_curve.csv")
offset <- rd("study_offset_mtd_byN.csv")
elim <- rd("study_elim_rule.csv")
scenario_name <- tapply(sim$sc_name, sim$sc, function(x) unique(x)[1])
name_for <- function(sc) unname(scenario_name[as.character(sc)])

primary <- c("adaptive_iso", "boin", "aboin", "crm", "mtpi2", "gboins")
monotone <- c(1:6, 9)

p <- fig_pcs(sim)
b <- ggplot2::ggplot_build(p)$data[[2]]
check_values("fig_pcs PCS", b$x, sim$pcs[sim$design %in% primary])

base <- sim[sim$design == "boin" & sim$sc %in% monotone, c("sc", "pcs")]
ai <- sim[sim$design == "adaptive_iso" & sim$sc %in% monotone, c("sc", "pcs")]
bb <- sim[sim$design == "boin_bern" & sim$sc %in% monotone, c("sc", "pcs")]
base <- base[order(base$sc), ]; ai <- ai[order(ai$sc), ]; bb <- bb[order(bb$sc), ]
p <- fig_decomposition(sim)
b <- ggplot2::ggplot_build(p)$data[[1]]
check_values("fig_decomposition differences", b$x,
             c(ai$pcs - base$pcs, bb$pcs - base$pcs))

ess <- sim[sim$design == "boin_bern", ]
p <- fig_ess(params = list(), df = sim)
b <- ggplot2::ggplot_build(p)$data
check_values("fig_ess medians", b[[1]]$y, ess$ess_med)
check_values("fig_ess tenth percentiles", b[[2]]$y, ess$ess_p10)
ess_scale <- max(ess$ess_med, na.rm = TRUE) / 100
check_values("fig_ess fallback", b[[4]]$y / ess_scale, ess$pct_fb)

ns_list <- lapply(seq_len(nrow(ns)), function(i) list(
  sc_idx = ns$sc[i], sc_name = name_for(ns$sc[i]), design = ns$design[i],
  N_val = ns$N[i], pcs = ns$pcs[i]
))
keep_ns <- ns$sc %in% c(6, 9, 1, 5)
p <- fig_nsens(ns_list)
b <- ggplot2::ggplot_build(p)$data[[2]]
check_xy("fig_nsensitivity", b, ns$N[keep_ns], ns$pcs[keep_ns])

loss_list <- lapply(seq_len(nrow(loss)), function(i) list(
  sc_idx = loss$sc[i], sc_name = name_for(loss$sc[i]),
  k_over = loss$k_over[i], pcs = loss$pcs[i],
  pct_over_sel = loss$pct_over[i]
))
p <- fig_loss(loss_list)
b <- ggplot2::ggplot_build(p)$data[[2]]
check_xy("fig_loss_tradeoff", b, loss$pct_over, loss$pcs)

p <- fig_steep(steep)
b <- ggplot2::ggplot_build(p)$data[[2]]
check_xy("fig_steepness", b, c(steep$gap, steep$gap),
         c(steep$pcs, steep$pod))

p <- fig_straddle(offset)
b <- ggplot2::ggplot_build(p)$data[[2]]
check_xy("fig_offset_mtd", b, c(offset$N, offset$N),
         c(offset$pcs, offset$pod))

elim90 <- elim[elim$N == 90, ]
p <- fig_elim(elim90)
b <- ggplot2::ggplot_build(p)$data[[1]]
check_values("fig_elim_rule", b$y,
             c(elim90$pcs, elim90$pct_mtd_elim))

message("All plotted numerical coordinates match the archived outputs.")
