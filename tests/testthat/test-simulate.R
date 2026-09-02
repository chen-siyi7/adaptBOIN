test_that("a small run returns coherent operating characteristics", {
  r <- run_scenario(6, "boin", n_sim = 200L)
  expect_true(r$pcs >= 0 && r$pcs <= 100)
  expect_equal(r$true_mtd, 3L)
  expect_equal(length(as.numeric(r$mtd_dist)), adapt_params()$J + 1L)
  expect_equal(sum(as.numeric(r$mtd_dist)), 200)
})

test_that("selection percentages sum to 100", {
  r  <- run_scenario(1, "adaptive_iso", n_sim = 200L)
  sd <- selection_dist(r)
  expect_equal(sum(sd$dose) + sd$no_rec, 100, tolerance = 1e-6)
})

test_that("results are reproducible for a fixed seed", {
  a <- run_scenario(2, "boin", n_sim = 100L)
  b <- run_scenario(2, "boin", n_sim = 100L)
  expect_equal(a$pcs, b$pcs)
})

test_that("the high-MTD diagnostic returns one row per curve and design", {
  z <- high_mtd_steepness_study(
    n_sim = 20L, designs = c("adaptive_iso", "boin"))
  expect_equal(nrow(z), 8L)
  expect_equal(sort(unique(z$mtd)), c(5L, 6L))
  expect_true(all(z$pcs >= 0 & z$pcs <= 100))
})

test_that("the engine reports elimination of the true MTD", {
  r <- run_scenario(1, "boin", n_sim = 200L)
  expect_false(is.null(r$pct_mtd_elim))
  expect_true(r$pct_mtd_elim >= 0 && r$pct_mtd_elim <= 100)
})

test_that("every design code runs", {
  for (dn in names(design_codes)) {
    r <- run_scenario(6, dn, n_sim = 50L)
    expect_true(is.finite(r$pcs), info = dn)
  }
})

test_that("the Bernstein prior mean curve matches its closed form", {
  K <- 5L; J <- 6L
  pm <- adaptBOIN:::prior_mean_curve(K, J)
  d  <- (seq_len(J) - 1) / (J - 1)
  expect_equal(as.numeric(pm$prior_pi), (K * d + 1) / (K + 2), tolerance = 1e-8)
})

test_that("mTPI-2 uses equal-width interval keys and isotonic selection", {
  p <- adapt_params()
  got <- adaptBOIN:::mtpi2_actions_cpp(p$N, p$m, p$phi1, p$phi2,
                                       p$a0, p$b0)
  ref_action <- function(s, n) {
    cuts <- unique(c(0, seq(p$phi1, 0, by = -(p$phi2 - p$phi1)),
                     p$phi1, p$phi2,
                     seq(p$phi2, 1, by = p$phi2 - p$phi1), 1))
    cuts <- sort(cuts[cuts >= 0 & cuts <= 1])
    mass <- diff(stats::pbeta(cuts, p$a0 + s, p$b0 + n - s))
    upm <- mass / diff(cuts)
    mid <- (cuts[-1] + cuts[-length(cuts)]) / 2
    key <- which.max(upm)
    if (mid[key] < p$phi1) 0L else if (mid[key] > p$phi2) 2L else 1L
  }
  for (ii in seq_len(p$N / p$m)) {
    n <- ii * p$m
    expect_equal(as.integer(got[ii, seq_len(n + 1L)]),
                 vapply(0:n, ref_action, integer(1), n = n))
  }

  tr <- adaptBOIN:::one_trial(scenarios[[6]]$pi, design_codes[["mtpi2"]],
                             adapt_table(p), p, 12345L)
  expect_true(tr$mtd_sel %in% 1:p$J || tr$mtd_sel == 0L)
})
