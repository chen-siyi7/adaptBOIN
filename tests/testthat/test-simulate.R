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
  pm <- prior_mean_curve(K, J)
  d  <- (seq_len(J) - 1) / (J - 1)
  expect_equal(as.numeric(pm$prior_pi), (K * d + 1) / (K + 2), tolerance = 1e-8)
})
