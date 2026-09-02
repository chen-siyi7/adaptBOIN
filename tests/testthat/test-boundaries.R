test_that("the decision table is a valid E-then-S-then-D threshold rule", {
  for (k in c(1, 1.5, 2)) {
    p <- adapt_params(); p$loss <- make_loss(k_over = k)
    tb <- adapt_table(p)
    for (nm in names(tb)) {
      a <- as.integer(tb[[nm]])
      expect_true(all(diff(a) >= 0),
                  info = paste("non-monotone actions at k =", k, "n =", nm))
    }
  }
})

test_that("the symmetric default reproduces the published boundary table", {
  b <- boundary_rates()
  expect_equal(b$adaptive_high[b$n == 30], 0.333, tolerance = 1e-3)
  expect_equal(b$adaptive_low[b$n == 30],  0.167, tolerance = 1e-3)
  expect_equal(b$adaptive_high[b$n == 3],  0.333, tolerance = 1e-3)
})

test_that("overdose aversion tightens the de-escalation boundary", {
  b <- boundaries_by_loss(k_overs = c(1, 1.5, 3), ns = 30)
  hi <- b$lambda_high[order(b$k_over)]
  expect_true(all(diff(hi) <= 0))
  expect_lt(hi[length(hi)], hi[1])
})

test_that("k_over = 3 is inadmissible because escalation stops at s = 0", {
  p <- adapt_params(); p$loss <- make_loss(k_over = 3)
  a3 <- as.integer(adapt_table(p)[["3"]])
  expect_gt(a3[1], 0)   # action at s = 0 is not Escalate
})

test_that("aBOIN follows the published no-history shrinking formula", {
  p <- adapt_params()
  b <- aboin_boundaries(c(3, 6, 9, 30), p)
  expect_equal(b$lambda_e[1:2], rep(p$lam1, 2))
  expect_equal(b$lambda_d[1:2], rep(p$lam2, 2))

  n <- 9
  phi1 <- p$phi_tgt - p$ab_delta1 / n^(p$ab_g1 / 2)
  phi2 <- p$phi_tgt + p$ab_delta2 / n^(p$ab_g2 / 2)
  le <- log((1 - phi1) / (1 - p$phi_tgt)) /
    log(p$phi_tgt * (1 - phi1) / ((1 - p$phi_tgt) * phi1))
  ld <- log((1 - p$phi_tgt) / (1 - phi2)) /
    log(phi2 * (1 - p$phi_tgt) / ((1 - phi2) * p$phi_tgt))
  expect_equal(b$lambda_e[b$n == n], le, tolerance = 1e-12)
  expect_equal(b$lambda_d[b$n == n], ld, tolerance = 1e-12)
  expect_true(all(diff(b$lambda_e) >= 0))
  expect_true(all(diff(b$lambda_d) <= 0))
})
