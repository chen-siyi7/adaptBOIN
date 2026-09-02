test_that("the tie convention is honoured in both directions", {
  pi <- c(.05, .10, .20, .30, .45, .60)   # d3 and d4 both 0.05 from 0.25
  expect_equal(mtd_of(pi, tie_high = TRUE), 4L)
  expect_equal(mtd_of(pi, tie_high = FALSE), 3L)
})

test_that("exactly two built-in scenarios contain exact ties", {
  ties <- which(vapply(scenarios, function(s) is_tie(s$pi), logical(1)))
  expect_equal(ties, c(1L, 4L))
})

test_that("stored MTDs match the higher-dose convention", {
  for (s in scenarios) expect_equal(mtd_of(s$pi), s$mtd)
})

test_that("acceptable sets contain the MTD and both tied doses", {
  expect_equal(acceptable_set(scenarios[[1]]$pi), c(3L, 4L))
  expect_equal(acceptable_set(scenarios[[4]]$pi), c(1L, 2L))
  expect_equal(acceptable_set(scenarios[[6]]$pi), 3L)
})

test_that("scenarios with an MTD above target are identified", {
  expect_equal(affected_scenarios(), c(1L, 4L))
})

test_that("the high-MTD diagnostic varies local steepness at d5 and d6", {
  expect_equal(vapply(high_mtd_scenarios, `[[`, integer(1), "mtd"),
               c(5L, 5L, 6L, 6L))
  expect_equal(vapply(high_mtd_scenarios, function(s) mtd_of(s$pi), integer(1)),
               c(5L, 5L, 6L, 6L))
  expect_equal(vapply(high_mtd_scenarios, `[[`, numeric(1), "gap"),
               c(.05, .10, .07, .13))
})

test_that("the high-MTD diagnostic covers both proposed configurations", {
  z <- high_mtd_steepness_study(n_sim = 2L)
  expect_setequal(unique(z$design),
                  c("adaptive_iso", "adaptive_bern", "boin_bern", "boin",
                    "aboin", "crm", "mtpi2", "gboins"))
  expect_equal(nrow(z), 32L)
})

test_that("gBOINS calibration sensitivity brackets the primary interpolation", {
  z <- gboins_calibration_sensitivity(scenarios_idx = 1:2, n_sim = 2L)
  expect_setequal(unique(z$calibration),
                  c("published_020", "midpoint_025", "published_030"))
  expect_equal(nrow(z), 6L)
  expect_equal(sort(unique(z$multiplier)), c(1.05, 1.075, 1.10))
  p <- adapt_params()
  expect_equal(p$gb_c1, z$c1[z$calibration == "midpoint_025"][1])
  expect_equal(p$gb_c2, z$c2[z$calibration == "midpoint_025"][1])
})

test_that("adapt_params_original restores the previous conventions", {
  p <- adapt_params_original()
  expect_false(p$tie_high)
  expect_equal(p$phi_elim, 0.35)
  expect_equal(p$elim_a0, 0.5)
})
