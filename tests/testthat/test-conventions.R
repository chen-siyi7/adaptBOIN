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

test_that("adapt_params_original restores the previous conventions", {
  p <- adapt_params_original()
  expect_false(p$tie_high)
  expect_equal(p$phi_elim, 0.35)
  expect_equal(p$elim_a0, 0.5)
})
