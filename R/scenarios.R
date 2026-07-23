#' Dose-toxicity scenarios
#'
#' Eight scenarios. The MTD is the dose closest to the target rate with ties
#' broken toward the higher dose. Scenarios 1 and 4 contain exact ties.
#' Scenarios 7 and 8 are non-monotone robustness checks.
#'
#' @format A list of lists, each with \code{name}, \code{pi} and \code{mtd}.
#' @export
scenarios <- list(
  list(name = "Standard sigmoidal",  pi = c(.05,.10,.20,.30,.45,.60), mtd = 4L),
  list(name = "Flat plateau",        pi = c(.05,.10,.22,.25,.28,.35), mtd = 4L),
  list(name = "Abrupt threshold",    pi = c(.05,.08,.12,.15,.25,.60), mtd = 5L),
  list(name = "Early toxicity",      pi = c(.20,.30,.45,.55,.65,.75), mtd = 2L),
  list(name = "Late toxicity",       pi = c(.02,.04,.06,.10,.18,.25), mtd = 6L),
  list(name = "Steep mid-range",     pi = c(.03,.07,.25,.50,.70,.85), mtd = 3L),
  list(name = "Mild non-monotone",   pi = c(.10,.15,.24,.27,.20,.30), mtd = 3L),
  list(name = "Strong non-monotone", pi = c(.08,.15,.25,.30,.18,.40), mtd = 3L)
)

#' Design codes understood by the C++ engine
#'
#' @format A named integer vector.
#' @export
design_codes <- c(adaptive_iso = 0L, adaptive_bern = 1L, boin_bern = 2L,
                  boin = 3L, crm = 4L, mtpi2 = 5L, gboins = 6L)

monotone_idx <- 1:6
primary_designs <- c("adaptive_iso", "boin", "crm", "mtpi2", "gboins")

#' True MTD under a stated tie convention
#'
#' @param pi True DLT probabilities.
#' @param phi Target toxicity rate.
#' @param tie_high Break exact ties toward the higher dose.
#' @return Integer dose index.
#' @examples
#' mtd_of(c(.05, .10, .20, .30, .45, .60))          # 4, higher-dose convention
#' mtd_of(c(.05, .10, .20, .30, .45, .60), tie_high = FALSE)  # 3
#' @export
mtd_of <- function(pi, phi = 0.25, tie_high = TRUE) {
  d <- abs(pi - phi); idx <- which(abs(d - min(d)) < 1e-9)
  if (tie_high) max(idx) else min(idx)
}

#' Doses at least as close to the target as the MTD
#'
#' @inheritParams mtd_of
#' @return Integer vector of dose indices.
#' @export
acceptable_set <- function(pi, phi = 0.25) {
  m <- mtd_of(pi, phi)
  which(abs(pi - phi) <= abs(pi[m] - phi) + 1e-9)
}

#' Does a scenario contain an exact tie?
#'
#' @inheritParams mtd_of
#' @return Logical.
#' @export
is_tie <- function(pi, phi = 0.25) {
  sum(abs(abs(pi - phi) - min(abs(pi - phi))) < 1e-9) > 1
}

#' Built-in scenarios whose true MTD lies above the target
#'
#' These are the scenarios affected by the elimination behaviour described in the
#' manuscript, since a target-based rule tests a hypothesis that is false at
#' such an MTD.
#'
#' @param phi Target toxicity rate.
#' @return Integer vector of scenario indices.
#' @export
affected_scenarios <- function(phi = 0.25) {
  which(vapply(seq_along(scenarios), function(i) {
    s <- scenarios[[i]]; s$pi[s$mtd] > phi + 1e-9
  }, logical(1)))
}
