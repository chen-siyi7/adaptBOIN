#' Build the adaptive decision table
#'
#' @param params Parameter list from \code{\link{adapt_params}}.
#' @param check If true, warn when the loss matrix violates the no-DLT
#'   escalation condition.
#' @return A list of action codes indexed by cohort size; see
#'   the internal C++ table builder.
#' @examples
#' adapt_table()[["30"]]
#' @export
adapt_table <- function(params = adapt_params(), check = TRUE) {
  tb <- build_adap_tbl(params$phi1, params$phi2, params$a0, params$b0,
                       params$N, params$m, params$loss)
  if (isTRUE(check)) {
    # The no-DLT condition of the count-threshold representation: at zero events
    # the Bayes action must be to escalate. Checked on the built table, so it
    # holds for any loss matrix, not only the k_over parameterisation.
    bad <- names(tb)[vapply(tb, function(a) as.integer(a)[1] != 0L, logical(1))]
    if (length(bad))
      warning("The loss matrix gives a design that does not escalate after a ",
              "cohort with no events, at n = ", paste(bad, collapse = ", "),
              ". This violates the no-DLT condition and is not a usable design. ",
              "See kappa_max() for the largest admissible overdose aversion.",
              call. = FALSE)
  }
  tb
}

#' Effective rate boundaries of the adaptive table
#'
#' Converts the integer count thresholds into rates, alongside BOIN's fixed
#' cutoffs. A boundary of \code{NA} means no stay region exists at that cohort
#' size.
#'
#' @param params Parameter list.
#' @param ns Cohort sizes to report.
#' @return A data frame.
#' @examples
#' boundary_rates()
#' @export
boundary_rates <- function(params = adapt_params(),
                           ns = c(3, 6, 9, 12, 15, 30)) {
  atbl <- adapt_table(params)
  do.call(rbind, lapply(ns, function(n) {
    key <- as.character(n)
    if (is.null(atbl[[key]])) return(NULL)
    a <- as.integer(atbl[[key]])
    s_lo <- which(a >= 1)[1] - 1L      # first non-escalate
    hi   <- which(a <= 1)              # last non-de-escalate
    s_hi <- if (length(hi)) max(hi) - 1L else NA_integer_
    data.frame(n = n,
               adaptive_low  = if (is.na(s_lo)) NA_real_ else round(s_lo / n, 3),
               adaptive_high = if (is.na(s_hi)) NA_real_ else round(s_hi / n, 3),
               boin_low = params$lam1, boin_high = params$lam2)
  }))
}

#' Effective rate boundaries across overdose-aversion values
#'
#' @param k_overs Overdose-aversion multipliers.
#' @param ns Cohort sizes to report.
#' @param params Parameter list.
#' @return A data frame with one row per \code{(k_over, n)} pair.
#' @examples
#' boundaries_by_loss(k_overs = c(1, 2))
#' @export
boundaries_by_loss <- function(k_overs = c(1.0, 1.5, 2.0, 3.0),
                               ns = c(3, 6, 9, 12, 15, 30),
                               params = adapt_params()) {
  do.call(rbind, lapply(k_overs, function(k) {
    pp <- params; pp$loss <- make_loss(k_over = k)
    b <- boundary_rates(pp, ns)
    data.frame(k_over = k, n = b$n,
               lambda_low = b$adaptive_low, lambda_high = b$adaptive_high)
  }))
}

#' aBOIN adaptive boundaries
#'
#' Computes the no-history adaptive BOIN boundaries of Li and Pan (2020).
#' The first \code{ab_N0} patients at a dose use the fixed BOIN boundaries;
#' thereafter the two point alternatives shrink toward the target according to
#' the published acceleration factors.
#'
#' @param ns Cumulative sample sizes at a dose.
#' @param params Parameter list from \code{\link{adapt_params}}.
#' @return A data frame with the escalation and de-escalation boundaries.
#' @export
aboin_boundaries <- function(ns = c(3, 6, 9, 12, 15, 30),
                             params = adapt_params()) {
  do.call(rbind, lapply(ns, function(n) {
    b <- aboin_bounds_cpp(params$phi_tgt, as.integer(n),
                          params$ab_delta1, params$ab_delta2,
                          params$ab_g1, params$ab_g2, params$ab_N0,
                          params$lam1, params$lam2)
    data.frame(n = n, lambda_e = unname(b[1]), lambda_d = unname(b[2]))
  }))
}
