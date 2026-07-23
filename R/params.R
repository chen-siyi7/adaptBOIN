#' Construct a three-state loss matrix
#'
#' Rows are actions (escalate, stay, de-escalate); columns are states
#' (under-dosing, target, over-toxic). The diagonal is zero.
#'
#' @param mild Penalty for an adjacent misclassification.
#' @param sev Penalty for the most discordant decisions.
#' @param k_over Overdose-aversion multiplier applied to the two entries that
#'   penalise failing to de-escalate from an overtoxic dose. \code{k_over = 1}
#'   gives the symmetric matrix; larger values encode the judgement that
#'   overdosing is the more serious error. Values of 3 or above are inadmissible,
#'   because the rule then fails to escalate after a cohort with no events.
#' @return A 3 by 3 matrix with dimnames.
#' @examples
#' make_loss()
#' make_loss(k_over = 1.5)
#' @export
make_loss <- function(mild = 0.5, sev = 1.5, k_over = 1.0) {
  matrix(c(
    0,    mild, k_over * sev,
    mild, 0,    k_over * mild,
    sev,  mild, 0
  ), nrow = 3, byrow = TRUE,
  dimnames = list(c("E", "S", "D"), c("under", "target", "over")))
}

#' Default design parameters
#'
#' @return A list of design settings. Fields of particular note are
#'   \code{tie_high} (higher-dose tie convention), \code{phi_elim},
#'   \code{elim_a0} and \code{elim_b0} (elimination rule, defaulting to the
#'   standard BOIN rule), \code{loss} (the decision loss), and \code{gb_c1},
#'   \code{gb_c2}, \code{gb_eps}, \code{gb_N0} (gBOINS calibration).
#' @examples
#' p <- adapt_params()
#' p$phi_elim
#' @export
adapt_params <- function() {
  list(
    phi_tgt = 0.25, phi1 = 0.15, phi2 = 0.35, rho = 0.95,
    a0 = 0.5, b0 = 0.5, N = 30L, m = 3L, J = 6L,
    lam1 = 0.197, lam2 = 0.298,
    K = 5L, alpha_d = 1.0, M_IS = 2000L, ESS_thr = 100.0,
    crm_skel = c(0.05, 0.12, 0.25, 0.40, 0.55, 0.70),
    crm_sd = 1.0, crm_ng = 150L, crm_thr = 0.35,
    tie_high = TRUE,
    phi_elim = 0.25, elim_a0 = 1.0, elim_b0 = 1.0,
    gb_c1 = log(1.075), gb_c2 = log(1.075) / 3, gb_eps = 0.5, gb_N0 = 6L,
    loss = make_loss(k_over = 1.0)
  )
}

#' Parameters reproducing an earlier version of this work
#'
#' Restores the lower-dose tie convention and the equivalence-bound elimination
#' rule under the Jeffreys prior. Provided so that the corrections described in
#' the accompanying manuscript can be verified directly; it should not be used
#' for new work.
#'
#' @return A parameter list; see \code{\link{adapt_params}}.
#' @export
adapt_params_original <- function() {
  p <- adapt_params()
  p$tie_high <- FALSE
  p$phi_elim <- p$phi2
  p$elim_a0  <- p$a0
  p$elim_b0  <- p$b0
  p
}
