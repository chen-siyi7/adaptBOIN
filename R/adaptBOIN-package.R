#' adaptBOIN: loss-calibrated escalation boundaries for Phase I dose finding
#'
#' Escalation boundaries derived by minimising posterior risk under an explicit
#' three-state loss, with an overdose-aversion multiplier that tightens the
#' de-escalation boundary; comparator designs (BOIN, gBOINS, mTPI-2, CRM); a
#' monotone Bernstein polynomial end-of-trial estimator; and diagnostics for the
#' interaction between the dose-elimination rule and the closest-dose definition
#' of the MTD.
#'
#' Start with \code{\link{adapt_params}} and \code{\link{run_all}}, or call
#' \code{\link{reproduce_paper}} to regenerate every table and figure.
#'
#' @section Conventions that matter:
#' Two settings materially change results and are therefore explicit parameters
#' rather than hard-coded choices. \code{tie_high} determines whether exact ties
#' in \eqn{|\pi_j - \phi_{tgt}|} are broken toward the higher dose, and is
#' applied consistently to the definition of the true MTD and to every
#' end-of-trial estimator. \code{phi_elim} together with \code{elim_a0} and
#' \code{elim_b0} determines the dose-elimination rule; the default is the
#' standard BOIN rule. See \code{\link{adapt_params_original}} for the settings
#' used in an earlier version of this work.
#'
#' @keywords internal
"_PACKAGE"

## quiet R CMD check for the internal data objects
utils::globalVariables(c("scenarios", "design_codes"))
