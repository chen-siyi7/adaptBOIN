#' Scenario definitions table
#'
#' @return A data frame with the true toxicity probabilities, the MTD and an
#'   exact-tie flag.
#' @export
table_scenarios <- function() {
  do.call(rbind, lapply(seq_along(scenarios), function(i) {
    s <- scenarios[[i]]
    data.frame(sc = i, description = s$name,
               t(stats::setNames(s$pi, paste0("pi", seq_along(s$pi)))),
               MTD = paste0("d", s$mtd),
               exact_tie = is_tie(s$pi), stringsAsFactors = FALSE)
  }))
}

#' Dose-selection distribution table
#'
#' @param results List from \code{\link{run_all}}.
#' @param sc_idx Scenario index.
#' @param designs Design names.
#' @return A data frame of selection percentages by dose.
#' @export
table_selection <- function(results, sc_idx = 1, designs = primary_designs) {
  rows <- lapply(designs, function(dn) {
    r <- Find(function(x) x$sc_idx == sc_idx && x$design == dn, results)
    if (is.null(r)) return(NULL)
    sd <- selection_dist(r)
    J <- length(sd$dose)
    stats::setNames(
      data.frame(design = dn, t(round(c(sd$dose, sd$no_rec), 1)),
                 pod = round(r$pod, 1)),
      c("design", paste0("d", 1:J), "no_rec", "POD"))
  })
  do.call(rbind, rows)
}

#' Strict and acceptable-selection accuracy for a tie scenario
#'
#' @param results List from \code{\link{run_all}}.
#' @param sc_idx Scenario index.
#' @param designs Design names.
#' @return A data frame with both criteria.
#' @export
table_tie_pcs <- function(results, sc_idx = 1, designs = primary_designs) {
  acc <- acceptable_set(scenarios[[sc_idx]]$pi)
  rows <- lapply(designs, function(dn) {
    r <- Find(function(x) x$sc_idx == sc_idx && x$design == dn, results)
    if (is.null(r)) return(NULL)
    sd <- selection_dist(r)
    data.frame(design = dn, PCS_strict = round(r$pcs, 1),
               PCS_acc = round(sum(sd$dose[acc]), 1), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Safety operating characteristics
#'
#' @param df Data frame from \code{\link{summarise_results}}.
#' @param designs Design names.
#' @return A list of three wide data frames: overdose selection, mean DLTs, and
#'   the percentage of patients treated above the MTD.
#' @export
table_safety <- function(df, designs = primary_designs) {
  sub <- df[df$design %in% designs, ]
  wide <- function(metric) {
    m <- stats::reshape(sub[, c("sc", "design", metric)],
                        idvar = "sc", timevar = "design", direction = "wide")
    names(m) <- sub("\\.", "_", names(m)); m
  }
  list(over = wide("pct_over"), dlts = wide("mean_dlts"), pod = wide("pod"))
}

#' CRM skeleton sensitivity table
#'
#' @param crm_res List from \code{\link{crm_sensitivity}}.
#' @return A data frame.
#' @export
table_crm <- function(crm_res) {
  do.call(rbind, lapply(crm_res, function(r)
    data.frame(sc = r$sc_idx, MTD = paste0("d", r$true_mtd),
               skeleton = r$skel, pcs = round(r$pcs, 1),
               stringsAsFactors = FALSE)))
}
