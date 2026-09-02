#' Per-dose allocation diagnostic
#'
#' Loops the internal single-trial simulator to recover where patients actually
#' went and how selection related to allocation, which the aggregate simulator does not
#' expose. This separates two explanations for poor selection: escalation
#' drifting off the correct dose, versus correct allocation followed by a failure
#' at the end-of-trial step.
#'
#' @param pi True DLT probabilities.
#' @param design_name Design name.
#' @param N Maximum sample size.
#' @param n_sim Number of simulated trials.
#' @param params Parameter list.
#' @param base_seed Base seed.
#' @return A data frame with one row per dose, including mean allocation,
#'   observed rate, selection percentage, and accuracy conditional on the correct
#'   dose having received the most patients.
#' @export
diagnose_alloc <- function(pi, design_name, N, n_sim = 2000L,
                           params = adapt_params(), base_seed = 900000L) {
  pp <- params; pp$N <- as.integer(N)
  atbl <- adapt_table(pp)
  J <- pp$J
  mtd <- mtd_of(pi)
  nmat <- matrix(0L, n_sim, J); smat <- matrix(0L, n_sim, J)
  sel  <- integer(n_sim)
  for (i in seq_len(n_sim)) {
    tr <- one_trial(pi, design_codes[[design_name]], atbl, pp,
                    as.integer(base_seed + i))
    nmat[i, ] <- as.integer(tr$n_vec)
    smat[i, ] <- as.integer(tr$s_vec)
    sel[i]    <- as.integer(tr$mtd_sel)   # 0 = no recommendation
  }
  argmax_n <- apply(nmat, 1, which.max)
  correct_alloc <- argmax_n == mtd
  data.frame(
    design = design_name, N = N, mtd = mtd, dose = 1:J,
    mean_n = round(colMeans(nmat), 2),
    mean_rate = round(colSums(smat) / pmax(colSums(nmat), 1), 3),
    pct_visited = round(100 * colMeans(nmat > 0), 1),
    pct_selected = round(100 * vapply(1:J, function(j) mean(sel == j), 0), 1),
    pct_norec = round(100 * mean(sel == 0), 1),
    pct_alloc_max = round(100 * vapply(1:J, function(j) mean(argmax_n == j), 0), 1),
    pcs_given_correct_alloc = round(100 * mean(sel[correct_alloc] == mtd), 1),
    pct_correct_alloc = round(100 * mean(correct_alloc), 1),
    stringsAsFactors = FALSE)
}

#' Allocation diagnostic for the offset-MTD scenario
#'
#' Runs \code{\link{diagnose_alloc}} across sample sizes and designs for the
#' scenario whose MTD lies above the target, where accuracy is non-monotone in
#' sample size under the standard elimination rule.
#'
#' @param n_sim Number of simulated trials.
#' @param params Parameter list.
#' @param Nvals Maximum sample sizes.
#' @param designs Design names.
#' @return A data frame.
#' @export
diagnose_offset_anomaly <- function(n_sim = 2000L, params = adapt_params(),
                                    Nvals = c(30L, 60L, 90L),
                                    designs = c("adaptive_iso", "boin",
                                                "gboins")) {
  pi <- offset_scenarios[[2]]$pi
  out <- list()
  for (N in Nvals) for (dn in designs) {
    message(sprintf("  diagnosing %s at N=%d", dn, N))
    # Match the appended offset-scenario seed in elim_rule_study() so that the
    # allocation diagnostic and aggregate PCS table describe the same trials.
    elim_index <- length(affected_scenarios()) + 1L
    out[[length(out) + 1]] <- diagnose_alloc(pi, dn, N, n_sim, params,
      base_seed = 950000L + 1000L * elim_index + N)
  }
  do.call(rbind, out)
}

# Bernstein basis evaluated at the equally spaced standardized dose labels.
.bern_basis_matrix <- function(K, J) {
  d <- (seq_len(J) - 1) / (J - 1)
  t(vapply(d, function(x) stats::dbinom(0:K, K, x), numeric(K + 1L)))
}

.bern_curves <- function(v, B) {
  w <- t(apply(v[, -ncol(v), drop = FALSE], 1, cumsum))
  p <- w %*% t(B)
  p[p < 1e-12] <- 1e-12
  p[p > 1 - 1e-12] <- 1 - 1e-12
  p
}

.bern_is <- function(n, s, params, M = params$M_IS, seed = 41001L) {
  set.seed(seed)
  q <- params$K + 2L
  v <- matrix(stats::rgamma(M * q, params$alpha_d), M, q)
  v <- v / rowSums(v)
  p <- .bern_curves(v, .bern_basis_matrix(params$K, length(n)))
  ll <- rowSums(sweep(log(p), 2, s, `*`) +
                sweep(log1p(-p), 2, n - s, `*`))
  ww <- exp(ll - max(ll)); ww <- ww / sum(ww)
  list(mean = colSums(p * ww), ess = 1 / sum(ww^2))
}

.bern_logpost_z <- function(z, n, s, alpha, B) {
  ez <- exp(c(z, 0) - max(c(z, 0)))
  v <- ez / sum(ez)
  w <- cumsum(v)[-length(v)]
  p <- pmin(1 - 1e-12, pmax(1e-12, as.vector(B %*% w)))
  # Under the additive-log-ratio parameterisation the Jacobian contributes
  # prod(v), so Dirichlet(alpha) gives alpha * sum(log(v)).
  sum(s * log(p) + (n - s) * log1p(-p)) + alpha * sum(log(v))
}

.bern_mcmc_chain <- function(n, s, params, n_iter, burn, thin, step, seed) {
  set.seed(seed)
  q <- params$K + 2L; B <- .bern_basis_matrix(params$K, length(n))
  z <- stats::rnorm(q - 1L, 0, 0.5)
  lp <- .bern_logpost_z(z, n, s, params$alpha_d, B)
  keep <- matrix(NA_real_, (n_iter - burn) %/% thin, length(n))
  out_i <- 0L; accepted <- 0L; window <- 0L
  for (i in seq_len(n_iter)) {
    zp <- z + stats::rnorm(q - 1L, 0, step)
    lpp <- .bern_logpost_z(zp, n, s, params$alpha_d, B)
    if (log(stats::runif(1)) < lpp - lp) {
      z <- zp; lp <- lpp; accepted <- accepted + 1L; window <- window + 1L
    }
    if (i <= burn && i %% 100L == 0L) {
      ar <- window / 100
      if (ar < 0.20) step <- step * 0.8
      if (ar > 0.45) step <- step * 1.2
      window <- 0L
    }
    if (i > burn && (i - burn) %% thin == 0L) {
      out_i <- out_i + 1L
      ez <- exp(c(z, 0) - max(c(z, 0))); v <- ez / sum(ez)
      keep[out_i, ] <- as.vector(B %*% cumsum(v)[-q])
    }
  }
  list(draws = keep, accept = accepted / n_iter, step = step)
}

#' Compare Bernstein importance sampling with an independent MCMC posterior
#'
#' The reference sampler uses four independently initialized random-walk
#' Metropolis chains on additive-log-ratio coordinates for the Dirichlet
#' increments. It targets the same posterior as the operational importance
#' sampler but shares neither its prior draws nor its weighting calculation.
#'
#' @param n,s Per-dose patient and DLT counts.
#' @param params Parameter list.
#' @param M Number of operational importance draws.
#' @param n_iter,burn,thin MCMC controls.
#' @param chains Number of independent chains.
#' @param seed Reproducibility seed.
#' @return A list with one-row summary and per-dose posterior means.
#' @export
bern_reference_check <- function(n, s, params = adapt_params(),
                                 M = params$M_IS, n_iter = 50000L,
                                 burn = 15000L, thin = 10L, chains = 4L,
                                 seed = 41001L) {
  stopifnot(length(n) == length(s), all(s <= n), chains >= 2L,
            n_iter > burn, (n_iter - burn) %% thin == 0L)
  if (!requireNamespace("coda", quietly = TRUE))
    stop("Package 'coda' is required for the reference-posterior diagnostic.")

  is <- .bern_is(n, s, params, M, seed)
  ch <- lapply(seq_len(chains), function(k)
    .bern_mcmc_chain(n, s, params, n_iter, burn, thin, 0.45,
                     seed + 1000L * k))
  ml <- coda::mcmc.list(lapply(ch, function(x) coda::mcmc(x$draws)))
  mcmc_mean <- colMeans(do.call(rbind, lapply(ch, `[[`, "draws")))
  rhat <- coda::gelman.diag(ml, autoburnin = FALSE,
                            multivariate = FALSE)$psrf[, 1]
  mc_ess <- coda::effectiveSize(ml)
  is_sel <- mtd_of(is$mean, params$phi_tgt, params$tie_high)
  ref_sel <- mtd_of(mcmc_mean, params$phi_tgt, params$tie_high)
  dose <- data.frame(dose = seq_along(n), n = n, s = s,
                     is_mean = is$mean, mcmc_mean = mcmc_mean,
                     abs_diff = abs(is$mean - mcmc_mean),
                     rhat = rhat, mcmc_ess = mc_ess)
  summary <- data.frame(
    is_ess = is$ess, max_abs_diff = max(dose$abs_diff),
    is_selection = is_sel, mcmc_selection = ref_sel,
    selection_agrees = is_sel == ref_sel,
    max_rhat = max(rhat), min_mcmc_ess = min(mc_ess),
    mean_accept = mean(vapply(ch, `[[`, numeric(1), "accept")))
  list(summary = summary, dose = dose)
}

#' Reference-posterior study for representative simulated trials
#'
#' Selects trials near the 10th and 50th percentiles of the operational
#' importance-sampling ESS in four scenarios, then compares each posterior with
#' independent MCMC using \code{\link{bern_reference_check}}.
#'
#' @param scenarios_idx Scenario indices to diagnose.
#' @param candidate_trials Number of simulated trials used to select each pair.
#' @inheritParams bern_reference_check
#' @return A list with summary and dose-level data frames.
#' @export
bern_reference_study <- function(scenarios_idx = c(1L, 2L, 5L, 10L),
                                 candidate_trials = 100L,
                                 params = adapt_params(), M = params$M_IS,
                                 n_iter = 50000L, burn = 15000L, thin = 10L,
                                 chains = 4L, seed = 51001L) {
  atbl <- adapt_table(params); ss <- list(); dd <- list()
  for (sc in scenarios_idx) {
    tr <- lapply(seq_len(candidate_trials), function(i)
      one_trial(scenarios[[sc]]$pi, design_codes[["adaptive_bern"]], atbl,
                params, as.integer(seed + 10000L * sc + i)))
    ev <- vapply(tr, `[[`, numeric(1), "ess")
    for (prob in c(0.10, 0.50)) {
      target <- as.numeric(stats::quantile(ev, prob, names = FALSE))
      ii <- which.min(abs(ev - target)); x <- tr[[ii]]
      ck <- bern_reference_check(x$n_vec, x$s_vec, params, M, n_iter,
                                 burn, thin, chains,
                                 seed + 100000L * sc + round(100 * prob))
      tag <- sprintf("Sc%d_q%02d", sc, round(100 * prob))
      ss[[tag]] <- cbind(dataset = tag, scenario = sc, ess_quantile = prob,
                         candidate_ess = ev[ii], ck$summary)
      dd[[tag]] <- cbind(dataset = tag, scenario = sc, ess_quantile = prob,
                         ck$dose)
    }
  }
  list(summary = do.call(rbind, ss), dose = do.call(rbind, dd))
}
