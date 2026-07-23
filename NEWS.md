# adaptBOIN 0.3.0

## Corrections

Two defects in earlier versions changed simulation results, and both are fixed.
Users who ran an earlier version should re-run.

* **Tie-breaking was inconsistent with the documentation.** The end-of-trial
  estimators and the internal definition of the true MTD resolved exact ties in
  `|pi_j - phi_tgt|` with a strict inequality, breaking ties toward the *lower*
  dose, whereas the documentation specified the higher dose. In scenarios
  containing an exact tie, accuracy was therefore scored against a different dose
  from the one designated as the MTD. Scenarios 1 and 4 are affected. The
  convention is now the parameter `tie_high`, applied consistently to the true
  MTD and to every estimator.

  This reverses the sign of two results reported previously. In Scenario 1 the
  boundary contribution changes from +3.0 to -3.5 percentage points and the
  Bernstein contribution from +25.8 to -30.8.

* **The elimination rule was not the published BOIN rule.** It tested
  `Pr(pi_j > phi_2)` under the Jeffreys prior rather than
  `Pr(pi_j > phi_tgt)` under Beta(1, 1). The rule is now parameterised by
  `phi_elim`, `elim_a0` and `elim_b0`, and defaults to the standard BOIN rule so
  that the `boin` arm is the published design.

`adapt_params_original()` restores both previous settings, so the change can be
verified directly.

## New features

* `make_loss()` gains `k_over`, an overdose-aversion multiplier that tightens the
  de-escalation boundary. The admissible range is `k_over < 3`; at 3 the rule no
  longer escalates after a cohort with no events, violating the nondegeneracy
  condition required for the count-threshold representation.
* `gboins` added as a design, implementing the shrinkage boundaries of Mu, Hu, Xu
  and Pan (2021) for a binary endpoint. Validated against their published
  boundary table at target rates of 0.20 and 0.30.
* The engine reports `pct_mtd_elim`, the percentage of trials in which the true
  MTD was eliminated, and `one_trial()` returns `elim_vec`.
* `elim_rule_study()` compares target-based against equivalence-bound
  elimination for scenarios whose MTD lies above the target.
* `steep_curve_study()` and `straddle_N_study()` characterise the regimes where
  the fixed equivalence band outperforms shrinking boundaries.
* `diagnose_alloc()` and `diagnose_offset_anomaly()` recover per-dose allocation,
  separating escalation drifting off the correct dose from a failure at the
  end-of-trial step.
* `loss_sensitivity()` sweeps the overdose-aversion frontier.
* `reproduce_paper()` regenerates every table and figure in one call.

## Other changes

* Sensitivity analyses default to 2000 replicates, matching the primary analysis.
  Earlier 200-replicate results were too noisy to interpret.
* The targeted studies share one base seed across designs within a scenario,
  giving a partially paired comparison that reduces the variance of design
  differences.
* Documentation now states the Dirichlet concentration correctly: values below
  one concentrate mass at the simplex vertices, and values above one concentrate
  toward the centroid. Earlier text had this reversed.
