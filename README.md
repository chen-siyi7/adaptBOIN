# adaptBOIN

Loss-calibrated escalation boundaries and dose-elimination diagnostics for
Phase I dose-finding trials.

This package accompanies the manuscript *A Practical Evaluation of Finite-Sample
BOIN Calibration and Bayesian Smoothing in Phase I Dose Finding* and reproduces
every table and figure in it.

## What it does

The escalation boundaries are derived by minimising posterior risk under an
explicit three-state loss over the actions escalate, stay and de-escalate. The
resulting rule always has a count-threshold form, so it deploys as a pre-computed
lookup table of the same shape as BOIN's. An overdose-aversion multiplier
`k_over` tightens the de-escalation boundary, making the tradeoff between
overdosing and underdosing an explicit design input rather than an implicit
consequence of fixed cutoffs.

The package also provides the comparator designs used in the paper (BOIN, gBOINS,
mTPI-2 and the continual reassessment method), a monotone Bernstein polynomial
end-of-trial estimator, and tools for examining how the dose-elimination rule
interacts with the closest-dose definition of the MTD.

## Installation

```r
# install.packages("remotes")
remotes::install_github("chen-siyi7/adaptBOIN")
```

A C++ compiler is required, since the simulation engine is written in C++ and
compiled on installation. On macOS run `xcode-select --install` first. On Windows
install Rtools.

## Quick start

```r
library(adaptBOIN)

# the decision table at the symmetric default
adapt_table()[["30"]]

# effective rate boundaries against BOIN's fixed cutoffs
boundary_rates()

# how overdose aversion tightens the upper boundary
boundaries_by_loss(k_overs = c(1, 1.5, 2))

# one scenario, one design
run_scenario(6, "adaptive_iso", n_sim = 500L)$pcs

# a grid, summarised
res <- run_all(n_sim = 500L, designs = c("boin", "adaptive_iso", "gboins"))
summarise_results(res)
```

## Reproducing the manuscript

```r
reproduce_paper(n_sim = 200L)                  # fast check, a few minutes
reproduce_paper(n_sim = 2000L, n_cores = 6L)   # full run
```

or from a shell,

```bash
Rscript -e 'adaptBOIN::reproduce_paper(n_sim = 2000L, n_cores = 6L)'
```

The bundled script offers the same thing with flags:

```bash
Rscript $(Rscript -e 'cat(system.file("scripts/reproduce_paper.R", package="adaptBOIN"))') --quick
```

Output lands in `adaptboin_output/` as CSV tables and, if **ggplot2** is
installed, PDF figures. Start with a small `n_sim`: the two Bernstein
configurations draw 2000 importance samples per trial and dominate the runtime.

## Two conventions that change results

Both are explicit parameters rather than hard-coded choices, because both
materially affect operating characteristics.

**Tie-breaking.** `tie_high` decides whether exact ties in
`|pi_j - phi_tgt|` are broken toward the higher dose, and is applied consistently
to the definition of the true MTD and to every end-of-trial estimator. This
matters more than it may appear. Because end-of-trial isotonic estimates are
discrete and are pooled by the pool-adjacent-violators algorithm, exact ties
arise frequently at `N = 30` even in scenarios whose true toxicity probabilities
contain no tie.

**Dose elimination.** `phi_elim`, `elim_a0` and `elim_b0` set the elimination
rule. The default is the standard BOIN rule, which tests whether toxicity exceeds
the target rate under a Beta(1, 1) prior. Because the MTD is defined as the dose
*closest* to the target rather than a dose *at* the target, the true MTD often has
toxicity above the target on a discrete grid, and a target-based rule then
eliminates it with probability that grows with sample size. Setting
`phi_elim = 0.35` tests the upper equivalence bound instead and does not have this
behaviour. `elim_rule_study()` quantifies the difference.

An earlier version of this work used the lower-dose tie convention and the
equivalence-bound rule. `adapt_params_original()` restores both, so that the
corrections reported in the manuscript can be verified directly:

```r
run_all(params = adapt_params_original(), n_sim = 2000L)
```

## Designs

| Name | Escalation | End-of-trial selection |
|---|---|---|
| `adaptive_iso` | posterior-risk boundaries | isotonic regression |
| `adaptive_bern` | posterior-risk boundaries | Bernstein posterior |
| `boin_bern` | standard BOIN | Bernstein posterior |
| `boin` | standard BOIN | isotonic regression |
| `crm` | continual reassessment method | posterior mean |
| `mtpi2` | mTPI-2 unit probability mass | isotonic regression |
| `gboins` | gBOINS shrinkage boundaries | isotonic regression |

The gBOINS implementation was validated against the boundary table published by
Mu, Hu, Xu and Pan (2021) at target rates of 0.20 and 0.30 before use.

## Targeted studies

```r
steep_curve_study(n_sim = 2000L)   # MTD at target, gap above it increasing
straddle_N_study(n_sim = 2000L)    # no dose at the target, swept over N
elim_rule_study(n_sim = 2000L)     # target-based against equivalence-bound
diagnose_offset_anomaly(n_sim = 2000L)  # per-dose allocation detail
```

`diagnose_alloc()` is the one to reach for when selection looks wrong and the
cause is unclear: it separates escalation drifting off the correct dose from
correct allocation followed by a failure at the end-of-trial step.

## Reproducibility notes

Seeds are deterministic. `run_scenario()` derives its seed from the scenario and
design, so results are stable across sessions. The targeted studies deliberately
share one base seed across designs within a scenario, which gives a partially
paired comparison: trials start from a common random stream and diverge only once
decisions differ. This reduces the variance of design differences appreciably. It
is not full common random numbers, which would require pre-generating patient
outcomes, and remains worth doing.

Note that `std::binomial_distribution` is not bit-portable across C++ standard
library versions, so absolute values may shift very slightly on a different
platform. Comparisons between designs within a platform are unaffected.
