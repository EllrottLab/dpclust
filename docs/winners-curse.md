# Winner's Curse Correction

DPClust can automatically correct for winner's curse bias during clustering.

Winner's curse bias occurs when low-purity or low-CCF mutations are only detected if random sampling gives them unusually high mutant read support. In those cases, the observed CCF distribution is biased upward, and fitting the ordinary binomial read-count likelihood can place cluster peaks too high.

## What DPClust Does

By default, DPClust uses:

```text
winner_curse_correction = auto
winner_curse_threshold = 3
winner_curse_compare_plots = true
```

In `auto` mode, DPClust first evaluates each loaded sample for low expected mutant-read support. If the sample appears vulnerable to winner's curse bias, DPClust enables a corrected likelihood in the Gibbs sampler.

This is a model-level correction. It is not a post-processing shift of the final cluster table.

When enabled, the mutation read-count likelihood is conditioned on mutation detectability:

```text
P(read counts | CCF, mutation was detectable)
```

For the default threshold, detectability means:

```text
mutant reads >= 3
```

Cluster locations are then sampled under the corrected likelihood.

## Implementation Checklist

This section maps the implementation back to the original standalone winner's curse code.

The original code had this slightly confusing threshold convention:

```r
wcc <- function(..., T, ...) {
    T = T - 1
    ...
    for (j in 0:T) {
        sp = sp + dbinom(j, rd, tbp)
    }
    sf = 1 / (1 - sp)
}
```

With `T = 3`, that code subtracts one and sums the missed mutations with `0`, `1`, or `2` mutant reads. The correction is therefore conditioning on mutations having:

```text
mutant reads >= 3
```

DPClust implements the same threshold meaning in C++:

| Original WCC behavior | DPClust implementation |
|---|---|
| `T = T - 1` means default `T = 3` corrects for missing `0`, `1`, and `2` mutant-read observations | `src/dpclust_rcpp.cpp:37-68` computes missed detection as `P(X < threshold)`; `src/dpclust_rcpp.cpp:53-55` is the optimized `threshold == 3` branch that sums `0`, `1`, and `2` mutant reads |
| `sf = 1 / (1 - sp)` divides by detection probability | `src/dpclust_rcpp.cpp:75-79` subtracts `log_detection_prob_from_logq(...)`, which is the log-likelihood equivalent of dividing by detection probability |
| The correction applies to observed mutations that pass the detection threshold | `src/dpclust_rcpp.cpp:225` sets `row.observed_at_detection_threshold = row.mut >= winner_curse_threshold`; `src/dpclust_rcpp.cpp:83-85` only applies the conditional likelihood for rows at or above that threshold |
| The default detection threshold is `3` mutant reads | `R/DirichletProcessClustering.R:28-36` and `R/pipeline_cli.R:49-50` default `winner_curse_threshold = 3L` |

This addresses the important `T - 1` concern from the original script: default `T = 3` still means missed observations are `0`, `1`, and `2` mutant reads, and detected observations are `>= 3` mutant reads.

This is not a post-hoc shift of the final CCF table. The standalone script uses a separate one-dimensional Brent optimization after clusters have already been identified. DPClust instead samples cluster locations inside the Gibbs model under the conditional read-count likelihood. Therefore the corrected DPClust result should follow the same detection-threshold logic, but it is not expected to match the standalone script bit-for-bit or guarantee the exact same `1.0` cluster position in every case.

## Why This Is Implemented Inside DPClust

A post-hoc winner's curse script is useful for proving the effect, but it is a poor long-term place to apply the correction.

The bias affects the read-count evidence that DPClust uses to decide where clusters are. If the correction is applied only after fitting, the clustering step has already seen biased evidence. That means the final table can be adjusted, but the cluster discovery, mutation assignment probabilities, density plots, and cluster labels may still reflect the uncorrected likelihood.

Putting the correction inside the DPClust likelihood keeps the model internally consistent:

| Pipeline step | Why model-level correction is better |
|---|---|
| Cluster location fitting | Cluster positions are sampled using the corrected evidence, rather than moved after the fact |
| Mutation assignment | Assignment probabilities are calculated against the corrected cluster positions |
| Plots | The plotted posterior density reflects the same model that produced the cluster table |
| Outputs | Corrected cluster info, assignments, and QC plots come from one coherent fit |
| Reproducibility | The correction is controlled by normal DPClust parameters and logged in the run, not by a separate script with separate assumptions |

The tradeoff is that the corrected model is not a literal reproduction of the standalone script. The standalone script asks, "given an already chosen peak, where would this peak move after conditioning on detection?" DPClust asks, "where should the clusters be if detection bias is part of the observation model from the start?" The second question is the better software and statistical integration for routine pipeline use, even though the numbers are not guaranteed to be bit-for-bit identical to the standalone check.

The before/after cluster outputs are written here:

| Requested output | DPClust file |
|---|---|
| Uncorrected cluster positions | `*_winnerCurse_uncorrected_bestClusterInfo.txt` |
| Corrected cluster positions | `*_winnerCurse_corrected_bestClusterInfo.txt` |
| Side-by-side cluster position comparison | `*_winnerCurse_clusterComparison.txt` |
| 1D before/after density and cluster-position plot | `SAMPLE_winnerCurse_before_after_1D.png` |

The code that creates those files is in `R/DirichletProcessClustering.R`: `R/DirichletProcessClustering.R:677-711` runs the uncorrected comparator fit with `winner_curse_correction = FALSE`; `R/DirichletProcessClustering.R:717-747` runs the corrected fit with the effective winner's curse decision; `R/DirichletProcessClustering.R:332-346` defines `.write_winner_curse_cluster_info()`; `R/DirichletProcessClustering.R:748-762` writes the uncorrected and corrected `bestClusterInfo` tables; `R/DirichletProcessClustering.R:297-329` writes the side-by-side comparison table; and `R/DirichletProcessClustering.R:769-777` writes the 1D before/after plot.

## Auto-Detection Criteria

DPClust computes expected mutant reads per mutation using the loaded read counts, purity, copy-number scaling, and multiplicity adjustment.

The diagnostic high-risk window is:

```text
expected_alt_reads <= 2 * winner_curse_threshold
```

With the default threshold of `3`, this means:

```text
expected_alt_reads <= 6
```

In `auto` mode, correction is enabled if any sample satisfies one of these criteria:

| Criterion | Meaning |
|---|---|
| At least 10% of mutations have expected alt reads in the high-risk window | Many mutations are close to the detection boundary |
| Median expected alt reads is no more than `3 * threshold` and at least 10% of observed mutant counts are low | The sample has broadly weak mutant read support |
| Purity is at most 0.30 and at least 5% of mutations are in the high-risk window | Low purity makes even moderate low-support fractions risky |

These thresholds are intended to avoid applying the heavier corrected model to well-powered samples.

## Logging

Each run logs the winner's curse decision before fitting:

```text
Winner's curse model assumptions:
  policy                     : auto
  effective enabled           : TRUE
  decision reason             : auto: enabled because expected mutant-read support is near the detection boundary
  detection threshold         : >= 3 mutant reads
  likelihood adjustment       : condition SNV/indel read-count likelihood on mutation being detectable
  observation model           : binomial mutant reads with per-locus purity/copy-number/multiplicity scaling
  diagnostic high-risk window : expected_alt_reads <= 2 * threshold
  refit behavior              : cluster locations are sampled under the corrected likelihood, not post-hoc shifted
```

DPClust also logs per-sample diagnostics, for example:

```text
Winner's curse diagnostics (SAMPLE) sample 1: purity=0.126, expected_alt_reads q05=4.00 q50=6.00 q95=11.00, 57.2% expected <= 6 reads, 65.6% observed <= 6 reads
```

The key field is `effective enabled`. This is the actual boolean passed to the C++ Gibbs sampler.

## User Controls

Normal runs should use the default:

```text
--winner_curse_correction auto
```

For validation or comparison, the behavior can be forced:

```text
--winner_curse_correction true
--winner_curse_correction false
```

For transparent QC, DPClust runs an uncorrected comparator fit and writes before/after outputs by default when correction is effectively enabled:

```text
--winner_curse_compare_plots true
```

This does not replace the primary output. The primary output still follows `--winner_curse_correction` (`auto`, `true`, or `false`). When comparison plots are enabled and correction is effectively enabled, DPClust first runs an uncorrected comparator, preserves those plots with a `winnerCurse_uncorrected` filename tag, then runs the corrected fit normally.

For a uniform correction/reporting pass across all samples, use:

```text
--winner_curse_correction true
--winner_curse_compare_plots true
```

To suppress the extra comparator run:

```text
--winner_curse_compare_plots false
```

For one-dimensional runs, this also writes:

```text
SAMPLE_winnerCurse_before_after_1D.png
```

which overlays the uncorrected and corrected density curves and cluster locations. All runs with comparison enabled write:

```text
*_winnerCurse_clusterComparison.txt
*_winnerCurse_uncorrected_bestClusterInfo.txt
*_winnerCurse_corrected_bestClusterInfo.txt
```

These files record the before/after cluster positions directly. The comparison table aligns cluster positions before and after the winner's curse-aware fit by cluster rank.

The detection threshold can also be changed:

```text
--winner_curse_threshold 4
```

Additional sampler controls are available:

```text
--winner_curse_mh_sd 0.12
--winner_curse_mh_steps 8
```

These control the Metropolis updates used for cluster locations when the corrected likelihood is active. Most users should leave them at their defaults.

## Runtime Notes

The corrected model is slower than the original DPClust likelihood. It evaluates a truncated binomial observation model and uses non-conjugate cluster-location updates.

For exploratory runs on low-purity samples, prefer fewer iterations first:

```text
--iterations 5000 --burnin 1000
```

Then compare against a longer run if the cluster locations or plots are unstable.

## Interpretation

When correction is enabled, the primary cluster outputs and plots already reflect the winner's curse-aware fit. There is no separate post-hoc corrected cluster table unless `--winner_curse_compare_plots true` is enabled, in which case the extra files are QC comparison artifacts from an additional uncorrected fit.

If `effective enabled` is `FALSE`, DPClust uses the ordinary likelihood because the sample did not appear to need the correction under the auto-detection criteria, or because correction was forced off.
