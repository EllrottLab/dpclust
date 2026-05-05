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
```

which compares cluster positions before and after the winner's curse-aware fit by cluster rank.

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
