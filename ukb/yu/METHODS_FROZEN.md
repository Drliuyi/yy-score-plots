# Frozen methods

## Design

- Baseline UKB-PPP Olink Explore 3072 data.
- Participants with any of 14 baseline cardiovascular diseases are excluded.
- Proteins with missingness above 30% are excluded.
- Random two-thirds derivation and one-third hold-out split; local split seed
  `20260715` because the publication did not report split EIDs.
- Full-panel Cox models use the reported adjustment set for all 14 outcomes.
- Derivation-only Bonferroni-positive proteins enter endpoint-specific
  preliminary LightGBM models. Proteins accounting for the first 30% cumulative
  gain are unioned across outcomes.
- Final predictor sets are `SCORE2`, `Protein`, and `Protein+SCORE2`.
- LightGBM uses 500 estimators, depth 15, 10 leaves, row fraction 0.7, learning
  rate 0.01 and feature fraction 0.7, with ten-fold derivation diagnostics.
- SCORE2 is calculated for the UK low-risk region and recalibrated by isotonic
  regression using derivation data only.
- Hold-out evaluation uses 1,000 participant-level bootstrap samples.

## Local decisions

- Follow-up cutoff: `2023-09-30`.
- Protein missing values remain missing for LightGBM native handling.
- `subsample_freq=1` makes the reported row fraction active.
- Classification thresholds are selected from derivation out-of-fold
  predictions by the Youden criterion and then frozen for the hold-out set.
- Local counts and selected proteins are reported as observed and are not forced
  to match the publication.

## Boundary

The authors did not release scripts, split EIDs or the LightGBM tuning search
space. This project is a source-locked independent implementation of the
reported method, not an execution of original author code. Published,
inferred and local decisions are recorded in `f/config/method_provenance.csv`.
The accepted workflow contains no YYScore extension.
