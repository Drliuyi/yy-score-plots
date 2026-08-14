# Frozen release contract

Project: **Yin–Yang Score Plots**

- Release: `v1.0.0`
- Locked date: `2026-08-15`
- Primary endpoint preset: CAD
- Public interface: `yy score` followed by `yy plot`
- Primary figure preset: `yy plot --main`

The release freezes method names, score estimands, common-cohort counts,
five-year horizon, outer-fold count, model hyperparameters, plot defaults, and
input-path contracts. Exact code identity is the Git commit recorded by the
release tag/branch.

The source tree contains no UKB participant data or EIDs. In particular, the
locked Yin/Yang fold manifests are protected runtime inputs validated by the
hashes recorded in `score/config/fair.json`; they are not Git files.

The two strict models and two fair models remain separate named estimands. A
future change to cohort membership, endpoint definition, source panel,
hyperparameters, or fold manifests requires a new release version.
