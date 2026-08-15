# YY score methods

## Main CAD reproduction

Huang-lab reproduction is exposed as two commands:

```bash
yy score --main
yy plot --main
```

The first command runs the four locked CAD methods (`pradeep-strict`,
`yu-strict`, `pradeep-fair`, and `yu-fair`) and then creates the two strict
common Yin/Yang projections. The fair computations already write plot-ready
individual scores. The second command draws the locked baseline-centred raw
four-score trajectory and ROC figure. Interrupted work can be resumed with the
same first command.

Each formal method computation holds an output lock until its final child
process exits. A second invocation of the same method/disease fails early,
preventing duplicate folds from writing to one project.

The strict input resolver searches the Huang `data.BIG`, phenotype and `ppp`
trees. Both strict methods reuse `prot.tab.gz`; Yu does not require a second
1.1-GB decompressed protein table. Its prediction-only preflight does not
require the article PDF/workbook. When the processing-date table is missing,
the command retrieves official UKB Resource 1019 and validates the locked
SHA-256 before use.

This project exposes four distinct protein-score methods without using a
`legacy`/`unified` version switch:

| Method ID | Design | Native target |
|---|---|---|
| `pradeep-strict` | Pradeep/Schuermans 1.5k-panel reproduction | selected incident cardiac outcome in the frozen 80/20 split |
| `yu-strict` | Yu/Chen source-locked reproduction | selected incident outcome in the frozen derivation/holdout split |
| `pradeep-fair` | all-2,910-protein LASSO-logistic | 5-year incident CAD, locked common five-fold cohort |
| `yu-fair` | all-2,910-protein LightGBM | 5-year incident CAD, locked common five-fold cohort |

The strict methods can either delegate to their complete GitHub workflows under
`<SCRIPT_ROOT>/ukb/pradeep` and `<SCRIPT_ROOT>/ukb/yu`, or—normally—connect the
already completed derived projects under
`<YY_OUTDIR>/score/source-projects` with `--project`. Projection reads the
frozen Pradeep coefficients or Yu LightGBM booster and does not refit either
model. The fair methods first
rebuild their common input from `<PHEDIR>/Rdata/all.rds`,
`<PHEDIR>/Rdata/prot.rds`. When participant-level fold manifests are absent
under `<YY_OUTDIR>/reference/cad_fivefold_v1`, the workflow deterministically
regenerates them from those RDS sources and verifies normalized content hashes;
the fold files are never stored in Git. Models fit only within the locked
outer-training folds. That preparation
also writes the participant tables, 2,910-feature
float matrices and locked Yin/Yang target RDS files under
`<YY_OUTDIR>/score/common-fair-inputs`; after its `COMPLETE` marker exists,
`yy plot` can use this public-path rebuild directly.  The two native strict
AUCs and the two fair five-year IPCW AUCs are
different estimands and must not be ranked as if they were one matched test.

## Huang-lab defaults

```text
DIR0=/mnt/d
PHEDIR=/mnt/d/data/ukb/phe
SCRIPT_ROOT=/mnt/d/scripts
ANALYSIS_ROOT=/mnt/d/analysis
YY_OUTDIR=/mnt/d/analysis/yy
YY_SCORE_FOLD_ROOT=/mnt/d/analysis/yy/reference/cad_fivefold_v1
```

All roots accept environment overrides.  No entry point writes into input or
shared-helper directories.

## Public commands

```bash
yy score --h
yy score --status
yy score pradeep-strict --preflight
yy score yu-strict --preflight
yy score pradeep-strict --disease afib --preflight
yy score yu-strict --disease heart_failure --preflight
yy score pradeep-strict --project
yy score yu-strict --project
yy score pradeep-fair --project
yy score yu-fair --project
yy score pradeep-fair --preflight
yy score yu-fair --preflight
```

Formal computation is explicit:

```bash
yy score pradeep-strict --disease afib --compute --confirm-compute --workers=10 --resume
yy score yu-strict --disease heart_failure --compute --confirm-compute --workers=10 --resume
yy score pradeep-fair --compute --confirm-compute --workers=10 --resume
yy score yu-fair --compute --confirm-compute --workers=10 --resume
```

After a successful fair computation, the individual Yin OOF and Yang projected
scores are already plot-ready under `<YY_OUTDIR>/score/pradeep-fair` and
`<YY_OUTDIR>/score/yu-fair`; no additional `--project` call is needed. Use
`--project` in place of `--compute` only to connect an already completed common
cohort source project.

The default disease is `cad`. Pradeep supports `cad`, `afib`, `hfail` and
`ao_sten`. Yu supports the 14 IDs listed by `yy score --h`. Common aliases are
translated to each upstream project's native ID. Native outputs are isolated
under `<ANALYSIS_ROOT>/pradeep/<disease>` and
`<ANALYSIS_ROOT>/yu/<disease>`. Pradeep always constructs the four paper
outcome controls needed for its cohort rules and retains the original
four-outcome Bonferroni family while fitting only the requested outcome.

The two fair methods and `--project` remain CAD-only. Their common five-year
cohort, Yin/Yang participants and endpoint contract are CAD-specific; a
non-CAD native model is therefore never silently attached to a CAD plot.

The default command never fits a model.  `--compute` without
`--confirm-compute` stops before model fitting.

After both native projections complete, figures can be forced to use them
instead of a historical fallback cache:

```bash
yy plot --baseline --yy --traj --roc --require-score \
  --proteins=GDF15,PCSK9,NTPROBNP \
  --score pradeep-strict yu-strict
```
