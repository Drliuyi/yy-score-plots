# Inputs and end-to-end boundaries

All logical defaults follow the Huang-lab D-drive contract:

```text
DIR0=/mnt/d
PHEDIR=/mnt/d/data/ukb/phe
SCRIPT_ROOT=/mnt/d/scripts
ANALYSIS_ROOT=/mnt/d/analysis
YY_OUTDIR=/mnt/d/analysis/yy
YY_SCORE_FOLD_ROOT=/mnt/d/analysis/yy/reference/cad_fivefold_v1
```

Inspect the effective paths without fitting anything:

```bash
yy score --inputs
yy score pradeep-strict --inputs
yy score yu-strict --inputs
yy score pradeep-fair --inputs
yy score yu-fair --inputs
```

## Pradeep strict

The dispatcher delegates to `<SCRIPT_ROOT>/ukb/pradeep/pradeep.sh` with the
selected `--outcome` (`cad`, `afib`, `hfail`, or `ao_sten`) and runs the 1.5k-panel source workflow. Its primary raw inputs are
`<PHEDIR>/Rdata/all.rds`, the unimputed `prot.tab.gz`, the 1.5k Olink map and
the PPP protein BED.  It rebuilds the cohort and panel, runs the association
steps, and fits the native frozen 80/20 LASSO-logistic prediction analysis.
Outputs are isolated as a disease child project such as
`<ANALYSIS_ROOT>/pradeep/cad` or `<ANALYSIS_ROOT>/pradeep/afib`.
The other cardiac endpoints remain paper-defined cohort controls and
time-dependent covariates. The original four-outcome multiplicity family is
retained even when only one outcome is fitted.

The normal no-refit route is `yy score pradeep-strict --project`. It resolves
an already completed derived model project under
`<YY_OUTDIR>/score/source-projects/pradeep-strict-cad` (or an explicit
`--source-root`) and reads:

```text
outputs/ukbppp_cardiac_analysis_base.rds
outputs/lasso/cad_coefficients.csv
outputs/lasso/cad_predictions.csv
```

It bridges the common 2,910-protein matrix to the frozen Pradeep preprocessing
scale using outcome-blind overlapping measurements, verifies the replay
against stored held-out scores, and writes individual common Yin/Yang scores
under `<YY_OUTDIR>/score/pradeep-strict`.

The 1.5k raw table and assay map are needed only for a new strict refit. They
must be exposed through reviewed public data paths or explicit environment
overrides. They are not needed when projecting an already completed derived
model project.

## Yu strict

The dispatcher delegates to `<SCRIPT_ROOT>/ukb/yu/yu.sh` with the selected
native `--disease` ID. Its primary raw
inputs are `<PHEDIR>/raw/prot_full_unimputed.tsv`,
`<PHEDIR>/Rdata/all.rds`, `<PHEDIR>/pheno.tsv.gz` and the reviewed 3k assay
map. Olink processing dates are passed explicitly from
`<DIR0>/files/yu-protein-analysis/references/raw` when rebuilding the native
cohort. The public supplementary workbook and methods PDF are consulted only
by the full source-audit stage; they are not numerical inputs to a frozen CAD
score or to the fair comparison. The dispatcher never accepts a saved legacy
path configuration silently.
Steps 1-4 rebuild the input QC, selected incident-disease cohort, full-panel
Cox screen, local protein selection and native LightGBM held-out evaluation.
Outputs are isolated under a disease child such as
`<ANALYSIS_ROOT>/yu/cad` or `<ANALYSIS_ROOT>/yu/heart_failure`.

The no-refit common Yin/Yang projection described below is currently available
only for CAD, because its participant and endpoint contract is CAD-specific.

The normal no-refit route is `yy score yu-strict --project`. It reads the
existing derived project under `<YY_OUTDIR>/score/source-projects/yu-strict-cad`,
including its `09_models/cad__Protein.txt` booster and
`09_models/test_predictions.csv.gz`, applies the booster to the original
unimputed protein table (preserving LightGBM native missing-value handling),
verifies replay of stored held-out predictions, and writes individual common
Yin/Yang scores under `<YY_OUTDIR>/score/yu-strict`.

## Common fair inputs

Both fair methods independently rebuild a locked comparison input from
`<PHEDIR>/Rdata/all.rds` and `<PHEDIR>/Rdata/prot.rds`. If Yin and Yang fold
manifests are absent under `<YY_SCORE_FOLD_ROOT>`, the workflow regenerates
them from those two source files with the frozen event/duration-stratified
rules and seeds, then verifies their counts and published hashes before
installation. These participant-level outputs contain EIDs and are
intentionally excluded from GitHub. The preparation recomputes the CAD endpoint, checks
every EID and event against the frozen contract, then writes participant
tables, ordered 2,910-protein float matrices and target RDS files under
`<YY_OUTDIR>/score/common-fair-inputs`.

Pradeep fair and Yu fair then use the same 37,127-person Yin cohort, five
outer folds, five-year target and 2,910-feature order.  Yang participants are
never used in fitting; they are projected only after each outer model has
been trained.  Outputs are isolated under
`<YY_OUTDIR>/score/pradeep-fair` and `<YY_OUTDIR>/score/yu-fair`.

An existing locked fair-comparison project can likewise be connected with
`yy score pradeep-fair --project` and `yy score yu-fair --project`. That route
copies only the method-specific individual score rows and metrics into the
public score contract; it does not fit or tune a model.

No formal fit starts without both `--compute` and `--confirm-compute`.
`--project` is separate and never fits or tunes a model. If common inputs do
not yet exist, a strict projection may build them once from `all.rds` and
`prot.rds`; subsequent projections and plots reuse the `COMPLETE` marker.
