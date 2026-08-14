# Yin–Yang Score Plots

This project draws baseline- or diagnosis-centred Yin–Yang trajectories for
proteins and participant-level prediction scores. It also provides reproducible
CAD score workflows for the Pradeep/Schuermans and Yu/Chen methods.

The public repository contains **code, small method configuration files, tests,
and documentation only**. It does not contain UK Biobank rows, EIDs, phenotype
tables, protein matrices, fold assignments, fitted models, predictions, or
analysis results.

## Install the `yy` command

From any WSL checkout of this standalone repository:

```bash
bash yy/install.sh
yy --h
```

The installer creates `~/.local/bin/yy` as a symbolic link to
the current checkout's `yy/cli.sh`. The command resolves `SCRIPT_ROOT` from the
checkout automatically; setting `SCRIPT_ROOT` remains available as an override.
Ensure `~/.local/bin` is on `PATH`.

## Main CAD reproduction

The complete score-to-figure workflow is intentionally two commands:

```bash
yy score --compute --confirm-compute --workers=10 --resume
yy plot --main
```

The first command computes four named methods and creates plot-ready scores:

- `pradeep-strict`: source-method 1.5k-panel LASSO-logistic reproduction;
- `yu-strict`: source-method Yu LightGBM reproduction;
- `pradeep-fair`: all-2,910-protein LASSO-logistic on the locked common
  five-year/five-fold CAD cohort;
- `yu-fair`: all-2,910-protein LightGBM on the same common cohort.

The second command draws one baseline-centred combined trajectory/ROC figure
with participant-count bars. Interrupted computations may be resumed with the
same score command.

Strict and fair AUCs are different estimands. Strict curves retain the native
paper-style holdout definition; fair curves use common five-year out-of-fold
predictions. Their juxtaposition is descriptive and is not a paired comparison.

## Inspect before computing

```bash
yy score --inputs
yy score --preflight --workers=10
yy score --status
yy plot --status
```

No formal model fit starts unless both `--compute` and `--confirm-compute` are
present.

## Custom figures

```bash
yy plot --baseline --yy --traj --roc --bar \
  --proteins=GDF15,PCSK9,NTPROBNP \
  --score pradeep-strict yu-strict pradeep-fair yu-fair

yy plot --diagnosis --yy --traj --roc \
  --adj=age+sex+le8+medications \
  --proteins=GDF15,PCSK9 \
  --score pradeep-fair yu-fair
```

Run `yy plot --h` for side, anchor, adjustment, mean/SD, bar, protein, score,
and ROC options.

## Huang-lab path contract

```text
DIR0=/mnt/d
PHEDIR=/mnt/d/data/ukb/phe
SCRIPT_ROOT=/mnt/d/scripts
ANALYSIS_ROOT=/mnt/d/analysis
YY_OUTDIR=/mnt/d/analysis/yy
```

All roots accept environment overrides. Required inputs and their publication
policy are listed in [INPUTS.md](INPUTS.md). Frozen method details are in
[score/METHODS.md](score/METHODS.md).
