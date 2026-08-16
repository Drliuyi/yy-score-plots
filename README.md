# Yin–Yang Score Plots

Code-only workflow for drawing CAD Yin–Yang protein trajectories together with
Pradeep-style and Yu-style participant prediction scores and ROC curves.

The repository contains source code, method configuration, tests, and
documentation only. It does **not** contain UK Biobank participant data, EIDs,
phenotypes, protein matrices, fold assignments, fitted models, predictions, or
analysis results.

## Install

```bash
git clone https://github.com/Drliuyi/yy-score-plots.git
cd yy-score-plots
bash yy/install.sh
yy --h
```

The entrypoint resolves the repository checkout as `SCRIPT_ROOT`. Huang-lab
data and output defaults remain `/mnt/d/data/ukb/phe`, `/mnt/d/analysis`, and
`/mnt/d/analysis/yy`; a checkout located at `/mnt/d/scripts` keeps the original
Huang script-root layout. Every root can be overridden by an environment
variable.

## Main CAD workflow

```bash
yy score --inputs
yy score --main
yy plot --main
```

`yy score --main` is the locked, resumable CAD preset: it runs the four score
methods and produces their plot-ready individual scores. Use
`yy score --preflight --workers=10` first when only checking a new machine.
Per-method process locks reject a second concurrent computation instead of
silently launching duplicate folds into the same output directory.

The score command keeps four methods separate:

- `pradeep-strict`: original 1.5k-panel LASSO-logistic reproduction;
- `yu-strict`: original Yu LightGBM reproduction;
- `pradeep-fair`: all-2,910-protein LASSO-logistic on the common cohort;
- `yu-fair`: all-2,910-protein LightGBM on the same common cohort.

Strict and fair ROC curves have different estimands and are not a paired AUC
comparison. See [yy/README.md](yy/README.md) for plotting examples and
[yy/INPUTS.md](yy/INPUTS.md) for the full external-input contract.

The same locked `yy score` common-cohort products can also be passed to
`yy/R/area_auc_from_score.R` for cross-fitted protein-level Yin/Yang 0–5-year
trajectory-area versus five-year IPCW-AUC analysis; see
[yy/README.md](yy/README.md#cross-fitted-trajectory-area-and-auc5).

## External references not included in Git

The following small, non-participant metadata files are deliberately not
bundled in this code-only repository. The command searches the reviewed
`data.BIG`, phenotype and `ppp` locations before reporting a missing path:

```text
/mnt/d/data.BIG/gwas/ppp/map.raw/olink_protein_map_1.5k_v1.tsv
/mnt/d/data.BIG/gwas/ppp/ppp_3k_b38.bed
/mnt/d/data.BIG/gwas/ppp/olink_protein_map_3k_v1.tsv
```

The underlying UK Biobank assay catalogue is public Resource 1013:
https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=1013. The two local mapping
tables and BED contain assay/gene metadata rather than participant data.

Yu preferentially reuses the existing `prot.tab.gz`; no decompressed
`prot_full_unimputed.tsv` duplicate is required. If the Olink processing-date
table is absent, Yu preflight/compute downloads official UKB Resource 1019 and
accepts it only after its locked SHA-256 checksum passes:
https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=1019.

Participant-level five-fold manifests are generated under
`/mnt/d/analysis/yy/reference/cad_fivefold_v1` from Huang's `all.rds/prot.rds`
when absent. They are never uploaded; only their deterministic rules and
validation hashes are stored in the repository.
