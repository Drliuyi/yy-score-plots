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
yy score --preflight --workers=10
yy score --compute --confirm-compute --workers=10 --resume
yy plot --main
```

The score command keeps four methods separate:

- `pradeep-strict`: original 1.5k-panel LASSO-logistic reproduction;
- `yu-strict`: original Yu LightGBM reproduction;
- `pradeep-fair`: all-2,910-protein LASSO-logistic on the common cohort;
- `yu-fair`: all-2,910-protein LightGBM on the same common cohort.

Strict and fair ROC curves have different estimands and are not a paired AUC
comparison. See [yy/README.md](yy/README.md) for plotting examples and
[yy/INPUTS.md](yy/INPUTS.md) for the full external-input contract.

## External references not included in Git

The following small, non-participant metadata files are deliberately not
bundled in this code-only repository. They are already present on the audited
WinPC Huang paths; a fresh installation may rebuild them from the UK Biobank
Olink assay resource and reviewed gene coordinates:

```text
/mnt/d/data.BIG/gwas/ppp/map.raw/olink_protein_map_1.5k_v1.tsv
/mnt/d/data.BIG/gwas/ppp/ppp_3k_b38.bed
/mnt/d/data.BIG/gwas/ppp/olink_protein_map_3k_v1.tsv
```

The underlying UK Biobank assay catalogue is public Resource 1013:
https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=1013. The two local mapping
tables and BED contain assay/gene metadata rather than participant data.

Participant-level five-fold manifests are generated under
`/mnt/d/analysis/yy/reference/cad_fivefold_v1` from Huang's `all.rds/prot.rds`
when absent. They are never uploaded; only their deterministic rules and
validation hashes are stored in the repository.
