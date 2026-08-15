# Input contract and publication boundary

## Publication rule

GitHub contains source code and small non-participant configuration files only.
Every file listed below is an **external runtime input** and must remain outside
the repository. The workflow prints resolved paths and fails early when a
required file is unavailable.

## Core Huang-lab inputs

| Role | Default WSL path | Used by | Repository policy |
|---|---|---|---|
| Baseline phenotype RDS | `<PHEDIR>/Rdata/all.rds` | all methods and plots | restricted; never upload |
| Standardized protein RDS | `<PHEDIR>/Rdata/prot.rds` | fair methods and protein curves | restricted; never upload |
| Pradeep raw 1.5k protein table | `<DIR0>/data.BIG/gwas/ppp/prot.tab.gz` | `pradeep-strict` | restricted; never upload |
| Pradeep 1.5k assay map | `<DIR0>/data.BIG/gwas/ppp/map.raw/olink_protein_map_1.5k_v1.tsv` | `pradeep-strict` | external reference; do not upload |
| PPP protein BED | `<DIR0>/data.BIG/gwas/ppp/ppp_3k_b38.bed` | Pradeep preflight/figures | external reference; do not upload |
| Yu unimputed protein table | `<PHEDIR>/raw/prot_full_unimputed.tsv` | `yu-strict` | restricted; never upload |
| Yu raw phenotype table | `<PHEDIR>/pheno.tsv.gz` | `yu-strict` | restricted; never upload |
| Yu 3k assay map | `<DIR0>/data.BIG/gwas/ppp/olink_protein_map_3k_v1.tsv` | `yu-strict` | external reference; do not upload |
| Yu supplementary workbook | `<DIR0>/files/yu-protein-analysis/references/raw/pwaf072_supplementary_table_1.xlsx` | `yu-strict` source lock | external paper file; do not upload |
| Yu supplementary methods | `<DIR0>/files/yu-protein-analysis/references/raw/pwaf072_supplementary_figure_1.pdf` | `yu-strict` source lock | external paper file; do not upload |
| Olink processing dates | `<DIR0>/files/yu-protein-analysis/references/raw/olink_processing_start_date.dat` | `yu-strict` | controlled UKB resource; never upload |

The 1.5k/3k maps and PPP BED are small, but they are deliberately treated as
external references rather than bundled data. If Huang's machine lacks them,
place reviewed copies at the paths above or set the documented environment
overrides; do not add them to Git.

## Deterministically generated fold outputs

The common fair comparison uses participant-level fold manifests generated
from `<PHEDIR>/Rdata/all.rds` and `<PHEDIR>/Rdata/prot.rds`. Because they
contain EIDs, they are written under the analysis tree rather than the source
tree:

```text
<YY_OUTDIR>/reference/cad_fivefold_v1/fold_assignment_yin.csv
<YY_OUTDIR>/reference/cad_fivefold_v1/fold_assignment_yang.csv
```

Override the directory with `YY_SCORE_FOLD_ROOT`. When both files are absent,
fair preflight reconstructs them in memory and fair compute writes them using
the frozen event-stratified rules (`seed=2026` for Yin; duration-stratified
`seed=2027` for Yang). Existing files are never overwritten. A partial pair is
an error. The code validates counts, event labels, fold structure, participant
separation, and frozen canonical content identity. Canonical hashes ignore
Windows/WSL line-ending differences while remaining sensitive to every EID,
event and fold assignment. Only hashes are public:

```text
Yin SHA-256:  3c325decc66fb6ae9a3a24190ff55b3481bf23e38add6ef9913e589f436f380d
Yang SHA-256: bf65dc6d34043c67471dab29d8a866fb6be6ec2a94d5b8b5627f5be0bcb3fdbc
```

## Runtime dependencies

- R is resolved from `RSCRIPT` (default `/opt/R/4.3.2/bin/Rscript`).
- Formal Yu computation requires Python 3.9 with the packages pinned in
  `ukb/yu/f/config/requirements-py39.txt`.
- Set `YU_PYTHON` explicitly, put `python3.9` on `PATH`, or install the Yu
  environment under `<DIR0>/software`.
- PLINK2 and external LD/genotype inputs are required only for optional
  Pradeep MR/colocalization steps, not for the four-score CAD figure. Keep
  genotype data on its reviewed external/Zspace root.

## Derived outputs

All fitted models, individual scores, matrices, plots, logs, and manifests are
written below `<ANALYSIS_ROOT>` or `<YY_OUTDIR>`. They are runtime products and
are not committed to GitHub.
