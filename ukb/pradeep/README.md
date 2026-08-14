# Pradeep/Schuermans UKB-PPP reproduction

Command-driven reproduction of the Pradeep/Schuermans UKB-PPP cardiac
proteomics workflow, including paper-style figures and optional local
cis-pQTL MR/coloc.

## Quick start

```bash
bash pradeep.sh --h
```

From WinPC PowerShell, run the project in WSL:

```powershell
wsl bash /mnt/d/scripts/ukb/pradeep/pradeep.sh --h
```

From an existing WSL terminal, use Bash continuation (`\\`), not the
PowerShell backtick:

```bash
cd /mnt/d/scripts/ukb/pradeep
./pradeep.sh \
  --step core \
  --profile winpc \
  --panel 1.5k \
  --allow-kinship-fallback \
  --workers 12
```

`--workers 12` remains the overall worker request. The memory-heavy Cox and
sex-interaction stages automatically cap themselves at four workers and reuse
completed chunk files, preventing WSL out-of-memory failures without changing
the tested proteins or statistical models.

## Protein panels

`--panel 1.5k` is the Pradeep paper-panel mode. It reads the un-imputed raw
protein table, resolves the 1,459 paper assays, and creates dedicated
`all.rds` and `prot.rds` files under `<analysis>/inputs/pradeep_15k/` before
Step 1. The generated `all.rds` is derived from the supplied phenotype
`all.rds` after the paper-required race and PC availability gate; it is not
the existing 3k analysis input. The builder blocks pre-imputed raw protein
input and incomplete assay mapping.

`--panel 3k` uses the existing `all.rds` and `prot.rds` supplied by the local
phenotype project and does not apply the 1.5k assay map. It is a full-panel
extension, not a strict paper replication.

The default result directory is always `pradeep`. Select either `1.5k` or `3k`
for a project run; they are alternative inputs, not two analyses that must both
be performed. The selected panel is locked in the project audit manifest. A
later command requesting the other panel fails before writing outputs. To keep
both analyses intentionally, provide different `--analysis-project` names.

```bash
# Pradeep 1.5k panel
./pradeep.sh --step core --panel 1.5k --workers 12

# Existing full 3k panel
./pradeep.sh --step core --panel 3k --workers 12
```

All analysis and figure code is under `f/`. Run `--h` for steps, inputs,
outputs, path overrides, and examples.

## Workflow

1. Build the cohort, outcomes, protein QC and analysis base.
2. Run protein Cox associations for CAD, AF, HF and aortic stenosis.
3. Run sex-stratified and interaction analyses.
4. Fit clinical, protein and combined logistic LASSO models.
5. Generate the complete paper-style figure set.
6. Run GO enrichment and generate enrichment figures.
7. Standardize FinnGen GWAS files.
8. Run local cis-pQTL MR/coloc and generate the combined figures.

## Validated environment

| Program | Version |
|---|---|
| WSL Ubuntu | 22.04.5 LTS |
| Bash | 5.1.16 |
| R | 4.3.2 |
| PLINK2 | 2.0.0-a.6.9LM (29 Jan 2025; Step 8 only) |

| R package | Version | R package | Version |
|---|---:|---|---:|
| data.table | 1.18.2.1 | dplyr | 1.2.0 |
| survival | 3.8-6 | impute | 1.76.0 |
| glmnet | 4.1-10 | pROC | 1.19.0.1 |
| ggplot2 | 4.0.2 | ggrepel | 0.9.8 |
| patchwork | 1.3.2 | clusterProfiler | 4.10.1 |
| org.Hs.eg.db | 3.22.0 | AnnotationDbi | 1.72.0 |
| ggtext | 0.1.2 | ggpubr | 0.6.3 |
| gtools | 3.9.5 | openxlsx | 4.2.8.1 |
| tidyr | 1.3.2 | TwoSampleMR | 0.7.4 |
| MendelianRandomization | 0.10.0 | coloc | 5.2.3 |

Install missing packages in WSL R:

```r
install.packages(c(
  "data.table", "dplyr", "survival", "glmnet", "pROC", "ggplot2",
  "ggrepel", "patchwork", "ggtext", "ggpubr", "gtools", "openxlsx",
  "tidyr", "MendelianRandomization", "coloc", "BiocManager"
), repos = "https://cloud.r-project.org")

BiocManager::install(c(
  "impute", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi"
), ask = FALSE, update = FALSE)

install.packages(
  "TwoSampleMR",
  repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org")
)
```

Preflight checks package availability for the selected steps. The versions
above document the tested environment; the workflow never silently installs
or upgrades packages.

Default Huang-lab paths follow `D:/data/ukb/phe`, `D:/scripts` and
`D:/analysis/ukb`. Use `--profile winpc` for the legacy `D:/UKB_data` tree.
Every path can be overridden; missing inputs stop before computation.

Key outputs are saved under the selected analysis directory:

- `audit/command_manifest.csv`: resolved inputs and parameters.
- `audit/run_status.csv`: per-step status.
- `audit/output_manifest.csv`: output inventory and checksums.
- `audit/figure_completeness.csv`: required-figure acceptance gate.
- `outputs/figures/`: complete figure set.

The repository contains code only. UK Biobank phenotype, proteomic, genotype,
pQTL, GWAS, and LD-reference data are not included.
