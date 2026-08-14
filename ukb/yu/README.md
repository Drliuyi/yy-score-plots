# Yu Protein Analysis

Reproducible implementation of the Yu/Chen 2025 UKB-PPP proteomic analysis.
The project exposes a single public entry point:

```bash
./yu.sh --help
```

## Quick start

```bash
./yu.sh setup
./yu.sh doctor
./yu.sh all --confirm-heavy --resume
```

Use the same entry point through WSL from WinPC PowerShell:

```powershell
wsl bash -lc "cd /mnt/d/scripts/ukb/yu && ./yu.sh --help"
wsl bash -lc "cd /mnt/d/scripts/ukb/yu && ./yu.sh all --confirm-heavy --resume"
```

`setup` validates and saves paths without starting an analysis. `doctor` checks
the data, software, and code environment. `all` runs the complete workflow.
Common utility commands are:

```bash
./yu.sh status
./yu.sh figures --resume
./yu.sh finalize --resume
./yu.sh package
```

Steps can be run individually, as a list, or as a range:

```bash
./yu.sh 3
./yu.sh 3,4
./yu.sh 1-4 --resume
```

| Step | Description |
|---:|---|
| 1 | Literature sources, input validation, and panel QC |
| 2 | Incident cohort without baseline CVD |
| 3 | Full-panel Cox analysis for 14 outcomes |
| 4 | Training-set protein selection, LightGBM, and held-out evaluation |
| 5 | Cardiac MRI associations |
| 6 | MR |
| 7 | Mediation candidate analysis |
| 8 | CMAverse mediation analysis |
| 9 | PRS-protein analysis for 13 outcomes |
| 10 | Enrichment, transcription-factor, and STRING-PPI analyses |
| 11 | Figures 1-6 and the results report |

Steps 8 and 9 are computationally intensive and require `--confirm-heavy`.
`--resume` reuses only stages and shards that pass integrity checks.

## Paths

The default interface follows the Huang-lab D-drive layout:
`D:/data/ukb/phe`, `D:/scripts/ukb/yu`, `D:/analysis`, and
`D:/data.BIG/gwas/ppp`. The legacy `D:/UKB_data` tree is used only when
explicitly supplied with `--dir0 D:/UKB_data` or saved through `setup`.
Genotypes remain a separate root and are read directly from
`Z:/projects/genotype_pc_nas/imputed_pgen_autosomes` by default.

If a required path does not exist, an interactive terminal prompts for the
correct location and saves the confirmed value in the current user's
configuration directory. Paths can also be overridden explicitly:

```bash
./yu.sh 1-4 \
  --analysis-project yu_proteomic_repo_v3 \
  --raw-protein-file D:/data/ukb/phe/raw/prot_full_unimputed.tsv \
  --phenotype-rds D:/data/ukb/phe/Rdata/all.rds \
  --resume
```

Legacy WinPC override, when intentionally required:

```bash
./yu.sh 1-4 --dir0 D:/UKB_data \
  --raw-protein-file D:/UKB_data/phe/raw/prot_full_unimputed.tsv \
  --phenotype-rds D:/UKB_data/phe/Rdata/all.rds --resume
```

Use a new analysis project when changing the disease or model proteins:

```bash
./yu.sh 2-4 \
  --analysis-project yu_avs_custom_v1 \
  --disease aortic_valve_stenosis \
  --protein-panel custom \
  --model-proteins GDF15,NPPB,ADM,CST3
```

## Repository structure

```text
yu.sh          Single public command
f/entry        Low-level R entry points
f/R            Analysis and figure code
f/python       Modeling and PRS code
f/tools        Internal orchestration and installation tools
f/config       Frozen parameters and mappings
f/tests        Automated checks
references     Public supplementary files and source manifests
```

The public workflow is frozen as the Yu/Chen reproduction. It does not run
YYScore or read FairK, ProtWAS, Top-K, or other legacy-project results. See
[METHODS_FROZEN.md](METHODS_FROZEN.md) for the statistical boundaries and
[QC_REVIEW_CHECKLIST.md](QC_REVIEW_CHECKLIST.md) for delivery checks.
