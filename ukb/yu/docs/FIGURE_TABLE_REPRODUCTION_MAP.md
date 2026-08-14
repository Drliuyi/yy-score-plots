# Figure and table reproduction map

| Published item | Published input | Local stage | Local output |
|---|---|---|---|
| Table S1 | Baseline cohort | `cohort` | `04_cohort/table_s1_characteristics.csv` |
| Table S2 | Cox incident associations | `cox` | `05_cox/table_s2_incident_associations.csv.gz` |
| Table S3 | Outcome-specific mortality | external/local extension | `06_sensitivity/` |
| Table S4 | Fine-Gray sensitivity | external/local extension | `06_sensitivity/` |
| Tables S5-S7 | Age, sex, extended adjustment | external/local extension | `06_sensitivity/` |
| Tables S8-S9 | 19 cardiac MRI phenotypes and protein linear models | `cmr_parallel` | `07_cmr/cmr_associations.csv.gz`, `cmr_association_summary.csv` |
| Tables S10-S12 | Prediction metrics, NRI/IDI, importance | `evaluate` | `10_evaluation/` |
| Tables S13-S14 | Forward/reverse MR | external summary statistics required | `11_mr/` |
| Table S15 | Drug-target enrichment | external database export required | `14_enrichment/` |
| Table S16 | Mediation | participant risk-factor inputs required | `12_mediation/` |
| Table S17 | PRS-protein associations | genotype/PRS input required | `13_prs/` |
| Tables S18-S22 | Enrichment, TF, PPI modules | external database export required | `14_enrichment/` |
| Tables S23-S26 | Protein, endpoint, risk factor, GWAS dictionaries | `sources` | `03_source_tables/` |
| Figure 1 | Study-flow summary | `figures` | `15_figures/figure1_workflow.*` |
| Figure 2 | Significant counts and 14 volcano panels | `figures` | PDF/TIFF/PNG plus `figure2_source_data_*.csv*` |
| Figure 3 | CMR Manhattan and top-15 heat map | `cmr_parallel` then `figures` | `15_figures/figure3_local_cmr.*`; official S9 is fallback only |
| Figure 4A | Published-model hold-out AUC bars | `figures` | `figure4a_prediction_auc.*` and source data |
| Figure 4B | Recurrent top-15 protein gain stacked by outcome | `figures` | `figure4b_recurrent_importance.*` and source data |
| Figure 5 | MR Manhattan and mediation panels | gated/reference | `15_figures/reference_figure5*` until local MR exists |
| Figure 6A | PRS-protein heat map | local `prs` stage | `13_prs/prs_protein_associations.csv.gz` and local Figure 6A |
| Figure 6B-D | Enrichment, PPI and TF panels | planned local downstream stage | use locally significant PRS proteins; official S18 is fallback only |
| Figures S1-S5 | Event pins, Cox breadth, sensitivity, PPI | `figures`/gated | `16_supplementary_figures/` |

External-data stages (`mr`, `prs`, `enrichment`) require explicit source files.
They fail closed when inputs are absent and never substitute old FairK results. For
local PRS downstream analyses, the foreground is frozen as proteins significant in
at least one local outcome and the background is all proteins tested in the local
PRS association matrix.
