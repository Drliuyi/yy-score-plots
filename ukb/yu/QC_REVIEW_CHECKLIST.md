# QC review checklist

- [ ] Public supplement workbook and methods PDF match recorded SHA256 hashes.
- [ ] Raw panel has 2,923 assays before local QC and 2,920 after the frozen 30% rule, or every difference is documented.
- [ ] Cohort removes all 14 baseline CVD outcomes before incident follow-up.
- [ ] Derivation/test EIDs and ten-fold IDs are frozen and hashed.
- [ ] Protein selection uses derivation data only; hold-out outcomes are never read.
- [ ] The local 30% cumulative-gain union and its protein count are recorded.
- [ ] SCORE2 calibration is fitted on derivation data only.
- [ ] Final models are exactly SCORE2, Protein and Protein+SCORE2 with frozen LightGBM 3.3.2 parameters.
- [ ] Cox, CMR, MR, mediation, PRS and systems stages label local versus reference evidence.
- [ ] Figure 1-6 each have non-empty PDF, PNG and TIFF plus source data.
- [ ] No FairK, ProtWAS, Top-K, YYScore or unrelated model result enters this project.
- [ ] `RESULTS_AND_QC.md`, stage markers and dependency manifests are present.
