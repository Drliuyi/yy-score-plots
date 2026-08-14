# Frozen method definitions

## Strict reproductions

`pradeep-strict` and `yu-strict` are disease-selectable article-method reproductions.
Their cohort definitions, protein panels, screening, train/test split,
hyperparameters and native binary AUC are owned by the corresponding upstream
projects. The Pradeep run retains the original four-outcome cohort controls
and multiplicity family while fitting the requested one of `cad`, `afib`,
`hfail`, or `ao_sten`; the Yu run uses one of its 14 upstream disease IDs.
This dispatcher does not substitute another endpoint.

`--project` is a frozen-parameter connector, not a fifth model. For Pradeep it
reads the saved protein-model coefficients and restores their native scaled
input space; for Yu it reads the saved booster and retains native missing-value
handling from the unimputed protein table. Both connectors must reproduce the
stored native held-out scores before writing common Yin/Yang individual scores.
The connector is currently CAD-only. Non-CAD models remain valid native
reproductions but require a separately defined disease-specific Yin/Yang
participant and endpoint contract before projection.

## Fair comparison

The fair methods share the following fixed elements:

- 37,127 baseline-CAD-free participants and 3,442 incident CAD events;
- 1,766 baseline prevalent-CAD participants for score projection only;
- the same 2,910 proteins in the same order from `prot.rds`;
- the same five outer folds and the same five-year endpoint;
- training-fold-only imputation/preprocessing and no use of Yang labels when
  fitting either incident prediction model.

`pradeep-fair` uses a protein-only LASSO-logistic model.  Within each outer
training fold, participants with known five-year status are used for training,
proteins are median-imputed and standardized from that training subset, and
10-fold stratified inner CV selects `lambda.1se` using AUC.

`yu-fair` uses all 2,910 proteins in a LightGBM binary classifier with the
publication-aligned fixed parameters in `config/fair.json`.  It uses the same
known-five-year training subset and outer test fold as `pradeep-fair`.

Yang participants never enter fitting.  Each fitted outer-fold model is only
projected onto Yang.  The five Yang projections are averaged for trajectory
display.
