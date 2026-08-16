#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLI="${PROJECT_ROOT}/score.sh"

help_out="$(bash "$CLI" --h)"
grep -q 'CAD fair end-to-end' <<<"$help_out"
grep -q 'do not run --project afterwards' <<<"$help_out"
grep -q 'Main CAD reproduction in two commands' <<<"$help_out"
grep -q 'yy score --main' <<<"$help_out"
grep -q 'yy score --aa' <<<"$help_out"

main_plan="$(YY_SCORE_MAIN_PLAN_ONLY=1 bash "$CLI" --compute --confirm-compute --workers=10 --resume)"
[[ "$(grep -c '^MAIN_CAD_STEP ' <<<"$main_plan")" -eq 6 ]]
grep -q 'pradeep-strict --compute' <<<"$main_plan"
grep -q 'yu-strict .*--project' <<<"$main_plan"
grep -q 'MAIN_CAD_SCORE_PLAN_COMPLETE' <<<"$main_plan"

main_alias_plan="$(YY_SCORE_MAIN_PLAN_ONLY=1 bash "$CLI" --main)"
[[ "$(grep -c '^MAIN_CAD_STEP ' <<<"$main_alias_plan")" -eq 6 ]]
grep -q 'yy score action: compute' <<<"$main_alias_plan"
grep -q 'pradeep-strict --compute --confirm-compute --workers 10 --disease cad --resume' <<<"$main_alias_plan"
grep -q 'MAIN_CAD_SCORE_PLAN_COMPLETE' <<<"$main_alias_plan"

area_auc_plan="$(YY_SCORE_AA_PLAN_ONLY=1 bash "$CLI" --aa --workers=7)"
grep -q 'yy score preset: area-auc' <<<"$area_auc_plan"
grep -q 'yy score action: compute' <<<"$area_auc_plan"
grep -q -- '--stage=compute' <<<"$area_auc_plan"
grep -q -- '--workers=7' <<<"$area_auc_plan"
grep -q -- '--resume' <<<"$area_auc_plan"
grep -q -- '--output-root=/mnt/d/analysis/yy/area-auc' <<<"$area_auc_plan"

area_auc_preflight_plan="$(YY_SCORE_AA_PLAN_ONLY=1 bash "$CLI" --aa --preflight)"
grep -q 'yy score action: preflight' <<<"$area_auc_preflight_plan"
grep -q -- '--stage=preflight' <<<"$area_auc_preflight_plan"
if grep -q -- '--resume' <<<"$area_auc_preflight_plan"; then
  echo 'Area-AUC preflight unexpectedly enabled resume.' >&2
  exit 1
fi

area_auc_status_plan="$(YY_SCORE_AA_PLAN_ONLY=1 bash "$CLI" --aa --status)"
grep -q 'yy score action: status' <<<"$area_auc_status_plan"
grep -q -- '--stage=status' <<<"$area_auc_status_plan"

if YY_SCORE_AA_PLAN_ONLY=1 bash "$CLI" --main --aa >/dev/null 2>&1; then
  echo 'Mutually exclusive --main and --aa were accepted.' >&2
  exit 1
fi

YY_ROOT="$(cd "${PROJECT_ROOT}/.." && pwd)"
SCRIPT_ROOT_FOR_TEST="$(cd "${YY_ROOT}/.." && pwd)"
TRUE_BIN="/usr/bin/true"
plot_plan="$(RSCRIPT="$TRUE_BIN" SCRIPT_ROOT="$SCRIPT_ROOT_FOR_TEST" bash "${YY_ROOT}/cli.sh" plot --main)"
grep -q 'yy plot preset: main-cad-four-score' <<<"$plot_plan"
grep -q 'pradeep-strict yu-strict pradeep-fair yu-fair' <<<"$plot_plan"

pradeep_cad_out="$(bash "$CLI" pradeep-strict --disease cad --inputs)"
grep -q 'YY_PROJECTION_OUTPUT=/mnt/d/analysis/yy/score/pradeep-strict' <<<"$pradeep_cad_out"
grep -q 'DERIVED_COMMON_INPUT=/mnt/d/analysis/yy/score/common-fair-inputs' <<<"$pradeep_cad_out"

yu_cad_out="$(bash "$CLI" yu-strict --disease cad --inputs)"
grep -q 'YY_PROJECTION_OUTPUT=/mnt/d/analysis/yy/score/yu-strict' <<<"$yu_cad_out"

pradeep_out="$(bash "$CLI" pradeep-strict --disease atrial_fibrillation --inputs)"
grep -q 'yy score disease: afib' <<<"$pradeep_out"
grep -q 'NEW_COMPUTE_OUTPUT=/mnt/d/analysis/pradeep/afib' <<<"$pradeep_out"
grep -q 'YY_PROJECTION_OUTPUT=NOT_APPLICABLE_CAD_ONLY' <<<"$pradeep_out"

yu_out="$(bash "$CLI" yu-strict --disease hfail --inputs)"
grep -q 'yy score disease: heart_failure' <<<"$yu_out"
grep -q 'NEW_COMPUTE_OUTPUT=/mnt/d/analysis/yu/heart_failure' <<<"$yu_out"
grep -q 'YY_PROJECTION_OUTPUT=NOT_APPLICABLE_CAD_ONLY' <<<"$yu_out"

if bash "$CLI" pradeep-strict --disease ischemic_stroke --inputs >/dev/null 2>&1; then
  echo 'Unsupported Pradeep disease was accepted.' >&2
  exit 1
fi

if bash "$CLI" yu-fair --disease heart_failure --inputs >/dev/null 2>&1; then
  echo 'Non-CAD fair comparison was accepted.' >&2
  exit 1
fi

if bash "$CLI" yu-strict --disease heart_failure --project >/dev/null 2>&1; then
  echo 'Non-CAD common Yin/Yang projection was accepted.' >&2
  exit 1
fi

echo 'test_disease_cli.sh: PASS'
