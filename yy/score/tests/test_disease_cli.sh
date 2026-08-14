#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLI="${PROJECT_ROOT}/score.sh"

help_out="$(bash "$CLI" --h)"
grep -q 'CAD fair end-to-end' <<<"$help_out"
grep -q 'do not run --project afterwards' <<<"$help_out"
grep -q 'Main CAD reproduction in two commands' <<<"$help_out"

main_plan="$(YY_SCORE_MAIN_PLAN_ONLY=1 bash "$CLI" --compute --confirm-compute --workers=10 --resume)"
[[ "$(grep -c '^MAIN_CAD_STEP ' <<<"$main_plan")" -eq 6 ]]
grep -q 'pradeep-strict --compute' <<<"$main_plan"
grep -q 'yu-strict .*--project' <<<"$main_plan"
grep -q 'MAIN_CAD_SCORE_PLAN_COMPLETE' <<<"$main_plan"

YY_ROOT="$(cd "${PROJECT_ROOT}/.." && pwd)"
SCRIPT_ROOT_FOR_TEST="$(cd "${YY_ROOT}/.." && pwd)"
TRUE_BIN="/usr/bin/true"
plot_plan="$(RSCRIPT="$TRUE_BIN" SCRIPT_ROOT="$SCRIPT_ROOT_FOR_TEST" bash "${YY_ROOT}/cli.sh" plot --main)"
grep -q 'yy plot preset: main-cad-four-score' <<<"$plot_plan"
grep -q 'pradeep-strict yu-strict pradeep-fair yu-fair' <<<"$plot_plan"

pradeep_cad_out="$(bash "$CLI" pradeep-strict --disease cad --inputs)"
grep -q 'YY_PROJECTION_OUTPUT=/mnt/d/analysis/yy/score/pradeep-strict' <<<"$pradeep_cad_out"

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
