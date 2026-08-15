#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIR0="${DIR0:-/mnt/d}"
PHEDIR="${PHEDIR:-${DIR0}/data/ukb/phe}"
SCRIPT_ROOT="${SCRIPT_ROOT:-${DIR0}/scripts}"
ANALYSIS_ROOT="${ANALYSIS_ROOT:-${DIR0}/analysis}"
YY_OUTDIR="${YY_OUTDIR:-${ANALYSIS_ROOT}/yy}"
COMMON_ROOT="${YY_SCORE_COMMON_ROOT:-${YY_OUTDIR}/score/common-fair-inputs}"
FOLD_ROOT="${YY_SCORE_FOLD_ROOT:-${YY_OUTDIR}/reference/cad_fivefold_v1}"
SOURCE_PROJECT_ROOT="${YY_SCORE_SOURCE_PROJECT_ROOT:-${YY_OUTDIR}/score/source-projects}"
RSCRIPT="${RSCRIPT:-/opt/R/4.3.2/bin/Rscript}"

first_existing_or_default() {
  local fallback="$1"
  shift
  local candidate
  for candidate in "$@"; do
    if [[ -f "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  printf '%s\n' "$fallback"
}

file_sha256() {
  local path="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$path" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$path" | awk '{print $1}'
  else
    echo "No SHA-256 utility found (sha256sum or shasum)." >&2
    return 2
  fi
}

PRADEEP_NATIVE_ROOT="${PRADEEP_NATIVE_ROOT:-${ANALYSIS_ROOT}/pradeep}"
PRADEEP_STRICT_PROJECT="cad"
if [[ -z "${PRADEEP_STRICT_RAW_PROTEIN:-}" ]]; then
  PRADEEP_STRICT_RAW_PROTEIN="$(first_existing_or_default \
    "${DIR0}/data.BIG/gwas/ppp/prot.tab.gz" \
    "${DIR0}/data.BIG/gwas/ppp/prot.tab.gz" \
    "${PHEDIR}/rap/raw/prot.tab.gz" \
    "${PHEDIR}/raw/prot.tab.gz" \
    "${DIR0}/ppp/prot.tab.gz")"
fi
if [[ -z "${PRADEEP_STRICT_PROTEIN_MAP:-}" ]]; then
  PRADEEP_STRICT_PROTEIN_MAP="$(first_existing_or_default \
    "${DIR0}/data.BIG/gwas/ppp/map.raw/olink_protein_map_1.5k_v1.tsv" \
    "${DIR0}/data.BIG/gwas/ppp/map.raw/olink_protein_map_1.5k_v1.tsv" \
    "${DIR0}/ppp/map.raw/olink_protein_map_1.5k_v1.tsv" \
    "${DIR0}/ppp/olink_protein_map_1.5k_v1.tsv")"
fi
if [[ -z "${PRADEEP_STRICT_BED:-}" ]]; then
  PRADEEP_STRICT_BED="$(first_existing_or_default \
    "${DIR0}/data.BIG/gwas/ppp/ppp_3k_b38.bed" \
    "${DIR0}/data.BIG/gwas/ppp/ppp_3k_b38.bed" \
    "${DIR0}/data.BIG/gwas/ppp/ppp.b38.bed" \
    "${DIR0}/ppp/ppp_3k_b38.bed" \
    "${DIR0}/ppp/ppp.b38.bed")"
fi
YU_NATIVE_ROOT="${YU_NATIVE_ROOT:-${ANALYSIS_ROOT}/yu}"
YU_STRICT_PROJECT="cad"
if [[ -z "${YU_STRICT_RAW_PROTEIN:-}" ]]; then
  YU_STRICT_RAW_PROTEIN="$(first_existing_or_default \
    "${PRADEEP_STRICT_RAW_PROTEIN}" \
    "${PHEDIR}/raw/prot_full_unimputed.tsv" \
    "${PHEDIR}/raw/prot_full_unimputed.tsv.gz" \
    "${PRADEEP_STRICT_RAW_PROTEIN}" \
    "${PHEDIR}/rap/raw/prot.tab.gz" \
    "${PHEDIR}/raw/prot.tab.gz" \
    "${DIR0}/ppp/prot.tab.gz")"
fi
if [[ -z "${YU_STRICT_RAW_PHENOTYPE:-}" ]]; then
  YU_STRICT_RAW_PHENOTYPE="$(first_existing_or_default \
    "${PHEDIR}/pheno.tsv.gz" \
    "${PHEDIR}/pheno.tsv.gz" \
    "${PHEDIR}/rap/pheno.tsv.gz" \
    "${PHEDIR}/rap/raw/pheno.tsv.gz" \
    "${PHEDIR}/common/pheno.tsv.gz")"
fi
if [[ -z "${YU_STRICT_PANEL_MAP:-}" ]]; then
  YU_STRICT_PANEL_MAP="$(first_existing_or_default \
    "${DIR0}/data.BIG/gwas/ppp/olink_protein_map_3k_v1.tsv" \
    "${DIR0}/data.BIG/gwas/ppp/olink_protein_map_3k_v1.tsv" \
    "${DIR0}/data.BIG/gwas/ppp/map.raw/olink_protein_map_3k_v1.tsv" \
    "${DIR0}/ppp/olink_protein_map_3k_v1.tsv" \
    "${DIR0}/ppp/map.raw/olink_protein_map_3k_v1.tsv")"
fi
YU_STRICT_REFERENCE_ROOT="${YU_STRICT_REFERENCE_ROOT:-${DIR0}/files/yu-protein-analysis/references/raw}"
YU_STRICT_SUPPLEMENT_WORKBOOK="${YU_STRICT_SUPPLEMENT_WORKBOOK:-${YU_STRICT_REFERENCE_ROOT}/pwaf072_supplementary_table_1.xlsx}"
YU_STRICT_SUPPLEMENT_METHODS="${YU_STRICT_SUPPLEMENT_METHODS:-${YU_STRICT_REFERENCE_ROOT}/pwaf072_supplementary_figure_1.pdf}"
if [[ -z "${YU_STRICT_OLINK_DATES:-}" ]]; then
  YU_STRICT_OLINK_DATES="$(first_existing_or_default \
    "${YU_STRICT_REFERENCE_ROOT}/olink_processing_start_date.dat" \
    "${YU_STRICT_REFERENCE_ROOT}/olink_processing_start_date.dat" \
    "${SCRIPT_ROOT}/ukb/yu/references/raw/olink_processing_start_date.dat" \
    "${PHEDIR}/rap/raw/olink_processing_start_date.dat" \
    "${DIR0}/ppp/olink_processing_start_date.dat")"
fi
YU_OLINK_DATES_URL="${YU_OLINK_DATES_URL:-https://biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/olink_processing_start_date.dat}"
YU_OLINK_DATES_SHA256="${YU_OLINK_DATES_SHA256:-249d5400603ea57c647ef812ff3ab6cb4ef4990bf0d4aeb39bb5693e573f4380}"
YU_PYTHON="${YU_PYTHON:-}"

ensure_yu_olink_dates() {
  local target="$YU_STRICT_OLINK_DATES"
  local observed=""
  if [[ -f "$target" ]]; then
    observed="$(file_sha256 "$target")" || return $?
    if [[ "$observed" != "$YU_OLINK_DATES_SHA256" ]]; then
      echo "UKB Resource 1019 checksum mismatch: ${target}" >&2
      echo "Expected: ${YU_OLINK_DATES_SHA256}" >&2
      echo "Observed: ${observed}" >&2
      return 2
    fi
    echo "UKB Resource 1019: VERIFIED ${target}"
    return 0
  fi

  mkdir -p "$(dirname "$target")"
  local temporary
  temporary="$(mktemp "${target}.download.XXXXXX")"
  if command -v curl >/dev/null 2>&1; then
    if ! curl --fail --location --silent --show-error "$YU_OLINK_DATES_URL" --output "$temporary"; then
      unlink "$temporary"
      echo "Failed to download UKB Resource 1019 from ${YU_OLINK_DATES_URL}" >&2
      return 2
    fi
  elif command -v wget >/dev/null 2>&1; then
    if ! wget -qO "$temporary" "$YU_OLINK_DATES_URL"; then
      unlink "$temporary"
      echo "Failed to download UKB Resource 1019 from ${YU_OLINK_DATES_URL}" >&2
      return 2
    fi
  else
    unlink "$temporary"
    echo "Missing downloader: install curl or wget, then rerun the same command." >&2
    return 2
  fi
  observed="$(file_sha256 "$temporary")" || { unlink "$temporary"; return 2; }
  if [[ "$observed" != "$YU_OLINK_DATES_SHA256" ]]; then
    unlink "$temporary"
    echo "Downloaded UKB Resource 1019 failed checksum validation." >&2
    echo "Expected: ${YU_OLINK_DATES_SHA256}" >&2
    echo "Observed: ${observed}" >&2
    return 2
  fi
  mv "$temporary" "$target"
  echo "UKB Resource 1019: DOWNLOADED_AND_VERIFIED ${target}"
  echo "Source: ${YU_OLINK_DATES_URL}"
}

acquire_compute_lock() {
  local method_name="$1"
  local disease_name="$2"
  command -v flock >/dev/null 2>&1 || {
    echo "Formal computation requires 'flock' to prevent duplicate model runs." >&2
    return 2
  }
  local lock_root="${YY_OUTDIR}/score/.locks"
  local lock_file="${lock_root}/${method_name}-${disease_name}.lock"
  local legacy_marker=""
  local legacy_output=""
  case "$method_name" in
    pradeep-fair)
      legacy_marker="02_pradeep_fair.R"
      legacy_output="--output-root=${YY_OUTDIR}/score/pradeep-fair"
      ;;
    yu-fair)
      legacy_marker="03_yu_fair.py"
      legacy_output="--output-root=${YY_OUTDIR}/score/yu-fair"
      ;;
  esac
  if [[ -n "$legacy_marker" ]]; then
    local process_path process_pid process_arg process_command
    for process_path in /proc/[0-9]*; do
      [[ -r "${process_path}/cmdline" ]] || continue
      process_pid="${process_path##*/}"
      [[ "$process_pid" == "$$" ]] && continue
      process_command=""
      while IFS= read -r -d '' process_arg; do
        process_command+="${process_arg} "
      done <"${process_path}/cmdline" || true
      if [[ "$process_command" == *"$legacy_marker"* && "$process_command" == *"$legacy_output"* ]]; then
        echo "A legacy ${method_name} process is already writing to the same output (PID ${process_pid})." >&2
        echo "Wait for it to finish; the workflow will not launch another fold set." >&2
        return 3
      fi
    done
  fi
  mkdir -p "$lock_root"
  exec 9>"$lock_file"
  if ! flock -n 9; then
    echo "A ${method_name} (${disease_name}) computation is already running." >&2
    echo "Lock: ${lock_file}" >&2
    return 3
  fi
  printf 'pid=%s\nstarted_at=%s\nmethod=%s\ndisease=%s\n' \
    "$$" "$(date -Iseconds)" "$method_name" "$disease_name" 1>&9
}

resolve_yu_python() {
  local candidate version
  local candidates=()
  [[ -n "$YU_PYTHON" ]] && candidates+=("$YU_PYTHON")
  candidates+=(
    "${DIR0}/software/conda/envs/yu_proteomic_repo_py39/bin/python"
    "${DIR0}/software/conda/envs/yu_proteomic_repo_py39/python.exe"
    "${DIR0}/software/python/yu_proteomic_repo_py39/python.exe"
  )
  for candidate in /mnt/c/Users/*/anaconda3/envs/yu_proteomic_repo_py39/python.exe \
                   /mnt/c/Users/*/miniconda3/envs/yu_proteomic_repo_py39/python.exe; do
    [[ -e "$candidate" ]] && candidates+=("$candidate")
  done
  command -v python3.9 >/dev/null 2>&1 && candidates+=("$(command -v python3.9)")
  for candidate in "${candidates[@]}"; do
    # Windows executables mounted through drvfs may not expose a Unix execute
    # bit even though WSL can invoke them directly.
    [[ -f "$candidate" ]] || continue
    version="$($candidate -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")' 2>/dev/null || true)"
    version="${version//$'\r'/}"
    if [[ "$version" == 3.9 ]]; then printf '%s\n' "$candidate"; return 0; fi
  done
  echo "Yu Python 3.9 runtime not found. Set YU_PYTHON or install the environment under <DIR0>/software." >&2
  return 1
}

help_text() {
  cat <<'EOF'
Usage
  yy score --h
  yy score --main [--workers=N]
  yy score --main --status
  yy score --inputs
  yy score --status
  yy score --preflight [--workers=N]
  yy score --project
  yy score --compute --confirm-compute [--workers=N] [--resume]
  yy score METHOD --inputs
  yy score METHOD --status
  yy score METHOD --preflight [--disease=NAME] [--workers=N]
  yy score METHOD --project [--source-root=PATH]
  yy score METHOD --compute --confirm-compute [--disease=NAME] [--workers=N] [--resume]

Methods
  pradeep-strict   Pradeep 1.5k article-method reproduction
  yu-strict        Yu article-method reproduction
  pradeep-fair     All-2,910 LASSO-logistic on the common five-year cohort
  yu-fair          All-2,910 LightGBM on the same common five-year cohort

Main CAD reproduction
  `yy score --main` computes or resumes all four CAD scores in this order:
  pradeep-strict, yu-strict, pradeep-fair, yu-fair.
  A successful main run also projects both strict models onto the
  common Yin/Yang participants, so all four outputs are immediately plot-ready.

Disease selection
  Default: cad
  Pradeep: cad, afib, hfail, ao_sten
  Yu: abdominal_aneurysm, atrial_fibrillation, aortic_valve_stenosis,
      cad, cardiomyopathy, deep_vein_thrombosis, heart_failure,
      intracerebral_hemorrhage, ischemic_stroke,
      peripheral_arterial_disease, pulmonary_embolism,
      subarachnoid_hemorrhage, thoracic_aneurysm,
      transient_ischemic_attack
  Common aliases are mapped automatically: afib/atrial_fibrillation,
  hfail/heart_failure and ao_sten/aortic_valve_stenosis.
  The two fair methods and --project are currently CAD-only because their
  common Yin/Yang cohort and endpoint contract are CAD-specific.

Huang-lab defaults
  DIR0=/mnt/d
  PHEDIR=/mnt/d/data/ukb/phe
  SCRIPT_ROOT=/mnt/d/scripts
  ANALYSIS_ROOT=/mnt/d/analysis
  YY_OUTDIR=/mnt/d/analysis/yy
  YY_SCORE_FOLD_ROOT=/mnt/d/analysis/yy/reference/cad_fivefold_v1

Strict-method raw inputs
  Pradeep raw protein: auto-detected prot.tab.gz under data.BIG, phe or ppp
  Pradeep 1.5k map:   <DIR0>/data.BIG/gwas/ppp/map.raw/olink_protein_map_1.5k_v1.tsv
  Yu raw protein:     reuses the same prot.tab.gz (a decompressed copy is not required)
  Yu raw phenotype:   auto-detected pheno.tsv.gz under PHEDIR
  Yu 3k assay map:    <DIR0>/data.BIG/gwas/ppp/olink_protein_map_3k_v1.tsv
  Yu processing date: official UKB Resource 1019, downloaded only if absent
                      and accepted only with the locked SHA-256 checksum
  Generated fold files: <YY_SCORE_FOLD_ROOT>/fold_assignment_{yin,yang}.csv
    If absent, fair preflight reconstructs them in memory from all.rds/prot.rds;
    fair compute installs them only after frozen count, balance and hash checks.

Privacy
  GitHub contains code, configuration and tests only. UKB rows, EIDs, protein
  matrices, fold manifests and fitted results remain outside the repository.

Derived-project boundary
  --project reads a completed model project under <ANALYSIS_ROOT> or
  <YY_OUTDIR>/score/source-projects. Coefficients, boosters and predictions
  are derived outputs; they are never searched for under PHEDIR or data.BIG.

Safety
  --project never refits a model. It reads the completed native project and
  projects its frozen parameters onto the common Yin/Yang participants.
  No model is fitted unless both --compute and --confirm-compute are present.
  --resume preserves completed folds and refuses partial fold products.

Examples
  yy score pradeep-strict --preflight
  yy score pradeep-strict --disease afib --preflight
  yy score pradeep-strict --compute --confirm-compute --workers=10 --resume
  yy score yu-strict --preflight
  yy score yu-strict --disease heart_failure --preflight
  yy score yu-strict --compute --confirm-compute --workers=10 --resume
  yy score pradeep-strict --project
  yy score yu-strict --project
  yy score pradeep-fair --compute --confirm-compute --workers=10 --resume
  yy score yu-fair --compute --confirm-compute --workers=10 --resume

CAD strict end-to-end
  yy score pradeep-strict --disease cad --compute --confirm-compute --workers=10 --resume
  yy score yu-strict --disease cad --compute --confirm-compute --workers=10 --resume
  yy score pradeep-strict --disease cad --project
  yy score yu-strict --disease cad --project
  yy plot --baseline --yy --traj --roc --bar --score pradeep-strict yu-strict --require-score

CAD fair end-to-end
  yy score pradeep-fair --preflight
  yy score yu-fair --preflight
  yy score pradeep-fair --compute --confirm-compute --workers=10 --resume
  yy score yu-fair --compute --confirm-compute --workers=10 --resume
  yy plot --baseline --yy --traj --roc --bar --score pradeep-fair yu-fair --require-score

  A completed fair --compute already writes plot-ready individual scores;
  do not run --project afterwards. Use --project instead of --compute only
  when connecting a previously completed common-cohort source project.

Main CAD reproduction in two commands
  yy score --main
  yy plot --main
EOF
}

is_pradeep_source() {
  local root="$1"
  [[ -f "${root}/outputs/ukbppp_cardiac_analysis_base.rds" && \
     -f "${root}/outputs/lasso/${PRADEEP_STRICT_PROJECT}_coefficients.csv" && \
     -f "${root}/outputs/lasso/${PRADEEP_STRICT_PROJECT}_predictions.csv" ]]
}

is_yu_source() {
  local root="$1"
  [[ -f "${root}/09_models/${YU_STRICT_PROJECT}__Protein.txt" && \
     -f "${root}/09_models/test_predictions.csv.gz" ]]
}

resolve_pradeep_source() {
  local explicit="${1:-}"
  local candidates=()
  [[ -n "$explicit" ]] && candidates+=("$explicit")
  [[ -n "${PRADEEP_STRICT_MODEL_ROOT:-}" ]] && candidates+=("${PRADEEP_STRICT_MODEL_ROOT}")
  candidates+=(
    "${SOURCE_PROJECT_ROOT}/pradeep-strict-${PRADEEP_STRICT_PROJECT}"
    "${PRADEEP_NATIVE_ROOT}/${PRADEEP_STRICT_PROJECT}"
  )
  local root
  for root in "${candidates[@]}"; do
    if is_pradeep_source "$root"; then printf '%s\n' "$root"; return 0; fi
  done
  echo "No completed Pradeep derived project found under the public analysis interface. Checked: ${candidates[*]}" >&2
  echo "Run the GitHub workflow with --compute, or expose an existing completed derived project with --source-root." >&2
  return 1
}

resolve_yu_source() {
  local explicit="${1:-}"
  local candidates=()
  [[ -n "$explicit" ]] && candidates+=("$explicit")
  [[ -n "${YU_STRICT_MODEL_ROOT:-}" ]] && candidates+=("${YU_STRICT_MODEL_ROOT}")
  candidates+=(
    "${SOURCE_PROJECT_ROOT}/yu-strict-${YU_STRICT_PROJECT}"
    "${YU_NATIVE_ROOT}/${YU_STRICT_PROJECT}"
  )
  local root
  for root in "${candidates[@]}"; do
    if is_yu_source "$root"; then printf '%s\n' "$root"; return 0; fi
  done
  echo "No completed Yu derived project found under the public analysis interface. Checked: ${candidates[*]}" >&2
  echo "Run the GitHub workflow with --compute, or expose an existing completed derived project with --source-root." >&2
  return 1
}

is_fair_source() {
  local root="$1"
  [[ -f "${root}/03_source_data/unified_yin_oof_scores.csv.gz" && \
     -f "${root}/03_source_data/unified_yang_ensemble_scores.csv.gz" && \
     -f "${root}/03_source_data/panel_b_unified_model_metrics.csv" ]]
}

resolve_fair_source() {
  local explicit="${1:-}"
  local candidates=()
  [[ -n "$explicit" ]] && candidates+=("$explicit")
  [[ -n "${YY_FAIR_MODEL_ROOT:-}" ]] && candidates+=("${YY_FAIR_MODEL_ROOT}")
  candidates+=("${SOURCE_PROJECT_ROOT}/fair-unified")
  local root
  for root in "${candidates[@]}"; do
    if is_fair_source "$root"; then printf '%s\n' "$root"; return 0; fi
  done
  echo "No completed fair-comparison derived project found. Checked: ${candidates[*]}" >&2
  echo "Run the selected fair method with --compute, or pass --source-root to a completed unified project." >&2
  return 1
}

prepare_common_projection_input() {
  local yin_participant="" yang_participant=""
  [[ -f "${COMMON_ROOT}/participants_yin.csv.gz" ]] && yin_participant="${COMMON_ROOT}/participants_yin.csv.gz"
  [[ -f "${COMMON_ROOT}/participants_yin.csv" ]] && yin_participant="${COMMON_ROOT}/participants_yin.csv"
  [[ -f "${COMMON_ROOT}/participants_yang.csv.gz" ]] && yang_participant="${COMMON_ROOT}/participants_yang.csv.gz"
  [[ -f "${COMMON_ROOT}/participants_yang.csv" ]] && yang_participant="${COMMON_ROOT}/participants_yang.csv"
  if [[ -n "$yin_participant" && -n "$yang_participant" && \
        -f "${COMMON_ROOT}/protein_features.csv" && \
        -f "${COMMON_ROOT}/protein_yin.f32" && \
        -f "${COMMON_ROOT}/protein_yang.f32" ]]; then
    echo "Existing common Yin/Yang projection input preserved: ${COMMON_ROOT}"
    return 0
  fi
  local prepare_args=(
    "--dir0=${DIR0}" "--phe-dir=${PHEDIR}" "--script-root=${SCRIPT_ROOT}"
    "--analysis-root=${ANALYSIS_ROOT}" "--yy-outdir=${YY_OUTDIR}"
    "--common-root=${COMMON_ROOT}" "--fold-root=${FOLD_ROOT}"
  )
  "$RSCRIPT" --vanilla "${PROJECT_ROOT}/R/01_prepare_fair_inputs.R" \
    --stage=prepare "${prepare_args[@]}"
}

normalize_method() {
  case "${1,,}" in
    pradeep-strict|pradeep_strict) echo pradeep-strict ;;
    yu-strict|yu_strict) echo yu-strict ;;
    pradeep-fair|pradeep_fair) echo pradeep-fair ;;
    yu-fair|yu_fair) echo yu-fair ;;
    *) return 1 ;;
  esac
}

normalize_disease_for_method() {
  local method="$1"
  local requested="${2,,}"
  requested="${requested//-/_}"
  case "$method" in
    pradeep-strict)
      case "$requested" in
        cad|coronary_artery_disease) echo cad ;;
        af|afib|atrial_fibrillation) echo afib ;;
        hf|hfail|heart_failure) echo hfail ;;
        as|ao_sten|aortic_stenosis|aortic_valve_stenosis) echo ao_sten ;;
        *)
          echo "Unsupported Pradeep disease: ${2}. Valid native IDs: cad, afib, hfail, ao_sten." >&2
          return 1
          ;;
      esac
      ;;
    yu-strict)
      case "$requested" in
        coronary_artery_disease) requested=cad ;;
        af|afib) requested=atrial_fibrillation ;;
        hf|hfail) requested=heart_failure ;;
        as|ao_sten|aortic_stenosis) requested=aortic_valve_stenosis ;;
      esac
      case "$requested" in
        abdominal_aneurysm|atrial_fibrillation|aortic_valve_stenosis|cad|cardiomyopathy|deep_vein_thrombosis|heart_failure|intracerebral_hemorrhage|ischemic_stroke|peripheral_arterial_disease|pulmonary_embolism|subarachnoid_hemorrhage|thoracic_aneurysm|transient_ischemic_attack)
          echo "$requested"
          ;;
        *)
          echo "Unsupported Yu disease: ${2}. Run 'yy score --h' for the 14 valid IDs." >&2
          return 1
          ;;
      esac
      ;;
    pradeep-fair|yu-fair)
      case "$requested" in
        cad|coronary_artery_disease) echo cad ;;
        *)
          echo "${method} is CAD-only because the locked fair-comparison cohort and endpoint are CAD-specific." >&2
          return 1
          ;;
      esac
      ;;
  esac
}

pradeep_strict_args() {
  printf '%s\0' \
    --profile huang \
    --dir0 "$DIR0" \
    --phe-dir "$PHEDIR" \
    --script-root "$SCRIPT_ROOT" \
    --helper-dir "${SCRIPT_ROOT}/0f" \
    --output-root "$PRADEEP_NATIVE_ROOT" \
    --analysis-project "$PRADEEP_STRICT_PROJECT" \
    --outcome "$PRADEEP_STRICT_PROJECT" \
    --all-rds "${PHEDIR}/Rdata/all.rds" \
    --raw-protein "$PRADEEP_STRICT_RAW_PROTEIN" \
    --protein-map "$PRADEEP_STRICT_PROTEIN_MAP" \
    --bed-file "$PRADEEP_STRICT_BED" \
    --allow-kinship-fallback
}

yu_strict_args() {
  printf '%s\0' \
    --path-config "${YY_OUTDIR}/score/config/yu_paths.json" \
    --reset-paths \
    --dir0 "$DIR0" \
    --phe-dir "$PHEDIR" \
    --analysis-root "$YU_NATIVE_ROOT" \
    --analysis-project "$YU_STRICT_PROJECT" \
    --disease "$YU_STRICT_PROJECT" \
    --raw-protein-file "$YU_STRICT_RAW_PROTEIN" \
    --phenotype-rds "${PHEDIR}/Rdata/all.rds" \
    --raw-phenotype-file "$YU_STRICT_RAW_PHENOTYPE" \
    --panel-mapping-file "$YU_STRICT_PANEL_MAP" \
    --supplement-workbook-file "$YU_STRICT_SUPPLEMENT_WORKBOOK" \
    --supplement-methods-file "$YU_STRICT_SUPPLEMENT_METHODS" \
    --olink-processing-start-date-file "$YU_STRICT_OLINK_DATES" \
    --rscript "$RSCRIPT" \
    --model-python "$YU_PYTHON" \
    --path-prompt off
}

method_status() {
  local method="$1"
  case "$method" in
    pradeep-strict)
      local source_root=""
      source_root="$(resolve_pradeep_source 2>/dev/null || true)"
      echo "METHOD ${method}"
      echo "  disease=${PRADEEP_STRICT_PROJECT}"
      echo "  native_project=$([[ -n "$source_root" ]] && echo AVAILABLE || echo MISSING)"
      [[ -n "$source_root" ]] && echo "  native_project_root=${source_root}"
      echo "  fitted_outcome=${PRADEEP_STRICT_PROJECT}"
      if [[ "$PRADEEP_STRICT_PROJECT" == cad ]]; then
        echo "  common_input=$([[ -f "${COMMON_ROOT}/COMPLETE" ]] && echo COMPLETE || echo MISSING)"
        echo "  projection=$([[ -f "${YY_OUTDIR}/score/${method}/COMPLETE" ]] && echo COMPLETE || echo MISSING)"
        if [[ -f "${YY_OUTDIR}/score/${method}/metrics.csv" ]]; then cat "${YY_OUTDIR}/score/${method}/metrics.csv"; fi
      else
        echo "  common_input=NOT_APPLICABLE_CAD_ONLY"
        echo "  projection=NOT_APPLICABLE_CAD_ONLY"
      fi
      ;;
    yu-strict)
      local source_root=""
      source_root="$(resolve_yu_source 2>/dev/null || true)"
      echo "METHOD ${method}"
      echo "  disease=${YU_STRICT_PROJECT}"
      echo "  native_project=$([[ -n "$source_root" ]] && echo AVAILABLE || echo MISSING)"
      [[ -n "$source_root" ]] && echo "  native_project_root=${source_root}"
      echo "  fitted_outcome=${YU_STRICT_PROJECT}"
      if [[ "$YU_STRICT_PROJECT" == cad ]]; then
        echo "  common_input=$([[ -f "${COMMON_ROOT}/COMPLETE" ]] && echo COMPLETE || echo MISSING)"
        echo "  projection=$([[ -f "${YY_OUTDIR}/score/${method}/COMPLETE" ]] && echo COMPLETE || echo MISSING)"
        if [[ -f "${YY_OUTDIR}/score/${method}/metrics.csv" ]]; then cat "${YY_OUTDIR}/score/${method}/metrics.csv"; fi
      else
        echo "  common_input=NOT_APPLICABLE_CAD_ONLY"
        echo "  projection=NOT_APPLICABLE_CAD_ONLY"
      fi
      ;;
    pradeep-fair|yu-fair)
      local root="${YY_OUTDIR}/score/${method}"
      local source_root=""
      source_root="$(resolve_fair_source 2>/dev/null || true)"
      echo "METHOD ${method}"
      echo "  derived_source_project=$([[ -n "$source_root" ]] && echo AVAILABLE || echo MISSING)"
      [[ -n "$source_root" ]] && echo "  derived_source_root=${source_root}"
      echo "  common_input=$([[ -f "${COMMON_ROOT}/COMPLETE" ]] && echo COMPLETE || echo MISSING)"
      echo "  final=$([[ -f "${root}/COMPLETE" ]] && echo COMPLETE || echo MISSING)"
      if [[ -f "${root}/COMPLETE" ]] && grep -q 'action=project_existing_scores_without_refit' "${root}/COMPLETE"; then
        echo "  existing_project_projection=COMPLETE"
      else
        for fold in 1 2 3 4 5; do
          printf '  fold%02d=%s\n' "$fold" "$([[ -f "${root}/models/fold$(printf '%02d' "$fold")_COMPLETE" ]] && echo COMPLETE || echo MISSING)"
        done
      fi
      if [[ -f "${root}/metrics.csv" ]]; then
        cat "${root}/metrics.csv"
      fi
      ;;
  esac
}

method_inputs() {
  local method="$1"
  echo "METHOD ${method}"
  case "$method" in
    pradeep-strict)
      echo "  GITHUB_METHOD_CODE=${SCRIPT_ROOT}/ukb/pradeep/pradeep.sh"
      echo "  FITTED_OUTCOME=${PRADEEP_STRICT_PROJECT}"
      echo "  RAW_PHENOTYPE=${PHEDIR}/Rdata/all.rds"
      echo "  RAW_PROTEIN=${PRADEEP_STRICT_RAW_PROTEIN}"
      echo "  RAW_PROTEIN_MAP=${PRADEEP_STRICT_PROTEIN_MAP}"
      echo "  RAW_PROTEIN_BED=${PRADEEP_STRICT_BED}"
      echo "  DERIVED_NATIVE_PROJECT=$(resolve_pradeep_source 2>/dev/null || echo MISSING)"
      echo "  NATIVE_PARENT_PROJECT=${PRADEEP_NATIVE_ROOT}"
      echo "  NEW_COMPUTE_OUTPUT=${PRADEEP_NATIVE_ROOT}/${PRADEEP_STRICT_PROJECT}"
      echo "  DERIVED_COMMON_INPUT=${COMMON_ROOT}"
      if [[ "$PRADEEP_STRICT_PROJECT" == cad ]]; then
        echo "  YY_PROJECTION_OUTPUT=${YY_OUTDIR}/score/${method}"
      else
        echo "  YY_PROJECTION_OUTPUT=NOT_APPLICABLE_CAD_ONLY"
      fi
      ;;
    yu-strict)
      echo "  GITHUB_METHOD_CODE=${SCRIPT_ROOT}/ukb/yu/yu.sh"
      echo "  FITTED_OUTCOME=${YU_STRICT_PROJECT}"
      echo "  RAW_PHENOTYPE=${PHEDIR}/Rdata/all.rds"
      echo "  RAW_PROTEIN=${YU_STRICT_RAW_PROTEIN}"
      echo "  RAW_PHENOTYPE_TABLE=${YU_STRICT_RAW_PHENOTYPE}"
      echo "  RAW_PROTEIN_MAP=${YU_STRICT_PANEL_MAP}"
      echo "  SOURCE_AUDIT_SUPPLEMENT_WORKBOOK=${YU_STRICT_SUPPLEMENT_WORKBOOK}"
      echo "  SOURCE_AUDIT_SUPPLEMENT_METHODS=${YU_STRICT_SUPPLEMENT_METHODS}"
      echo "  NATIVE_REBUILD_OLINK_DATES=${YU_STRICT_OLINK_DATES}"
      echo "  NATIVE_REBUILD_OLINK_DATES_STATUS=$([[ -f "$YU_STRICT_OLINK_DATES" ]] && echo AVAILABLE || echo AUTO_DOWNLOAD_ON_PREFLIGHT_OR_COMPUTE)"
      echo "  NATIVE_REBUILD_OLINK_DATES_URL=${YU_OLINK_DATES_URL}"
      echo "  NATIVE_REBUILD_OLINK_DATES_SHA256=${YU_OLINK_DATES_SHA256}"
      echo "  R_RUNTIME=${RSCRIPT}"
      echo "  MODEL_PYTHON_RUNTIME=${YU_PYTHON:-AUTO_RESOLVE_PYTHON_3_9}"
      echo "  DERIVED_NATIVE_PROJECT=$(resolve_yu_source 2>/dev/null || echo MISSING)"
      echo "  NATIVE_PARENT_PROJECT=${YU_NATIVE_ROOT}"
      echo "  NEW_COMPUTE_OUTPUT=${YU_NATIVE_ROOT}/${YU_STRICT_PROJECT}"
      echo "  DERIVED_COMMON_INPUT=${COMMON_ROOT}"
      if [[ "$YU_STRICT_PROJECT" == cad ]]; then
        echo "  YY_PROJECTION_OUTPUT=${YY_OUTDIR}/score/${method}"
      else
        echo "  YY_PROJECTION_OUTPUT=NOT_APPLICABLE_CAD_ONLY"
      fi
      ;;
    pradeep-fair|yu-fair)
      echo "  GITHUB_METHOD_CODE=${PROJECT_ROOT}"
      echo "  RAW_PHENOTYPE=${PHEDIR}/Rdata/all.rds"
      echo "  RAW_PROTEIN=${PHEDIR}/Rdata/prot.rds"
      echo "  GENERATED_OR_VALIDATED_FOLD_ROOT=${FOLD_ROOT}"
      echo "  DERIVED_SOURCE_PROJECT=$(resolve_fair_source 2>/dev/null || echo MISSING)"
      echo "  DERIVED_COMMON_INPUT=${COMMON_ROOT}"
      echo "  YY_PROJECTION_OUTPUT=${YY_OUTDIR}/score/${method}"
      ;;
  esac
}

run_score_child() {
  local args=("$@")
  printf 'MAIN_CAD_STEP'
  printf ' %q' "${args[@]}"
  printf '\n'
  if [[ "${YY_SCORE_MAIN_PLAN_ONLY:-0}" == 1 ]]; then return 0; fi
  bash "${PROJECT_ROOT}/score.sh" "${args[@]}"
}

run_main_cad_score() {
  local action="$1" workers="$2" confirm="$3" resume="$4"
  local methods=(pradeep-strict yu-strict pradeep-fair yu-fair)
  local method args
  case "$action" in
    status)
      for method in "${methods[@]}"; do method_status "$method"; done
      ;;
    inputs)
      for method in "${methods[@]}"; do method_inputs "$method"; done
      ;;
    preflight)
      for method in "${methods[@]}"; do
        args=("$method" --preflight --workers "$workers")
        [[ "$method" == *-strict ]] && args+=(--disease cad)
        run_score_child "${args[@]}"
      done
      ;;
    project)
      for method in "${methods[@]}"; do
        args=("$method" --project)
        [[ "$method" == *-strict ]] && args+=(--disease cad)
        run_score_child "${args[@]}"
      done
      ;;
    compute)
      ((confirm == 1)) || {
        echo "Refusing all-method fit: add --confirm-compute." >&2
        return 2
      }
      for method in "${methods[@]}"; do
        args=("$method" --compute --confirm-compute --workers "$workers")
        [[ "$method" == *-strict ]] && args+=(--disease cad)
        ((resume == 1)) && args+=(--resume)
        run_score_child "${args[@]}"
      done
      for method in pradeep-strict yu-strict; do
        run_score_child "$method" --disease cad --project
      done
      if [[ "${YY_SCORE_MAIN_PLAN_ONLY:-0}" == 1 ]]; then
        echo "MAIN_CAD_SCORE_PLAN_COMPLETE"
        return 0
      fi
      local missing=()
      for method in "${methods[@]}"; do
        [[ -f "${YY_OUTDIR}/score/${method}/COMPLETE" ]] || missing+=("$method")
      done
      ((${#missing[@]} == 0)) || {
        echo "Main CAD score contract failed; missing COMPLETE: ${missing[*]}" >&2
        return 1
      }
      mkdir -p "${YY_OUTDIR}/score"
      printf 'status\tmethods\tcompleted_at\nCOMPLETE\tpradeep-strict|yu-strict|pradeep-fair|yu-fair\t%s\n' \
        "$(date -Iseconds)" >"${YY_OUTDIR}/score/MAIN_CAD_COMPLETE.tsv"
      echo "MAIN_CAD_SCORE_COMPLETE"
      echo "NEXT: yy plot --main"
      ;;
    *)
      echo "Unsupported main CAD score action: ${action}" >&2
      return 2
      ;;
  esac
}

if (($# == 0)); then help_text; exit 0; fi
case "$1" in --h|-h|--help|help) help_text; exit 0 ;; esac

if [[ "$1" == --* ]]; then
  main_action=""
  main_workers="${YY_SCORE_WORKERS:-10}"
  main_confirm=0
  main_resume=0
  main_preset=0
  main_disease="cad"
  while (($#)); do
    case "$1" in
      --main) main_preset=1; shift ;;
      --status) main_action=status; shift ;;
      --inputs) main_action=inputs; shift ;;
      --preflight) main_action=preflight; shift ;;
      --project) main_action=project; shift ;;
      --compute) main_action=compute; shift ;;
      --confirm-compute) main_confirm=1; shift ;;
      --resume) main_resume=1; shift ;;
      --workers=*) main_workers="${1#*=}"; shift ;;
      --workers)
        [[ $# -ge 2 ]] || { echo "--workers requires a value." >&2; exit 2; }
        main_workers="$2"; shift 2 ;;
      --disease=*) main_disease="${1#*=}"; shift ;;
      --disease)
        [[ $# -ge 2 ]] || { echo "--disease requires a value." >&2; exit 2; }
        main_disease="$2"; shift 2 ;;
      --h|-h|--help|help) help_text; exit 0 ;;
      *) echo "Unknown main CAD yy score option: $1" >&2; exit 2 ;;
    esac
  done
  if ((main_preset == 1)) && [[ -z "$main_action" ]]; then
    main_action=compute
    main_confirm=1
    main_resume=1
  fi
  [[ -n "$main_action" ]] || { help_text; exit 0; }
  [[ "$main_workers" =~ ^[1-9][0-9]*$ ]] || { echo "--workers must be a positive integer." >&2; exit 2; }
  case "${main_disease,,}" in
    cad|coronary_artery_disease) ;;
    *) echo "The no-METHOD yy score workflow is the locked four-score CAD reproduction. Select METHOD explicitly for another disease." >&2; exit 2 ;;
  esac
  echo "yy command: score"
  echo "yy score preset: main-cad-four-score"
  echo "yy score action: ${main_action}"
  echo "yy score methods: pradeep-strict yu-strict pradeep-fair yu-fair"
  run_main_cad_score "$main_action" "$main_workers" "$main_confirm" "$main_resume"
  exit $?
fi

method="$(normalize_method "$1")" || {
  echo "Unknown score method: $1" >&2
  help_text >&2
  exit 2
}
shift

action="status"
workers="${YY_SCORE_WORKERS:-8}"
confirm=0
resume=0
source_root=""
disease_requested="cad"
while (($#)); do
  case "$1" in
    --status) action=status; shift ;;
    --inputs) action=inputs; shift ;;
    --preflight) action=preflight; shift ;;
    --project) action=project; shift ;;
    --compute) action=compute; shift ;;
    --confirm-compute) confirm=1; shift ;;
    --resume) resume=1; shift ;;
    --disease=*) disease_requested="${1#*=}"; shift ;;
    --disease)
      [[ $# -ge 2 ]] || { echo "--disease requires a value." >&2; exit 2; }
      disease_requested="$2"; shift 2 ;;
    --workers=*) workers="${1#*=}"; shift ;;
    --workers)
      [[ $# -ge 2 ]] || { echo "--workers requires a value." >&2; exit 2; }
      workers="$2"; shift 2 ;;
    --source-root=*) source_root="${1#*=}"; shift ;;
    --source-root)
      [[ $# -ge 2 ]] || { echo "--source-root requires a path." >&2; exit 2; }
      source_root="$2"; shift 2 ;;
    --h|-h|--help) help_text; exit 0 ;;
    *) echo "Unknown yy score option: $1" >&2; exit 2 ;;
  esac
done
[[ "$workers" =~ ^[1-9][0-9]*$ ]] || { echo "--workers must be a positive integer." >&2; exit 2; }
disease="$(normalize_disease_for_method "$method" "$disease_requested")" || exit 2
case "$method" in
  pradeep-strict) PRADEEP_STRICT_PROJECT="$disease" ;;
  yu-strict) YU_STRICT_PROJECT="$disease" ;;
esac

echo "yy command: score"
echo "yy score method: ${method}"
echo "yy score action: ${action}"
echo "yy score disease: ${disease}"
echo "DIR0=${DIR0}"
echo "PHEDIR=${PHEDIR}"
echo "SCRIPT_ROOT=${SCRIPT_ROOT}"
echo "ANALYSIS_ROOT=${ANALYSIS_ROOT}"
echo "YY_OUTDIR=${YY_OUTDIR}"

if [[ "$action" == status ]]; then method_status "$method"; exit $?; fi
if [[ "$action" == inputs ]]; then method_inputs "$method"; exit 0; fi
if [[ "$action" == project && "$disease" != cad ]]; then
  echo "--project is currently CAD-only: the common Yin/Yang participant and endpoint contract is CAD-specific." >&2
  echo "The native ${disease} reproduction remains available through --preflight and --compute." >&2
  exit 2
fi
if [[ "$action" == compute ]]; then
  acquire_compute_lock "$method" "$disease" || exit $?
fi
[[ -x "$RSCRIPT" ]] || { echo "Missing Rscript: ${RSCRIPT}" >&2; exit 2; }
"$RSCRIPT" --vanilla "${PROJECT_ROOT}/tests/test_contracts.R"
"$RSCRIPT" --vanilla "${PROJECT_ROOT}/tests/test_fold_generation.R"

case "$method" in
  pradeep-strict)
    if [[ "$action" == project ]]; then
      resolved_source="$(resolve_pradeep_source "$source_root")"
      echo "STRICT_SOURCE_ROOT=${resolved_source}"
      prepare_common_projection_input
      exec "$RSCRIPT" --vanilla "${PROJECT_ROOT}/R/05_project_pradeep_strict.R" \
        --stage=project --source-root="${resolved_source}" --common-root="${COMMON_ROOT}" \
        --dir0="${DIR0}" --phe-dir="${PHEDIR}" --script-root="${SCRIPT_ROOT}" \
        --analysis-root="${ANALYSIS_ROOT}" --yy-outdir="${YY_OUTDIR}" \
        --output-root="${YY_OUTDIR}/score/pradeep-strict"
    fi
    command_path="${SCRIPT_ROOT}/ukb/pradeep/pradeep.sh"
    [[ -f "$command_path" ]] || { echo "Missing Pradeep entrypoint: ${command_path}" >&2; exit 2; }
    strict_args=()
    while IFS= read -r -d '' value; do strict_args+=("$value"); done < <(pradeep_strict_args)
    if [[ "$action" == preflight ]]; then
      exec bash "$command_path" --step preflight --panel 1.5k --workers "$workers" "${strict_args[@]}"
    fi
    ((confirm == 1)) || { echo "Refusing fit: add --confirm-compute." >&2; exit 2; }
    args=(--step 1,2,3,4 --panel 1.5k --workers "$workers" "${strict_args[@]}")
    ((resume == 1)) && args+=(--resume)
    exec bash "$command_path" "${args[@]}"
    ;;
  yu-strict)
    if [[ "$action" == project ]]; then
      resolved_source="$(resolve_yu_source "$source_root")"
      echo "STRICT_SOURCE_ROOT=${resolved_source}"
      prepare_common_projection_input
      YU_PYTHON="$(resolve_yu_python)" || exit 2
      python_script="${PROJECT_ROOT}/python/05_project_yu_strict.py"
      if [[ "$YU_PYTHON" == /mnt/c/* ]]; then
        to_win() { wslpath -m "$1"; }
        exec "$YU_PYTHON" "$(to_win "$python_script")" --stage project \
          --source-root "$(to_win "$resolved_source")" --common-root "$(to_win "$COMMON_ROOT")" \
          --output-root "$(to_win "${YY_OUTDIR}/score/yu-strict")" \
          --raw-protein-file "$(to_win "$YU_STRICT_RAW_PROTEIN")"
      fi
      exec "$YU_PYTHON" "$python_script" --stage project \
        --source-root "$resolved_source" --common-root "$COMMON_ROOT" \
        --output-root "${YY_OUTDIR}/score/yu-strict" \
        --raw-protein-file "$YU_STRICT_RAW_PROTEIN"
    fi
    command_path="${SCRIPT_ROOT}/ukb/yu/yu.sh"
    [[ -f "$command_path" ]] || { echo "Missing Yu entrypoint: ${command_path}" >&2; exit 2; }
    ensure_yu_olink_dates || exit $?
    if [[ -z "$YU_PYTHON" ]]; then YU_PYTHON="$(resolve_yu_python)" || exit 2; fi
    strict_args=()
    while IFS= read -r -d '' value; do strict_args+=("$value"); done < <(yu_strict_args)
    if [[ "$action" == preflight ]]; then
      YU_DOCTOR_SCOPE=prediction exec bash "$command_path" doctor "${strict_args[@]}"
    fi
    ((confirm == 1)) || { echo "Refusing fit: add --confirm-compute." >&2; exit 2; }
    args=(1-4 --workers "$workers" "${strict_args[@]}")
    ((resume == 1)) && args+=(--resume)
    exec bash "$command_path" "${args[@]}"
    ;;
  pradeep-fair|yu-fair)
    if [[ "$action" == project ]]; then
      resolved_source="$(resolve_fair_source "$source_root")"
      echo "FAIR_SOURCE_ROOT=${resolved_source}"
      exec "$RSCRIPT" --vanilla "${PROJECT_ROOT}/R/05_project_fair_existing.R" \
        --method="${method}" --source-root="${resolved_source}" \
        --output-root="${YY_OUTDIR}/score/${method}"
    fi
    prepare_args=(
      "--dir0=${DIR0}" "--phe-dir=${PHEDIR}" "--script-root=${SCRIPT_ROOT}"
      "--analysis-root=${ANALYSIS_ROOT}" "--yy-outdir=${YY_OUTDIR}"
      "--common-root=${COMMON_ROOT}" "--fold-root=${FOLD_ROOT}"
    )
    if [[ "$action" == preflight ]]; then
      exec "$RSCRIPT" --vanilla "${PROJECT_ROOT}/R/01_prepare_fair_inputs.R" \
        --stage=preflight "${prepare_args[@]}"
    fi
    ((confirm == 1)) || { echo "Refusing fit: add --confirm-compute." >&2; exit 2; }
    output_root="${YY_OUTDIR}/score/${method}"
    if [[ -f "${output_root}/COMPLETE" ]] &&
       grep -q 'action=project_existing_scores_without_refit' "${output_root}/COMPLETE"; then
      echo "Replacing projected ${method} outputs in place with newly computed outputs."
      for name in \
        scores_yin_oof.csv.gz scores_yang.csv.gz metrics.csv fold_metrics.csv \
        roc_curves.csv.gz input_manifest.csv COMPLETE; do
        [[ -f "${output_root}/${name}" ]] && unlink "${output_root}/${name}"
      done
    fi
    "$RSCRIPT" --vanilla "${PROJECT_ROOT}/R/01_prepare_fair_inputs.R" \
      --stage=prepare "${prepare_args[@]}"
    mkdir -p "${output_root}/logs"
    if [[ "$method" == pradeep-fair ]]; then
      outer_workers="${YY_SCORE_OUTER_WORKERS:-5}"
      ((outer_workers > workers)) && outer_workers="$workers"
      ((outer_workers > 5)) && outer_workers=5
      active=()
      active_folds=()
      failed=0
      for fold in 1 2 3 4 5; do
        OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
          "$RSCRIPT" --vanilla "${PROJECT_ROOT}/R/02_pradeep_fair.R" \
          --stage=fold --fold="$fold" --dir0="${DIR0}" --phe-dir="${PHEDIR}" \
          --script-root="${SCRIPT_ROOT}" --analysis-root="${ANALYSIS_ROOT}" \
          --yy-outdir="${YY_OUTDIR}" --common-root="${COMMON_ROOT}" \
          --output-root="${output_root}" >"${output_root}/logs/fold${fold}.log" 2>&1 &
        active+=("$!"); active_folds+=("$fold")
        if ((${#active[@]} == outer_workers || fold == 5)); then
          progress_seconds="${YY_SCORE_PROGRESS_SECONDS:-10}"
          [[ "$progress_seconds" =~ ^[1-9][0-9]*$ ]] || {
            echo "YY_SCORE_PROGRESS_SECONDS must be a positive integer." >&2
            exit 2
          }
          while :; do
            running=0
            status_line="$(date '+%Y-%m-%d %H:%M:%S') |"
            for index in "${!active[@]}"; do
              fold_id="${active_folds[$index]}"
              marker="${output_root}/models/fold$(printf '%02d' "$fold_id")_COMPLETE"
              if [[ -f "$marker" ]]; then
                fold_state="COMPLETE"
              elif kill -0 "${active[$index]}" 2>/dev/null; then
                fold_state="RUNNING"
                running=$((running + 1))
              else
                fold_state="EXITED"
              fi
              status_line+=" fold${fold_id}=${fold_state}"
            done
            echo "$status_line"
            ((running == 0)) && break
            sleep "$progress_seconds"
          done
          for index in "${!active[@]}"; do
            if ! wait "${active[$index]}"; then
              echo "Pradeep fair fold ${active_folds[$index]} failed; see ${output_root}/logs/fold${active_folds[$index]}.log" >&2
              failed=1
            fi
          done
          active=(); active_folds=()
        fi
      done
      ((failed == 0)) || exit 1
    else
      yu_python="$(resolve_yu_python)" || exit 2
      python_script="${PROJECT_ROOT}/python/03_yu_fair.py"
      if [[ "$yu_python" == /mnt/c/* ]]; then
        to_win() { wslpath -m "$1"; }
        for fold in 1 2 3 4 5; do
          "$yu_python" "$(to_win "$python_script")" --stage fold --fold "$fold" \
            --common-root "$(to_win "$COMMON_ROOT")" --output-root "$(to_win "$output_root")" \
            --config "$(to_win "${PROJECT_ROOT}/config/fair.json")" --workers "$workers" \
            >"${output_root}/logs/fold${fold}.log" 2>&1
        done
      else
        for fold in 1 2 3 4 5; do
          "$yu_python" "$python_script" --stage fold --fold "$fold" \
            --common-root "$COMMON_ROOT" --output-root "$output_root" \
            --config "${PROJECT_ROOT}/config/fair.json" --workers "$workers" \
            >"${output_root}/logs/fold${fold}.log" 2>&1
        done
      fi
    fi
    "$RSCRIPT" --vanilla "${PROJECT_ROOT}/R/04_finalize_fair.R" \
      --stage=finalize --method="${method}" --dir0="${DIR0}" --phe-dir="${PHEDIR}" \
      --script-root="${SCRIPT_ROOT}" --analysis-root="${ANALYSIS_ROOT}" \
      --yy-outdir="${YY_OUTDIR}" --common-root="${COMMON_ROOT}" \
      --output-root="${output_root}"
    ;;
esac
