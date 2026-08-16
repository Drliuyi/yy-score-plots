#!/usr/bin/env bash
set -euo pipefail

global_help() {
  cat <<'EOF'
Yin–Yang Score Plots

Usage
  yy --h
  yy plot --h
  yy plot --main
  yy score --h
  yy score --aa

Commands
  plot    Draw Yin/Yang trajectories and optional ROC panels.
  score   Compute four Pradeep/Yu scores and downstream area-AUC analysis.
EOF
}

plot_help() {
  cat <<'EOF'
Usage
  yy plot --status
  yy plot --main
  yy plot [--yin|--yang|--yy] [--baseline|--diagnosis]
          [--traj] [--roc] [--bar] [--mean] [--sd] [--sd_line]
          [--proteins=P1,P2,...]
          [--score METHOD...]
          [--adj=TOKEN+TOKEN+...]
          [--require-score]
          [--recompute]

Score methods
  pradeep-strict   Pradeep article-method reproduction
  yu-strict        Yu article-method reproduction
  pradeep-fair     Pradeep LASSO on the common 5-year/five-fold cohort
  yu-fair          Yu LightGBM on the same common 5-year/five-fold cohort

  `pradeep` and `yu` remain short aliases for the fair methods only.
  There is no plot version switch. Any four methods may be drawn together.

Main CAD preset
  --main draws the locked baseline-centred Yin/Yang trajectory plus ROC figure
  for pradeep-strict, yu-strict, pradeep-fair and yu-fair. It requires the
  plot-ready products created by the main `yy score` command.

Adjustment
  --adj=raw
  --adj=age+sex
  --adj=age+sex+pc+center+tdi
  --adj=age+sex+le8+medications

  Tokens are independent. For example, age+sex+le8+medications does not add
  PC, centre or TDI. Available tokens are:
  age, sex, pc (PC1+PC2), center, tdi, le8, medications.

Side and anchor
  --yin / --yang / --yy     default: --yy
  --baseline                Yang left of zero, Yin right of zero; default
  --diagnosis               Yin left of zero, Yang right of zero

Layers
  --traj                     trajectory panel
  --roc                      incident CAD ROC panel; invalid with --yang alone
  --bar                      participant-count bars
  --mean                     point-wise mean labels
  --sd                       observed point-wise SD labels; raw only
  --sd_line                  observed mean +/- SD lines; raw only
  --recompute                replay frozen models as a parameter audit; never refits
  --require-score            require completed yy score outputs (the default public-path contract)

ROC interpretation
  Requested proteins and fair scores use the common five-fold 5-year IPCW AUC.
  Strict scores retain their own native held-out eventual-event AUC. If strict
  and fair curves are juxtaposed, the figure is explicitly marked descriptive;
  those AUCs are not a paired comparison.

Examples
  yy plot --main
  yy plot --status
  yy plot --baseline --yy --traj --bar \
    --proteins=GDF15,PCSK9,NTPROBNP \
    --score pradeep-strict yu-strict pradeep-fair yu-fair
  yy plot --diagnosis --yy --traj --roc \
    --proteins=GDF15,PCSK9 --score pradeep-fair yu-fair
  yy plot --baseline --yy --traj --roc \
    --adj=age+sex+le8+medications \
    --proteins=GDF15,PCSK9 --score pradeep-fair yu-fair

Output
  D:/analysis/yy/plot/custom/<configuration-specific folder>
  One command writes one combined figure directory.
EOF
}

if (($# == 0)); then global_help; exit 0; fi
case "$1" in --h|-h|--help|help) global_help; exit 0 ;; esac

command_name="$1"
shift
DIR0="${DIR0:-/mnt/d}"
PHEDIR="${PHEDIR:-${DIR0}/data/ukb/phe}"
ENTRY_PATH="${BASH_SOURCE[0]}"
if [[ -L "${ENTRY_PATH}" ]]; then
  LINK_TARGET="$(readlink "${ENTRY_PATH}")"
  if [[ "${LINK_TARGET}" == /* ]]; then
    ENTRY_PATH="${LINK_TARGET}"
  else
    ENTRY_PATH="$(dirname "${ENTRY_PATH}")/${LINK_TARGET}"
  fi
fi
CHECKOUT_ROOT="$(cd -L "$(dirname "${ENTRY_PATH}")/.." && pwd -L)"
SCRIPT_ROOT="${SCRIPT_ROOT:-${CHECKOUT_ROOT}}"
ANALYSIS_ROOT="${ANALYSIS_ROOT:-${DIR0}/analysis}"
YY_OUTDIR="${YY_OUTDIR:-${ANALYSIS_ROOT}/yy}"

case "$command_name" in
  score)
    score_entry="${SCRIPT_ROOT}/yy/score/score.sh"
    [[ -f "$score_entry" ]] || { echo "Missing score entrypoint: ${score_entry}" >&2; exit 2; }
    DIR0="$DIR0" PHEDIR="$PHEDIR" SCRIPT_ROOT="$SCRIPT_ROOT" ANALYSIS_ROOT="$ANALYSIS_ROOT" \
      YY_OUTDIR="$YY_OUTDIR" exec bash "$score_entry" "$@"
    ;;
  plot)
    for arg in "$@"; do case "$arg" in --h|-h|--help|help) plot_help; exit 0 ;; esac; done
    side=yy; side_count=0; anchor=baseline; anchor_count=0
    status=0; main=0; traj=0; roc=0; bar=0; mean=0; sd=0; sd_line=0; recompute=0; require_score=0
    proteins=""; adjustment=raw; scores=()
    while (($#)); do
      case "$1" in
        --status) status=1; shift ;;
        --main) main=1; shift ;;
        --yin) side=yin; side_count=$((side_count + 1)); shift ;;
        --yang) side=yang; side_count=$((side_count + 1)); shift ;;
        --yy) side=yy; side_count=$((side_count + 1)); shift ;;
        --baseline) anchor=baseline; anchor_count=$((anchor_count + 1)); shift ;;
        --diagnosis) anchor=diagnosis; anchor_count=$((anchor_count + 1)); shift ;;
        --traj) traj=1; shift ;;
        --roc) roc=1; shift ;;
        --bar) bar=1; shift ;;
        --mean) mean=1; shift ;;
        --sd) sd=1; shift ;;
        --sd_line|--sd-line) sd_line=1; shift ;;
        --recompute) recompute=1; shift ;;
        --require-score) require_score=1; shift ;;
        --score=*) scores+=("${1#*=}"); shift ;;
        --score)
          shift; before=${#scores[@]}
          while (($#)) && [[ "$1" != --* ]]; do scores+=("$1"); shift; done
          ((${#scores[@]} > before)) || { echo "--score requires at least one method." >&2; exit 2; }
          ;;
        --proteins=*) proteins="${1#*=}"; shift ;;
        --proteins)
          [[ $# -ge 2 ]] || { echo "--proteins requires a comma-separated list." >&2; exit 2; }
          proteins="$2"; shift 2 ;;
        -proteins=*|-proteins)
          echo "Use --proteins with two leading hyphens." >&2; exit 2 ;;
        --adj=*) adjustment="${1#*=}"; shift ;;
        --adj)
          [[ $# -ge 2 ]] || { echo "--adj requires a value such as age+sex." >&2; exit 2; }
          adjustment="$2"; shift 2 ;;
        --v|--v=*)
          echo "--v was retired. Choose explicit methods with --score pradeep-strict/yu-strict/pradeep-fair/yu-fair." >&2
          exit 2 ;;
        *) echo "Unknown yy plot option: $1" >&2; echo "Run: yy plot --h" >&2; exit 2 ;;
      esac
    done
    if ((main == 1)); then
      ((status == 0)) || { echo "--main cannot be combined with --status." >&2; exit 2; }
      ((side_count == 0 && anchor_count == 0)) || { echo "--main locks --yy and --baseline; omit explicit side/anchor options." >&2; exit 2; }
      [[ -z "$proteins" && ${#scores[@]} -eq 0 ]] || { echo "--main locks the four score methods; omit --proteins and --score." >&2; exit 2; }
      [[ "$adjustment" == raw ]] || { echo "--main is the locked raw figure; use the full yy plot syntax for adjusted figures." >&2; exit 2; }
      ((mean == 0 && sd == 0 && sd_line == 0 && recompute == 0)) || { echo "--main cannot be combined with mean/SD/recompute modifiers." >&2; exit 2; }
      side=yy
      anchor=baseline
      traj=1
      roc=1
      bar=1
      require_score=1
      scores=(pradeep-strict yu-strict pradeep-fair yu-fair)
    fi
    ((side_count <= 1)) || { echo "Choose only one of --yin, --yang, or --yy." >&2; exit 2; }
    ((anchor_count <= 1)) || { echo "Choose only one of --baseline or --diagnosis." >&2; exit 2; }
    renderer="${SCRIPT_ROOT}/yy/R/render_plot.R"
    rscript="${RSCRIPT:-/opt/R/4.3.2/bin/Rscript}"
    [[ -f "$renderer" ]] || { echo "Missing plot renderer: ${renderer}" >&2; exit 2; }
    [[ -x "$rscript" ]] || { echo "Missing Rscript: ${rscript}" >&2; exit 2; }
    if ((status == 1)); then
      DIR0="$DIR0" SCRIPT_ROOT="$SCRIPT_ROOT" ANALYSIS_ROOT="$ANALYSIS_ROOT" YY_OUTDIR="$YY_OUTDIR" \
        exec "$rscript" --vanilla "$renderer" --status
    fi
    ((traj == 1 || roc == 1)) || { echo "Choose --traj and/or --roc." >&2; exit 2; }
    [[ -n "$proteins" || ${#scores[@]} -gt 0 ]] || { echo "Choose --proteins and/or --score." >&2; exit 2; }
    [[ "$side" != "yang" || "$roc" -eq 0 ]] || { echo "--roc is invalid with --yang alone." >&2; exit 2; }
    output_root="${YY_PLOT_OUTPUT_ROOT:-${YY_OUTDIR}/plot}"
    call_args=(
      "--side=${side}" "--anchor=${anchor}" "--adj=${adjustment}"
      "--traj=${traj}" "--roc=${roc}" "--bar=${bar}" "--mean=${mean}"
      "--sd=${sd}" "--sd-line=${sd_line}" "--recompute=${recompute}"
      "--require-score=${require_score}"
      "--output-root=${output_root}"
    )
    [[ -n "$proteins" ]] && call_args+=("--proteins=${proteins}")
    if ((${#scores[@]})); then score_csv="$(IFS=,; echo "${scores[*]}")"; call_args+=("--scores=${score_csv}"); fi
    echo "yy command: plot"
    ((main == 1)) && echo "yy plot preset: main-cad-four-score"
    echo "yy plot side: ${side}"
    echo "yy plot anchor: ${anchor}"
    echo "yy plot adjustment: ${adjustment}"
    echo "yy plot score methods: ${scores[*]:-none}"
    DIR0="$DIR0" PHEDIR="$PHEDIR" SCRIPT_ROOT="$SCRIPT_ROOT" ANALYSIS_ROOT="$ANALYSIS_ROOT" \
      YY_OUTDIR="$YY_OUTDIR" exec "$rscript" --vanilla "$renderer" "${call_args[@]}"
    ;;
  *) echo "Unknown yy command: ${command_name}" >&2; global_help >&2; exit 2 ;;
esac
