#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLI="${PROJECT_ROOT}/score.sh"
TEST_ROOT="$(mktemp -d)"
trap 'rm -rf "$TEST_ROOT"' EXIT

PHE_ROOT="${TEST_ROOT}/data/ukb/phe"
SCRIPT_ROOT_TEST="${TEST_ROOT}/scripts"
mkdir -p \
  "${PHE_ROOT}/Rdata" \
  "${PHE_ROOT}/rap" \
  "${TEST_ROOT}/ppp/map.raw" \
  "${SCRIPT_ROOT_TEST}/ukb/pradeep" \
  "${SCRIPT_ROOT_TEST}/ukb/yu" \
  "${TEST_ROOT}/bin"

touch \
  "${PHE_ROOT}/Rdata/all.rds" \
  "${PHE_ROOT}/Rdata/prot.rds" \
  "${PHE_ROOT}/rap/pheno.tsv.gz" \
  "${TEST_ROOT}/ppp/prot.tab.gz" \
  "${TEST_ROOT}/ppp/map.raw/olink_protein_map_1.5k_v1.tsv" \
  "${TEST_ROOT}/ppp/map.raw/olink_protein_map_3k_v1.tsv" \
  "${TEST_ROOT}/ppp/ppp.b38.bed"

inputs="$({
  DIR0="$TEST_ROOT" \
  PHEDIR="$PHE_ROOT" \
  SCRIPT_ROOT="$SCRIPT_ROOT_TEST" \
  ANALYSIS_ROOT="${TEST_ROOT}/analysis" \
  bash "$CLI" --inputs
})"

grep -q "RAW_PROTEIN=${TEST_ROOT}/ppp/prot.tab.gz" <<<"$inputs"
grep -q "RAW_PHENOTYPE_TABLE=${PHE_ROOT}/rap/pheno.tsv.gz" <<<"$inputs"
grep -q "RAW_PROTEIN_MAP=${TEST_ROOT}/ppp/map.raw/olink_protein_map_1.5k_v1.tsv" <<<"$inputs"
grep -q "RAW_PROTEIN_MAP=${TEST_ROOT}/ppp/map.raw/olink_protein_map_3k_v1.tsv" <<<"$inputs"
grep -q "RAW_PROTEIN_BED=${TEST_ROOT}/ppp/ppp.b38.bed" <<<"$inputs"
grep -q 'NATIVE_REBUILD_OLINK_DATES_STATUS=AUTO_DOWNLOAD_ON_PREFLIGHT_OR_COMPUTE' <<<"$inputs"

printf 'PlateID\tPanel\tProcessing_StartDate\nTEST\tCARDIOMETABOLIC\t2020-01-01\n' \
  >"${TEST_ROOT}/resource1019.dat"
if command -v sha256sum >/dev/null 2>&1; then
  resource_sha="$(sha256sum "${TEST_ROOT}/resource1019.dat" | awk '{print $1}')"
else
  resource_sha="$(shasum -a 256 "${TEST_ROOT}/resource1019.dat" | awk '{print $1}')"
fi

printf '%s\n' \
  '#!/usr/bin/env bash' \
  'if [[ "${1:-}" == "-c" ]]; then echo 3.9; fi' \
  'exit 0' >"${TEST_ROOT}/bin/fake-python"
printf '%s\n' '#!/usr/bin/env bash' 'exit 0' >"${SCRIPT_ROOT_TEST}/ukb/yu/yu.sh"
printf '%s\n' '#!/usr/bin/env bash' 'sleep 2' >"${SCRIPT_ROOT_TEST}/ukb/pradeep/pradeep.sh"
chmod +x "${TEST_ROOT}/bin/fake-python" "${SCRIPT_ROOT_TEST}/ukb/yu/yu.sh"
chmod +x "${SCRIPT_ROOT_TEST}/ukb/pradeep/pradeep.sh"

target="${TEST_ROOT}/files/yu-protein-analysis/references/raw/olink_processing_start_date.dat"
preflight="$({
  DIR0="$TEST_ROOT" \
  PHEDIR="$PHE_ROOT" \
  SCRIPT_ROOT="$SCRIPT_ROOT_TEST" \
  ANALYSIS_ROOT="${TEST_ROOT}/analysis" \
  RSCRIPT=/usr/bin/true \
  YU_PYTHON="${TEST_ROOT}/bin/fake-python" \
  YU_OLINK_DATES_URL="file://${TEST_ROOT}/resource1019.dat" \
  YU_OLINK_DATES_SHA256="$resource_sha" \
  bash "$CLI" yu-strict --preflight
})"

grep -q 'UKB Resource 1019: DOWNLOADED_AND_VERIFIED' <<<"$preflight"
[[ -f "$target" ]]
cmp -s "$target" "${TEST_ROOT}/resource1019.dat"

if command -v flock >/dev/null 2>&1; then
  DIR0="$TEST_ROOT" \
  PHEDIR="$PHE_ROOT" \
  SCRIPT_ROOT="$SCRIPT_ROOT_TEST" \
  ANALYSIS_ROOT="${TEST_ROOT}/analysis" \
  RSCRIPT=/usr/bin/true \
  bash "$CLI" pradeep-strict --compute --confirm-compute \
    >"${TEST_ROOT}/first-compute.log" 2>&1 &
  first_pid=$!
  for _ in 1 2 3 4 5 6 7 8 9 10; do
    [[ -s "${TEST_ROOT}/analysis/yy/score/.locks/pradeep-strict-cad.lock" ]] && break
    sleep 0.1
  done
  if DIR0="$TEST_ROOT" \
    PHEDIR="$PHE_ROOT" \
    SCRIPT_ROOT="$SCRIPT_ROOT_TEST" \
    ANALYSIS_ROOT="${TEST_ROOT}/analysis" \
    RSCRIPT=/usr/bin/true \
    bash "$CLI" pradeep-strict --compute --confirm-compute \
      >"${TEST_ROOT}/second-compute.log" 2>&1; then
    echo 'Duplicate compute lock did not reject a concurrent run.' >&2
    exit 1
  fi
  grep -q 'computation is already running' "${TEST_ROOT}/second-compute.log"
  wait "$first_pid"
fi

echo 'test_input_resolution.sh: PASS'
