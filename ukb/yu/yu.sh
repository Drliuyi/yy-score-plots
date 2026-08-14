#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
runner="${root}/f/tools/yu_runner.py"

for python in python3 python; do
  if command -v "${python}" >/dev/null 2>&1; then
    exec "${python}" "${runner}" "$@"
  fi
done

printf 'ERROR: Python 3 is required to start Yu Protein Analysis.\n' >&2
exit 127
