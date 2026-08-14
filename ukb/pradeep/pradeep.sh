#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PRADEEP_PROJECT_DIR="$PROJECT_DIR"

find_rscript() {
  if [[ -n "${RSCRIPT:-}" ]]; then
    if [[ -x "$RSCRIPT" ]] || command -v "$RSCRIPT" >/dev/null 2>&1; then
      printf '%s\n' "$RSCRIPT"
      return 0
    fi
    printf 'RSCRIPT was set but is not executable: %s\n' "$RSCRIPT" >&2
    return 1
  fi

  if command -v Rscript >/dev/null 2>&1; then
    command -v Rscript
    return 0
  fi

  local candidate
  for candidate in \
    "/c/Program Files/R/R-4.3.2/bin/x64/Rscript.exe" \
    "/c/Program Files/R/R-4.5.1/bin/x64/Rscript.exe"; do
    if [[ -x "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done

  printf '%s\n' "Rscript was not found." >&2
  printf '%s\n' "Install R or set RSCRIPT to the executable, for example:" >&2
  printf '%s\n' "  RSCRIPT=/path/to/Rscript bash pradeep.sh --h" >&2
  return 1
}

RSCRIPT_BIN="$(find_rscript)"
exec "$RSCRIPT_BIN" --vanilla "$PROJECT_DIR/f/cli.R" "$@"
