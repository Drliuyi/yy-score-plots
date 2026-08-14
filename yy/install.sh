#!/usr/bin/env bash
set -euo pipefail

CHECKOUT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_ROOT="${SCRIPT_ROOT:-${CHECKOUT_ROOT}}"
entry="${SCRIPT_ROOT}/yy/cli.sh"
target="${HOME}/.local/bin/yy"

[[ -f "$entry" ]] || { echo "Missing yy entrypoint: ${entry}" >&2; exit 2; }
mkdir -p "$(dirname "$target")"
ln -sfn "$entry" "$target"
echo "Installed: ${target} -> ${entry}"
case ":${PATH}:" in
  *":$(dirname "$target"):"*) ;;
  *) echo "Add $(dirname "$target") to PATH before running yy." ;;
esac
