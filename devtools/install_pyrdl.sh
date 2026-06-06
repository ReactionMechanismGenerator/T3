#!/usr/bin/env bash
set -eo pipefail

# Thin wrapper: installs PyRDL via ARC's install_pyrdl.sh.
# Default target environment is t3_env (overridable via $1).

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
T3_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
ARC_ROOT="${ARC_ROOT:-$T3_ROOT/../ARC}"

if [[ ! -f "$ARC_ROOT/devtools/install_pyrdl.sh" ]]; then
    echo "Error: ARC not found at $ARC_ROOT"
    echo "Clone ARC alongside T3, or set ARC_ROOT to point to it."
    exit 1
fi

exec bash "$ARC_ROOT/devtools/install_pyrdl.sh" "${1:-t3_env}"
