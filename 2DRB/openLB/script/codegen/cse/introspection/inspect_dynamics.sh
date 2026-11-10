#!/usr/bin/env bash
#
# This script generates CSE optimized code for a specific dynamics class.
#
# Usage: ./inspect_dynamics.sh [--disable-automatic-cse] "dynamics::Tuple<...>"
#
# The script outputs the generated header content to stdout, unless --install is specified.

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." >/dev/null 2>&1 && pwd )"
OLB_ROOT="$( cd "$SCRIPT_DIR/../../.." >/dev/null 2>&1 && pwd )"

# Parse arguments
PRECISION=float
DISABLE_CSE=""
DYNAMICS_RAW=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --disable-automatic-cse)
    	DISABLE_CSE="--disable_cse"
        shift
        ;;
    --precision)
	PRECISION=$2
	shift
	;;
    *)
      DYNAMICS_RAW="$1"
      shift
      ;;
  esac
done

if [ -z "$DYNAMICS_RAW" ]; then
    echo "Usage: $0 [--disable-automatic-cse] \"dynamics::Tuple<...>\"" >&2
    exit 1
fi

# Preprocessing
DYNAMICS_FORMATTED=$(python3 "$SCRIPT_DIR/source/format_common.py" "$DYNAMICS_RAW" "dummy")

if [ -z "$DYNAMICS_FORMATTED" ]; then
    echo "Error: Failed to format dynamics string." >&2
    exit 1
fi

# Compute hash for unique identification
DYNAMICS_HASH=$(echo -n "$DYNAMICS_FORMATTED" | sha256sum | awk '{print $1}')

# Setup temporary workspace
WORK_DIR=$(mktemp -d -t cse_optimize_dynamics.XXXXXX)

# Cleanup function
cleanup() {
    rm -rf "$WORK_DIR"
}
trap cleanup EXIT

# Extraction
cat <<EOF > "$WORK_DIR/Makefile"
EXAMPLE = $DYNAMICS_HASH
OLB_ROOT := $OLB_ROOT
include \$(OLB_ROOT)/default.mk
EOF

# Create a dummy input file to satisfy extraction.py
touch "$WORK_DIR/$DYNAMICS_HASH.txt"

# Run extraction.py
python3 "$SCRIPT_DIR/source/extraction.py" \
    "$DYNAMICS_FORMATTED" \
    "$SCRIPT_DIR/source/templates/dynamics/inspect.cpp.template" \
    "$WORK_DIR/$DYNAMICS_HASH.txt" \
    "$WORK_DIR/" \
    "$WORK_DIR/" \
    "dynamics" \
    --precision $PRECISION $DISABLE_CSE

