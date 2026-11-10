#!/usr/bin/env bash
#
# This script generates CSE optimized code for a specific dynamics class.
#
# Usage: ./compile_dynamics.sh [--precision float/double] "dynamics::Tuple<...>"
#
# The script outputs the generated header content to stdout, unless --install is specified.

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." >/dev/null 2>&1 && pwd )"
OLB_ROOT="$( cd "$SCRIPT_DIR/../../.." >/dev/null 2>&1 && pwd )"

# Parse arguments
DYNAMICS_RAW=""
PRECISION=float

echo "RUNNING"

while [[ $# -gt 0 ]]; do
  case $1 in
      --precision)
	PRECISION=$2
	shift
	;;
    *)
      DYNAMICS_RAW="$1"
      ;;
  esac
  shift
done

if [ -z "$DYNAMICS_RAW" ]; then
    echo "Usage: $0 [--precision] \"dynamics::Tuple<...>\"" >&2
    exit 1
fi

echo $PRECISION
echo $DYNAMICS_RAW

# Preprocessing
DYNAMICS_FORMATTED=$(python3 "$SCRIPT_DIR/source/format_common.py" "$DYNAMICS_RAW" "dummy")

if [ -z "$DYNAMICS_FORMATTED" ]; then
    echo "Error: Failed to format dynamics string." >&2
    exit 1
fi

# Compute hash for unique identification
DYNAMICS_HASH=$(echo -n "$DYNAMICS_FORMATTED" | sha256sum | awk '{print $1}')
echo $DYNAMICS_HASH
# Setup temporary workspace
mkdir -p $SCRIPT_DIR/benchmark/tmp
mkdir -p $SCRIPT_DIR/benchmark/tmp/$DYNAMICS_HASH
WORK_DIR=$SCRIPT_DIR/benchmark/tmp/$DYNAMICS_HASH

# Extraction
cat <<EOF > "$WORK_DIR/Makefile"
EXAMPLE = $DYNAMICS_HASH
OLB_ROOT := $OLB_ROOT
include \$(OLB_ROOT)/default.mk
EOF

# Run compile.py
python3 "$SCRIPT_DIR/source/compile.py" \
    "$DYNAMICS_FORMATTED" \
    "$DYNAMICS_HASH" \
    "$SCRIPT_DIR/source/templates/dynamics/benchmark.cpp.template" \
    "$WORK_DIR/" \
    "dynamics" \
    --precision $PRECISION

