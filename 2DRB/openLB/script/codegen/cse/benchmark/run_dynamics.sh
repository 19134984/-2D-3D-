#!/usr/bin/env bash
#
# This script generates CSE optimized code for a specific dynamics class.
#
# Usage: ./benchmark_dynamics.sh [--runs] [--size (million cells)] "dynamics::Tuple<...>"
#
# The script outputs the generated header content to stdout, unless --install is specified.

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." >/dev/null 2>&1 && pwd )"
OLB_ROOT="$( cd "$SCRIPT_DIR/../../.." >/dev/null 2>&1 && pwd )"

# Parse arguments
DYNAMICS_RAW=""
NCELLS=1 # million cells
NRUNS=5

while [[ $# -gt 0 ]]; do
  case $1 in
     --size)
     	NCELLS=$2
	    shift
     	;;
      --runs)
        NRUNS=$2
        shift
        ;;
    *)
      DYNAMICS_RAW="$1"
      shift
      ;;
  esac
done

if [ -z "$DYNAMICS_RAW" ]; then
    echo "Usage: $0 [--size] [--runs] \"dynamics::Tuple<...>\"" >&2
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

WORK_DIR=$SCRIPT_DIR/benchmark/tmp/$DYNAMICS_HASH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WORK_DIR

# Extraction
cat <<EOF > "$WORK_DIR/Makefile"
EXAMPLE = $DYNAMICS_HASH
OLB_ROOT := $OLB_ROOT
include \$(OLB_ROOT)/default.mk
EOF

# Run extraction.py
python3 "$SCRIPT_DIR/source/runBenchmark.py" \
    "$DYNAMICS_HASH" \
    "$WORK_DIR/" \
    "dynamics" \
    --size $NCELLS --runs $NRUNS
