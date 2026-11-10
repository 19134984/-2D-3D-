#!/usr/bin/env bash
#
# This script generates CSE optimized code for a specific dynamics class.
#
# Usage: ./optimize_dynamics.sh [--deterministic] [--install] [--skip-header-update] [--keep-tmp] "dynamics::Tuple<...>"
#
# The script outputs the generated header content to stdout, unless --install is specified.

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
OLB_ROOT="$( cd "$SCRIPT_DIR/../../.." >/dev/null 2>&1 && pwd )"

# Parse arguments
DETERMINISTIC_FLAG=""
INSTALL_FLAG=0
SKIP_HEADER_UPDATE=0
KEEP_TMP_FLAG=0
DYNAMICS_RAW=""
TARGET="Generic"

while [[ $# -gt 0 ]]; do
  case $1 in
    --deterministic)
      DETERMINISTIC_FLAG="--deterministic"
      shift
      ;;
    --install)
      INSTALL_FLAG=1
      shift
      ;;
    --skip-header-update)
      SKIP_HEADER_UPDATE=1
      shift
      ;;
    --keep-tmp)
      KEEP_TMP_FLAG=1
      shift
      ;;
    --target)
      TARGET="$2"
      shift 2
      ;;
    *)
      DYNAMICS_RAW="$1"
      shift
      ;;
  esac
done

# CodegenTarget controls which collision body is generated (Generic, StoredMoments).
# StoredMoments variants are installed with a .sm suffix and the hash is unchanged.
case "$TARGET" in
  Generic)       TARGET_SUFFIX="" ;;
  StoredMoments) TARGET_SUFFIX=".sm" ;;
  *) echo "Error: unknown --target '$TARGET' (expected Generic or StoredMoments)." >&2; exit 1 ;;
esac

if [ -z "$DYNAMICS_RAW" ]; then
    echo "Usage: $0 [--deterministic] [--install] [--skip-header-update] \"dynamics::Tuple<...>\"" >&2
    exit 1
fi

# Preprocessing
DYNAMICS_FORMATTED=$(python3 "$SCRIPT_DIR/source/format_common.py" "$DYNAMICS_RAW")

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
    if [ "$KEEP_TMP_FLAG" -eq 0 ]; then
        rm -rf "$WORK_DIR"
    else
        echo "Keepting tmp directory: $WORK_DIR"
    fi
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
    "$SCRIPT_DIR/source/templates/dynamics/extract.cpp.template" \
    "$WORK_DIR/$DYNAMICS_HASH.txt" \
    "$WORK_DIR/" \
    "$WORK_DIR/" \
    "dynamics" \
    --codegen_target "$TARGET" >&2

if [ ! -f "$WORK_DIR/$DYNAMICS_HASH.out" ]; then
    echo "Error: Extraction failed, .out file not found." >&2
    exit 1
fi

# Optimization
# Run optimize.py
python3 "$SCRIPT_DIR/source/optimize.py" \
    "$WORK_DIR/$DYNAMICS_HASH.out" \
    "$SCRIPT_DIR/source/templates/dynamics/dynamics.cse.template" \
    "$WORK_DIR/$DYNAMICS_HASH.cse.h" \
    $DETERMINISTIC_FLAG \
    --codegen_target "$TARGET" >&2

if [ ! -f "$WORK_DIR/$DYNAMICS_HASH.cse.h" ]; then
    echo "Error: Optimization failed, .cse.h file not found." >&2
    exit 1
fi

# Append Generation Info
if git rev-parse --git-dir > /dev/null 2>&1; then
    COMMIT_HASH=$(git rev-parse HEAD)
    echo "" >> "$WORK_DIR/$DYNAMICS_HASH.cse.h"
    echo "// Generation Info: commit=$COMMIT_HASH" >> "$WORK_DIR/$DYNAMICS_HASH.cse.h"
fi

# Output or Install
if [ "$INSTALL_FLAG" -eq 1 ]; then
    TARGET_DIR="$OLB_ROOT/src/cse/dynamics"
    INSTALL_NAME="${DYNAMICS_HASH}${TARGET_SUFFIX}"
    TARGET_FILE="$TARGET_DIR/$INSTALL_NAME.cse.h"
    HEADER_FILE="$TARGET_DIR/generated_cse.h"

    # Ensure target directory exists
    mkdir -p "$TARGET_DIR"

    # Install file
    cp "$WORK_DIR/$DYNAMICS_HASH.cse.h" "$TARGET_FILE"
    echo "Installed generated code to $TARGET_FILE" >&2

    if [ "$SKIP_HEADER_UPDATE" -eq 0 ]; then
        python3 "$SCRIPT_DIR/source/update_cse_includes.py" dynamics "$OLB_ROOT" >&2
    fi
else
    cat "$WORK_DIR/$DYNAMICS_HASH.cse.h"
fi
