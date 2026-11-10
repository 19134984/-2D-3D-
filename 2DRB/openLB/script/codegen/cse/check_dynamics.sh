#!/usr/bin/env bash
# Check code generation for the test dynamics against the committed references, for
# each CodegenTarget: Generic (<hash>.cse.h) and StoredMoments (<hash>.sm.cse.h,
# skipped if absent). Run under the pinned codegen toolchain (nix develop .#codegen).

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
OLB_ROOT="$( cd "$SCRIPT_DIR/../../.." >/dev/null 2>&1 && pwd )"
TEST_FILE="$SCRIPT_DIR/tests/test_dynamics.txt"
REF_DIR="$OLB_ROOT/src/cse/dynamics"

if [ ! -f "$TEST_FILE" ]; then
    echo "Error: Test definition file not found: $TEST_FILE"
    exit 1
fi

echo "Running Code Generation Checks..."

FAILED=0
PASS_COUNT=0
SKIP_COUNT=0

# check_target <dynamics> <hash> <target> <suffix>
# Returns 0 on pass, 1 on fail, 2 on skip (StoredMoments reference absent).
check_target() {
    local DYNAMICS="$1" HASH="$2" TARGET="$3" SUFFIX="$4"
    local REF_FILE="$REF_DIR/$HASH$SUFFIX.cse.h"

    if [ ! -f "$REF_FILE" ]; then
        if [ "$TARGET" = "StoredMoments" ]; then
            echo "  [$TARGET] SKIP (no $HASH$SUFFIX.cse.h reference)"
            return 2
        fi
        echo "  [$TARGET] FAIL: reference $REF_FILE missing"
        return 1
    fi

    local GEN_FILE
    GEN_FILE=$(mktemp)
    if ! "$SCRIPT_DIR/optimize_dynamics.sh" --deterministic --target "$TARGET" "$DYNAMICS" < /dev/null > "$GEN_FILE" 2>/dev/null; then
        echo "  [$TARGET] FAIL: generation failed"
        rm -f "$GEN_FILE"
        return 1
    fi

    # Compare, ignoring the generation info line which contains the commit hash.
    if diff -q -I "^// Generation Info:" "$REF_FILE" "$GEN_FILE" >/dev/null; then
        echo "  [$TARGET] PASS"
        rm -f "$GEN_FILE"
        return 0
    fi

    echo "  [$TARGET] FAIL: output differs from reference"
    diff -I "^// Generation Info:" "$REF_FILE" "$GEN_FILE" | head -n 10
    rm -f "$GEN_FILE"
    return 1
}

mapfile -t DYNAMICS_LIST < "$TEST_FILE"

INDEX=0
for DYNAMICS in "${DYNAMICS_LIST[@]}"; do
    if [ -z "$DYNAMICS" ]; then
        continue
    fi

    ((++INDEX))

    DYNAMICS_FORMATTED=$(python3 "$SCRIPT_DIR/source/format_common.py" "$DYNAMICS")
    HASH=$(echo -n "$DYNAMICS_FORMATTED" | sha256sum | awk '{print $1}')

    echo "[$INDEX] $HASH"

    for SPEC in "Generic:" "StoredMoments:.sm"; do
        TARGET="${SPEC%%:*}"
        SUFFIX="${SPEC#*:}"
        set +e
        check_target "$DYNAMICS" "$HASH" "$TARGET" "$SUFFIX"
        RC=$?
        set -e
        case $RC in
            0) ((++PASS_COUNT)) ;;
            2) ((++SKIP_COUNT)) ;;
            *) FAILED=1 ;;
        esac
    done
done

echo
echo "Passed: $PASS_COUNT, Skipped: $SKIP_COUNT"
if [ $FAILED -eq 0 ]; then
    echo "All code generation checks PASSED."
    exit 0
else
    echo "Some code generation checks FAILED."
    exit 1
fi
