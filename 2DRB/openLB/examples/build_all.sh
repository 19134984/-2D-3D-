#!/usr/bin/env bash

set -uo pipefail

MODE="build"
ROOT="."

if [[ $# -gt 0 ]]; then
    case "$1" in
        clean)
            MODE="clean"
            ;;
        build)
            MODE="build"
            ;;
        *)
            echo "Usage: $0 [build|clean]"
            exit 1
            ;;
    esac
fi

LOGFILE="make_${MODE}_$(date +%Y%m%d_%H%M%S).log"

passed=()
failed=()

echo "===================================" | tee "$LOGFILE"
echo "Mode : $MODE" | tee -a "$LOGFILE"
echo "Root : $(realpath "$ROOT")" | tee -a "$LOGFILE"
echo "Time : $(date)" | tee -a "$LOGFILE"
echo "===================================" | tee -a "$LOGFILE"

# Find all directories containing a Makefile
mapfile -t dirs < <(
    find "$ROOT" \( -name Makefile -o -name makefile \) -printf '%h\n' | sort -u
)

echo "Found ${#dirs[@]} make directories."
echo

for dir in "${dirs[@]}"; do
    echo "-----------------------------------" | tee -a "$LOGFILE"
    echo "$MODE: $dir" | tee -a "$LOGFILE"

    if [[ "$MODE" == "build" ]]; then
        (
            cd "$dir" && make -j 6
        ) >>"$LOGFILE" 2>&1
    else
        (
            cd "$dir" && make clean
        ) >>"$LOGFILE" 2>&1
    fi

    if [[ $? -eq 0 ]]; then
        passed+=("$dir")
        echo "PASS" | tee -a "$LOGFILE"
    else
        failed+=("$dir")
        echo "FAIL" | tee -a "$LOGFILE"
    fi
done

echo
echo "==================================="
echo "Summary"
echo "==================================="

echo "Passed: ${#passed[@]}"
for d in "${passed[@]}"; do
    echo "  ✓ $d"
done

echo
echo "Failed: ${#failed[@]}"
for d in "${failed[@]}"; do
    echo "  ✗ $d"
done

echo
echo "Detailed log written to:"
echo "  $LOGFILE"

if [[ ${#failed[@]} -gt 0 ]]; then
    exit 1
fi

exit 0