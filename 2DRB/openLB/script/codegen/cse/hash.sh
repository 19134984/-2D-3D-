DYNAMICS_RAW=$1

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
OLB_ROOT="$( cd "$SCRIPT_DIR/../../.." >/dev/null 2>&1 && pwd )"


DYNAMICS_FORMATTED=$(python3 "$SCRIPT_DIR/source/format_common.py" "$DYNAMICS_RAW")

if [ -z "$DYNAMICS_FORMATTED" ]; then
    echo "Error: Failed to format dynamics string." >&2
    exit 1
fi

# Compute hash for unique identification
DYNAMICS_HASH=$(echo -n "$DYNAMICS_FORMATTED" | sha256sum | awk '{print $1}')

echo $DYNAMICS_HASH

