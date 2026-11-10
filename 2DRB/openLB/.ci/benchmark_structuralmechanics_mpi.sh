#!/usr/bin/env nix
#!nix --extra-experimental-features nix-command --extra-experimental-features flakes develop .#gcc-mpi --command bash
export FLOATING_POINT_TYPE=double

set -euo pipefail

export OMPI_MCA_pml="ob1"

echo "Running in MPI mode"

# Compile the root project
make clean || true
make -j4

ROOT=$PWD

# Helper function to compare eoc.csv files with a 5% margin
compare_eoc_csv() {
    python3 -c '
import sys, csv, math

file_out, file_ref = sys.argv[1], sys.argv[2]

def read_eoc_csv(filepath):
    data = {}
    with open(filepath, mode="r", newline="") as f:
        reader = csv.reader(f)
        header = [col.strip() for col in next(reader)]
        for row in reader:
            if not row: continue
            metric_label = row[0].strip()
            row_vals = {}
            for col_name, val_str in zip(header[1:], row[1:]):
                try:
                    row_vals[col_name.strip()] = float(val_str)
                except ValueError:
                    pass
            data[metric_label] = row_vals
    return data

data_out = read_eoc_csv(file_out)
data_ref = read_eoc_csv(file_ref)

# Compare keys/rows between generated and reference
for row_key, ref_metrics in data_ref.items():
    if row_key not in data_out:
        print(f"Missing row {row_key} in output EOC file")
        sys.exit(1)
    
    out_metrics = data_out[row_key]
    for metric_name, ref_val in ref_metrics.items():
        if metric_name not in out_metrics:
            print(f"Missing metric {metric_name} under {row_key} in output")
            sys.exit(1)
        
        out_val = out_metrics[metric_name]
        
        # 5% relative tolerance comparison
        if not math.isclose(out_val, ref_val, rel_tol=0.05, abs_tol=1e-4):
            diff_pct = abs(out_val - ref_val) / abs(ref_val) * 100 if ref_val != 0 else 0
            print(f"EOC mismatch in {row_key} [{metric_name}]: {out_val:.4f} (out) vs {ref_val:.4f} (ref) | Diff: {diff_pct:.2f}%")
            sys.exit(1)

sys.exit(0)
' "$1" "$2"
}

EXAMPLES_BASE="examples/structuralmechanics"
REF_BASE=".ci/reference/structuralmechanics"

# Load makefiles into an array to avoid stdin conflicts
readarray -t makefiles < <(find "$EXAMPLES_BASE" -type f -name Makefile | sort)

for makefile_path in "${makefiles[@]}"; do
    [ -z "$makefile_path" ] && continue
    
    dir_path=$(dirname "$makefile_path")
    case_name=$(basename "$dir_path")
    
    if [ ! -f "$ROOT/$dir_path/run.sh" ]; then
        continue
    fi

    echo "==========================================================="
    echo " Building and Testing (MPI): $case_name"
    echo " Directory: $dir_path"
    echo "==========================================================="
    
    cd "$ROOT/$dir_path"
    
    make clean || true
    make -j4
    
    export EXEC_CMD="mpirun"
    bash run.sh </dev/null
    
    rel_path=${dir_path#"$EXAMPLES_BASE/"}
    ref_dir="$ROOT/$REF_BASE/$rel_path"
    out_dir="$ROOT/$dir_path/tmp/gnuplotData/data"
    out_file="$out_dir/eoc.csv"
    ref_file="$ref_dir/eoc.csv"

    if [ ! -f "$out_file" ]; then
        echo "ERROR: $out_file was not generated in $dir_path"
        exit 1
    fi

    if [ ! -f "$ref_file" ]; then
        echo "ERROR: Reference file $ref_file does not exist."
        exit 1
    fi

    # Check exact match first for speed
    if diff -q "$out_file" "$ref_file" > /dev/null; then
        echo "  [OK] eoc.csv matches reference perfectly."
    else
        # Fallback to EOC tolerance comparison (±5%)
        if error_msg=$(compare_eoc_csv "$out_file" "$ref_file"); then
            echo "  [OK] eoc.csv matches reference within 5% tolerance."
        else
            echo "  [FAIL] Case $case_name failed."
            echo "  [FAIL] EOC mismatch in eoc.csv"
            echo "         $error_msg"
            exit 1
        fi
    fi

    echo "  -> $case_name passed!"
    echo ""
    
    cd "$ROOT"

done

echo "All structural mechanics tests (MPI) completed successfully!"
