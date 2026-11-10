#!/usr/bin/env bash

rm -rf tmp/gnuplotData/data/* 2>/dev/null || true

EXEC_CMD=${EXEC_CMD:-mpirun}
# Automatically grab the executable name based on the current directory
EXE_NAME=$(basename "$PWD")

for N in 30 40 50 60 80 100 120; do
    $EXEC_CMD ./$EXE_NAME -- --RESOLUTION $N --VTK_ENABLED 0
done

python3 - << 'EOF'
import csv
import glob
import math
import os

data_dir = "tmp/gnuplotData/data"
output_file = os.path.join(data_dir, "eoc.csv")
metrics = ["disp", "stress"]

# Storage for results: {metric: (avg_l2, avg_linf)}
avg_results = {}

for metric in metrics:
    # Match files like disp_err.csv or disp.csv
    pattern = os.path.join(data_dir, f"*{metric}*.csv")
    matching_files = glob.glob(pattern)

    if not matching_files:
        print(f"Warning: No CSV file found for {metric} in {data_dir}")
        avg_results[metric] = ("", "")
        continue

    file_path = matching_files[0]
    
    rows = []
    with open(file_path, mode='r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                'N': float(row['N']),
                'L2': float(row['L2']),
                'LInf': float(row['LInf'])
            })

    # Sort rows by resolution N
    rows.sort(key=lambda r: r['N'])

    eoc_l2_list = []
    eoc_linf_list = []

    # Calculate step-by-step EOC: log(err1/err2) / log(N2/N1)
    for i in range(1, len(rows)):
        n1, n2 = rows[i-1]['N'], rows[i]['N']
        l2_1, l2_2 = rows[i-1]['L2'], rows[i]['L2']
        linf_1, linf_2 = rows[i-1]['LInf'], rows[i]['LInf']

        log_n_ratio = math.log(n2 / n1)
        
        eoc_l2 = math.log(l2_1 / l2_2) / log_n_ratio
        eoc_linf = math.log(linf_1 / linf_2) / log_n_ratio

        eoc_l2_list.append(eoc_l2)
        eoc_linf_list.append(eoc_linf)

    avg_l2 = sum(eoc_l2_list) / len(eoc_l2_list) if eoc_l2_list else ""
    avg_linf = sum(eoc_linf_list) / len(eoc_linf_list) if eoc_linf_list else ""

    avg_results[metric] = (avg_l2, avg_linf)

# Write results to eoc.csv
with open(output_file, mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Metric"] + metrics)
    writer.writerow(["Avg_EOC_L2"] + [avg_results[m][0] for m in metrics])
    writer.writerow(["Avg_EOC_LInf"] + [avg_results[m][1] for m in metrics])

print(f"Generated {output_file} successfully.")
EOF
