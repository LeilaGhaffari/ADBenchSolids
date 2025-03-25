#!/usr/bin/env python3

import re
import os
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

if len(sys.argv) != 2:
    print("Usage: ./parse_benchmark.py <benchmark_output_file>")
    sys.exit(1)

filename = sys.argv[1]

# Nested dictionaries:
residual_times = defaultdict(lambda: defaultdict(list))
jacobian_times = defaultdict(lambda: defaultdict(list))

with open(filename, "r") as f:
    lines = f.readlines()

    current_qpts = None
    for line in lines:
        qpts_match = re.search(r"Quadrature Points\s*=\s*(\d+)", line)
        if qpts_match:
            current_qpts = int(qpts_match.group(1))
            continue

        if (
            line.strip() == "" or
            line.startswith("Iteration") or
            line.startswith("AD Tool") or
            line.startswith("-")
        ):
            continue

        tokens = line.strip().split()
        if current_qpts is not None and len(tokens) >= 2:
            model = tokens[0]

            try:
                res_time = float(tokens[1])
                residual_times[model][current_qpts].append(res_time)
            except ValueError:
                continue

            if len(tokens) >= 4:
                try:
                    jac_time = float(tokens[3])
                    jacobian_times[model][current_qpts].append(jac_time)
                except ValueError:
                    pass

            current_qpts = None

# Begin plotting
plt.figure(figsize=(10, 6))

for model in residual_times:
    qpts_data = residual_times[model]
    qpts_sorted = sorted(qpts_data.keys())
    avg_res_times = [np.mean(qpts_data[q]) for q in qpts_sorted]
    throughput_res = [q / t for q, t in zip(qpts_sorted, avg_res_times)]
    plt.plot(qpts_sorted, throughput_res, marker='o', label=f"{model} (residual)")

for model in jacobian_times:
    qpts_data = jacobian_times[model]
    qpts_sorted = sorted(qpts_data.keys())
    avg_jac_times = [np.mean(qpts_data[q]) for q in qpts_sorted]
    throughput_jac = [q / t for q, t in zip(qpts_sorted, avg_jac_times)]
    plt.plot(qpts_sorted, throughput_jac, marker='x', linestyle='--', label=f"{model} (jacobian)")

# Plot total throughput (residual + jacobian) where both exist
for model in residual_times:
    if model not in jacobian_times:
        continue

    common_qpts = sorted(set(residual_times[model]) & set(jacobian_times[model]))
    total_times = [
        np.mean(residual_times[model][q]) + np.mean(jacobian_times[model][q])
        for q in common_qpts
    ]
    throughput_total = [q / t for q, t in zip(common_qpts, total_times)]
    plt.plot(common_qpts, throughput_total, marker='s', linestyle='-.', label=f"{model} (total)")

# Finalize plot
plt.xlabel("Quadrature Points")
plt.ylabel("Throughput (qpts / time)")
plt.title("Throughput vs Quadrature Points (Residual, Jacobian, Total)")
plt.xscale("log")
plt.grid(True)
plt.legend(title="Model")
plt.tight_layout()

# Save plot
basename = os.path.splitext(os.path.basename(filename))[0]
plot_filename = f"{basename}.png"
plt.savefig(plot_filename, dpi=300)
print(f"Plot saved as '{plot_filename}'")
