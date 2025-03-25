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

# Store times: model -> qpts -> list of times
residual_times = defaultdict(lambda: defaultdict(list))
jacobian_times = defaultdict(lambda: defaultdict(list))

with open(filename, "r") as f:
    lines = f.readlines()

    current_qpts = None
    for line in lines:
        # Match quadrature points
        qpts_match = re.search(r"Quadrature Points\s*=\s*(\d+)", line)
        if qpts_match:
            current_qpts = int(qpts_match.group(1))
            continue

        # Skip non-data lines
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

            # Residual time
            try:
                res_time = float(tokens[1])
                residual_times[model][current_qpts].append(res_time)
            except ValueError:
                pass

            # Jacobian time if available
            if len(tokens) >= 4:
                try:
                    jac_time = float(tokens[3])
                    jacobian_times[model][current_qpts].append(jac_time)
                except ValueError:
                    pass

# Plotting section
plt.figure(figsize=(10, 6))

# Union of all models seen
all_models = set(residual_times) | set(jacobian_times)

for model in sorted(all_models):
    # Residual throughput
    if model in residual_times:
        qpts_sorted = sorted(residual_times[model])
        avg_res_times = [np.mean(residual_times[model][q]) for q in qpts_sorted]
        throughput_res = [q / t for q, t in zip(qpts_sorted, avg_res_times)]
        if model == "stream":
            label = f"{model} (triad)"
        else:
            label = f"{model} (residual)"
        plt.plot(qpts_sorted, throughput_res, marker='o', label=label)


    # Jacobian throughput
    if model in jacobian_times:
        qpts_sorted = sorted(jacobian_times[model])
        avg_jac_times = [np.mean(jacobian_times[model][q]) for q in qpts_sorted]
        throughput_jac = [q / t for q, t in zip(qpts_sorted, avg_jac_times)]
        plt.plot(qpts_sorted, throughput_jac, marker='x', linestyle='--', label=f"{model} (jacobian)")

    # Total throughput where both residual and jacobian are available
    if model in residual_times and model in jacobian_times:
        common_qpts = sorted(set(residual_times[model]) & set(jacobian_times[model]))
        if common_qpts:
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
