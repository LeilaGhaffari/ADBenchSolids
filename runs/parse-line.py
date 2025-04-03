#!/usr/bin/env python3

import re
import os
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

# ------------------------------------------
# Parse command-line arguments
# ------------------------------------------

filename = None
mode = "bandwidth"  # default

for arg in sys.argv[1:]:
    if arg.startswith("--input="):
        filename = arg.split("=", 1)[1]
    elif arg.startswith("--mode="):
        mode_value = arg.split("=", 1)[1]
        if mode_value in ("bandwidth", "throughput"):
            mode = mode_value
        else:
            print("Invalid mode. Use --mode=bandwidth or --mode=throughput.")
            sys.exit(1)
    else:
        print(f"Unknown argument: {arg}")
        print("Usage: ./parse_benchmark.py --input=<file> [--mode=bandwidth|throughput]")
        sys.exit(1)

if not filename:
    print("Missing required argument: --input=<file>")
    sys.exit(1)

# ------------------------------------------
# Parse data
# ------------------------------------------

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
                pass

            if len(tokens) >= 4:
                try:
                    jac_time = float(tokens[3])
                    jacobian_times[model][current_qpts].append(jac_time)
                except ValueError:
                    pass

# ------------------------------------------
# Constants for bandwidth mode
# ------------------------------------------

# Residual evaluation:
#   Load two 3×3 matrices = 2×9×8 = 144 B
#   Write one 3×3 matrix = 72 B
#   Write 16-element vector = 128 B
#   → Total = 344 B
BYTES_PER_QPT_RESIDUAL = 344

# Jacobian evaluation:
#   Load 3×3 matrix = 72 B
#   Load 16-vector = 128 B
#   Write 3×3 matrix = 72 B
#   → Total = 272 B
BYTES_PER_QPT_JACOBIAN = 272

# ------------------------------------------
# Plotting
# ------------------------------------------

plt.figure(figsize=(10, 6))
all_models = set(residual_times) | set(jacobian_times)

for model in sorted(all_models):
    # Residual
    if model in residual_times:
        qpts_sorted = sorted(residual_times[model])
        avg_res_times = [np.mean(residual_times[model][q]) for q in qpts_sorted]

        if mode == "bandwidth":
            y_res = [q * BYTES_PER_QPT_RESIDUAL / t / 1e9 for q, t in zip(qpts_sorted, avg_res_times)]
        else:
            y_res = [q / t for q, t in zip(qpts_sorted, avg_res_times)]

        if model == "stream":
            label = f"{model} (triad)"
        else:
            label = f"{model} (residual)"

        plt.plot(qpts_sorted, y_res, marker='o', label=label)

    # Jacobian
    if model in jacobian_times:
        qpts_sorted = sorted(jacobian_times[model])
        avg_jac_times = [np.mean(jacobian_times[model][q]) for q in qpts_sorted]

        if mode == "bandwidth":
            y_jac = [q * BYTES_PER_QPT_JACOBIAN / t / 1e9 for q, t in zip(qpts_sorted, avg_jac_times)]
        else:
            y_jac = [q / t for q, t in zip(qpts_sorted, avg_jac_times)]
        label = f"{model} (jacobian)"

        plt.plot(qpts_sorted, y_jac, marker='x', linestyle='--', label=label)

# ------------------------------------------
# Finalize plot
# ------------------------------------------

plt.xlabel("Quadrature Points")
plt.ylabel("Bandwidth (GB/s)" if mode == "bandwidth" else "Throughput (qpts/s)")
plt.title(f"{mode.capitalize()} vs Quadrature Points")
plt.xscale("log")
plt.grid(True)
plt.legend(title="Model")
plt.tight_layout()

basename = os.path.splitext(os.path.basename(filename))[0]
plot_filename = f"{basename}-{mode}.png"
plt.savefig(plot_filename, dpi=300)
print(f"Plot saved as '{plot_filename}'")
