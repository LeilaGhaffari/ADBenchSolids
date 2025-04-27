#!/usr/bin/env python3

import re
import os
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np


def compute_bytes_per_qpts(model, op):
    STORED_VALUES = 0 # stream-triad
    if (model == "stream-residual"):
        STORED_VALUES = 16
    elif (model == "analytic-c"):
        STORED_VALUES = 16
    elif (model == "analytic-rust"):
        STORED_VALUES = 15
    elif (model == "enzyme-c"):
        STORED_VALUES = 15
    elif (model == "enzyme-rust"):
        STORED_VALUES = 15
    elif (model == "tapenade"):
        STORED_VALUES = 15
    elif (model == "adolc"):
        STORED_VALUES = 21

    BYTES_PER_QPTS = 8 * STORED_VALUES
    if (op == "residual"):
        BYTES_PER_QPTS = BYTES_PER_QPTS + 216 # 3 matrices (3*9*8)
    elif (op == "jacobian"):
        BYTES_PER_QPTS = BYTES_PER_QPTS + 144 # 2 matrices (2*9*8)

    return BYTES_PER_QPTS

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
# Plotting setup
# ------------------------------------------
plt.figure(figsize=(12, 7))

# Predefined unique markers and colors
markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'h', 'H', '+', 'x', 'd', '|', '_']
colors = [
    '#1f77b4',  # blue
    '#ff7f0e',  # orange
    '#2ca02c',  # green
    '#d62728',  # red
    '#9467bd',  # purple
    '#8c564b',  # brown
    '#e377c2',  # pink
    '#7f7f7f',  # gray
    '#bcbd22',  # olive
    '#17becf',  # cyan
    '#000000',  # black
    '#ffff00',  # yellow
    '#ff00ff',  # magenta
    '#00ffff',  # aqua
    '#800000',  # dark red
    '#008000',  # dark green
    '#000080',  # navy blue
    '#ffa500',  # strong orange
]

all_models = sorted(set(residual_times.keys()) | set(jacobian_times.keys()))

# ------------------------------------------
# Plot each model
# ------------------------------------------
for idx, model in enumerate(all_models):
    color = colors[idx % len(colors)]
    marker = markers[idx % len(markers)]

    # Residual
    if model in residual_times:
        qpts_sorted = sorted(residual_times[model])
        avg_res_times = [np.mean(residual_times[model][q]) for q in qpts_sorted]

        BYTES_PER_QPTS_RES = compute_bytes_per_qpts(model, "residual")

        if mode == "bandwidth":
            y_res = [q * BYTES_PER_QPTS_RES / t / 1e9 for q, t in zip(qpts_sorted, avg_res_times)]
        else:
            y_res = [q / t for q, t in zip(qpts_sorted, avg_res_times)]

        label = f"{model} (residual)" if "stream" not in model else model
        total_size = [q * BYTES_PER_QPTS_RES / 1e6 for q in qpts_sorted]

        plt.plot(total_size, y_res, marker=marker, color=color, label=label)

    # Jacobian
    if model in jacobian_times:
        qpts_sorted = sorted(jacobian_times[model])
        avg_jac_times = [np.mean(jacobian_times[model][q]) for q in qpts_sorted]

        BYTES_PER_QPTS_JAC = compute_bytes_per_qpts(model, "jacobian")
        if mode == "bandwidth":
            y_jac = [q * BYTES_PER_QPTS_JAC / t / 1e9 for q, t in zip(qpts_sorted, avg_jac_times)]
        else:
            y_jac = [q / t for q, t in zip(qpts_sorted, avg_jac_times)]

        label = f"{model} (jacobian)"
        total_size = [q * BYTES_PER_QPTS_JAC / 1e6 for q in qpts_sorted]

        plt.plot(total_size, y_jac, marker=marker, linestyle='--', color=color, label=label)

# ------------------------------------------
# Finalize plot
# ------------------------------------------
plt.xlabel("Input data (MB)")
plt.ylabel("Bandwidth (GB/s)" if mode == "bandwidth" else "Throughput (qpts/s)")
plt.title(f"{mode.capitalize()} vs Input Data")
plt.xscale("log")
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend(title="Model")
plt.tight_layout()

basename = os.path.splitext(os.path.basename(filename))[0]
plot_filename = f"{basename}-{mode}.png"
plt.savefig(plot_filename, dpi=300)
print(f"Plot saved as '{plot_filename}'")
