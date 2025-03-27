#!/usr/bin/env python3

import re
import os
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

first_qpts = None

if len(sys.argv) != 2:
    print("Usage: ./parse_total_times.py <benchmark_output_file>")
    sys.exit(1)

filename = sys.argv[1]

# model -> list of residual and jacobian times
residual_times = defaultdict(list)
jacobian_times = defaultdict(list)

with open(filename, "r") as f:
    lines = f.readlines()

    current_qpts = None
    for line in lines:
        # Match "Quadrature Points = ####"
        if "Quadrature Points" in line:
            match = re.search(r"Quadrature Points\s*=\s*(\d+)", line)
            if match:
                current_qpts = int(match.group(1))
                if first_qpts is None:
                    first_qpts = current_qpts
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

        # Skip if line doesn't have enough tokens
        if len(tokens) < 4:
            continue

        # First token(s) = model name; rest = numeric values
        numeric_start = next((i for i, t in enumerate(tokens) if re.match(r"^[\deE\.\+\-]+$", t)), None)
        if numeric_start is None or numeric_start + 3 >= len(tokens):
            continue

        model = " ".join(tokens[:numeric_start])
        try:
            res_time = float(tokens[numeric_start])
            jac_time = float(tokens[numeric_start + 2])
        except ValueError:
            continue

        residual_times[model].append(res_time)
        jacobian_times[model].append(jac_time)
        current_qpts = None

# Compute average total time per model
models = sorted(set(residual_times.keys()) | set(jacobian_times.keys()))
total_times = []

for model in models:
    res_avg = np.mean(residual_times[model]) if model in residual_times else 0.0
    jac_avg = np.mean(jacobian_times[model]) if model in jacobian_times else 0.0
    total_times.append(res_avg + jac_avg)

# Plot bar chart
plt.figure(figsize=(10, 6))
plt.bar(models, total_times, color='skyblue')
plt.yscale("log")
plt.ylabel("Time (s)")
plt.title(f"Time per Model ({first_qpts:,} data points)")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

# Save output
basename = os.path.splitext(os.path.basename(filename))[0]
plot_filename = f"{basename}.png"
plt.savefig(plot_filename, dpi=300)
print(f"Bar plot saved as '{plot_filename}'")
