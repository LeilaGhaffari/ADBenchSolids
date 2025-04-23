#!/usr/bin/env python3

import re
import os
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

# Residual evaluation:
#   Load two 3×3 matrices = 2×9×8 = 144 B
#   Write one 3×3 matrix = 72 B
#   Write 15-element vector = 120 B
#   → Total = 336 B
BYTES_PER_QPT_RESIDUAL = 336

# Jacobian evaluation:
#   Load 3×3 matrix = 72 B
#   Load 15-vector = 120 B
#   Write 3×3 matrix = 72 B
#   → Total = 264 B
BYTES_PER_QPT_JACOBIAN = 264

if len(sys.argv) != 2:
    print("Usage: ./parse_total_times.py <benchmark_output_file>")
    sys.exit(1)

filename = sys.argv[1]

residual_times = defaultdict(list)
jacobian_times = defaultdict(list)
first_qpts = None

with open(filename, "r") as f:
    lines = f.readlines()

    current_qpts = None
    for line in lines:
        if "Quadrature Points" in line:
            match = re.search(r"Quadrature Points\s*=\s*(\d+)", line)
            if match:
                current_qpts = int(match.group(1))
                if first_qpts is None:
                    first_qpts = current_qpts
            continue

        if (
            line.strip() == "" or
            line.startswith("Iteration") or
            line.startswith("AD Tool") or
            line.startswith("-")
        ):
            continue

        tokens = line.strip().split()
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

models = sorted(set(residual_times.keys()) | set(jacobian_times.keys()))
x_pos = np.arange(len(models))

# Compute average throughput (GB/s)
res_throughputs = [
    (first_qpts * BYTES_PER_QPT_RESIDUAL) / (1e6 * np.mean(residual_times[m])) for m in models
]
jac_throughputs = [
    (first_qpts * BYTES_PER_QPT_JACOBIAN) / (1e6 * np.mean(jacobian_times[m])) for m in models
]

# Plot
plt.figure(figsize=(10, 6))
plt.scatter(x_pos, res_throughputs, color='blue', label='Residual', s=50)
plt.scatter(x_pos, jac_throughputs, color='orange', label='Jacobian', s=50, marker='x')

plt.yscale('log')
plt.ylabel("Throughput (MB/s) [log scale]")
plt.title(f"Residual and Jacobian Throughput")
plt.xticks(x_pos, models, rotation=45, ha='right')
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()

# Save output
basename = os.path.splitext(os.path.basename(filename))[0]
plot_filename = f"{basename}_throughput.png"
plt.savefig(plot_filename, dpi=300)
print(f"Dot plot saved as '{plot_filename}'")
