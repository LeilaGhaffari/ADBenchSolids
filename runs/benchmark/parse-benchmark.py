import re
import os
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

if len(sys.argv) != 2:
    print("Usage: python parse_benchmark.py <benchmark_output_file>")
    sys.exit(1)

filename = sys.argv[1]

# Nested dictionary: model -> qpts -> list of residual times
model_qpts_times = defaultdict(lambda: defaultdict(list))

with open(filename, "r") as f:
    lines = f.readlines()

    current_qpts = None
    for line in lines:
        # Match "Quadrature Points = ####"
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

        # Match line with model name and residual time
        model_match = re.match(r"\s*(\w[\w\-]*)\s+([\deE\.\+-]+)", line)
        if model_match and current_qpts is not None:
            model = model_match.group(1)
            residual_time = float(model_match.group(2))
            model_qpts_times[model][current_qpts].append(residual_time)


# Plot throughput for each model
plt.figure(figsize=(9, 6))

for model, qpts_data in model_qpts_times.items():
    qpts_sorted = sorted(qpts_data.keys())
    avg_times = [np.mean(qpts_data[q]) for q in qpts_sorted]
    throughput = [q / t for q, t in zip(qpts_sorted, avg_times)]
    plt.plot(qpts_sorted, throughput, marker='o', label=model)

plt.xlabel("Quadrature Points")
plt.ylabel("Throughput (qpts / residual time)")
plt.title("Throughput vs Quadrature Points")
plt.xscale("log")
plt.grid(True)
plt.legend(title="Model")
plt.tight_layout()

# Save plot
basename = os.path.splitext(os.path.basename(filename))[0]
plot_filename = f"{basename}.png"
plt.savefig(plot_filename, dpi=300)
print(f"Plot saved as '{plot_filename}'")
