import re
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

if len(sys.argv) != 2:
    print("Usage: python parse_benchmark.py <benchmark_output_file>")
    sys.exit(1)

filename = sys.argv[1]

# Dictionary to collect times per qpts
qpts_times = defaultdict(list)

with open(filename, "r") as f:
    lines = f.readlines()

    current_qpts = None
    for line in lines:
        # Match qpts line
        qpts_match = re.search(r"Quadrature Points\s*=\s*(\d+)", line)
        if qpts_match:
            current_qpts = int(qpts_match.group(1))
            continue

        # Match model + residual time line
        model_match = re.match(r"\s*\w+\s+([\deE\.\+-]+)", line)
        if model_match and current_qpts is not None:
            residual_time = float(model_match.group(1))
            qpts_times[current_qpts].append(residual_time)
            current_qpts = None

# Sort qptss and compute averages and derived metric (qpts / time)
qpts_sorted = sorted(qpts_times.keys())
avg_times = [np.mean(qpts_times[qpts]) for qpts in qpts_sorted]
throughput = [qpts / t for qpts, t in zip(qpts_sorted, avg_times)]

# Plot
plt.figure(figsize=(8, 5))
plt.plot(qpts_sorted, throughput, marker='x', linestyle='--', label='qpts / Time (Throughput)')
plt.xlabel("Quadrature Points")
plt.ylabel("Throughput (#quadrature points/time)")
plt.xscale("log")
#plt.yscale("log")
plt.title("Throughput vs Quadrature Points")
plt.grid(True)
plt.legend()
plt.tight_layout()

# Save the plot
plot_filename = "throughput.png"
plt.savefig(plot_filename, dpi=300)
print(f"Plot saved as '{plot_filename}'")
