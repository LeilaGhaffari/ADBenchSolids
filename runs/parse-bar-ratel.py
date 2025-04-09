#!/usr/bin/env python3
import os
import sys
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <log_directory>")
    sys.exit(1)

log_dir = sys.argv[1]
log_files = [f for f in os.listdir(log_dir) if f.endswith(".log")]

data = []

for filename in log_files:
    path = os.path.join(log_dir, filename)
    with open(path, "r") as f:
        content = f.read()

        # Determine method
        if "adolc" in filename.lower():
            method = "ADOL-C"
        elif "enzyme" in filename.lower():
            method = "Enzyme"
        else:
            method = "Analytic"

        # Get config from inside the file
        model_match = re.search(r"Model:\s+(.+?)\s*(?:\n|$)", content)
        model_desc = model_match.group(1).strip().lower() if model_match else "unknown"

        if "current" in model_desc:
            config = "Current"
        elif "initial" in model_desc:
            config = "Initial"
        else:
            config = "Unknown"

        # Get solve time
        time_match = re.search(r"Ratel Solve:\s+([0-9.eE+-]+)", content)
        if time_match:
            time = float(time_match.group(1))
            data.append({
                "Method": method,
                "Config": config,
                "Time": time
            })

# Prepare DataFrame
df = pd.DataFrame(data)

# Enforce method order
method_order = ["Analytic", "Enzyme", "ADOL-C"]
df["Method"] = pd.Categorical(df["Method"], categories=method_order, ordered=True)
df = df.sort_values(["Method", "Config"])

# Plot setup
fig, ax = plt.subplots(figsize=(8, 6))

bar_width = 0.35
x = np.arange(len(method_order))  # the label locations

initial_times = []
current_times = []

for method in method_order:
    times = df[df["Method"] == method].set_index("Config")["Time"]
    initial_times.append(times.get("Initial", np.nan))
    current_times.append(times.get("Current", np.nan))

# Bar plots
bars1 = ax.bar(x - bar_width/2, initial_times, bar_width, label="Initial", color="#1f77b4")
bars2 = ax.bar(x + bar_width/2, current_times, bar_width, label="Current", color="#ff7f0e")

# Labels and formatting
ax.set_ylabel("Time (s)")
#ax.set_title("Ratel Matrix-Free Solve Time by for 86490 Global DoFs")
ax.set_title("Ratel Matrix-Based Solve Time by for 26460 Global DoFs")
ax.set_xticks(x)
ax.set_xticklabels(method_order)
ax.set_yscale("log")
ax.legend(title="Config")
plt.tight_layout()
plt.savefig("ratel_solve_barplot.png")
print("Saved plot as 'ratel_solve_barplot.png'")
