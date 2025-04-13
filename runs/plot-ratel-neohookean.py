#!/usr/bin/env python3
import os
import sys
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ------------------- CONFIGURATION -------------------
METHOD_MAP = {
    "adolc": "ADOL-C",
    "enzyme": "Enzyme",
}
CONFIG_KEYS = ["initial", "current"]
MATRIX_KEYS = {"free": "Matrix-Free", "sparse": "Sparse-Matrix"}

STYLE_MAP = {
    "Analytic": {"color": "red", "marker": "s"},
    "Enzyme": {"color": "green", "marker": "^"},
    "ADOL-C": {"color": "black", "marker": "o"},
}

METHOD_ORDER = ["Analytic", "Enzyme", "ADOL-C"]

GROUP_POSITIONS = {
    "Sparse-Matrix / Initial": 0,
    "Sparse-Matrix / Current": 1,
    "Matrix-Free / Initial": 3,
    "Matrix-Free / Current": 4,
}

# ------------------- HELPER FUNCTION -------------------
def parse_log_file(filepath):
    fname = os.path.basename(filepath).lower()

    config = next((c.title() for c in CONFIG_KEYS if c in fname), "Unknown")
    matrix_type = next((MATRIX_KEYS[k] for k in MATRIX_KEYS if k in fname), "Unknown")
    method = next((v for k, v in METHOD_MAP.items() if k in fname), "Analytic")

    with open(filepath, "r") as f:
        content = f.read()
        dofs_match = re.search(r"Global DoFs:\s+(\d+)", content)
        time_match = re.search(r"Ratel Solve:\s+([0-9.eE+-]+)", content)

        if dofs_match and time_match:
            dofs = int(dofs_match.group(1))
            time = float(time_match.group(1))
            throughput = dofs / time if time > 0 else np.nan
            return {
                "Method": method,
                "Config": config,
                "MatrixType": matrix_type,
                "Throughput": throughput
            }
    return None

# ------------------- MAIN SCRIPT -------------------
if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <log_directory>")
    sys.exit(1)

log_dir = sys.argv[1]
log_files = [f for f in os.listdir(log_dir) if f.endswith(".log")]

data = [parse_log_file(os.path.join(log_dir, f)) for f in log_files]
data = [d for d in data if d is not None]
df = pd.DataFrame(data)

# ------------------- PLOTTING -------------------
fig, ax = plt.subplots(figsize=(10, 6))

for method in METHOD_ORDER:
    for category, x in GROUP_POSITIONS.items():
        matrix_type, config = category.split(" / ")
        subset = df[(df["Method"] == method) & (df["Config"] == config) & (df["MatrixType"] == matrix_type)]
        if not subset.empty:
            y = subset["Throughput"].values[0]
            ax.scatter(
                x, y,
                color=STYLE_MAP[method]["color"],
                marker=STYLE_MAP[method]["marker"],
                s=60,
                label=method if x == 0 else None
            )

# Axes & Labels
ax.set_yscale("log")
ax.set_ylabel("Throughput (DoFs/s)")
ax.set_title("Ratel Solve Throughput (Grouped by Matrix Type)", pad=25)
ax.set_xticks(list(GROUP_POSITIONS.values()))
ax.set_xticklabels(["Initial", "Current", "Initial", "Current"])

# Group-level labels
ax.text(0.25, -0.15, "Sparse-Matrix", ha="center", va="center", transform=ax.transAxes, fontsize=12)
ax.text(0.75, -0.15, "Matrix-Free", ha="center", va="center", transform=ax.transAxes, fontsize=12)
ax.axvline(2, color='gray', linestyle='--', linewidth=1)
ax.grid(True, which='both', linestyle='--', linewidth=0.6, alpha=0.6)

# Legend (auto-positioned)
handles = [
    plt.Line2D([0], [0], marker=STYLE_MAP[m]["marker"], color='w',
               label=m, markerfacecolor=STYLE_MAP[m]["color"],
               markeredgecolor='k', markersize=8)
    for m in METHOD_ORDER
]
ax.legend(handles=handles, title="Method", loc="best")

plt.tight_layout(rect=[0, 0, 1, 0.9])
plt.savefig("ratel_neohookean.png")
print("✅ Saved plot as 'ratel_neohookean.png'")
