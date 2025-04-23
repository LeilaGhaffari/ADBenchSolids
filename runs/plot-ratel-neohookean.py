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
    "Analytic": {"color": "red"},
    "Enzyme": {"color": "green"},
    "ADOL-C": {"color": "black"},
}

# P for tetrahedra: triangles; Q for hexahedra: squares
ELEMENT_STYLE = {
    "P1": {"marker": "^", "size": 40},
    "P2": {"marker": "^", "size": 160},
    "Q1": {"marker": "s", "size": 40},
    "Q2": {"marker": "s", "size": 160},
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

    # Extract p/q element type
    element_match = re.search(r"-(p\d|q\d)-", fname)
    element_type = element_match.group(1).upper() if element_match else "Unknown"

    # Extract mesh size (optional)
    size_match = re.search(r"size(\d+)", fname)
    mesh_size = int(size_match.group(1)) if size_match else None

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
                "ElementType": element_type,
                "MeshSize": mesh_size,
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

ELEMENT_OFFSETS = {
    "P1": -0.25,
    "P2": -0.10,
    "Q1": 0.10,
    "Q2": 0.25,
}
METHOD_OFFSETS = {
    "Analytic": -0.05,
    "Enzyme": 0.05,
    "ADOL-C": 0.0,
}

for method in METHOD_ORDER:
    method_subset = df[df["Method"] == method]
    for category, base_x in GROUP_POSITIONS.items():
        matrix_type, config = category.split(" / ")
        group = method_subset[
            (method_subset["Config"] == config)
            & (method_subset["MatrixType"] == matrix_type)
        ]
        for _, row in group.iterrows():
            element = row["ElementType"]
            y = row["Throughput"]

            element_offset = ELEMENT_OFFSETS.get(element, 0.0)
            method_offset = METHOD_OFFSETS.get(method, 0.0)
            x = base_x + element_offset + method_offset

            marker = ELEMENT_STYLE.get(element, {"marker": "x", "size": 70})
            ax.scatter(
                x, y,
                color=STYLE_MAP[method]["color"],
                marker=marker["marker"],
                s=marker["size"],
                alpha=0.9,
                label=None
            )

# Axes & Labels
ax.set_yscale("log")
ax.set_ylabel("Throughput (DoFs/s)")
ax.set_title("Ratel Solve Throughput", pad=25)
ax.set_xticks(list(GROUP_POSITIONS.values()))
ax.set_xticklabels(["Initial", "Current", "Initial", "Current"])

# Group-level labels
ax.text(0.25, -0.15, "Sparse-Matrix", ha="center", va="center", transform=ax.transAxes, fontsize=12)
ax.text(0.75, -0.15, "Matrix-Free", ha="center", va="center", transform=ax.transAxes, fontsize=12)
ax.axvline(2, color='gray', linestyle='--', linewidth=1)
ax.grid(True, which='both', linestyle='--', linewidth=0.6, alpha=0.6)

# Legend 1: Method
method_handles = [
    plt.Line2D([0], [0], marker='o', color='w',
               label=m, markerfacecolor=STYLE_MAP[m]["color"],
               markeredgecolor='k', markersize=8)
    for m in METHOD_ORDER
]
legend1 = ax.legend(handles=method_handles, title="Method", loc="upper left")

# Legend 2: Element Type
element_order = ["Q1", "Q2", "P1", "P2"]
used_elements = [e for e in element_order if e in df["ElementType"].unique()]

element_handles = [
    plt.Line2D(
        [0], [0],
        marker=ELEMENT_STYLE[e]["marker"],
        color='w',
        label=e,
        markerfacecolor='gray',
        markeredgecolor='k',
        markersize=ELEMENT_STYLE[e]["size"] ** 0.5
    )
    for e in used_elements
]

legend2 = ax.legend(
    handles=element_handles,
    title="Element Type",
    loc="center right"
)


ax.add_artist(legend1)

plt.tight_layout(rect=[0, 0, 1, 0.9])
plt.savefig("ratel_neohookean.png")
print("✅ Saved plot as 'ratel_neohookean.png'")
