import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 12,
    'figure.dpi': 300,
    'axes.grid': True,
    'grid.linestyle': '--',
    'grid.linewidth': 0.5,
    'grid.color': '0.85',
    'grid.alpha': 0.7,
    'axes.axisbelow': True,
})

PROCESSED_DIR = Path(os.environ["PROCESSED_DIR"])
PLOT_DIR = Path(os.environ["PLOT_DIR"])

AFFINITY_COL = "confidence_score"

files = {
    "C16 (PISA)": PROCESSED_DIR / "pisa_c16_human.csv",
    "C16 (Random)": PROCESSED_DIR / "split_predictions" / "confidence_predictions_random_c16.csv",
    "C16 (Effectors)": PROCESSED_DIR / "confidence_predictions_ceramide_effectors" / "CER_d18_1_16_0_C16.csv",

    "C16 DH (PISA)": PROCESSED_DIR / "pisa_c16dihydro.csv",
    "C16 DH (Random)": PROCESSED_DIR / "split_predictions" / "confidence_predictions_random_c16dihydro.csv",
    "C16 DH (Effectors)": PROCESSED_DIR / "confidence_predictions_ceramide_effectors" / "CER_d18_0_16_0_C16.csv",
}

dfs = []

for label, filepath in files.items():
    df = pd.read_csv(filepath)
    df = df[["protein", AFFINITY_COL]].copy()
    df["Group"] = label

    # Assign ligand
    if "dihydro" in label.lower() or "dh" in label.lower():
        df["Ligand"] = "C16 DH"
    else:
        df["Ligand"] = "C16"

    # Assign dataset
    if "PISA" in label:
        df["Dataset"] = "PISA"
    elif "Random" in label:
        df["Dataset"] = "Random"
    elif "Effectors" in label:
        df["Dataset"] = "Effectors"

    dfs.append(df)

plot_data = pd.concat(dfs, ignore_index=True)

# ---------------------------------------------------------------------
# Keep ONLY proteins present for BOTH C16 and C16 DH within each Dataset
# ---------------------------------------------------------------------

protein_counts = (
    plot_data
    .groupby(["Dataset", "protein"])["Ligand"]
    .nunique()
)

shared_proteins = (
    protein_counts[protein_counts == 2]
    .reset_index()[["Dataset", "protein"]]
)

plot_data = plot_data.merge(
    shared_proteins,
    on=["Dataset", "protein"],
    how="inner"
)

# ---------------------------------------------------------------------
# Calculate n = number of UNIQUE proteins shared between C16 and C16 DH
# for each dataset
# ---------------------------------------------------------------------

counts = (
    plot_data
    .groupby("Dataset")["protein"]
    .nunique()
)

# Optional: print the number of shared proteins
print("\nNumber of proteins common to C16 and C16 DH:")
for dataset in ["PISA", "Random", "Effectors"]:
    print(f"{dataset}: n={counts.get(dataset, 0)}")

# ---------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------

order = [
    "PISA",
    "Random",
    "Effectors",
]

FIGSIZE = (9, 7)
LEFT, RIGHT, TOP, BOTTOM = 0.15, 0.95, 0.90, 0.15

fig, ax = plt.subplots(figsize=FIGSIZE)

fig.subplots_adjust(
    left=LEFT,
    right=RIGHT,
    top=TOP,
    bottom=BOTTOM
)

sns.violinplot(
    data=plot_data,
    x="Dataset",
    y=AFFINITY_COL,
    hue="Ligand",
    order=order,
    hue_order=["C16", "C16 DH"],
    palette={
        "C16": "#66C2A5",
        "C16 DH": "#FC8D62",
    },
    inner="box",
    density_norm="area",
    ax=ax,
)

ax.grid(axis='y')
ax.grid(axis='x', visible=False)
ax.set_axisbelow(True)

ax.set_ylabel("Boltz Confidence Score")
ax.set_xlabel("")

# Keep legend
ax.legend(
    title="Ligand",
    frameon=True,
    loc="upper left"
)

# ---------------------------------------------------------------------
# Sample size annotations
# n = unique proteins shared between C16 and C16 DH
# ---------------------------------------------------------------------

ymin, ymax = ax.get_ylim()
y_offset = ymin - 0.07 * (ymax - ymin)

for i, label in enumerate(order):
    n = counts.get(label, 0)

    ax.text(
        i,
        y_offset,
        f"n={n}",
        ha="center",
        va="top",
        fontsize=16,
        color="0.3"
    )

outfile = PLOT_DIR / "hue_confidence_c16_vs_dihydro_pisa_random_effectors_violinplot.svg"

plt.savefig(
    outfile,
    dpi=300,
    bbox_inches="tight"
)

plt.show()
