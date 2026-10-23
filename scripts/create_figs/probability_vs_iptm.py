#!/usr/bin/env python3
"""Plot predicted affinity probability against ipTM for one or more proteins.

Each input CSV should contain the same Boltz prediction columns. Files are
combined before plotting, so each ligand panel includes data from all inputs.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


plt.rcParams.update({
    "font.size": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 16,
    "legend.fontsize": 10,
    "figure.dpi": 300,
})


LIGANDS = [
    "CER_d18_1_16_0_C16",
]
LIGAND_LABELS = {
    "CER_d18_1_16_0_C16": "C16",
}
LIGAND_COLORS = {
    "CER_d18_1_16_0_C16": "#4C72B0",
}


def load_csv_files(csv_files):
    """Read and combine CSVs, recording each input filename as a fallback ID."""
    frames = []
    for csv_file in csv_files:
        csv_path = Path(csv_file)
        if not csv_path.is_file():
            raise FileNotFoundError(f"CSV file not found: {csv_path}")

        frame = pd.read_csv(csv_path)
        frame["_source_file"] = csv_path.stem

        # A CSV may contain one protein without a protein column. Use its
        # filename in that case; also fill blank protein values this way.
        if "protein" not in frame.columns:
            frame["protein"] = csv_path.stem
        else:
            blank_protein = frame["protein"].isna() | (
                frame["protein"].astype(str).str.strip() == ""
            )
            frame.loc[blank_protein, "protein"] = csv_path.stem

        frames.append(frame)

    if not frames:
        raise ValueError("Provide at least one CSV file.")
    return pd.concat(frames, ignore_index=True, sort=False)


def create_probability_vs_iptm_scatter(
    data,
    output_path,
    ligand_col="ligand",
    protein_col="protein",
    probability_col="affinity_probability_binary",
    interface_col="ligand_iptm",
    highlight=None,
):
    """Create the two-panel C16 probability-versus-ipTM plot."""
    required = [ligand_col, protein_col, probability_col, interface_col]
    missing = [column for column in required if column not in data.columns]
    if missing:
        raise ValueError(
            "Input CSVs are missing required column(s): " + ", ".join(missing)
        )

    plot_data = data[data[ligand_col].isin(LIGANDS)].copy()
    if plot_data.empty:
        raise ValueError(
            f"No rows found for {', '.join(LIGANDS)} in column {ligand_col!r}."
        )

    plot_data[probability_col] = pd.to_numeric(
        plot_data[probability_col], errors="coerce"
    )
    plot_data[interface_col] = pd.to_numeric(
        plot_data[interface_col], errors="coerce"
    )
    plot_data = plot_data.dropna(subset=[probability_col, interface_col])
    if plot_data.empty:
        raise ValueError(
            f"No numeric values found in {probability_col!r} and {interface_col!r}."
        )

        fig, ax = plt.subplots(figsize=(7, 6))

    top_right_mask = (
        (plot_data[probability_col] > 0.5)
        & (plot_data[interface_col] > 0.5)
    )

    # Save top right quadrant to CSV
    top_right_columns = [protein_col]
    if "gene" in plot_data.columns:
        top_right_columns.append("gene")
    top_right_columns += [probability_col, interface_col]

    top_right = (
        plot_data.loc[top_right_mask, top_right_columns]
        .sort_values(probability_col, ascending=False)
    )

    fig, ax = plt.subplots(figsize=(7, 6))

    ligand = LIGANDS[0]
    ligand_data = plot_data[plot_data[ligand_col] == ligand]
    n_points = len(ligand_data)

    sns.scatterplot(
        data=ligand_data,
        x=probability_col,
        y=interface_col,
        color=LIGAND_COLORS[ligand],
        s=35,
        alpha=0.45,
        edgecolor="none",
        ax=ax,
    )

    # Highlight matching block
    highlighted_indices = set()

    if highlight:
        targets = {
            str(value).strip().casefold()
            for value in highlight
            if str(value).strip()
        }

        highlight_mask = pd.Series(False, index=ligand_data.index)
        for column in (protein_col, "gene"):
            if column in ligand_data.columns:
                highlight_mask |= (
                    ligand_data[column]
                    .astype("string")
                    .str.strip()
                    .str.casefold()
                    .isin(targets)
                )

        highlight_data = ligand_data.loc[highlight_mask]
        highlighted_indices = set(highlight_data.index)

        if highlight_data.empty:
            print(f"No matches found for --highlight {', '.join(highlight)}.")

        ax.scatter(
            highlight_data[probability_col],
            highlight_data[interface_col],
            color="red",
            s=70,
            alpha=1,
            edgecolor="black",
            zorder=3,
        )

        for _, row in highlight_data.iterrows():
            gene = row.get("gene")
            label = (
                str(gene)
                if pd.notna(gene) and str(gene).strip()
                else str(row[protein_col])
            )
            ax.annotate(
                label,
                xy=(row[probability_col], row[interface_col]),
                xytext=(6, 6),
                textcoords="offset points",
                color="red",
                fontsize=10,
                fontweight="bold",
                ha="left",
                va="bottom",
            )

    ax.axvline(0.5, color="gray", linestyle="--", linewidth=1, alpha=0.5)
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=1, alpha=0.5)

    top5 = ligand_data.nlargest(5, probability_col)
    for index, row in top5.iterrows():
        if index in highlighted_indices:
            continue
        gene = row.get("gene")
        label = (
            str(row[protein_col])
            if pd.isna(gene) or str(gene).strip() == ""
            else str(gene)
        )
        ax.annotate(
            label,
            xy=(row[probability_col], row[interface_col]),
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=9,
            fontweight="bold",
            ha="left",
            va="bottom",
        )

    # Labels
    ax.set_title(LIGAND_LABELS[ligand], fontsize=18)
    ax.set_xlabel("Predicted Affinity Probability")
    ax.set_ylabel("Ligand Interface Confidence (Ligand ipTM)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(True, linestyle="--", linewidth=0.5, color="0.85")
    ax.set_axisbelow(True)
    ax.text(
        0.98,
        0.02,
        f"N = {n_points}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=11,
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none")
    )

    fig.suptitle(
        "Boltz-2 Predicted Affinity Probability vs Ligand Interface Confidence",
        fontsize=20,
        y=1.02,
    )
    fig.tight_layout()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Print the top-right quadrant to csv
    top_right_path = output_path.with_name(
        f"{output_path.stem}_top_right.csv"
    )
    top_right.to_csv(top_right_path, index=False)
    print(f"Top-right quadrant list saved to: {top_right_path}")

    # Save scatterplot to file
    fig.savefig(output_path, format=output_path.suffix.lstrip(".") or "svg",
                bbox_inches="tight")
    print(f"Scatterplot saved to: {output_path}")
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot affinity probability versus ipTM for C16 "
            "using one or more protein prediction CSVs."
        )
    )
    parser.add_argument(
        "csv_files",
        nargs="+",
        type=Path,
        help="One or more prediction CSV files (all are combined in the plot).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("probability_vs_ligand_iptm.svg"),
        help="Output plot path (default: probability_vs_ligand_iptm.svg).",
    )
    parser.add_argument(
        "--highlight",
        nargs="+",
        help="Protein names or gene symbols to highlight in red.",
    )
    args = parser.parse_args()

    data = load_csv_files(args.csv_files)
    create_probability_vs_iptm_scatter(data, args.output, highlight=args.highlight)


if __name__ == "__main__":
    main()
