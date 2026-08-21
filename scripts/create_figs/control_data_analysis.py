import argparse
import os
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from pyrolite.plot import pyroplot
from matplotlib.lines import Line2D

from generate_plots import make_plot_folder

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 16,
    'axes.labelsize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 18,
    'figure.dpi':300,
})

OUTPUT_DIR = Path(os.environ["OUTPUT_DIR"])
SOURCE_DIR = Path(os.environ["SOURCE_DIR"])
PROCESSED_DIR = Path(os.environ["PROCESSED_DIR"])


def scatter_plot_experimental_vs_predicted(data, output_file, log_scale=False):
    predicted_affinity_field = 'micromolar_affinity_pred_value'
    experimental_affinity_field = 'Experimental_Affinity'

    x = data[experimental_affinity_field]
    y = data[predicted_affinity_field]

    plt.scatter(x, y, alpha=0.5)

    if log_scale:
        log_x, log_y = np.log10(x), np.log10(y)
        
        result = linregress(log_x, log_y)
        m = result.slope
        b = result.intercept
        r2 = result.rvalue ** 2

        min_val = min(log_x.min(), log_y.min())
        max_val = max(log_x.max(), log_y.max())

        x_line = np.geomspace(x.min(), x.max(), 300)
        y_line = (10 ** b) * (x_line ** m)

        plt.plot(x_line, y_line, color='red', 
                 linewidth=2, 
                 label=f"Fit: y={10**b:.2f}x^{m:.2f}, R^2={r2:.3f}")

        plt.xscale('log')
        plt.yscale('log')
        plt.title("Experimental vs Predicted Affinity (µM, log-log scale)")

    else:
        result = linregress(x, y)
        m = result.slope
        b = result.intercept
        r2 = result.rvalue ** 2

        min_val = min(x.min(), y.min())
        max_val = max(x.max(), y.max())

        x_line = np.linspace(min_val, max_val, 300)
        y_line = m * x_line + b

        plt.plot(x_line, y_line, color='red', 
                 linewidth=2, 
                 label=f"Fit: y={m:.2f}x+{b:.2f}, R^2={r2:.3f}")
        
        plt.title("Experimental vs Predicted Affinity (µM)")

    # X = Y
    max_val = max(x.max(), y.max())
    plt.plot([0, max_val], [0, max_val], color='gray', linestyle='--', label="y=x")

    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend()
    plt.xlabel("Experimental Affinity (µM)")
    plt.ylabel("Predicted Affinity (µM)")
    plt.savefig(output_file.with_suffix('.svg'))
    plt.close()


def heat_scatter_plot_experimental_vs_predicted(data, output_file):
    predicted_affinity_field = 'neg_log_molar_affinity_pred_value'
    experimental_affinity_field = 'neg_log_experimental_affinity' 
    
    plot_data = data[[experimental_affinity_field, predicted_affinity_field]]

    fig = plt.figure(figsize=(7, 7))

    # Main plot position: [left, bottom, width, height]
    ax = fig.add_axes([0.15, 0.15, 0.80, 0.75])

    plot_data.pyroplot.heatscatter(s=10, ax=ax)

    x = plot_data[experimental_affinity_field]
    y = plot_data[predicted_affinity_field]

    ax.set_xlabel("Experimental Affinity (-log10(M))")
    ax.set_ylabel("Predicted Affinity (-log10(M))")

    m, b, r_value, _, _ = linregress(x, y)
    r2 = r_value ** 2

    min_val = min(x.min(), y.min())
    max_val = max(x.max(), y.max())

    x_line = np.linspace(min_val, max_val, 300)
    y_line = m * x_line + b

    fit_line, = ax.plot(x_line, y_line, color='red', 
            linewidth=2, 
            label=f"Linear regression"
            )

    r2_text = f"$R^2$ = {r2:.3f}\nn = 1549"
    r2_label = Line2D([], [], color='none', label=r2_text)

    # ax.text(
    #     0.05, 0.95,
    #     f"$R^2$ = {r2:.3f}",
    #     transform=ax.transAxes,
    #     va="top",
    #     ha="left",
    #     fontsize=14,
    #     bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="none")
    # )

    #plt.title("Heat Scatter Plot of Predicted vs Experimental Affinities")
    ax.legend(handles=[fit_line, r2_label], loc="upper left", frameon=True, bbox_to_anchor=(0.02, 0.98))
    ax.grid(True, linestyle='--', linewidth=0.5, color='0.85', alpha=0.7, zorder=0)
    ax.set_axisbelow(True)  # ensures grid stays behind the data points
    bbox=ax.get_position()
    print(bbox.x0, bbox.y0, bbox.width, bbox.height)
    plt.savefig(output_file.with_suffix('.svg'), bbox_inches=None)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    # ap.add_argument("--input_folder", type=str, required=True)
    ap.add_argument("--output_folder", type=str, default="control_tests")
    ap.add_argument("--ref_file", type=str, default="control_proteins.csv")
    args = ap.parse_args()

    ref_file_path = SOURCE_DIR / args.ref_file
    # output_folder_path = OUTPUT_DIR / args.output_folder

    reference_data = pd.read_csv(ref_file_path)
    confidence_df = pd.read_csv(PROCESSED_DIR / f"confidence_predictions_{args.output_folder}.csv")

    combined_data = pd.merge(
        confidence_df,
        reference_data,
        left_on=["protein", "ligand"],
        right_on=["PDB", "ligand_id"],
        how="left",
    )

    combined_data = combined_data.drop_duplicates(subset=["protein", "ligand", "micromolar_affinity_pred_value", "Experimental_Affinity"])

    # Add -log10(molar experimental affinity)
    molar_experimental_affinity = combined_data["Experimental_Affinity"] / 1e6
    combined_data["neg_log_experimental_affinity"] = -np.log10(molar_experimental_affinity)

    # print(f'Found {len(combined_data)} unique protein-ligand pairs')

    csv_file = PROCESSED_DIR / f"confidence_predictions_{args.output_folder}_combined.csv"

    combined_data.to_csv(csv_file, index=False)

    plot_folder = make_plot_folder(args.output_folder)
    
    scatter_plot_experimental_vs_predicted(combined_data, 
                                           plot_folder / f"scatter_plot_{args.output_folder}.png")
    
    scatter_plot_experimental_vs_predicted(combined_data, 
                                           plot_folder / f"log_scatter_plot_{args.output_folder}.png", 
                                           log_scale=True)

    heat_scatter_plot_experimental_vs_predicted(combined_data,
                                               plot_folder / f"heat_scatter_plot_{args.output_folder}.png")


if __name__ == "__main__":
    main()
