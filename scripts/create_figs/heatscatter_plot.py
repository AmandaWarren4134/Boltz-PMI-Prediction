import os
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy import stats
from pathlib import Path
from pyrolite.plot import pyroplot
from matplotlib.lines import Line2D

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 16,
    'axes.labelsize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 18,
    'figure.dpi': 300,
})


# PROJECT_ROOT = Path(os.environ["PROJECT_ROOT"])
# SCRIPT_DIR = Path(os.environ["SCRIPT_DIR"])
# YAML_DIR = Path(os.environ["YAML_DIR"])
# OUTPUT_DIR = Path(os.environ["OUTPUT_DIR"])
PROCESSED_DIR = Path(os.environ["PROCESSED_DIR"])
PLOT_DIR = Path(os.environ["PLOT_DIR"])

def load_csv_data(csv_file: str):
    csv_file_path = PROCESSED_DIR / csv_file
    
    if not csv_file_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_file}")
    
    data = pd.read_csv(csv_file_path)
    return data


def create_heat_scatter_plot(data: pd.DataFrame, ligand_name_1: str, ligand_name_2: str, affinity_col: str, regression_line=False, xy_line=False):
    
    if regression_line and xy_line:
        raise ValueError("Cannot have both regression_line and xy_line set to True. Please choose one.")
       
    x_col = f"{ligand_name_1}_{affinity_col}"
    y_col = f"{ligand_name_2}_{affinity_col}"

    x = data[x_col]
    y = data[y_col]

    data = data[[x_col, y_col]]

    fig, ax = plt.subplots(figsize=(7, 7))
    plt.subplots_adjust(left=0.15, right=0.95, top=0.90, bottom=0.15)

    data.pyroplot.heatscatter(s=10, ax=ax)
    
    ax.set_xlabel(f"C16 Predicted Affinity (-log10(M))")
    ax.set_ylabel(f"C2 Predicted Affinity (-log10(M))")
    ax.grid(True, linestyle='--', linewidth=0.5, color='0.85', alpha=0.7, zorder=0)
    ax.set_axisbelow(True)  # ensures grid stays behind the data points

    if xy_line:
        # Add X = Y line for reference
        min_val = min(x.min(), y.min())
        max_val = max(x.max(), y.max())
        line,  = ax.plot([min_val, max_val], [min_val, max_val], 
                         color='blue', linestyle='solid', label='X = Y')
        
        t_stat, p_val = perform_paired_t_test(x, y)
        t_text = f"t = {t_stat:.3f}"

        if p_val < 0.0001:
            p_text = "$p < 0.0001$"
        else:
            p_text = f"$p = {p_val:.3f}$"
        stat_label = Line2D([], [], color='none', label=f"{t_text}, {p_text}")

        ax.text(
            0.05, 0.85,
            f"Paired $t$-test\n{p_text}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=16,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=0.3)
        )
        ax.legend(handles=[line], loc="upper left", frameon=True, bbox_to_anchor=(0.02, 0.98))
        #ax.set_title("Comparison of Ceramide Predicted Affinities")
        plot_title = f"xy_line_heat_scatter_with_p_val_{ligand_name_1}_vs_{ligand_name_2}.svg"
        
    if regression_line:
        slope, intercept, r_value, p_value, std_err = linear_regression(x, y)
        m = slope
        b = intercept
        r2 = r_value ** 2

        min_val = min(x.min(), y.min())
        max_val = max(x.max(), y.max())
        x_line = np.linspace(min_val, max_val, 1000)
        y_line = m * x_line + b

        fit_line, = ax.plot(x_line, y_line, color='red', 
                 linestyle='solid', linewidth=2, 
                 label=f"Linear regression")
        
        r2_text = f"$R^2$ = {r2:.3f}"
        r2_label = Line2D([], [], color='none', label=r2_text)

        ax.legend(handles=[fit_line, r2_label], loc="upper left", frameon=True, bbox_to_anchor=(0.02, 0.98))
        ax.set_title("Regression of Ceramide Predicted Affinities")
        plot_title = f"regress_heat_scatter_{ligand_name_1}_vs_{ligand_name_2}.svg"

    plot_file_path = PLOT_DIR / plot_title
    plt.savefig(plot_file_path.with_suffix('.svg'))
    print(f"Heat scatter plot saved to: {plot_title}")
    bbox=ax.get_position()
    print(bbox.x0, bbox.y0, bbox.width, bbox.height)
    plt.show()


def perform_paired_t_test(x, y):
    t_statistic, p_value = stats.ttest_rel(x, y)
    return t_statistic, p_value


def linear_regression(x, y):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    
    return slope, intercept, r_value, p_value, std_err


def main():
    # This script makes a heat scatter plot showing the affinity predictions for two different
    # ligands on a set of proteins. Each point corresponds to a protein and the X/Y coordinates
    # correspond to predicted affinity

    # It assumes two CSV files, each run with a different ligand, but same set of proteins
    # The CSV file names are specified with --csvs
    # The --ligand_names argument can be used to label the ligands for the two CSVs
    
    ap = argparse.ArgumentParser()
    ap.add_argument("--csvs", type=str, nargs=2, 
                    default=["split_predictions/confidence_predictions_random_c16.csv", "split_predictions/confidence_predictions_random_c16dihydro.csv"])
    ap.add_argument("--ligand_names", type=str, nargs=2, 
                    default=["c16", "c16_dihydro"])
    args = ap.parse_args()

    file_name_1, file_name_2 = args.csvs
    ligand_name_1, ligand_name_2 = args.ligand_names

    data_1 = load_csv_data(file_name_1)
    data_2 = load_csv_data(file_name_2)

    affinity_col = "neg_log_molar_affinity_pred_value"
    target_cols = ["protein", affinity_col]

    data_1 = data_1[target_cols]
    data_2 = data_2[target_cols]

    data_1 = data_1.rename(columns={affinity_col: f"{ligand_name_1}_{affinity_col}"})
    data_2 = data_2.rename(columns={affinity_col: f"{ligand_name_2}_{affinity_col}"})

    merged_data = pd.merge(data_1, data_2, on="protein")
    merged_data = merged_data.drop("protein", axis=1)

    #create_heat_scatter_plot(merged_data, ligand_name_1, ligand_name_2, affinity_col, regression_line=True)
    create_heat_scatter_plot(merged_data, ligand_name_1, ligand_name_2, affinity_col, xy_line=True)


if __name__ == "__main__":
    main()
