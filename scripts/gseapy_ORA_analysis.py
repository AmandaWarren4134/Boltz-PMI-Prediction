#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import gseapy as gp
from mygene import MyGeneInfo
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

# -------------------- Helper Functions -------------------- #
def get_good_binders(df, affinity_col, method="percentile", value=0.1, absolute_cutoff=None):
    """
    Select top binders based on affinity score (-log10 already applied).
    Percentile method selects top x% strongest binders (higher score better).
    Absolute method selects binders above a threshold.
    Returns subset dataframe and cutoff value.
    """
    if method == "percentile":
        cutoff = df[affinity_col].quantile(1 - value)
        good = df[df[affinity_col] >= cutoff]
    elif method == "absolute":
        if absolute_cutoff is None:
            raise ValueError("absolute_cutoff must be provided for absolute method")
        cutoff = absolute_cutoff
        good = df[df[affinity_col] >= absolute_cutoff]
    else:
        raise ValueError("method must be 'percentile' or 'absolute'")
    return good, cutoff


def map_to_gene_symbols(accessions):
    mg = MyGeneInfo()
    results = mg.querymany(
        accessions,
        scopes="uniprot",
        fields="symbol",
        species=9606
    )
    mapping = {r["query"]: r["symbol"] for r in results if "symbol" in r and "query" in r}
    return mapping


def run_enrichment(protein_list, gene_sets, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    enr = gp.enrichr(
        gene_list=protein_list,
        gene_sets=gene_sets,
        organism="homo sapiens",
        outdir=str(outdir),
        cutoff=0.05,
        no_plot=True
    )
    return enr


def plot_top_enrichment(results, outdir, p, title, ligand=None, top_n=15, x_col="Combined Score"):
    """
    Plot top enriched pathways.
    """
    if results is None or results.empty:
        print("No enrichment results to plot.")
        return

    results_sorted = results.sort_values("Adjusted P-value").head(top_n)

    # For FDR, take -log10
    if x_col.lower() in ["adjusted p-value", "fdr"]:
        results_sorted["neg_log10_fdr"] = -np.log10(results_sorted["Adjusted P-value"])
        x_col = "neg_log10_fdr"

    results_sorted = results_sorted.sort_values(x_col, ascending=False).head(top_n)

    plt.figure(figsize=(10, 6))
    sns.barplot(
        data=results_sorted,
        y="Term",
        x=x_col,
        color="skyblue"
    )

    if x_col == "neg_log10_fdr":
        plt.axvline(x=1.3, color='red', linestyle='--', linewidth=1)
        plt.text(1.32, 18, "FDR = 0.05", color='red', fontsize=10)

    xlabel_map = {
        "Combined Score": "Combined Score",
        "neg_log10_fdr": "-log10(FDR)",
        "Enrichment Ratio": "Enrichment Ratio",
        "Overlap": "Gene Count"
    }

    plt.xlabel(xlabel_map.get(x_col, x_col))
    plt.ylabel("")
    title_prefix = f"{ligand} - " if ligand else ""
    plt.title(f"{title} - {title_prefix}Top {int(p*100)}% Predicted Binders")
    plt.tight_layout()
    filename = f"top_enrichment_barplot_{x_col.replace(' ', '_')}.png"
    plt.savefig(outdir / filename)
    plt.close()
    print(f"Saved barplot to {filename}")


def run_enrichment_for_dataframe(df_subset, outdir, args, ligand=None):
    """
    Process a subset of dataframe (all proteins or per-ligand)
    """
    for p in args.percentiles:
        print(f"\nProcessing top {int(p*100)}% binders...")
        good_df, cutoff = get_good_binders(
            df_subset,
            affinity_col=args.affinity_col,
            method=args.method,
            value=p,
            absolute_cutoff=args.absolute_cutoff
        )
        print(f"Cutoff used: {cutoff}")
        print(f"Number of good binders: {len(good_df)}")

        gene_list = good_df["gene_symbol"].dropna().unique().tolist()
        print(f"Mapped to {len(gene_list)} gene symbols")

        if len(gene_list) == 0:
            print("No valid gene symbols found. Skipping this percentile.")
            continue

        percentile_outdir = outdir / f"top_{int(p*100)}pct"
        percentile_outdir.mkdir(parents=True, exist_ok=True)
        print("Running enrichment...")
        enr = run_enrichment(gene_list, args.gene_sets, percentile_outdir)

        results = enr.results
        if results is not None and not results.empty:
            results_sorted = results.sort_values("Adjusted P-value")
            results_sorted.to_csv(percentile_outdir / f"enrichment_top_{int(p*100)}pct.csv", index=False)
            print(f"Saved results to {percentile_outdir / f'enrichment_top_{int(p*100)}pct.csv'}")
            plot_top_enrichment(results_sorted, percentile_outdir, p, title=args.title, ligand=ligand, top_n=20, x_col=args.x_col)


# -------------------- Main -------------------- #
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--affinity_col", default="neg_log_molar_affinity_pred_value")
    parser.add_argument("--method", choices=["percentile", "absolute"], default="percentile")
    parser.add_argument("--percentiles", nargs="+", type=float, default=[0.01, 0.05, 0.10])
    parser.add_argument("--absolute_cutoff", type=float, default=None)
    parser.add_argument("--x_col", default="Adjusted P-value")
    parser.add_argument("--gene_sets", nargs="+", default=[
        "GO_Biological_Process_2021",
        "GO_Molecular_Function_2021",
        "KEGG_2021_Human"
    ])
    parser.add_argument("--by_ligand", action="store_true", help="Run enrichment separately for each ligand")
    parser.add_argument("--ligand_column", default="ligand", help="Column name for ligands if --by_ligand is used")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory to place the enrichment analysis results")
    parser.add_argument("--title", default="Enrichment", help="Title text before the top %.")
    args = parser.parse_args()

    # Load data
    df = pd.read_csv(args.csv)

    # Map UniProt IDs to gene symbols
    accessions = df["protein"].dropna().unique().tolist()
    print(f"Mapping {len(accessions)} accessions to gene symbols...")
    mapping = map_to_gene_symbols(accessions)
    df["gene_symbol"] = df["protein"].map(mapping)

    # Create main output folder
    csv_name = Path(args.csv).stem
    folder_name = csv_name.split("predictions_")[-1]
    outdir = args.outdir / folder_name
    outdir.mkdir(parents=True, exist_ok=True)

    if args.by_ligand:
        ligands = df[args.ligand_column].dropna().unique()
        print(f"Found ligands: {ligands}")
        for ligand in ligands:
            ligand_df = df[df[args.ligand_column] == ligand]
            if ligand_df.empty:
                print(f"No data for ligand {ligand}, skipping.")
                continue
            ligand_outdir = outdir / f"ligand_{ligand}"
            ligand_outdir.mkdir(parents=True, exist_ok=True)
            run_enrichment_for_dataframe(ligand_df, ligand_outdir, args, ligand=ligand)
    else:
        run_enrichment_for_dataframe(df, outdir, args, ligand=None)


if __name__ == "__main__":
    main()