#!/usr/bin/env python3
"""Restore gene names in a results table using the original protein input."""

import argparse
from pathlib import Path

import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    # Read .tsv/.txt files as tab-separated; treat other extensions as CSV.
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Add the gene column from an original protein input file."
    )
    parser.add_argument("results", type=Path, help="Results CSV to update")
    parser.add_argument(
        "protein_input",
        type=Path,
        help="Original protein input CSV/TSV containing accession and gene columns",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output CSV (default: <results name>_with_gene.csv)",
    )
    args = parser.parse_args()

    results = read_table(args.results)
    proteins = read_table(args.protein_input)

    if "protein" not in results.columns:
        parser.error("results is missing required column: protein")
    if "accession" not in proteins.columns:
        parser.error("protein input is missing required column: accession")
    if "gene" not in proteins.columns:
        parser.error("protein input is missing required column: gene")

    # Each protein should map to one gene. Ignore repeated identical rows, but
    # stop if the source file assigns different genes to the same protein.
    mapping = proteins[["accession", "gene"]].drop_duplicates()
    conflicts = mapping.groupby("accession", dropna=False)["gene"].nunique(dropna=False)
    conflicts = conflicts[conflicts > 1]
    if not conflicts.empty:
        parser.error(
            "protein input maps some proteins to multiple genes; examples: "
            + ", ".join(map(str, conflicts.index[:5]))
        )
    mapping = mapping.drop_duplicates(subset="accession").rename(
        columns={"accession": "protein"}
    )

    # Remove a stale gene column, if present, then restore it from the source.
    results = results.drop(columns=["gene"], errors="ignore")
    output = results.merge(mapping, on="protein", how="left", validate="many_to_one")

    output_path = args.output or args.results.with_name(
        f"{args.results.stem}_with_gene.csv"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, index=False)

    unmatched = output["gene"].isna().sum()
    print(f"Saved: {output_path}")
    if unmatched:
        print(f"Warning: {unmatched} result row(s) had no matching protein gene.")


if __name__ == "__main__":
    main()
