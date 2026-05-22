#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
from collections import defaultdict


def split_by_column(
    input_csv: Path,
    output_root: Path,
    split_column: str,
) -> None:

    output_folder = output_root / input_csv.stem
    output_folder.mkdir(parents=True, exist_ok=True)

    writers = {}
    files = {}
    counts = defaultdict(int)

    try:
        with open(input_csv, "r", newline="", encoding="utf-8-sig") as infile:
            reader = csv.DictReader(infile)

            if reader.fieldnames is None:
                raise ValueError(f"No header found in input file: {input_csv}")

            if split_column not in reader.fieldnames:
                raise ValueError(
                    f"Column '{split_column}' not found.\n"
                    f"Found columns: {reader.fieldnames}"
                )

            fieldnames = reader.fieldnames

            for row in reader:
                value = (row.get(split_column) or "").strip()

                # fallback name for empty values
                if not value:
                    value = "EMPTY"

                # make filename-safe
                safe_value = "".join(
                    c if c.isalnum() or c in ("-", "_") else "_"
                    for c in value
                )

                if safe_value not in writers:
                    out_path = output_folder / f"{safe_value}.csv"

                    outfile = open(out_path, "w", newline="", encoding="utf-8")
                    writer = csv.DictWriter(outfile, fieldnames=fieldnames)

                    writer.writeheader()

                    files[safe_value] = outfile
                    writers[safe_value] = writer

                writers[safe_value].writerow(row)
                counts[safe_value] += 1

    finally:
        for f in files.values():
            f.close()

    print(f"Input file: {input_csv}")
    print(f"Output folder: {output_folder}")
    print("\nRows written per ligand:")

    for ligand, count in sorted(counts.items()):
        print(f"  {ligand}: {count}")


def main():
    parser = argparse.ArgumentParser(
        description="Split a CSV into separate files based on a column value."
    )

    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--output_root", default=".")
    parser.add_argument("--split_column", default="ligand")

    args = parser.parse_args()

    split_by_column(
        input_csv=Path(args.input_csv),
        output_root=Path(args.output_root),
        split_column=args.split_column,
    )


if __name__ == "__main__":
    main()