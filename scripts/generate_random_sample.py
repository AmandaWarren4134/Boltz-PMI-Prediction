import csv
import random
from pathlib import Path

import argparse


def read_csv_rows(path: Path):
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames
    return rows, fieldnames


def main():

    parser = argparse.ArgumentParser(description="Generate random subset of proteome CSV")
    parser.add_argument("--proteome_csv", required=True, help="Input proteome CSV file path")
    parser.add_argument("--output", required=True, help="Output CSV file path")
    parser.add_argument("--n_random", type=int, default=2000, help="Number of random rows to sample")
    parser.add_argument("--random_seed", type=int, default=42, help="Random seed for reproducibility")

    args = parser.parse_args()

    random.seed(args.random_seed)       

    proteome_rows, fieldnames = read_csv_rows(args.proteome_csv)

    if len(proteome_rows) < args.n_random:
        raise ValueError(
            f"Not enough rows to sample {args.n_random}. "
            f"Only found {len(proteome_rows)} rows."
        )

    sampled_rows = random.sample(proteome_rows, args.n_random)

    # Write output
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sampled_rows)

    print(f"Wrote {len(sampled_rows)} random rows to {output_path}")


if __name__ == "__main__":
    main()