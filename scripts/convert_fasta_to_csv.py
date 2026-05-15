#!/usr/bin/env python3
import csv
import re
from pathlib import Path
import argparse


def parse_uniprot_header(h: str):
    m = re.match(r"^(sp|tr)\|([^|]+)\|(\S+)\s+(.*)$", h)
    if m:
        accession = m.group(2)
        entry_name = m.group(3)
        rest = m.group(4)
    else:
        parts = h.split(maxsplit=1)
        accession = parts[0]
        entry_name = parts[0]
        rest = parts[1] if len(parts) > 1 else ""

    protein_name = rest.split(" OS=")[0].strip() if " OS=" in rest else rest.strip()

    organism = None
    mo = re.search(r"\bOS=([^=]+?)\sOX=", rest)
    if mo:
        organism = mo.group(1).strip()

    gene = None
    mg = re.search(r"\bGN=([^\s]+)", rest)
    if mg:
        gene = mg.group(1).strip()

    return accession, entry_name, protein_name, gene, organism


def fasta_to_csv(fasta_path: Path, csv_path: Path, limit: int = 0) -> int:
    rows_written = 0

    csv_path.parent.mkdir(parents=True, exist_ok=True)

    with open(csv_path, "w", encoding="utf-8", newline="") as out_f:
        writer = csv.writer(out_f)
        writer.writerow([
            "accession", "entry_name", "protein_name",
            "gene", "organism", "sequence", "seq_len"
        ])

        header = None
        seq_parts = []

        def flush():
            nonlocal header, seq_parts, rows_written
            if header is None:
                return
            sequence = "".join(seq_parts).replace(" ", "").strip()
            accession, entry_name, protein_name, gene, organism = parse_uniprot_header(header)
            writer.writerow([
                accession, entry_name, protein_name,
                gene, organism, sequence, len(sequence)
            ])
            rows_written += 1

        with open(fasta_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if header:
                        flush()
                        if limit and rows_written >= limit:
                            break
                    header = line[1:].strip()
                    seq_parts = []
                else:
                    seq_parts.append(line.strip())

            if header and not (limit and rows_written >= limit):
                flush()

    return rows_written


def main():
    parser = argparse.ArgumentParser(description="Convert proteome FASTA to CSV")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output CSV file")
    parser.add_argument("--limit", type=int, default=0, help="Limit number of proteins (0 = all)")

    args = parser.parse_args()

    fasta_path = Path(args.input)
    csv_path = Path(args.output)
    limit = args.limit

    if not fasta_path.exists():
        raise FileNotFoundError(f"Missing {fasta_path}")

    n = fasta_to_csv(fasta_path, csv_path, limit=limit)
    print(f"Wrote {n} proteins to {csv_path.resolve()}")


if __name__ == "__main__":
    main()