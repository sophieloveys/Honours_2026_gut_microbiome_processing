#!/usr/bin/env python3
"""
combine_abundance.py
====================
Merges per-sample CoverM coverage TSV files into a single vOTU abundance
matrix suitable for R analysis.

Output format:
  Rows    = vOTU representative contigs
  Columns = sample names (e.g. L1.1, L1.2, ... G2.1 ... D3.10)
  Values  = trimmed mean coverage

Usage:
    python3 combine_abundance.py
        --covs  <space-separated list of per-sample CoverM TSV files>
        --out   <output combined abundance matrix TSV>
"""

import os
import argparse
import logging
import sys
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Combine CoverM coverage files")
    parser.add_argument("--covs", nargs="+", required=True,
                        help="Per-sample CoverM coverage TSV files")
    parser.add_argument("--out",  required=True,
                        help="Output abundance matrix TSV")
    return parser.parse_args()

def main():
    args = parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stderr
    )

    dfs = {}

    for cov_path in args.covs:
        # Infer sample label from filename
        # Expected: {sample}_coverage.tsv  e.g. L1.4_coverage.tsv
        basename = os.path.basename(cov_path)
        label    = basename.replace("_coverage.tsv", "")

        try:
            df = pd.read_csv(cov_path, sep="\t")
        except Exception as e:
            logging.warning(f"Could not read {cov_path}: {e}")
            continue

        # CoverM output: first column is contig name,
        # coverage column contains "Trimmed Mean" or "Mean"
        coverage_cols = [c for c in df.columns
                         if "Trimmed Mean" in c or "Mean" in c]

        if not coverage_cols:
            logging.warning(f"No coverage column in {cov_path} — skipping")
            continue

        df = df.rename(columns={
            df.columns[0]:    "contig",
            coverage_cols[0]: label
        })[["contig", label]]

        dfs[label] = df.set_index("contig")
        logging.info(f"  Loaded: {label} ({len(df)} contigs)")

    if not dfs:
        logging.error("No coverage files could be read.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Merge all samples — fill missing with 0 (contig absent = 0 coverage)
    # ------------------------------------------------------------------
    combined = pd.concat(dfs.values(), axis=1).fillna(0)
    combined.index.name = "contig"

    # Sort columns by individual then timepoint for readability
    def sort_key(col):
        # e.g. L1.4 -> ("L1", 4), D3.10 -> ("D3", 10)
        parts = col.split(".")
        if len(parts) == 2:
            return (parts[0], int(parts[1]))
        return (col, 0)

    sorted_cols = sorted(combined.columns, key=sort_key)
    combined    = combined[sorted_cols]

    combined.to_csv(args.out, sep="\t")

    logging.info(f"Abundance matrix: {combined.shape[0]} vOTUs x "
                 f"{combined.shape[1]} samples")
    logging.info(f"Written to: {args.out}")
    logging.info("Done.")

if __name__ == "__main__":
    main()
