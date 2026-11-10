#!/usr/bin/env python3
"""
filter_checkv.py
================
Filters CheckV output retaining contigs that meet minimum completeness
and length thresholds.

Usage:
    python3 filter_checkv.py
        --quality    <CheckV quality_summary.tsv>
        --viruses    <CheckV viruses.fna>
        --proviruses <CheckV proviruses.fna>
        --out-fasta  <filtered contigs FASTA>
        --out-tsv    <filtered quality TSV>
        --min-completeness <float, default 50>
        --min-length       <int, default 5000>
"""

import argparse
import logging
import sys
import pandas as pd
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Filter CheckV output")
    parser.add_argument("--quality",           required=True)
    parser.add_argument("--viruses",           required=True)
    parser.add_argument("--proviruses",        required=True)
    parser.add_argument("--out-fasta",         required=True)
    parser.add_argument("--out-tsv",           required=True)
    parser.add_argument("--min-completeness",  type=float, default=50)
    parser.add_argument("--min-length",        type=int,   default=5000)
    return parser.parse_args()

def main():
    args = parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stderr
    )

    logging.info(f"Filtering: completeness >= {args.min_completeness}%, "
                 f"length >= {args.min_length} bp")

    # ------------------------------------------------------------------
    # 1. Load quality summary
    # ------------------------------------------------------------------
    df = pd.read_csv(args.quality, sep="\t")
    logging.info(f"  Total contigs in CheckV output: {len(df)}")

    keep_qualities = {"Complete", "High-quality", "Medium-quality"}

    def should_keep(row):
        if pd.notna(row.get("completeness")):
            if float(row["completeness"]) >= args.min_completeness:
                return True
        if row.get("checkv_quality") in keep_qualities:
            return True
        return False

    df["keep"] = df.apply(should_keep, axis=1)

    if "contig_length" in df.columns:
        df["keep"] = df["keep"] & (df["contig_length"] >= args.min_length)

    filtered_df  = df[df["keep"]].drop(columns=["keep"])
    filtered_ids = set(filtered_df["contig_id"].tolist())
    logging.info(f"  Contigs passing filter: {len(filtered_ids)}")

    # ------------------------------------------------------------------
    # 2. Write filtered quality TSV
    # ------------------------------------------------------------------
    filtered_df.to_csv(args.out_tsv, sep="\t", index=False)

    # ------------------------------------------------------------------
    # 3. Extract matching sequences from viruses.fna + proviruses.fna
    # ------------------------------------------------------------------
    written = 0
    with open(args.out_fasta, "w") as out_fasta:
        for fasta_file in [args.viruses, args.proviruses]:
            for rec in SeqIO.parse(fasta_file, "fasta"):
                base_id = rec.id.split("/")[0]
                if base_id in filtered_ids or rec.id in filtered_ids:
                    SeqIO.write(rec, out_fasta, "fasta")
                    written += 1

    logging.info(f"  Sequences written: {written}")
    logging.info("Done.")

if __name__ == "__main__":
    main()
