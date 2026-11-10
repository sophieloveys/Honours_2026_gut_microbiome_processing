#!/usr/bin/env python3
"""
pool_contigs.py
===============
Pools all CheckV-filtered contigs from all samples into a single FASTA
ready for vOTU clustering. Consolidates provenance tables into one master
provenance file.

Usage:
    python3 pool_contigs.py
        --fastas        <space-separated list of per-sample filtered FASTAs>
        --provenances   <space-separated list of per-sample provenance TSVs>
        --out-fasta     <pooled output FASTA>
        --out-provenance <master provenance TSV>
"""

import argparse
import logging
import sys
import pandas as pd
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Pool contigs across all samples")
    parser.add_argument("--fastas",          nargs="+", required=True,
                        help="Per-sample filtered_contigs.fasta files")
    parser.add_argument("--provenances",     nargs="+", required=True,
                        help="Per-sample provenance TSV files")
    parser.add_argument("--out-fasta",       required=True,
                        help="Pooled output FASTA")
    parser.add_argument("--out-provenance",  required=True,
                        help="Master provenance TSV")
    return parser.parse_args()

def main():
    args = parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stderr
    )

    # ------------------------------------------------------------------
    # 1. Pool FASTA files — deduplicate by header
    # ------------------------------------------------------------------
    seen_ids      = set()
    total_written = 0

    with open(args.out_fasta, "w") as out_fasta:
        for fasta_path in args.fastas:
            try:
                for rec in SeqIO.parse(fasta_path, "fasta"):
                    if rec.id not in seen_ids:
                        seen_ids.add(rec.id)
                        SeqIO.write(rec, out_fasta, "fasta")
                        total_written += 1
            except Exception as e:
                logging.warning(f"  Could not read {fasta_path}: {e}")

    logging.info(f"Total contigs pooled: {total_written}")

    # ------------------------------------------------------------------
    # 2. Concatenate provenance tables
    # ------------------------------------------------------------------
    prov_dfs = []
    for prov_path in args.provenances:
        try:
            prov_dfs.append(pd.read_csv(prov_path, sep="\t"))
        except Exception as e:
            logging.warning(f"  Could not read {prov_path}: {e}")

    if prov_dfs:
        master = pd.concat(prov_dfs, ignore_index=True)
        master.to_csv(args.out_provenance, sep="\t", index=False)
        logging.info(f"Provenance records: {len(master)}")
    else:
        open(args.out_provenance, "w").close()
        logging.warning("No provenance files read — empty provenance written.")

    logging.info("Done.")

if __name__ == "__main__":
    main()
