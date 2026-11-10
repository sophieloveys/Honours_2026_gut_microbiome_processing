#!/usr/bin/env python3
"""
union_viral_contigs.py
======================
Takes the union of VIBRANT and geNomad viral contig calls for a single
sample. Records which tool(s) called each contig for downstream confidence
tier assignment.

Usage:
    python3 union_viral_contigs.py
        --vibrant   <VIBRANT phage FASTA>
        --genomad   <geNomad virus summary TSV>
        --assembly  <original assembly FASTA>
        --sample    <sample name e.g. L1.4>
        --out-fasta <output union FASTA>
        --out-provenance <output provenance TSV>
"""

import argparse
import logging
import sys
import pandas as pd
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Union of VIBRANT and geNomad viral calls")
    parser.add_argument("--vibrant",         required=True,  help="VIBRANT phage FASTA")
    parser.add_argument("--genomad",         required=True,  help="geNomad virus summary TSV")
    parser.add_argument("--assembly",        required=True,  help="Original assembly FASTA")
    parser.add_argument("--sample",          required=True,  help="Sample name e.g. L1.4")
    parser.add_argument("--out-fasta",       required=True,  help="Output union FASTA")
    parser.add_argument("--out-provenance",  required=True,  help="Output provenance TSV")
    return parser.parse_args()

def main():
    args = parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stderr
    )

    sample = args.sample
    logging.info(f"Processing sample: {sample}")

    # ------------------------------------------------------------------
    # 1. Load assembly contigs
    # ------------------------------------------------------------------
    assembly = {rec.id: rec for rec in SeqIO.parse(args.assembly, "fasta")}
    logging.info(f"  Assembly contigs loaded: {len(assembly)}")

    # ------------------------------------------------------------------
    # 2. VIBRANT calls — contig IDs from FASTA headers
    # ------------------------------------------------------------------
    vibrant_ids = set()
    try:
        for rec in SeqIO.parse(args.vibrant, "fasta"):
            base_id = rec.id.split("_fragment_")[0]
            vibrant_ids.add(base_id)
        logging.info(f"  VIBRANT calls: {len(vibrant_ids)}")
    except Exception as e:
        logging.warning(f"  Could not parse VIBRANT output: {e}")

    # ------------------------------------------------------------------
    # 3. geNomad calls — parse virus summary TSV
    # ------------------------------------------------------------------
    genomad_ids = set()
    try:
        gdf = pd.read_csv(args.genomad, sep="\t")
        if "seq_name" in gdf.columns:
            genomad_ids = set(
                gdf["seq_name"].str.split("|").str[0].tolist()
            )
        logging.info(f"  geNomad calls: {len(genomad_ids)}")
    except Exception as e:
        logging.warning(f"  Could not parse geNomad output: {e}")

    # ------------------------------------------------------------------
    # 4. Union
    # ------------------------------------------------------------------
    all_viral_ids = vibrant_ids | genomad_ids
    logging.info(f"  Union total: {len(all_viral_ids)}")

    # ------------------------------------------------------------------
    # 5. Write union FASTA and provenance TSV
    # ------------------------------------------------------------------
    provenance_rows = []
    written = 0

    with open(args.out_fasta, "w") as fasta_out:
        for contig_id in sorted(all_viral_ids):
            if contig_id not in assembly:
                logging.warning(f"  Contig {contig_id} not in assembly — skipping")
                continue

            rec = assembly[contig_id]
            new_id = f"{sample}__{contig_id}"
            rec.id          = new_id
            rec.description = ""
            SeqIO.write(rec, fasta_out, "fasta")
            written += 1

            provenance_rows.append({
                "contig_id":      new_id,
                "original_id":    contig_id,
                "sample":         sample,
                "called_vibrant": contig_id in vibrant_ids,
                "called_genomad": contig_id in genomad_ids,
            })

    logging.info(f"  Written to union FASTA: {written}")

    prov_df = pd.DataFrame(provenance_rows)
    prov_df.to_csv(args.out_provenance, sep="\t", index=False)
    logging.info(f"  Provenance written: {args.out_provenance}")
    logging.info("Done.")

if __name__ == "__main__":
    main()
