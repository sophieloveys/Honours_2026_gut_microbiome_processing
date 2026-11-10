#!/usr/bin/env python3
"""
assign_confidence_tiers.py
==========================
Assigns a confidence tier to each vOTU representative based on:
  - How many tools called it (VIBRANT, geNomad, or both)
  - How many samples it appears in (used as proxy for robustness)

Tier definitions (bulk metagenomics only — no VLP fraction):
  1 — Both tools, appears in multiple samples     (highest confidence)
  2 — Both tools, appears in single sample        (high confidence)
  3 — One tool,  appears in multiple samples      (medium-high confidence)
  4 — One tool,  appears in single sample         (medium confidence)
  5 — Unknown tool support                        (low — scrutinise)

Usage:
    python3 assign_confidence_tiers.py
        --reps        <vOTU representative FASTA>
        --cluster     <CD-HIT .clstr file>
        --provenance  <master provenance TSV>
        --out         <output confidence tiers TSV>
"""

import re
import argparse
import logging
import sys
import pandas as pd
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Assign vOTU confidence tiers")
    parser.add_argument("--reps",        required=True, help="vOTU representative FASTA")
    parser.add_argument("--cluster",     required=True, help="CD-HIT .clstr file")
    parser.add_argument("--provenance",  required=True, help="Master provenance TSV")
    parser.add_argument("--out",         required=True, help="Output confidence tiers TSV")
    return parser.parse_args()

def parse_clstr(clstr_path):
    """Parse CD-HIT .clstr file into {rep_id: [member_ids]} dict."""
    cluster_map = {}
    rep_members = {}
    current_rep = None

    with open(clstr_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                current_rep = None
                continue
            match = re.search(r">([^\s.]+)", line)
            if match:
                contig_id = match.group(1).rstrip(".")
                if line.endswith("*"):
                    current_rep = contig_id
                    rep_members[current_rep] = []
                if current_rep:
                    cluster_map[contig_id] = current_rep
                    rep_members[current_rep].append(contig_id)

    return rep_members

def main():
    args = parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stderr
    )

    # ------------------------------------------------------------------
    # 1. Load provenance table
    # ------------------------------------------------------------------
    prov = pd.read_csv(args.provenance, sep="\t")
    prov_dict = prov.set_index("contig_id").to_dict(orient="index")
    logging.info(f"Provenance records: {len(prov_dict)}")

    # ------------------------------------------------------------------
    # 2. Parse CD-HIT cluster file
    # ------------------------------------------------------------------
    rep_members = parse_clstr(args.cluster)
    logging.info(f"vOTU clusters: {len(rep_members)}")

    # ------------------------------------------------------------------
    # 3. Assign tiers
    # ------------------------------------------------------------------
    rows = []

    for votu_idx, (rep_id, members) in enumerate(rep_members.items(), 1):
        called_vibrant = False
        called_genomad = False
        samples_seen   = set()

        for member_id in members:
            p = prov_dict.get(member_id, {})
            if p.get("called_vibrant"):
                called_vibrant = True
            if p.get("called_genomad"):
                called_genomad = True
            if p.get("sample"):
                samples_seen.add(p["sample"])

        n_samples = len(samples_seen)

        if called_vibrant and called_genomad:
            tool_support = "both"
        elif called_vibrant:
            tool_support = "VIBRANT_only"
        elif called_genomad:
            tool_support = "geNomad_only"
        else:
            tool_support = "unknown"

        multi_sample = n_samples > 1

        if tool_support == "both" and multi_sample:
            tier = 1
            note = "Highest confidence — both tools, multiple samples"
        elif tool_support == "both" and not multi_sample:
            tier = 2
            note = "High confidence — both tools, single sample"
        elif tool_support != "both" and tool_support != "unknown" and multi_sample:
            tier = 3
            note = "Medium-high confidence — one tool, multiple samples"
        elif tool_support != "both" and tool_support != "unknown" and not multi_sample:
            tier = 4
            note = "Medium confidence — one tool, single sample"
        else:
            tier = 5
            note = "Low confidence — scrutinise before interpreting"

        rows.append({
            "votu_id":               f"vOTU_{votu_idx:05d}",
            "representative_contig": rep_id,
            "tier":                  tier,
            "tool_support":          tool_support,
            "n_samples":             n_samples,
            "samples":               ";".join(sorted(samples_seen)),
            "n_member_contigs":      len(members),
            "note":                  note,
        })

    tiers_df = pd.DataFrame(rows).sort_values("tier")
    tiers_df.to_csv(args.out, sep="\t", index=False)

    for t in range(1, 6):
        count = (tiers_df["tier"] == t).sum()
        logging.info(f"  Tier {t}: {count} vOTUs")

    logging.info(f"Confidence tiers written: {args.out}")
    logging.info("Done.")

if __name__ == "__main__":
    main()
