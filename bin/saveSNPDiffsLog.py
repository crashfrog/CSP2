#!/usr/bin/env python3
"""saveSNPDiffsLog.py — write the per-comparison log file for a Phraya-produced .snpdiffs.

This replaces the logging portion of compileMUMmer.py (which parsed mummer output
to extract identity/length stats). Now that phraya filter --format snpdiffs emits
the full `#\t<query metadata>\t<reference metadata>\t<summary>` header directly,
we just read those key:value pairs and write a human-readable log.

Usage:
    python saveSNPDiffsLog.py --out pair.snpdiffs --log_file report.log \
        --query QUERY_NAME --reference REF_NAME
"""

import argparse
import os
import sys
from datetime import datetime


def parse_snpdiffs_header(snpdiffs_path):
    """Parse the first line of a .snpdiffs file (starts with '#').

    Returns a dict of key:value pairs from all three sections (query, reference,
    summary), with the section prefix stripped so keys are unique within their
    namespace (e.g. 'Query_ID', 'Reference_ID', 'SNPs', etc.).
    """
    with open(snpdiffs_path, "r") as f:
        first_line = f.readline().rstrip("\n")

    if not first_line.startswith("#"):
        raise ValueError(f"Expected '#' header line in {snpdiffs_path}, got: {first_line[:50]}")

    parts = first_line.lstrip("#").strip().split("\t")
    fields = {}
    for part in parts:
        if ":" in part:
            key, val = part.split(":", 1)
            fields[key] = val
    return fields


def count_variants(snpdiffs_path):
    """Count SNP and Indel categories from the data rows (non-comment, non-BED lines)."""
    snps = 0
    indels = 0
    with open(snpdiffs_path, "r") as f:
        for line in f:
            stripped = line.rstrip("\n")
            if not stripped or stripped.startswith("#"):
                continue
            if stripped.startswith("##"):
                continue  # BED row
            fields = stripped.split("\t")
            cat = fields[-1]  # Cat is the last column
            if cat == "SNP":
                snps += 1
            elif cat == "Indel":
                indels += 1
    return snps, indels


def main():
    parser = argparse.ArgumentParser(
        description="Write per-comparison log for a Phraya-produced .snpdiffs file."
    )
    parser.add_argument("--out", required=True, help="Path to the .snpdiffs output file")
    parser.add_argument("--log_file", required=True, help="Path to the log file to write")
    parser.add_argument("--query", required=True, help="Query isolate name")
    parser.add_argument("--reference", required=True, help="Reference isolate name")
    args = parser.parse_args()

    if not os.path.exists(args.out):
        print(f"Error: .snpdiffs file does not exist: {args.out}", file=sys.stderr)
        sys.exit(1)

    # Parse the header
    try:
        fields = parse_snpdiffs_header(args.out)
    except Exception as e:
        print(f"Error parsing .snpdiffs header: {e}", file=sys.stderr)
        sys.exit(1)

    # Count variants from data rows
    snp_count, indel_count = count_variants(args.out)

    # Write log
    os.makedirs(os.path.dirname(os.path.abspath(args.log_file)), exist_ok=True)
    with open(args.log_file, "w") as log:
        log.write(f"CSP2 MUMmer Compilation Log: {args.query}__vs__{args.reference}\n")
        log.write("-------------------------------------------------------\n\n")

        log.write("Phraya alignment → .snpdiffs conversion\n")
        log.write(f"  snpdiffs File: {os.path.abspath(args.out)}\n")
        log.write(f"  Query: {args.query}\n")
        log.write(f"  Reference: {args.reference}\n")

        # Query metadata
        log.write(f"  Query_Assembly: {fields.get('Query_Assembly', 'NA')}\n")
        log.write(f"  Query_Contig_Count: {fields.get('Query_Contig_Count', 'NA')}\n")
        log.write(f"  Query_Assembly_Bases: {fields.get('Query_Assembly_Bases', 'NA')}\n")
        log.write(f"  Query_N50: {fields.get('Query_N50', 'NA')}\n")
        log.write(f"  Query_N90: {fields.get('Query_N90', 'NA')}\n")
        log.write(f"  Query_L50: {fields.get('Query_L50', 'NA')}\n")
        log.write(f"  Query_L90: {fields.get('Query_L90', 'NA')}\n")

        # Reference metadata
        log.write(f"  Reference_Assembly: {fields.get('Reference_Assembly', 'NA')}\n")
        log.write(f"  Reference_Contig_Count: {fields.get('Reference_Contig_Count', 'NA')}\n")
        log.write(f"  Reference_Assembly_Bases: {fields.get('Reference_Assembly_Bases', 'NA')}\n")
        log.write(f"  Reference_N50: {fields.get('Reference_N50', 'NA')}\n")
        log.write(f"  Reference_N90: {fields.get('Reference_N90', 'NA')}\n")
        log.write(f"  Reference_L50: {fields.get('Reference_L50', 'NA')}\n")
        log.write(f"  Reference_L90: {fields.get('Reference_L90', 'NA')}\n")

        # Summary metrics
        log.write("\n--- Summary Metrics ---\n")
        log.write(f"  SNPs: {snp_count}\n")
        log.write(f"  Indels: {indel_count}\n")
        log.write(f"  Reference_Percent_Aligned: {fields.get('Reference_Percent_Aligned', 'NA')}\n")
        log.write(f"  Query_Percent_Aligned: {fields.get('Query_Percent_Aligned', 'NA')}\n")
        log.write(f"  Median_Percent_Identity: {fields.get('Median_Percent_Identity', 'NA')}\n")
        log.write(f"  Median_Alignment_Length: {fields.get('Median_Alignment_Length', 'NA')}\n")
        log.write(f"  Kmer_Similarity: {fields.get('Kmer_Similarity', 'NA')}\n")
        log.write(f"  Shared_Kmers: {fields.get('Shared_Kmers', 'NA')}\n")
        log.write(f"  Reference_Unique_Kmers: {fields.get('Reference_Unique_Kmers', 'NA')}\n")
        log.write(f"  Query_Unique_Kmers: {fields.get('Query_Unique_Kmers', 'NA')}\n")
        log.write(f"  gSNPs: {fields.get('gSNPs', 'NA')}\n")
        log.write(f"  gIndels: {fields.get('gIndels', 'NA')}\n")
        log.write("\n")
        log.write(f"Log generated: {datetime.now().isoformat()}\n")

    # Print the canonical output line for Nextflow channel compatibility
    print(f"{args.query},{args.reference},{os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
