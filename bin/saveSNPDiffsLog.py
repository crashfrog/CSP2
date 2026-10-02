#!/usr/bin/env python3
"""saveSNPDiffsLog.py — merge per-contig Phraya snpdiffs outputs into one .snpdiffs
file and write a per-comparison log.

Replaces the logging portion of compileMUMmer.py. In Phase A, runMUMmer calls
phraya plan + phraya align (reference-palette mode), which writes one .phraya
file per reference contig. Each .phraya is filtered to .snpdiffs separately
(phraya filter --format snpdiffs is single-input — it processes inputs[0] only
when given multiple .phraya files, see docs/phraya-snpdiffs-multi-contig-blocker.md),
and this script merges the per-contig results.

CRITICAL: a homology/coverage gate filters out per-contig alignments where
Ref_Aligned / Ref_Length < MIN_COVERAGE (default 0.01). Without this gate,
non-homologous contig pairs in palette mode produce spurious SNPs (validated:
mummer = 33,390 SNPs vs phraya-no-gate = 173,840 SNPs on SRR30874442 vs
SRR30874443).

Usage:
    python saveSNPDiffsLog.py --snpdiffs_dir per_contig/ \
        --out pair.snpdiffs --log_file report.log \
        --query QUERY_NAME --reference REF_NAME
"""

import argparse
import glob
import os
import sys
from datetime import datetime

# Minimum fraction of reference contig covered to retain per-contig variants.
# Below this, the alignment is treated as spurious (cross-contig noise in
# palette mode). 0.01 = 1% coverage floor.
MIN_COVERAGE = 0.01


def parse_snpdiffs_file(snpdiffs_path):
    """Parse a .snpdiffs file into (header_line, bed_rows, snp_rows).

    Returns:
        header_line: str (the # line, without leading #)
        bed_rows: list of str (##-prefixed lines, stripped of #)
        snp_rows: list of str (variant data lines)
    """
    with open(snpdiffs_path, "r") as f:
        lines = f.readlines()

    header_line = ""
    bed_rows = []
    snp_rows = []

    for line in lines:
        stripped = line.rstrip("\n")
        if not stripped:
            continue
        if stripped.startswith("##"):
            bed_rows.append(stripped.lstrip("#").strip())
        elif stripped.startswith("#"):
            header_line = stripped.lstrip("#").strip()
        else:
            snp_rows.append(stripped)

    return header_line, bed_rows, snp_rows


def parse_header_fields(header_line):
    """Parse tab-separated key:value fields from a header line into a dict."""
    fields = {}
    parts = header_line.split("\t")
    for part in parts:
        if ":" in part:
            key, val = part.split(":", 1)
            fields[key] = val
    return fields


def passes_coverage_gate(bed_row):
    """Check if a single-contig alignment passes the coverage threshold.

    BED row columns (after stripping ##):
      Ref_Contig  Ref_Start  Ref_End  Ref_Length  Ref_Aligned
      Query_Contig  Query_Start  Query_End  Query_Length  Query_Aligned  Perc_Iden
    """
    parts = bed_row.split("\t")
    if len(parts) < 11:
        return True  # can't evaluate — don't drop

    ref_length = int(parts[3]) if parts[3].isdigit() else 0   # Ref_Length
    ref_aligned = int(parts[4]) if parts[4].isdigit() else 0  # Ref_Aligned

    if ref_length == 0:
        return True  # can't evaluate — don't drop
    return (ref_aligned / ref_length) >= MIN_COVERAGE


def merge_snpdiffs(snpdiffs_dir, out_path, query, reference):
    """Merge all .snpdiffs files in snpdiffs_dir into a single output file.

    Strategy:
    - Parse each per-contig .snpdiffs file.
    - Apply coverage gate: skip files where Ref_Aligned / Ref_Length < MIN_COVERAGE
      (filters spurious cross-contig alignments from palette mode).
    - Use the first passing file's header as the base (Query/Reference metadata is
      shared across all files from the same FASTA inputs).
    - Concatenate BED rows and SNP/Indel rows from passing files.
    - Recompute summary metrics: SNP/Indel counts aggregated; alignment stats
      aggregated via BED-row-derived per-contig coverage weighted by contig length.
    """
    files = sorted(glob.glob(os.path.join(snpdiffs_dir, "*.snpdiffs")))
    if not files:
        print(f"Error: no .snpdiffs files found in {snpdiffs_dir}", file=sys.stderr)
        sys.exit(1)

    # Parse all files, apply coverage gate, collect passing data
    all_bed = []
    all_snp_rows = []
    all_headers = []
    skipped = []

    for f in files:
        header_line, bed_rows, snp_rows = parse_snpdiffs_file(f)
        passed = False
        for bed in bed_rows:
            if passes_coverage_gate(bed):
                all_bed.append(bed)
                passed = True
            else:
                skipped.append(f"{os.path.basename(f)}: coverage={bed.split(chr(9))[4]}/{bed.split(chr(9))[3]}")
        if passed:
            all_snp_rows.extend(snp_rows)
            all_headers.append(parse_header_fields(header_line))

    if not all_headers:
        print(f"Error: all {len(files)} per-contig files failed the coverage gate", file=sys.stderr)
        sys.exit(1)

    header_fields = all_headers[0].copy()  # base header from first passing file

    # Sort SNP rows by Ref_Pos (1-indexed, column 3) then Cat
    all_snp_rows.sort(key=lambda r: (
        int(r.split("\t")[2]) if len(r.split("\t")) > 2 and r.split("\t")[2].isdigit() else 0,
        r.split("\t")[-1] if r.split("\t") else "",
    ))

    # Aggregate alignment statistics across all passing reference contigs.
    # Per-contig Ref_Length / Ref_Aligned come from BED rows (per reference space).
    # Query_Length is constant across all files (same query FASTA).
    total_ref_len = 0
    total_ref_aligned = 0
    total_query_aligned = 0
    weighted_identity = 0.0
    total_weight_for_identity = 0
    query_length = None  # constant across all per-contig files

    for bed in all_bed:
        parts = bed.split("\t")
        if len(parts) < 11:
            continue
        ref_length = int(parts[3]) if parts[3].isdigit() else 0    # Ref_Length
        ref_aligned = int(parts[4]) if parts[4].isdigit() else 0   # Ref_Aligned
        q_aligned = int(parts[9]) if parts[9].isdigit() else 0     # Query_Aligned
        if query_length is None:
            query_length = int(parts[8]) if parts[8].isdigit() else 0  # Query_Length
        bed_perc_iden = 0.0
        try:
            bed_perc_iden = float(parts[10])
        except (ValueError, TypeError):
            pass

        total_ref_len += ref_length
        total_ref_aligned += ref_aligned
        total_query_aligned += q_aligned
        weighted_identity += ref_aligned * bed_perc_iden
        total_weight_for_identity += ref_aligned

    # Reference_Percent_Aligned: covered ref bases / total ref contig lengths
    if total_ref_len > 0:
        header_fields["Reference_Percent_Aligned"] = f"{total_ref_aligned / total_ref_len * 100:.2f}"
    # Query_Percent_Aligned: covered query bases / query length (constant).
    # Clamp to 100 (query may have multiple contigs each contributing aligned positions).
    if query_length is not None and query_length > 0:
        q_pct = total_query_aligned / query_length * 100.0
        header_fields["Query_Percent_Aligned"] = f"{min(q_pct, 100.0):.2f}"
    if total_weight_for_identity > 0:
        header_fields["Median_Percent_Identity"] = f"{weighted_identity / total_weight_for_identity:.2f}"

    # Recompute variant counts from merged rows
    snp_count = sum(1 for r in all_snp_rows if r.endswith("SNP"))
    indel_count = sum(1 for r in all_snp_rows if r.endswith("Indel"))
    header_fields["SNPs"] = str(snp_count)
    header_fields["Indels"] = str(indel_count)

    # Rebuild header line, preserving original field order from the first file.
    first_header = parse_snpdiffs_file(files[0])[0]  # raw header (no leading #)
    parts = first_header.split("\t")
    updated_parts = []
    for part in parts:
        if ":" in part:
            key, _ = part.split(":", 1)
            if key in header_fields:
                updated_parts.append(f"{key}:{header_fields[key]}")
            else:
                updated_parts.append(part)
        else:
            updated_parts.append(part)

    # Write merged output
    with open(out_path, "w") as f:
        f.write("#\t" + "\t".join(updated_parts) + "\n")
        for bed in all_bed:
            f.write("##\t" + bed + "\n")
        for row in all_snp_rows:
            f.write(row + "\n")

    return snp_count, indel_count, header_fields, skipped


def main():
    parser = argparse.ArgumentParser(
        description='Merge per-contig Phraya snpdiffs outputs into one .snpdiffs file + log.'
    )
    parser.add_argument("--snpdiffs_dir", required=True,
                        help="Directory of per-contig .snpdiffs files")
    parser.add_argument("--out", required=True, help="Output merged .snpdiffs file path")
    parser.add_argument("--log_file", required=True, help="Log file path to write")
    parser.add_argument("--query", required=True, help="Query isolate name")
    parser.add_argument("--reference", required=True, help="Reference isolate name")
    args = parser.parse_args()

    snp_count, indel_count, fields, skipped = merge_snpdiffs(
        args.snpdiffs_dir, args.out, args.query, args.reference
    )

    # Write log
    os.makedirs(os.path.dirname(os.path.abspath(args.log_file)), exist_ok=True)
    with open(args.log_file, "w") as log:
        log.write(f"Phraya → .snpdiffs merge log: {args.query}__vs__{args.reference}\n")
        log.write("-------------------------------------------------------\n\n")
        log.write(f"  snpdiffs dir: {args.snpdiffs_dir}\n")
        log.write(f"  merged output: {os.path.abspath(args.out)}\n")
        log.write(f"  Query: {args.query}\n")
        log.write(f"  Reference: {args.reference}\n")
        log.write(f"  Coverage gate: MIN_COVERAGE={MIN_COVERAGE} (skipped {len(skipped)} per-contig files)\n")
        log.write(f"  SNPs: {snp_count}\n")
        log.write(f"  Indels: {indel_count}\n")
        log.write("\n--- Header Metadata ---\n")
        for key in ["Query_ID", "Query_Assembly", "Query_Contig_Count", "Query_Assembly_Bases",
                     "Query_N50", "Reference_Assembly", "Reference_Contig_Count",
                     "Reference_Assembly_Bases", "Reference_N50",
                     "Reference_Percent_Aligned", "Query_Percent_Aligned",
                     "Median_Percent_Identity", "Kmer_Similarity"]:
            log.write(f"  {key}: {fields.get(key, 'NA')}\n")
        log.write("\n")
        log.write(f"Log generated: {datetime.now().isoformat()}\n")

    # Print the canonical output line for Nextflow channel compatibility
    print(f"{args.query},{args.reference},{os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
