#!/usr/bin/env python3
"""Validate a .snpdiffs file against the CSP2 schema (parseSNPDiffs compatibility).

Usage: python validate_snpdiffs.py <snpdiffs_file>
"""
import sys

EXPECTED_HEADER_KEYS = [
    "Query_ID", "Reference_ID", "SNPs", "Indels",
    "Reference_Percent_Aligned", "Query_Percent_Aligned",
    "Median_Percent_Identity", "Kmer_Similarity",
    "Shared_Kmers", "Reference_Unique_Kmers", "Query_Unique_Kmers",
    "gSNPs", "gIndels",
]

EXPECTED_BED_COLS = [
    "Ref_Contig", "Ref_Start", "Ref_End", "Ref_Length", "Ref_Aligned",
    "Query_Contig", "Query_Start", "Query_End", "Query_Length", "Query_Aligned", "Perc_Iden",
]

EXPECTED_SNP_COLS = [
    "Ref_Contig", "Start_Ref", "Ref_Pos",
    "Query_Contig", "Start_Query", "Query_Pos",
    "Ref_Loc", "Query_Loc",
    "Ref_Start", "Ref_End",
    "Query_Start", "Query_End",
    "Ref_Base", "Query_Base",
    "Dist_to_Ref_End", "Dist_to_Query_End",
    "Ref_Aligned", "Query_Aligned",
    "Query_Direction", "Perc_Iden", "Cat",
]


def main():
    f = sys.argv[1]
    with open(f) as fh:
        lines = fh.readlines()

    errors = []
    warnings = []

    if not lines:
        errors.append("File is empty")
        print("FAIL: " + "; ".join(errors))
        sys.exit(1)

    # Parse header
    header_line = lines[0].strip()
    if not header_line.startswith("#"):
        errors.append(f"Line 1 does not start with #: {header_line[:50]}")
    else:
        fields = header_line.lstrip("#").strip().split("\t")
        header_keys = set()
        for field in fields:
            if ":" in field:
                key = field.split(":", 1)[0]
                header_keys.add(key)
        for key in EXPECTED_HEADER_KEYS:
            if key not in header_keys:
                errors.append(f"Missing header key: {key}")

    # Parse BED rows (##)
    bed_count = 0
    for line in lines[1:]:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("##"):
            bed_cols = stripped.lstrip("#").strip().split("\t")
            if len(bed_cols) != len(EXPECTED_BED_COLS):
                errors.append(f"BED row has {len(bed_cols)} columns, expected {len(EXPECTED_BED_COLS)}")
            bed_count += 1

    # Parse SNP/Indel rows
    snp_count = 0
    indel_count = 0
    for line in lines[1:]:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            continue
        cols = stripped.split("\t")
        if len(cols) != len(EXPECTED_SNP_COLS):
            errors.append(f"SNP row has {len(cols)} columns, expected {len(EXPECTED_SNP_COLS)}: {stripped[:50]}")
            continue
        cat = cols[-1]
        # Cat must be SNP or Indel
        if cat not in ("SNP", "Indel"):
            errors.append(f"Invalid Cat value: {cat}")
        # Dist_to_Ref_End and Dist_to_Query_End must be numeric
        try:
            int(cols[-5])  # Dist_to_Ref_End
            int(cols[-4])  # Dist_to_Query_End
        except ValueError:
            errors.append(f"Non-numeric Dist_to_Ref_End/Query_End in row: {stripped[:50]}")
        # Perc_Iden must be numeric
        try:
            float(cols[-2])  # Perc_Iden
        except ValueError:
            errors.append(f"Non-numeric Perc_Iden in row: {stripped[:50]}")
        if cat == "SNP":
            snp_count += 1
        elif cat == "Indel":
            indel_count += 1

    if errors:
        print("FAIL: " + "; ".join(errors[:10]))
        sys.exit(1)
    else:
        print(f"PASS: {len(header_line)} header, {bed_count} BED rows, {snp_count} SNPs, {indel_count} Indels")
        sys.exit(0)


if __name__ == "__main__":
    main()
