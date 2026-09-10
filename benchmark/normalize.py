#!/usr/bin/env python3
"""
Canonicalize a CSP2 output directory into a deterministic equivalence manifest.

Why: CSP2's outputs carry run-to-run noise that is NOT scientifically meaningful --
absolute paths, `list(set(...))` row ordering, float formatting, timestamps. To assert
that a Nextflow/Python refactor "remains equivalent" we must compare the *meaningful*
content, not the bytes. This module reduces each scientific output to a canonical form
and hashes it, producing a manifest that is stable across equivalent runs but changes
the moment a SNP call, distance, alignment column, or QC verdict changes.

It is deliberately conservative: anything it does not recognize as noise is preserved.
Timing artifacts, logs, raw MUMmer deltas, mash sketches, and stochastic trees are
excluded (compared separately / not at all).

Usage:
    python normalize.py <csp2_output_dir> -o manifest.json
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import re
import sys
from pathlib import Path

FLOAT_DECIMALS = 6  # rounding headroom: real SNP/distance changes are >> 1e-6

# Path segments / globs whose contents are noise or non-deterministic -> excluded.
EXCLUDE_DIR_PARTS = {"work", ".nextflow", "MUMmer_Output", "sketch_dir", "MUMmer_Logs",
                     "Screening_Logs", "SNP_Logs"}
EXCLUDE_SUFFIXES = (".log", ".msh", ".html", ".delta", ".mdelta", ".1delta",
                    ".coords", ".mcoords", ".1coords", ".snps", ".rdiff", ".qdiff",
                    ".report", ".treefile", ".tre", ".iqtree", ".bionj", ".mldist",
                    ".ckp.gz", ".contree", ".nex", ".gz")
EXCLUDE_NAMES = {"CSP2_Params.txt", "trace.txt", "report.html", "timeline.html",
                 "All_SNPDiffs.txt", "SNP_Dirs.txt", "SNPDiffs.txt",  # intermediate path lists
                 "Mash_Sketches.txt", "Mash_Triangle"}

# Recognize absolute-ish filesystem paths embedded in fields -> reduce to basename so
# that run location / temp dirs / work paths don't register as differences.
PATH_RE = re.compile(r"(/[^\s,;:]+/)+([^\s,;:/]+)")


def _basename_paths(text: str) -> str:
    """Replace embedded absolute paths with their basename."""
    return PATH_RE.sub(lambda m: m.group(2), text)


def _round_token(tok: str) -> str:
    """If a token is a float literal, round it; otherwise return unchanged."""
    try:
        f = float(tok)
    except (ValueError, TypeError):
        return tok
    if tok.strip().lstrip("-").isdigit():   # keep pure ints exact
        return tok
    if f != f:                              # NaN
        return "NaN"
    return f"{round(f, FLOAT_DECIMALS):.{FLOAT_DECIMALS}f}"


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


# ---- per-file-class canonicalizers -----------------------------------------

def canon_snpdiffs(raw: str) -> str:
    """Header: normalize path-valued tokens to basename. BED + SNP rows: sort each
    block. SHA256 tokens (content hashes) are preserved verbatim."""
    header, bed, snp = [], [], []
    for line in raw.splitlines():
        if line.startswith("#\t") or (line.startswith("#") and not line.startswith("##")):
            toks = line.lstrip("#").strip().split("\t")
            header = [_basename_paths(t) for t in toks]
        elif line.startswith("##\t"):
            bed.append(_basename_paths(line))
        elif line.strip():
            snp.append(line)
    out = ["#\t" + "\t".join(header)]
    out += sorted(bed)
    out += sorted(snp)
    return "\n".join(out)


def canon_fasta(raw: str) -> str:
    """Alignment FASTA -> records sorted by id as 'id\\tseq'. Sequence (column) order is
    preserved (it is the biological content); only record order is canonicalized."""
    recs, cur_id, cur_seq = {}, None, []
    for line in raw.splitlines():
        if line.startswith(">"):
            if cur_id is not None:
                recs[cur_id] = "".join(cur_seq)
            cur_id, cur_seq = line[1:].strip(), []
        else:
            cur_seq.append(line.strip())
    if cur_id is not None:
        recs[cur_id] = "".join(cur_seq)
    return "\n".join(f"{i}\t{recs[i]}" for i in sorted(recs))


def canon_matrix(raw: str) -> str:
    """Symmetric distance matrix -> sorted (row_id, col_id, value) triples. Order-
    independent, so row/column permutations don't register."""
    lines = [l for l in raw.splitlines() if l.strip()]
    if not lines:
        return ""
    header = lines[0].split("\t")
    col_ids = header[1:]
    triples = []
    for row in lines[1:]:
        cells = row.split("\t")
        rid = cells[0]
        for cid, val in zip(col_ids, cells[1:]):
            triples.append((rid, cid, _round_token(val)))
    return "\n".join(f"{a}\t{b}\t{v}" for a, b, v in sorted(triples))


def canon_pairwise(raw: str) -> str:
    """Pairwise distance table -> unordered isolate pair, rounded numerics, sorted rows.
    First two columns are treated as the (unordered) pair; the rest as values."""
    lines = raw.splitlines()
    if not lines:
        return ""
    header = lines[0]
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        c = line.split("\t")
        if len(c) >= 2:
            pair = sorted([c[0], c[1]])
            rest = [_round_token(x) for x in c[2:]]
            rows.append("\t".join(pair + rest))
        else:
            rows.append(line)
    return header + "\n" + "\n".join(sorted(rows))


def canon_tsv(raw: str) -> str:
    """Generic TSV -> path fields to basename, floats rounded, DATA rows sorted (header
    kept first)."""
    lines = raw.splitlines()
    if not lines:
        return ""
    header = lines[0]
    out_rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        toks = [_round_token(_basename_paths(t)) for t in line.split("\t")]
        out_rows.append("\t".join(toks))
    return header + "\n" + "\n".join(sorted(out_rows))


def canon_lines(raw: str) -> str:
    """List-style content (snplist, Reference_IDs): a SET of lines — normalize paths and
    SORT, since row order here is a `list(set(...))` artifact, never meaningful. Column
    order of the alignment is handled separately (canon_alignment), so sorting snplist is
    safe."""
    return "\n".join(sorted(_basename_paths(l) for l in raw.splitlines() if l.strip()))


def canon_alignment(fasta_raw: str, snplist_raw: str) -> str:
    """Core-SNP alignment made COLUMN-ORDER-INDEPENDENT. snplist line j labels alignment
    column j; the column order is a nondeterministic `set` artifact but the biology (and
    every downstream consumer: distances, IQ-TREE) is invariant to it. We rebuild each
    locus's column vector and sort by locus, so a permutation is equal but a changed base,
    added/dropped locus, or changed isolate is NOT."""
    recs, cur, seq = {}, None, []
    for line in fasta_raw.splitlines():
        if line.startswith(">"):
            if cur is not None:
                recs[cur] = "".join(seq)
            cur, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if cur is not None:
        recs[cur] = "".join(seq)
    loci = [l.strip() for l in snplist_raw.splitlines() if l.strip()]
    ids = sorted(recs)
    # defensive: if lengths disagree, fall back to record-sorted fasta
    if any(len(recs[i]) != len(loci) for i in ids):
        return canon_fasta(fasta_raw)
    rows = []
    for j, locus in enumerate(loci):
        col = ",".join(f"{i}:{recs[i][j]}" for i in ids)
        rows.append(f"{locus}\t{col}")
    return "\n".join(sorted(rows))


def classify(rel: str):
    """Return (kind, canonicalizer) for a relative path, or (None, None) to skip."""
    name = os.path.basename(rel)
    parts = set(Path(rel).parts)
    if parts & EXCLUDE_DIR_PARTS:
        return None, None
    if name in EXCLUDE_NAMES:
        return None, None
    if rel.endswith(EXCLUDE_SUFFIXES):
        return None, None
    if rel.endswith(".snpdiffs"):
        return "snpdiffs", canon_snpdiffs
    if rel.endswith(".fasta") or rel.endswith(".fna") or rel.endswith(".fa"):
        # alignment outputs only; input assemblies never live in the output dir
        return "fasta", canon_fasta
    if "matrix" in name and rel.endswith(".tsv"):
        return "matrix", canon_matrix
    if "pairwise" in name and rel.endswith(".tsv"):
        return "pairwise", canon_pairwise
    if rel.endswith(".tsv") or rel.endswith(".csv"):
        return "tsv", canon_tsv
    if name.endswith(".txt"):
        return "lines", canon_lines
    return None, None


def build_manifest(root: Path) -> dict:
    entries = {}
    skipped = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        rel = str(path.relative_to(root))
        kind, fn = classify(rel)
        if fn is None:
            skipped.append(rel)
            continue
        try:
            raw = path.read_text(errors="replace")
        except Exception as e:            # noqa: BLE001
            entries[rel] = {"kind": kind, "error": str(e)}
            continue
        # Alignments are canonicalized against their snplist for column-order independence.
        if os.path.basename(rel).startswith("snpma") and rel.endswith(".fasta"):
            sib = os.path.basename(rel).replace("snpma", "snplist").replace(".fasta", ".txt")
            sibpath = path.parent / sib
            if sibpath.exists():
                canon, kind = canon_alignment(raw, sibpath.read_text()), "alignment"
            else:
                canon = fn(raw)
        else:
            canon = fn(raw)
        entries[rel] = {"kind": kind, "sha256": _sha(canon), "lines": canon.count("\n") + 1}
    return {"root": str(root), "entries": entries, "skipped": sorted(skipped)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("output_dir", type=Path)
    ap.add_argument("-o", "--out", type=Path, default=None,
                    help="write manifest JSON here (default: stdout)")
    args = ap.parse_args()
    if not args.output_dir.is_dir():
        sys.exit(f"not a directory: {args.output_dir}")
    manifest = build_manifest(args.output_dir)
    text = json.dumps(manifest, indent=2, sort_keys=True)
    if args.out:
        args.out.write_text(text)
        n = len(manifest["entries"])
        print(f"manifest: {n} canonicalized files -> {args.out} "
              f"({len(manifest['skipped'])} skipped)")
    else:
        print(text)


if __name__ == "__main__":
    main()
