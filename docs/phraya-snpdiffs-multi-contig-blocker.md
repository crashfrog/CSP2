# Blocker: Phraya `phraya filter --format snpdiffs` does NOT merge multiple .phraya inputs

**Severity**: Correctness blocker for Phase A
**For**: Justin (Phraya maintainer)

## Problem

CSP2's `runMUMmer` receives multi-contig reference FASTAs (confirmed: real SRR assemblies
at `/isilon/reggen/RAL_Pipelines/CSP2_Test/` produce 87+ contigs). Phraya's reference-palette
mode (`phraya align --reference multi.fasta`) writes **one `<contig_name>.phraya` per
reference contig** (ADR-0011). To produce CSP2's single `.snpdiffs` file per pair, all these
per-contig `.phraya` files must be merged.

However, `phraya filter` does NOT merge for `--format snpdiffs`:

```rust
// phraya-cli/src/main.rs:2218-2221
// Read all input .phraya files. For single-input formats (vcf, tsv, phraya,
// snpdiffs), only inputs[0] is used. For multi-sample formats (core-snp,
// distance-matrix), all inputs are read and each header.sample_id is used
// as the sample label.

// line 2448:
"Snpdiffs" => {
    let file = &phraya_files[0];  // <-- only first file used, rest ignored
```

Verified empirically on reedling2:
- `query1.phraya` → 1 SNP
- `query2.phraya` → 33 SNPs
- `phraya filter query1.phraya query2.phraya --format snpdiffs` → **1 SNP** (same as first
  file only; the 33 SNPs from query2.phraya are silently dropped)

## What CSP2 needs

A way to merge N per-contig `.phraya` files into one `.snpdiffs` file (all variant rows from
all reference contigs, single `#` header, single `##` BED block or one per contig).

## Proposed fixes (either is acceptable)

**Option A (Phraya-side, preferred):** Make `--format snpdiffs` aggregate observations from
ALL input `.phraya` files when ≥2 are provided (merge variant rows, sum SNP/indel counts,
extend the `##` BED block to one row per contig). This follows the existing multi-file
pattern used by `core-snp` and `distance-matrix` formats.

**Option B (CSP2-side workaround):** Call `phraya filter --format snpdiffs` on each `.phraya`
file individually, then concatenate the SNP/Indel data rows (skipping duplicate headers).
This is a lossy workaround — k-mer/identity stats would be per-contig, not aggregate — but
produces complete variant coverage.

**Option C (CSP2-side, Phraya plan tweak):** If Phraya gains a "single-reference mode" that
accepts a multi-FASTA as one target (treating all contigs as a single reference space), that
would sidestep the palette split entirely. ADR-0011 currently says "never collapsed to
whatever sat at index 0" — but for CSP2's use case (whole-assembly vs whole-assembly SNP
calling), one-space mode would be the natural fit.

## Recommendation

Option A or C. Option B is a stopgap that loses aggregate header metadata accuracy.

## Reproduction

```bash
phraya plan --reference ref_multi_contig.fasta --inputs query.fasta --output plan.phrayaplan
phraya align --reference ref_multi_contig.fasta plan.phrayaplan --output out/ --strategy balanced
ls out/      # → Contig_1.phraya Contig_2.phraya ... cross_space.phraya.queries
phraya filter out/*.phraya --format snpdiffs --output pair.snpdiffs  # BUG: only processes out/<first>.phraya
```
