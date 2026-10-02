# Feature Handoff: CSP2-compatible `.snpdiffs` output for `phraya align`

**For:** Justin (Phraya maintainer)
**From:** Phraya × CSP2 refactor planning
**Context:** CSP2 (bacterial-genomics SNP pipeline, `~/projects/CSP2`) is being refactored
to replace `dnadiff`/`nucmer` (mummer) with `phraya plan` + `phraya align --reference`.
The blocker is an output-format gap described below.

---

## Gap: Phraya needs a CSP2 `.snpdiffs` output format

CSP2's downstream pipeline consumes a per-pair `.snpdiffs` file produced today by
`bin/compileMUMmer.py`, which parses mummer's `.snps`, `.coords`, `.delta`, and `.1coords`
files into a single TSV. The schema is:

```
#\t<query metadata key:value>... \t <reference metadata> \t <summary metrics>
##\t<Ref_Contig> <Ref_Start> <Ref_End> <Ref_Length> <Ref_Aligned>
    <Query_Contig> <Query_Start> <Query_End> <Query_Length> <Query_Aligned> <Perc_Iden>
<Ref_Contig> <Start_Ref> <Ref_Pos> <Query_Contig> <Start_Query> <Query_Pos>
    <Ref_Loc> <Query_Loc> <Ref_Start> <Ref_End> <Query_Start> <Query_End>
    <Ref_Base> <Query_Base> <Dist_to_Ref_End> <Dist_to_Query_End>
    <Ref_Aligned> <Query_Aligned> <Query_Direction> <Perc_Iden> <Cat>
```

Where:
- `Ref_Pos` is 1-indexed (mummer convention).
- `Start_Ref` / `Start_Query` are 0-indexed starts.
- `Cat` is `SNP`, `Indel`, or `Invalid`.
- `Ref_Loc` / `Query_Loc` are `contig/position` strings.
- `Dist_to_Ref_End` / `Dist_to_Query_End` are `min(pos, length−pos−1)`.
- `Perc_Iden` is `100 × (1 − edit_distance / query_aligned_length)`.
- `Query_Direction` is `F` or `R`.

CSP2's `csp2_common.py` parsers (`parseSNPDiffs`, `processSNPs`, `processBED`) consume
this exact schema. A Python adapter shim in `bin/` would need to deserialize the binary
`.phraya` format (MessagePack + zstd) — a layer that belongs in Phraya, not in CSP2.

## Request: `--format snpdiffs` output option

Add a `snpdiffs` output format to `phraya align` (and/or `phraya filter`), alongside the
existing `vcf`, `tsv`, and `phraya` formats:

```bash
phraya plan --reference ref.fasta --output plan.phrayaplan query.fasta
phraya align --plan plan.phrayaplan --output pair.snpdiffs --format snpdiffs \
    --query-id QUERY --reference-id REF --reference-fasta ref.fasta --query-fasta query.fasta
```

### What the format must satisfy

1. **Variant rows** (SNP + Indel):
   - `position` → 1-indexed `Ref_Pos`; `Start_Ref = position` (0-indexed).
   - `ref_base` → for SNPs: the reference base character; for deletions: the deleted base(s);
     for insertions: `.` (reference has no base).
   - `all_alleles` → for SNPs: the alternate allele with highest count; for insertions: the
     inserted base(s); for deletions: `.` (query has no base).
   - `variant_type` → `Snp` maps to `Cat=SNP`; `Insertion`/`Deletion` map to `Cat=Indel`.
   - `confidence` → `Perc_Iden = 100 × confidence` (Phraya stores 0–1; CSP2 expects 0–100).
   - `cigar` → parse via `CigarStats` to derive `Ref_Aligned` and `Query_Aligned`
     (target_aligned_len, query_aligned_len respectively). If CIGAR is absent, fall back
     to `reference_length` / query length from `--query-fasta`.
   - `query_position` → `Start_Query` (0-indexed); `Query_Pos = query_position + 1`.
   - `strand` → `Query_Direction`: `Forward`→`F`, `Reverse`→`R`.
   - `Dist_to_Ref_End = min(position, ref_len − position − 1)`, `Dist_to_Query_End` from
     `query_position` similarly.
   - `Ref_Loc`/`Query_Loc` = `"contig/pos"`.

2. **BED rows** (`##-prefixed`): the aligned blocks. For CSP2's single-reference-space
   use case, emit one row covering the full reference contig span with
   `Percent_Aligned = 100 × (covered_positions / reference_length)` derived from the
   `CoverageTrack`. If the reference FASTA has multiple contigs, emit one row per record
   that received any alignment.

3. **Header** (`#` line): query/reference metadata (contig count, assembly bases, N50/L90,
   SHA256) + summary metrics (SNP count, indel count, median identity, median alignment
   length, percent aligned). K-mer similarity can be `NA` (CSP2 computes it via
   `kmercountexact.sh` separately and does not depend on it for correctness).

### Why in Phraya, not a Python shim

- `.phraya` is a binary MessagePack + zstd format. Deserializing it in Python requires
  `msgpack` + `zstandard` + exact knowledge of the struct layout — coupling CSP2 to
  Phraya's internal serialization format, which can change.
- The `VariantObservation` → `.snpdiffs` field mapping is non-trivial (coordinate
  systems, Phraya's inverted CIGAR alphabet, indel base placement, k-mer similarity opt-out).
  Keeping this translation in Rust guarantees it's validated by Phraya's test suite
  and evolves with the format.
- `phraya filter` already has the VCF and TSV formatters in `phraya-filter/src/`; a
  `Snpdiffs` formatter follows the same pattern.

### Suggested Phraya-side placement

- New module `phraya-filter/src/snpdiffs.rs` with `fn format_snpdiffs(phraya_file, query_id, reference_id, query_fasta_len, reference_fasta_len) -> String`
- Register `snpdiffs` in `phraya-cli/src/main.rs`'s `Filter` `format` validation
  (`["vcf", "tsv", "phraya", "snpdiffs"]`) and dispatch to the new formatter.
- CLI additions on `phraya align`/`phraya filter`:
  - `--format snpdiffs`
  - `--query-id <name>` and `--reference-id <name>` (for the header metadata)
  - `--reference-fasta <path>` and `--query-fasta <path>` (for contig lengths in BED rows)
- Existing `VCF`/`TSV` tests in `phraya-filter/src/` provide the template.

### Acceptance criterion (for this handoff)

The following CLI invocation produces a file byte-compatible with what
`bin/compileMUMmer.py` emits, for the SNP/Indel/empty-alignment cases:

```bash
phraya align --plan plan.phrayaplan --output pair.snpdiffs \
    --format snpdiffs --strategy balanced \
    --reference-fasta ref.fasta --query-fasta query.fasta \
    --reference-id REF --query-id QUERY
# then in CSP2:
python -c "from csp2_common import parseSNPDiffs; bed, snps = parseSNPDiffs('pair.snpdiffs', 'ref_first'); print(len(snps))"
```

Once this ships, Phase A of the CSP2 refactor (replacing `dnadiff` with `phraya align`)
is a ~10-line Nextflow subworkflow edit, gated by `bin/phraya2snpdiffs.py` being deleted.

## Out of scope for this handoff

- Distance matrix / core-SNP FASTA / multi-sample VCF / tree inference — all remain in
  CSP2's Python pipeline (Phase B). Phraya has no equivalent and PRD.md schedules them
  for Phase 5.
- `--format mash-triangle` (replacing `mashSketch`/`mashTriangle`) — optional future swap,
  not required for Phase A.
