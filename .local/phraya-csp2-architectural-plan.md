# Phraya × CSP2 Refactor — Architectural Plan & Decision

## Decision: Proceed with Phraya, phased.

Phraya's aligner engine is production-ready (SIMD dispatch fixed 2026-06-08, CSP2 spike
findings B1-B4 all resolved in Unreleased CHANGELOG). It can replace `runMUMmer`'s dnadiff
role for variant calling. **However**, Phraya is NOT a drop-in mummer replacement — its
output is a binary `.phraya` format (MessagePack + zstd), and it has **zero** distance-matrix,
core-SNP-concatenation, multi-sample-VCF, or tree-inference code. The downstream SNP-matrix
and tree stages must remain in Python (or be newly written).

This is explicitly documented as planned future work in Phraya's own PRD (Phase 5) and
CLAUDE.md (deliberately excluded).

## Capability Matrix

| CSP2 process | Current tool | Phraya equivalent | Action |
|---|---|---|---|
| `runMUMmer` (dnadiff/nucmer) | mummer | ✅ `phraya plan` + `phraya align --reference` | Replace — adapter converts `.phraya`→`.snpdiffs` |
| `saveMUMmerLog` | python | retains | Retains python, reads adapter output |
| `screenSNPDiffs` | python + bedtools | ⚠️ Partial (filters) | Retains python; `phraya filter` QC fields now available |
| `runSnpPipeline` | python + bedtools | ❌ No (distance matrix/core-SNP FASTA) | Retains python entirely |
| `compileResults` | python | ❌ No | Dead code, untouched |
| `runiqtree` | iqtree | ❌ No | Retains iqtree |
| `chooseRefs` | python + mash | ✅ Partial (centroid selection) | Optional future swap (mash triangle → phraya plan centroid) |
| `mashSketch`/`mashTriangle` | mash | ✅ Equivalent (minimizer sketch Jaccard) | Optional future swap |
| `fetchReads`/`getSNPDiffsData` | python | N/A (data fetching) | Retains python |
| `skesaAssemble` | skesa | ❌ No | Retains skesa |

## The Scope Seam (clean, pre-decided)

**Phraya = alignment engine + variant emission + QC-field-based filtering.**
**Pipeline = distance matrix, core-SNP matrix, multi-sample VCF, tree inference, I/O orchestration.**

This is architectural, not accidental:
- Phraya writes one `.phraya` file per query-vs-reference pair (reference space model).
- CSP2's `runSNPPipeline.py` pivots all pairs into a per-Ref_Loc × per-Query_ID matrix,
  computes pairwise SNP distances, writes `snpma_preserved.fasta`, and (when `--notree` is off)
  feeds it to iqtree. None of that logic exists in Phraya.
- Grep proof: zero hits for `distance_matrix|pairwise_dist|snpma|core_snp|concat` in
  phraya source. Zero hits for `tree|phylog|iqtree` in phraya source.

## Proposed Implementation — 2 Phases

### Phase A (single PR): Replace runMUMmer with Phraya align

**Goal**: `runMUMmer` calls `phraya plan` + `phraya align --reference` instead of dnadiff;
a new `phraya2snpdiffs.py` adapter converts `.phraya` + `.phraya.queries` → CSP2's `.snpdiffs`
format so `screenSNPDiffs`/`runSnpPipeline`/`saveMUMmerLog` are byte-for-byte unchanged.

**New component: `bin/phraya2snpdiffs.py`**
- Input: `.phraya` file (PhrayaFile via read_phraya), reference FASTA path, query/ref IDs.
- Reads `VariantObservation` fields: position (0-based), ref_base, all_alleles,
  mapq, confidence, cigar, edit_distance, local_coverage, avg_base_quality,
  variant_type (Snp/Insertion/Deletion), kmer_uniqueness, strand, query_position,
  snp_density_15/125/1000.
- Maps to CSP2's `.snpdiffs` SNP/Indel rows: Ref_Contig, Start_Ref, Ref_Pos (1-based),
  Query_Contig, Start_Query, Query_Pos, Ref_Loc (contig/pos), Query_Loc,
  Ref_Start/Ref_End (from coords), Query_Start/Query_End, Ref_Base, Query_Base,
  Dist_to_Ref_End, Dist_to_Query_End, Ref_Aligned, Query_Aligned, Query_Direction,
  Perc_Iden (= confidence, since Phraya's confidence = 1 − edit_dist/aligned), Cat (SNP/Indel).
- Indels: Phraya encodes deletions/ref_base=deleted_bases, insertions/ref_base="."; map
  to CSP2's `Cat=Indel` with Query_Base="." or Ref_Base="." accordingly.
- BED rows (##-prefixed header): reference/query contig, start, end, length, aligned
  length, percent identity. Derived from CoverageTrack breadth + PhrayaFile header.
- `#`-prefixed metadata header (## is BED rows, # is summary): total_snp_count,
  percent_ref_aligned, percent_query_aligned, median_percent_identity,
  median_alignment_length, kmer_similarity, shared_kmers, etc.
- Output: `<query>__vs__<ref>.snpdiffs` (same filename pattern), printed as
  `query,reference,snpdiffs_file` for Nextflow channel compatibility.

**Nextflow wiring changes** (`subworkflows/alignData/main.nf` `runMUMmer`):
```groovy
process runMUMmer {
    label 'mummerMem'
    publishDir snpdiffs_directory, mode: 'copy', pattern: '*.snpdiffs'
    input:
        tuple val(query_name), val(query_fasta), val(ref_name), val(ref_fasta)
    output:
        tuple val(query_name), val(ref_name), path("${query_name}__vs__${ref_name}.snpdiffs")
    script:
        report_id = "${query_name}__vs__${ref_name}"
        mummer_log = file("${mummer_log_directory}/${report_id}.log")
        """
        phraya plan --reference ${ref_fasta} --output ${work}/plan.phrayaplan ${query_fasta}
        phraya align --plan ${work}/plan.phrayaplan --reference ${ref_fasta} --output ${work}/aligned.phraya --strategy balanced
        python ${projectDir}/bin/phraya2snpdiffs.py --phraya ${work}/aligned.phraya --query ${query_name} --query_fasta ${query_fasta} --reference ${ref_name} --reference_fasta ${ref_fasta} --out ${report_id}.snpdiffs --log ${mummer_log}
        """
}
```

**Container**: Add a `phraya` conda package entry to `conf/CSP2.yaml` (Python wheels
available on PyPI — `phraya-aligner`). Add `process.container` profile for Singularity/Docker
using the existing `docker/Dockerfile` as a base layer with phraya installed.

### Phase B (separate PR, after benchmarking): Replace mash reference selection

If Phraya's centroid selection (`phraya plan` Case 3 auto-centroid) produces equivalent
reference choices to CSP2's `chooseRefs.py`, swap `mashSketch`/`mashTriangle`/`chooseRefs`
for `phraya plan` centroid mode. Gate with the benchmark harness. Not in scope for Phase A.

## Risk Register

| Risk | Likelihood | Mitigation |
|---|---|---|
| `.snpdiffs` field mapping differs → screenSNPDiffs/runSnpPipeline break | Medium | Adapter tested against real `.phraya` from live HPC, gated by `benchmark/compare.py` |
| Indel encoding mismatch (Phraya `.`/`ref_base` vs CSP2 Indel Cat) | Medium | Explicit mapping in adapter + unit test per variant_type |
| K-mer similarity (`compare_kmers` in compileMUMmer.py) has no Phraya equivalent | High | Adapter calls `kmercountexact.sh` for kmer metrics only; or drop kmer_similarity column (check if downstream consumes it) |
| Phraya strategy default (balanced) ≠ mummer identity behavior | Medium | Gate equivalence: same input → same `.snpdiffs` content via benchmark harness |
| Container build fails (no phraya wheel on bioconda for CSP2 pin) | Low-Medium | Use pip-installed wheel; verify CSP2.yaml phraya package builds |

## Validation Plan (reuses existing harness)

- **Pre/post equivalence**: Same 4-assembly test set. Run pipeline with mummer (current tip)
  → `.snpdiffs` manifest. Then run with Phraya adapter → `.snpdiffs` manifest. Compare with
  `benchmark/compare.py` expecting EQUIVALENT.
- **Field-level unit tests** for `phraya2snpdiffs.py` against a fixed `.phraya` fixture,
  asserting exact column output for SNP, Insertion, Deletion, empty-alignment cases.
- **Local gate** (§8a): `standard` profile, conda env, 4 genomes — no SLURM needed.
- **HPC gate** (§8c): `csp2_reedling` profile on live Reedling SLURM.

## Out of Scope (not touched)
- `runSnpPipeline` distance matrix / core-SNP FASTA (Phase B candidate — Phraya has no equivalent)
- `runiqtree`/iqtree (no Phraya tree code exists)
- `compileResults` (dead code)
- `fetchReads`/`getSNPDiffsData`/`skesaAssemble` (I/O and assembly, not alignment)
