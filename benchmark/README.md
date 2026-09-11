# CSP2 Equivalence Benchmark

Purpose: let us refactor CSP2's Nextflow and Python **without changing results**. The
refactor proceeds in rings of increasing blast radius:

1. **Upstream-able Nextflow/Python fixes** (this harness's job): restructure orchestration
   and de-duplicate Python while keeping the MUMmer pipeline's output *identical*.
2. **Phraya rearchitecture** (later, separate effort): swap the aligner. Equivalence
   becomes "within documented tolerance" rather than exact.

This harness freezes a **baseline** from the current MUMmer pipeline, then mechanically
checks that a candidate (a PR branch) reproduces it. See `../docs/refactor-phraya.md` for
the larger plan.

## What "equivalent" means

CSP2 outputs carry run-to-run noise that is not scientific: absolute paths,
`list(set(...))` row ordering, float formatting, timestamps. `normalize.py` reduces each
scientific output to a **canonical** form and hashes it:

- `.snpdiffs` — header paths → basenames; BED and SNP blocks sorted; content SHA256s kept.
- `snpma*.fasta` — records sorted by isolate id (column/sequence order preserved).
- `*matrix*.tsv` — reduced to sorted `(row, col, value)` triples (permutation-independent).
- `*pairwise*.tsv` — unordered isolate pairs, rounded numerics, sorted rows.
- other `*.tsv` — path fields → basenames, floats rounded, data rows sorted.
- **excluded**: `work/`, `.nextflow/`, logs, `CSP2_Params.txt`, raw `MUMmer_Output/`, mash
  sketches, timing artifacts, and IQ-TREE outputs (stochastic — we run with `--notree skip`).

`compare.py` diffs two manifests: **EQUIVALENT** iff every canonical file matches. A real
change to a SNP call, distance, alignment column, or QC verdict flips a hash.

## One-time setup

```bash
benchmark/setup.sh          # builds the `csp2-bench` conda env + nextflow, verifies tools
```

Runtime notes for this machine (4 vCPU, ~3.8 GiB RAM, Zscaler proxy):
- Bio tools + python come from conda (`environment.yml`, pinned to CSP2's pandas-1.x era).
- `nextflow_bench.config` caps task memory (stock labels ask 4–12 GB, which won't schedule
  locally here) — it changes **scheduling only, never algorithm parameters**.
- Conda downloads use the CA bundle at `~/.certs/ca-bundle-with-zscaler.crt`.

## Datasets (git-ignored, fetched from CFSAN-Biostatistics/data-commons)

- `data/refsel_full/`  — 10 real related draft assemblies (~4.7 MB each). Authoritative.
- `data/refsel_small/` — 4-assembly subset for the fast dev loop.

Re-fetch with the sparse-clone shown in the run log, or drop your own FASTAs into a dir.

## Freeze the baseline (current MUMmer pipeline)

```bash
# fast subset first (smoke), then the authoritative full run
benchmark/run.sh baseline_small snp benchmark/data/refsel_small --n_ref 1
benchmark/run.sh baseline_full  snp benchmark/data/refsel_full  --n_ref 3
```

Each writes `benchmark/runs/<name>/{out/,manifest.json,trace.txt,wall_time.txt}`.
Commit the `manifest.json` (small JSON) as the frozen reference; the `out/` and `work/`
trees are git-ignored.

## Validate a refactor PR

```bash
# on the PR branch, re-run with the SAME name suffix into a candidate dir, then:
benchmark/run.sh candidate_full snp benchmark/data/refsel_full --n_ref 3
python benchmark/compare.py \
    benchmark/runs/baseline_full/manifest.json \
    benchmark/runs/candidate_full/manifest.json
```

Exit 0 / `EQUIVALENT` ⇒ the refactor preserved every scientific output. Otherwise the
report lists added/removed/changed files; `--show <substr>` dumps the differing entries.

### `screen` mode is now equivalence-gateable

`snp` mode is fully deterministic (RefChooser `random_state=0` + single-thread nucmer).
`screen` mode used to be nondeterministic: with no reference, CSP2.nf built the all-vs-all
pair list with a mutable Groovy accumulator (`seen_combinations`), so the query/reference
ORIENTATION of each pair varied across runs of the *same code* (`A__vs__B` vs `B__vs__A` —
different `.snpdiffs` filenames, swapped Query/Reference columns). `fix/deterministic-pair-
selection` replaced that accumulator with a deterministic operator (canonical orientation:
lexicographically smaller id = query), so `screen` mode is now safe to gate the same way as
`snp` mode.

## Reproducibility check

Because RefChooser (KMeans `random_state=0`) and the filters are deterministic given fixed
inputs, two baseline runs must be EQUIVALENT to each other. Running the baseline twice and
comparing is the harness's own self-test (any diff means we haven't neutralized a source of
nondeterminism and the manifest can't be trusted yet).

## Files

| File | Role |
|------|------|
| `environment.yml` | pinned conda env (tools + python libs + java) |
| `setup.sh` | provision env + nextflow, verify tools |
| `nextflow_bench.config` | resource caps for a modest box (scheduling only) |
| `run.sh` | run one pipeline invocation → timing + manifest |
| `normalize.py` | canonicalize an output dir → equivalence manifest |
| `compare.py` | diff two manifests → EQUIVALENT / DIFFERENT |
