# Frozen baselines (current MUMmer pipeline)

These manifests are the **reference point** for the Nextflow/Python refactor PRs. A PR that
only restructures orchestration or de-duplicates Python MUST reproduce these exactly:

```bash
benchmark/run.sh candidate_full snp benchmark/data/refsel_full --n_ref 3
python benchmark/compare.py \
    benchmark/baselines/baseline_full.manifest.json \
    benchmark/runs/candidate_full/manifest.json      # expect: EQUIVALENT
```

(`compare.py` runs under the `csp2-bench` env python, or any Python 3.)

## Provenance

| | |
|---|---|
| CSP2 commit | `d77ec10` (pre-refactor `main`) |
| Dataset | `CFSAN-Biostatistics/data-commons` → `test/reference_selection/assemblies` |
| Runtime | nextflow 24.10.6, openjdk 17 (conda `csp2-bench`) |
| Python stack | python 3.8, pandas 1.2.5, numpy 1.22.4, scipy 1.8.1, scikit-learn 1.2.2 |
| Aligner | MUMmer4 `dnadiff`, **nucmer forced `-t 1`** (see determinism note) |
| Trees | disabled (`--notree skip`) — IQ-TREE is stochastic |
| Host | 4 vCPU, ~3.8 GiB RAM; executor serialized (`queueSize=1`) |

## What was run

| Baseline | Command | Wall | Canonical files |
|---|---|---|---|
| `baseline_full`  | `run.sh baseline_full snp data/refsel_full --n_ref 3` | 856 s | 60 |
| `baseline_small` | `run.sh baseline_small snp data/refsel_small --n_ref 1` | ~180 s | 17 |

### `baseline_full` scientific summary (for humans)

RefChooser selected 3 references; per-reference core-SNP alignment:

| Reference | isolates | core SNPs | preserved (max_missing 50) |
|---|---|---|---|
| SRR498285_contigs | 6 | 311 | 308 |
| SRR498286_contigs | 6 | 307 | 307 |
| SRR498387_contigs | 4 | 190 | 190 |

## Determinism notes (why these are reproducible)

Two sources of run-to-run noise were neutralized so that a diff means a *real* change:

1. **nucmer multithreading.** nucmer4 defaults to 2 threads and its anchoring is
   nondeterministic — back-to-back `dnadiff` runs gave different SNP/indel/coord counts
   (232 vs 233 SNPs, 105 vs 117 indels). `setup.sh` shadows the env `nucmer` with a `-t 1`
   wrapper (dnadiff calls nucmer from its own bin dir, ignoring PATH). Verified bit-identical
   across runs. This affects only the benchmark env; baseline and candidate both use it.
2. **`list(set(...))` ordering.** CSP2 orders some intermediate lists and the preserved
   alignment's columns via set iteration. `normalize.py` makes the manifest order-independent
   (rows sorted; alignments canonicalized as sorted per-locus column vectors against their
   snplist; intermediate path lists excluded). The preserved *distance matrix* is byte-identical
   across runs regardless — confirming the biology is stable and only column order floated.

Self-test: two independent `refsel_small` runs → `EQUIVALENT` (17/17). This is the harness's
own guarantee that the manifest can be trusted.

## Refreshing

If the dataset or pinned env changes, regenerate and re-freeze:

```bash
benchmark/setup.sh
benchmark/run.sh baseline_full  snp benchmark/data/refsel_full  --n_ref 3
benchmark/run.sh baseline_small snp benchmark/data/refsel_small --n_ref 1
cp benchmark/runs/baseline_full/manifest.json  benchmark/baselines/baseline_full.manifest.json
cp benchmark/runs/baseline_small/manifest.json benchmark/baselines/baseline_small.manifest.json
```
