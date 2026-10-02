Phase A Wiring Spec: Replace `runMUMmer`'s dnadiff with `phraya plan` + `phraya align` + `phraya filter --format snpdiffs`

**Status**: Gated resolved — Phraya's `.snpdiffs` formatter ships as `phraya filter --format
snpdiffs` (Phraya-side design decision: snpdiffs is an emission/interpretation-layer format,
not an alignment-layer one). See `~/Documents/specs_and_handoffs/phraya-csp2-handoff-back.md`
</!-- (The .snpdiffs formatter lives in phraya filter, not phraya align — see §2 below.) -->
subworkflows/alignData/main.nf` changes from:
```groovy
script:
    report_id = "${query_name}__vs__${ref_name}"
    mummer_log = file("${mummer_log_directory}/${report_id}.log")
    """
    ( cd ${mummer_directory} && dnadiff -p ${report_id} ${ref_fasta} ${query_fasta} )
    python ${mummerScript} --query "${query_name}" --query_fasta "${query_fasta}" --reference "${ref_name}" --reference_fasta "${ref_fasta}" --mummer_dir "${mummer_directory}" --snpdiffs_dir "." --temp_dir "${temp_dir}" --log_file "${mummer_log}"
    """
```

to:
```groovy
script:
    report_id = "${query_name}__vs__${ref_name}"
    mummer_log = file("${mummer_log_directory}/${report_id}.log")
    """
    plan_file = "${work}/${report_id}.phrayaplan"
    phraya plan --reference ${ref_fasta} --inputs ${query_fasta} --output ${plan_file}

    # Balanced strategy: Myers primary (K=2), WFA fallback for divergent/secondary hits.
    # Score-ratio ≥ 0.95 threshold is Phraya's default for variant reporting.
    phraya align --reference ${ref_fasta} ${plan_file} \
        --output ${work}/${report_id} \
        --strategy balanced

    phraya filter ${work}/${report_id}/ref.phraya \
        --format snpdiffs \
        --output ${report_id}.snpdiffs \
        --reference-fasta ${ref_fasta} \
        --query-fasta ${query_fasta} \
        --reference-id "${ref_name}"

    python ${projectDir}/bin/saveSNPDiffsLog.py \
        --out ${report_id}.snpdiffs \
        --log_file "${mummer_log}" \
        --query "${query_name}" --reference "${ref_name}"
    """
```

### What stays unchanged
- **Process I/O contract**: same `input` (tuple of query/ref fastas), same
  `output` tuple, same `publishDir` pattern `*.snpdiffs`, same filename
  `${query}__vs__${ref}.snpdiffs`. Downstream `alignGenomes` workflow is untouched.
- **`saveMUMmerLog`**: still runs on the emitted snpdiffs paths — the log writer
  script changes from `compileMUMmer.py` (which parsed mummer output) to a new
  `saveSNPDiffsLog.py` (which reads the `#` header from the `.snpdiffs` file that
-  `phraya filter --format snpdiffs` now emits).
- **`mummer_directory`**: still used as a scratch space in the `work` dir, but no
  longer needs the `dnadiff` subshell.

### What changes
- `dnadiff` + `compileMUMmer.py` → `phraya plan` + `phraya align` + `phraya filter --format snpdiffs`.
- `compileMUMmer.py` is replaced by Phraya's built-in snpdiffs formatter (`phraya filter --format
  snpdiffs` in `phraya-filter/src/snpdiffs.rs`), which emits the header with query/reference
  metadata + summary metrics.
- A new `bin/saveSNPDiffsLog.py` (trivial, ~30 lines) replaces the logging portion
  of `compileMUMmer.py` — reads the `.snpdiffs` `#` header, writes the log file
  in the same format `saveMUMmerLog` expects. This is a pipeline-layer concern,
  not a Phraya feature.

### Equivalence validation (gated by `phraya` >=0.2.0)

When the feature ships, run:

```bash
# Current (mummer) baseline:
nextflow run CSP2.nf -profile csp2_reedling --runmode snp --fasta testdata/ --n_ref 1 --out mummer_baseline
python benchmark/normalize.py mummer_baseline/ -o mummer_manifest.json

# Phraya replacement:
nextflow run CSP2.nf -profile csp2_reedling --runmode snp --fasta testdata/ --n_ref 1 --out phraya_run
python benchmark/normalize.py phraya_run/ -o phraya_manifest.json

python benchmark/compare.py mummer_manifest.json phraya_manifest.json   # expect EQUIVALENT
```
