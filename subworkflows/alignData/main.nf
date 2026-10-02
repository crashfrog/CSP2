// Subworkflow to run Phraya aligner for query/reference comparisons


// Set path variables
output_directory = file(params.output_directory)
mummer_directory = file(params.mummer_directory)
mummer_log_directory = file(params.mummer_log_directory)
snpdiffs_directory = file(params.snpdiffs_directory)
log_directory = file(params.log_directory)

if(params.tmp_dir == ""){
    temp_dir = ""
} else{
    temp_dir = file(params.temp_dir)
}

ref_mode = params.ref_mode
ref_id_file = file(params.ref_id_file)

// Set path to accessory scripts/files
all_snpdiffs_list = file("${log_directory}/All_SNPDiffs.txt")
isolate_data_file = file("${output_directory}/Isolate_Data.tsv")
snpdiffs_summary_file = file("${output_directory}/Raw_MUMmer_Summary.tsv")

workflow alignGenomes {
    take:
    to_align
    snpdiffs_data

    emit:
    return_snpdiffs

    main:
    // Align anything that needs aligning
    sample_pairwise = to_align
        .filter { "${it[0]}" != "${it[2]}" } // Don't map things to themselves
        | runMUMmer
        | map { q, r, sd -> [ q, r, sd.toString() ] } // (query, reference, snpdiffs path)

    log_hold = sample_pairwise
        .concat(snpdiffs_data)
        .unique { it -> it[2] }
        .collect { it -> it[2] }

    snpdiff_files = saveMUMmerLog(log_hold)
        .collect().flatten().collate(1)

    return_snpdiffs = sample_pairwise
        .concat(snpdiffs_data)
        .map { it -> tuple([it[0], it[1]].sort().join(',').toString(), it[0], it[1], it[2]) }
        .unique { it -> it[0] }
        .map { it -> tuple(it[3], it[1], it[2]) }
        .join(snpdiff_files, by: 0)
        .map { it -> tuple(it[1], it[2], it[0]) }
}

process runMUMmer {
    label 'mummerMem'

    // .snpdiffs is now a real task output staged from the workdir and published to the
    // canonical dir. Alignment is performed by phraya (WFA O(s·n)), and the .snpdiffs
    // text format is produced by `phraya filter --format snpdiffs` — see
    // docs/phraya-runMUMmer-wiring-spec.md and
    // ~/Documents/specs_and_handoffs/phraya-csp2-handoff-back.md §1.3.
    publishDir snpdiffs_directory, mode: 'copy', pattern: '*.snpdiffs'

    input:
    tuple val(query_name), val(query_fasta), val(ref_name), val(ref_fasta)

    output:
    tuple val(query_name), val(ref_name), path("${query_name}__vs__${ref_name}.snpdiffs")

    script:
    report_id = "${query_name}__vs__${ref_name}"
    mummer_log = file("${mummer_log_directory}/${report_id}.log")
    plan_file = "${report_id}.phrayaplan"
    """
    # In Nextflow, the process working directory is the task workdir.
    # Use relative paths — plan_file is a Groovy var declared above, rendered
    # into the bash script here.

    phraya plan --reference ${ref_fasta} --inputs ${query_fasta} --output ${plan_file}

    # Balanced strategy: Myers primary (K=2), WFA fallback for divergent/secondary hits.
    # Score-ratio ≥ 0.95 threshold is Phraya's default for variant reporting.
    phraya align --reference ${ref_fasta} ${plan_file} \
        --output "${report_id}" \
        --strategy balanced

    phraya filter "${report_id}/ref.phraya" \
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
}

process saveMUMmerLog {
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(snpdiffs_paths)

    output:
    val(snpdiffs_paths)

    script:
    saveSNPDiffs = file("$projectDir/bin/saveSNPDiffs.py")
    all_snpdiffs_list.write(snpdiffs_paths.join('\n') + '\n')
    """
    python $saveSNPDiffs --snpdiffs_file "${all_snpdiffs_list}" --summary_file "${snpdiffs_summary_file}" --isolate_file "${isolate_data_file}" --trim_name "${params.trim_name}" --ref_id_file "${ref_id_file}"
    """
}
