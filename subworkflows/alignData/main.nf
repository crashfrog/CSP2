// Subworkflow to run MUMmer for query/reference comparisons

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
mummerScript = file("$projectDir/bin/compileMUMmer.py")

workflow {
    take:
    read_data
    snpdiffs_data

    main:
    // Align genomes
    snpdiffs = alignGenomes(to_align: read_data, snpdiffs_data: snpdiffs_data)

    emit:
    snpdiffs
}

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
        | splitCsv

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

    input:
    tuple val(query_name), val(query_fasta), val(ref_name), val(ref_fasta)

    output:
    stdout

    script:
    report_id = "${query_name}__vs__${ref_name}"
    mummer_log = file("${mummer_log_directory}/${report_id}.log")

    // Ensure MUMmer directories exist
    if (!mummer_directory.isDirectory()) {
        error "$mummer_directory does not exist..."
    } else {
        """
        $params.load_mummer_module
        $params.load_python_module
        $params.load_bedtools_module
        $params.load_bbtools_module

        cd ${mummer_directory}
        dnadiff -p ${report_id} ${ref_fasta} ${query_fasta}

        # rm -rf ${mummer_directory}/${report_id}.mdelta
        # rm -rf ${mummer_directory}/${report_id}.mcoords
        # rm -rf ${mummer_directory}/${report_id}.1delta
        # rm -rf ${mummer_directory}/${report_id}.delta

        python ${mummerScript} --query "${query_name}" --query_fasta "${query_fasta}" --reference "${ref_name}" --reference_fasta "${ref_fasta}" --mummer_dir "${mummer_directory}" --snpdiffs_dir "${snpdiffs_directory}" --temp_dir "${temp_dir}" --log_file "${mummer_log}"
        """
    }
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
    $params.load_python_module
    python $saveSNPDiffs --snpdiffs_file "${all_snpdiffs_list}" --summary_file "${snpdiffs_summary_file}" --isolate_file "${isolate_data_file}" --trim_name "${params.trim_name}" --ref_id_file "${ref_id_file}"
    """
}
