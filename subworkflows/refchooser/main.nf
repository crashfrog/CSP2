// Subworkflow to run RefChooser for list of queries

// Set directory structure
output_directory = file(params.output_directory)
log_directory = file(params.log_directory)
mash_directory = file(params.mash_directory)

workflow runRefChooser{
    take:
    query_data

    emit:
    reference_data

    main:

    // Sketch each unique query (1 CPU each). mashSketch now emits a typed
    // (name, .msh path) tuple staged through channels instead of writing to a shared
    // directory and echoing the path via stdout.
    sketches = query_data
    .unique{it -> it[1]}
    .map { [ it[0], it[1] ] }
    | mashSketch

    // Collect the sketch files and build the triangle. The .msh files are staged into
    // the triangle task's workdir; `ls *.msh` there reproduces the same lexical order the
    // old shared-directory glob produced, so the triangle (and RefChooser's order-
    // sensitive selection) is byte-identical.
    mash_refs = sketches
    .map { it -> it[1] }
    .collect()
    | mashTriangle
    | chooseRefs
    | splitCsv | collect | flatten | collate(1)

    reference_data = query_data
    .map{it -> tuple(it[1].toString(),it[0])}
    .join(mash_refs, by:0)
    .map{tuple(it[1],it[0])}
    .unique{it -> it[0]}.collect().flatten().collate(2)

    // Save reference data to file
    reference_data
    .collect{it -> it[0]}
    | saveRefIDs
}

process chooseRefs{

    executor = 'local'
    cpus = 1
    maxForks = 1

    publishDir mash_directory, mode: 'copy', pattern: 'CSP2_Ref_Selection.tsv'

    input:
    path(mash_triangle)

    output:
    stdout

    script:

    ref_count = params.n_ref.toInteger()
    ref_script = file("${projectDir}/bin/chooseRefs.py")
    """
    python $ref_script --ref_count $ref_count --mash_triangle_file $mash_triangle --trim_name "${params.trim_name}"
    """
}

process mashTriangle{

    publishDir mash_directory, mode: 'copy', pattern: 'Mash_{Triangle,Sketches.txt}'

    input:
    path(mash_sketches)

    output:
    path("Mash_Triangle")

    script:

    """
    ls *.msh > Mash_Sketches.txt
    mash triangle -p ${params.cores} -l Mash_Sketches.txt > Mash_Triangle
    """
}

process mashSketch{
    cpus = 1

    publishDir mash_directory, mode: 'copy', pattern: '*.msh'

    input:
    tuple val(query_name),val(query_fasta)

    output:
    tuple val(query_name), path("${query_name}.msh")

    script:

    """
    mash sketch -s 10000 -p 1 -o ${query_name} $query_fasta
    """
}

process saveRefIDs{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(ref_ids)

    script:
    ref_id_file = file(params.ref_id_file)
    ref_id_file.append(ref_ids.join('\n') + '\n')        
    """
    """
}
