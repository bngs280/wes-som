// QC check Paired End
process quality_check {
    tag "Quality Checking on ${sample_id}"
    publishDir "${params.outdir}/1.QC", mode: 'copy'

    input:
    tuple val(sample_id), path(read_files)

    output:
    path "${sample_id}", emit: fastqc_out

    script:
    """
    mkdir ${sample_id}
    fastqc -o ${sample_id} -f fastq -q ${read_files}
    """
}

// QC check Single End
process quality_checkONT {
    tag "Quality Checking on ${sample_id} (ONT)"
    publishDir "${params.outdir}/1.QC", mode: 'copy'

    input:
    tuple val(sample_id), path(read_files)

    output:
    path "${sample_id}", emit: fastqc_ont_out

    script:
    """
    mkdir ${sample_id}
    NanoPlot --fastq ${read_files} -o ${sample_id}
    """
}