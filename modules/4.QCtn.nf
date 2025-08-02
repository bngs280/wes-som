// QC of Tumor-Normal paired-end
process quality_checkTN {
    tag "QCTN on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/QCTN", mode: 'copy'
    input:
    tuple val(tumor_id), path(tumor_fastq1), path(tumor_fastq2),
          val(normal_id), path(normal_fastq1), path(normal_fastq2)
    
    output:
    path "${tumor_id}_${normal_id}", emit: fastqc_outTN

    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'
    
    script:
    """
    mkdir ${tumor_id}_${normal_id}
    fastqc -o ${tumor_id}_${normal_id} -f fastq -q ${tumor_fastq1} ${tumor_fastq2} ${normal_fastq1} ${normal_fastq2}
    """
}
