// Trimming of Tumor-Normal paired-end
process fastpTumorNormal {
    tag "FASTP on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/TrimmingTN", mode: 'copy'

    input:
    tuple val(tumor_id), path(tumor_fastq1), path(tumor_fastq2), 
          val(normal_id), path(normal_fastq1), path(normal_fastq2)

    output:
    tuple val(tumor_id), path("${tumor_id}_trimmed_R1.fastq.gz"), path("${tumor_id}_trimmed_R2.fastq.gz"),
          val(normal_id), path("${normal_id}_trimmed_R1.fastq.gz"), path("${normal_id}_trimmed_R2.fastq.gz"), path("${tumor_id}_fastp.html"), path("${normal_id}_fastp.html"), emit: trimmed_fastqs_tn

    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'

    script:
    """
    fastp -q 20 -i ${tumor_fastq1} -I ${tumor_fastq2} -o ${tumor_id}_trimmed_R1.fastq.gz -O ${tumor_id}_trimmed_R2.fastq.gz --html ${tumor_id}_fastp.html --report_title "Quality Control for ${tumor_id}"
    fastp -q 20 -i ${normal_fastq1} -I ${normal_fastq2} -o ${normal_id}_trimmed_R1.fastq.gz -O ${normal_id}_trimmed_R2.fastq.gz --html ${normal_id}_fastp.html --report_title "Quality Control for ${normal_id}"
    """
}