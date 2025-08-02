process align_tumor_normal {
    tag "Alignment of ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/AlignmentTN", mode: 'copy'
    cpus 8
    memory '16 GB'

    input:
    tuple val(tumor_id), path(tumor_fastq1), path(tumor_fastq2),
          val(normal_id), path(normal_fastq1), path(normal_fastq2)
    val(params.ref)

    output:
    tuple val(tumor_id), path("${tumor_id}.sorted.md.bam"),
          val(normal_id), path("${normal_id}.sorted.md.bam"), emit: sorted_bams_tn
    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'

    script:
    """
    # Align tumor FASTQ
    bwa mem -t ${task.cpus} -M \
            -R "@RG\\tID:${tumor_id}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${tumor_id}" \
            ${params.ref} ${tumor_fastq1} ${tumor_fastq2} | samtools view -@ ${task.cpus} -b - | samtools sort -o ${tumor_id}_sorted.bam
    sambamba markdup -r -t ${task.cpus} ${tumor_id}_sorted.bam ${tumor_id}.sorted.md.bam 
    samtools index ${tumor_id}.sorted.md.bam 

    # Align normal FASTQ
    bwa mem -t ${task.cpus} -M \
            -R "@RG\\tID:${normal_id}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${normal_id}" \
            ${params.ref} ${normal_fastq1} ${normal_fastq2} | samtools view -@ ${task.cpus} -b - | samtools sort -o ${normal_id}_sorted.bam
    sambamba markdup -r -t ${normal_id}_sorted.bam ${normal_id}.sorted.md.bam 
    samtools index ${normal_id}.sorted.md.bam 
    """
}