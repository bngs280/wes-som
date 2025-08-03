process align_paired {
    tag "Alignment on ${sample_id} (Paired-End)"
    publishDir "${params.outdir}/3.Alignment", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz")
    val(params.ref)

    output:
    tuple val(sample_id), path("${sample_id}_sorted_md.bam"), path("${sample_id}_sorted_md.bam.bai"), emit: valid_markdup_bams    
    when:
    params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI' && params.library == 'Paired'
    
    script:
    """
    bwa mem -t ${task.cpus} -M \
            -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${sample_id}" \
            ${params.ref} ${sample_id}_R1_trimmed.fastq.gz ${sample_id}_R2_trimmed.fastq.gz | samtools view -@ ${task.cpus} -b - | samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam
    sambamba markdup -r -t ${task.cpus} ${sample_id}_sorted.bam ${sample_id}_sorted.md.bam 
    samtools index ${sample_id}_sorted.md.bam  
    """
}

process align_single {
    tag "Alignment on ${sample_id} (Single-End)"
    publishDir "${params.outdir}/3.Alignment", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz")
    val(params.ref)

    output:
    tuple val(sample_id), path("${sample_id}_sorted_md.bam"), path("${sample_id}_sorted_md.bam.bai"), emit: valid_markdup_bams

    when:
    params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI' || params.platform == 'ThermoFisher' && params.library == 'Single'

    script:
    """
    bwa mem -t ${task.cpus} -M \
            -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${sample_id}" \
            ${params.ref} ${sample_id}_trimmed.fastq.gz | samtools view -@ ${task.cpus} -b - | samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam 
    sambamba markdup -r -t ${task.cpus} ${sample_id}_sorted.bam ${sample_id}_sorted.md.bam 
    samtools index ${sample_id}_sorted.md.bam 
    """
}
process align_ONT {
    tag "Alignment on ${sample_id} (ONT)"
    publishDir "${params.outdir}/3.Alignment", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz")
    val(params.ref)

    output:
    //tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: sorted_bams_ont
    tuple val(sample_id), path("${sample_id}_sorted_md.bam"), path("${sample_id}_sorted_md.bam.bai"), emit: sorted_bams_ont


    when:
    params.platform == 'Nanopore' && params.library == 'Single'

    script:
    """
    /usr/src/app/minimap2/./minimap2 -ax map-ont ${params.ref} ${sample_id}_trimmed.fastq.gz | samtools view -@ ${task.cpus} -b - | samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam
    sambamba markdup -r -t ${task.cpus} ${sample_id}_sorted.bam ${sample_id}_sorted.md.bam 
    samtools index ${sample_id}_sorted.md.bam
    """
}

process mdvalidate {
    tag "${sample_id}_mdvalidation"

    input:
    tuple val(sample_id), path("${sample_id}_sorted_md.bam")

    output:
    tuple val(sample_id), path("${sample_id}_sorted_md.bam"), emit: valid_markdup_bams

    script:
    """
    python3 /usr/src/app/ref17/Validation_script/picardOutPutvalidation.py --bam ${sample_id}_sorted_md.bam
    """
}