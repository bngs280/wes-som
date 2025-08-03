process CNV_somtn {
    tag "Somatic Copy Number on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/CNV_somatic", mode: 'copy'
    cpus 4

    input:
    tuple val(tumor_id), path(tumor_bam), path(tumor_bai),
          val(normal_id), path(normal_bam), path(normal_bai)
    val(params.ref)
    val(params.bed)

    output:
    path "${tumor_id}_vs_${normal_id}", emit: cnvSom

    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'
    
    script:
    """
    source activate p2Manta
    python /usr/src/app/CONTRA.v2.0.8/contra.py \
    --target ${params.bed} \
    --test ${tumor_bam} \
    --control ${normal_bam} \
    --fasta ${params.ref} \
    -p --sampleName ${tumor_id}_vs_${normal_id} -o ${tumor_id}_vs_${normal_id}
    """
}
