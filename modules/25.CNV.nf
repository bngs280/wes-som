process CNV_som {
    tag "Somatic Copy Number on ${sample_id}"
    publishDir "${params.outdir}/14.CNV_somatic", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(bam), path(bai)
    val(params.ref)
    val(params.baselineContra)
    val(params.bed)

    output:
    path "${sample_id}", emit: cnvSom
    
    script:
    """
    source activate p2Manta
    python /usr/src/app/CONTRA.v2.0.8/contra.py \
    --target ${params.bed} \
    --test ${bam} \
    --control ${params.baselineContra} \
    --fasta ${params.ref} \
    --largeDeletion --bed -p --sampleName ${sample_id} -o ${sample_id}
    """
}