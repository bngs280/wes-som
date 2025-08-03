process somVarCall {
    tag "Somatic Variant on ${sample_id}"
    publishDir "${params.outdir}/5.Variant", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(bam), path(bai)
    val(params.ref)
    val(params.bed)
    val(params.gnomad)
    val(params.db1000g)

    output:
    tuple val(sample_id), path("${sample_id}_raw.vcf.gz"), path("${sample_id}_raw.vcf.gz.tbi"), path("${sample_id}_raw.vcf.gz.stats"), emit: raw_vcfs

    script:
    """
    java -jar /usr/src/app/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar Mutect2 \
        -R ${params.ref} \
        --native-pair-hmm-threads ${task.cpus} \
        -L ${params.bed} \
        --germline-resource ${params.gnomad} \
        --panel-of-normals ${params.db1000g} \
        -I ${bam} \
        -O ${sample_id}_raw.vcf.gz
    """
}

process validatevcf {
    tag "${sample_id}_validation"

    input:
    tuple val(sample_id), path("${sample_id}_raw.vcf.gz")

    output:
    tuple val(sample_id), path("${sample_id}_raw.vcf.gz"), emit: valid_vcfs

    script:
    """
    python3 /usr/src/app/ref17/Validation_script/vcf_validation.py --vcf ${sample_id}_raw.vcf.gz
    """
}

process FilterMT {
    tag "Somatic Filter on ${sample_id}"
    publishDir "${params.outdir}/6.FilterSom", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(vcf), path(tbi), path(tsv)
    val(params.ref)
    
    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf.gz"), emit: filtered_vcfs

    script: 
    """
    java -jar /usr/src/app/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar FilterMutectCalls \
        -R ${params.ref} \
        -V ${vcf} \
        -O ${sample_id}_filtered.vcf.gz
    """
}

process KeepPASS {
    tag "Extract PASS variants on ${sample_id}"
    publishDir "${params.outdir}/7.FinalFilteredVCF", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(vcf)//worked
    
    output:
    tuple val(sample_id), path("${sample_id}_filtered_PASS.vcf"), emit: final_vcfs

    script: 
    """
    bcftools view -f PASS ${vcf} > ${sample_id}_filtered_PASS.vcf
    """
}
process somVarCall_ont {
    tag "Somatic Variant on ${sample_id} (ONT)"
    publishDir "${params.outdir}/5.Variant", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(bam), path(bai)
    val(params.ref)
    val(params.bed)
 
    output:
    //path "${sample_id}", emit: final_vcfs
    tuple val(sample_id), path("${sample_id}/${sample_id}.vcf"), emit: final_vcfs

    script:
    """
    source activate nanocaller_env 
    mkdir -p ${sample_id}
    /usr/src/app/NanoCaller/./NanoCaller \
        --bam ${bam} \
        --ref ${params.ref} \
        --bed ${params.bed} \
        --cpu 20 \
        --mode all \
        --sequencing ont \
        --output ${sample_id} \
        --prefix ${sample_id}
    """
}
process validateFinalvcf {
    tag "${sample_id}_validation"

    input:
    tuple val(sample_id), path("${sample_id}.vcf")

    output:
    tuple val(sample_id), path("${sample_id}.vcf"), emit: valid_final_vcfs

    script:
    """
    python3 /usr/src/app/ref17/Validation_script/vcf_Finalvalidation.py --vcf ${sample_id}.vcf
    """
}