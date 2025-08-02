// Tumor-Normal Variant calling
process somVarCall_tumor_normal {
    tag "Somatic Variant on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/VariantCallingTN", mode: 'copy'
    cpus 4

    input:
    tuple val(tumor_id), path(tumor_bam), path(tumor_bai),
          val(normal_id), path(normal_bam), path(normal_bai)
    val(params.ref)
    val(params.bed)
    val(params.gnomad)
    val(params.db1000g)

    output:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_raw.vcf.gz"), path("${tumor_id}_vs_${normal_id}_raw.vcf.gz.tbi"), path("${tumor_id}_vs_${normal_id}_raw.vcf.gz.stats"), emit: final_vcfs
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered.vcf.gz"), emit: filtered_vcfs
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered_PASS.vcf"), emit: final_vcfs
    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'

    script:
    """
    java -jar /usr/src/app/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar Mutect2 \
        -R ${params.ref} \
        -L ${params.bed} \
        -I ${tumor_bam} \
        -I ${normal_bam} \
        -normal ${normal_id} \
        --germline-resource ${params.gnomad} \
        --panel-of-normals ${params.db1000g} \
        -O ${tumor_id}_vs_${normal_id}_raw.vcf.gz

    # Validation
    python3 ../script/vcf_validation.py --vcf ${tumor_id}_vs_${normal_id}_raw.vcf.gz

    # Filter
    java -jar /usr/src/app/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar FilterMutectCalls \
        -R ${params.ref} \
        -V ${tumor_id}_vs_${normal_id}_raw.vcf.gz \
        -O ${tumor_id}_vs_${normal_id}_filtered.vcf.gz

    # Keep pass
    bcftools view -f PASS ${tumor_id}_vs_${normal_id}_filtered.vcf.gz > ${tumor_id}_vs_${normal_id}_filtered_PASS.vcf

    # Validation of final VCF
    python3 ../script/vcf_Finalvalidation.py --vcf ${tumor_id}_vs_${normal_id}_filtered_PASS.vcf
    """
}
