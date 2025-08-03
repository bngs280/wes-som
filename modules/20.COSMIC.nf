process cosmic_hg38 {
    tag "COSMIC annotation on ${sample_id}"
    publishDir "${params.outdir}/9.AdvanceAnnot", mode: 'copy'
    cpus 4

    input:
    val(params.cosmiccod)
    tuple val(sample_id), path("${sample_id}_amp_vcf_output.tsv")
 
    output:
    tuple val(sample_id), path("${sample_id}_amp_vcf_output_cosmic.tsv"), emit: vcf_vep_TSV1

    when:
    params.genome == 'hg38'
    script:
    """
    python /usr/src/app/ref17/Validation_script/cosmic_dataTSV_argument.py --cosmic_file ${params.cosmiccod} --input_file ${sample_id}_amp_vcf_output.tsv --output_file ${sample_id}_amp_vcf_output_cosmic.tsv
    """
}

process cosmic_hg19 {
    tag "COSMIC annotation on ${sample_id}"
    publishDir "${params.outdir}/9.AdvanceAnnot", mode: 'copy'
    cpus 4

    input:
    val(params.cosmiccod)
    tuple val(sample_id), path("${sample_id}_amp_vcf_output.tsv")
 
    output:
    tuple val(sample_id), path("${sample_id}_amp_vcf_output_cosmic.tsv"), emit: vcf_vep_TSV1

    when:
    params.genome == 'hg19'
    script:
    """
    python /usr/src/app/ref17/Validation_script/cosmic_dataTSV_argument.py --cosmic_file ${params.cosmiccod} --input_file ${sample_id}_amp_vcf_output.tsv --output_file ${sample_id}_amp_vcf_output_cosmic.tsv
    """
}