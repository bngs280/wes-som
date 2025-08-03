process vcf2maf {
    tag "VCF to MAF on ${sample_id}"
    publishDir "${params.outdir}/9.AdvanceAnnot", mode: 'copy'
    cpus 4

    input:
    //tuple val(sample_id), path(vcf)
    //val(params.assemblyAMP)
    tuple val(sample_id), path("${sample_id}_filtered_PASS.vcf")
    val(params.cachee)
    val(params.ref)
    val(params.assembly)

    output:
    tuple val(sample_id), path("${sample_id}_vcf2MAF.maf"), emit: vcf_maf

   
    script:
    """
    perl /usr/src/app/vcf2maf/vcf2maf.pl --vep-path /usr/src/app/ensembl-vep --vep-data ${params.cachee} --ref-fasta ${params.ref} --ncbi-build ${params.assembly} --input-vcf ${sample_id}_filtered_PASS.vcf --output-maf ${sample_id}_vcf2MAF.maf
    """
}

process oncokb {
    tag "Precison Analysis on ${sample_id}"
    publishDir "${params.outdir}/9.AdvanceAnnot", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path("${sample_id}_vcf2MAF.maf")
    val(params.assembly)
    val(params.tumerType)
    val(params.OKBAPI)
    
    output:
    tuple val(sample_id), path("${sample_id}_annotat_oncokb.tsv"), emit: oncokbA

  
    script:
    """
    python /usr/src/app/oncokb-annotator/MafAnnotator.py -i ${sample_id}_vcf2MAF.maf -o ${sample_id}_annotat_oncokb.tsv -b ${params.OKBAPI} -a -d -t ${params.tumerType} -r ${params.assembly}
    """
}

process precision_analysis {
    tag "Precison Analysis on ${sample_id}"
    publishDir "${params.outdir}/10.Precision", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path("${sample_id}_amp_vcf_output_cosmic.tsv")
    tuple val(sample_id), path("${sample_id}_annotat_oncokb.tsv")

    output:
    tuple val(sample_id), path("${sample_id}_precision_final_output.xlsx"), emit: precision_out

    script:
    """
    python /usr/src/app/ref17/Validation_script/cosmic_okb.py --cosmic ${sample_id}_amp_vcf_output_cosmic.tsv --okb ${sample_id}_annotat_oncokb.tsv --output ${sample_id}_precision_final_output.xlsx
    """
}