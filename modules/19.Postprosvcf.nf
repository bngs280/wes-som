process anno_ampMerge {
    tag "Merge annotation with AMP on ${sample_id}"
    publishDir "${params.outdir}/8.Annotation", mode: 'copy'
    cpus 4
    cache 'false'

    input:
    //tuple val(sample_id), path(vcf)
    //val(params.assemblyAMP)
    tuple val(sample_id), path("${sample_id}_filtered_PASS_annot_proccessed.tsv")
    tuple val(sample_id), path("${sample_id}.hg38_multianno.txt.cancervar")


    output:
    tuple val(sample_id), path("${sample_id}_amp_vep_output.tsv"), emit: vep_TSVf
    
    when:
    params.genome == 'hg38'

    script:
    """
    python /usr/src/app/ref17/Validation_script/merge_vep_amp.py --vep ${sample_id}_filtered_PASS_annot_proccessed.tsv --amp ${sample_id}.hg38_multianno.txt.cancervar --output ${sample_id}_amp_vep_output.tsv

    """
}

process anno_ampMerge19 {
    tag "Merge annotation with AMP on ${sample_id}"
    publishDir "${params.outdir}/8.Annotation", mode: 'copy'
    cpus 4

    input:
    //tuple val(sample_id), path(vcf)
    //val(params.assemblyAMP)
    tuple val(sample_id), path("${sample_id}_filtered_PASS_annot_proccessed.tsv")
    tuple val(sample_id), path("${sample_id}.hg19_multianno.txt.cancervar")


    output:
    tuple val(sample_id), path("${sample_id}_amp_vep_output.tsv"), emit: vep_TSVfhg19

    when:
    params.genome == 'hg19'
    script:
    """
    python /usr/src/app/ref17/Validation_script/merge_vep_amp.py --vep ${sample_id}_filtered_PASS_annot_proccessed.tsv --amp ${sample_id}.hg19_multianno.txt.cancervar --output ${sample_id}_amp_vep_output.tsv

    """
}

process vcf_annAMP {
    tag "VCF information into annotation on ${sample_id}"
    publishDir "${params.outdir}/8.Annotation", mode: 'copy'
    cpus 4

    input:
    //tuple val(sample_id), path(vcf)
    //val(params.assemblyAMP)
    tuple val(sample_id), path("${sample_id}_filtered_PASS.vcf")
    tuple val(sample_id), path("${sample_id}_amp_vep_output.tsv")


    output:
    tuple val(sample_id), path("${sample_id}_amp_vcf_output.tsv"), emit: vcf_vep_TSV

    when:
    params.genome == 'hg38'
    script:
    """
    python /usr/src/app/ref17/Validation_script/vcf_amp_sampleid.py --vcf ${sample_id}_filtered_PASS.vcf --tsv ${sample_id}_amp_vep_output.tsv --output ${sample_id}_amp_vcf_output.tsv

    """
}

process vcf_annAMPhg19 {
    tag "VCF information into annotation on ${sample_id}"
    publishDir "${params.outdir}/8.Annotation", mode: 'copy'
    cpus 4

    input:
    //tuple val(sample_id), path(vcf)
    //val(params.assemblyAMP)
    tuple val(sample_id), path("${sample_id}_filtered_PASS.vcf")
    tuple val(sample_id), path("${sample_id}_amp_vep_output.tsv")


    output:
    tuple val(sample_id), path("${sample_id}_amp_vcf_output.tsv"), emit: vcf_vep_TSV

    when:
    params.genome == 'hg19'
    script:
    """
    python /usr/src/app/ref17/Validation_script/vcf_amp_sampleid.py --vcf ${sample_id}_filtered_PASS.vcf --tsv ${sample_id}_amp_vep_output.tsv --output ${sample_id}_amp_vcf_output.tsv

    """
}