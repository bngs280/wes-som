process anno_ampMergetn {
    tag "Merge annotation with AMP on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/Annotation", mode: 'copy'
    cpus 4
    cache 'false'

    input:
    //tuple val(sample_id), path(vcf)
    //val(params.assemblyAMP)
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered_PASS_annot_proccessed.tsv")
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}.hg38_multianno.txt.cancervar")


    output:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_amp_vep_output.tsv"), emit: vep_TSVf
    
    when:
    params.genome == 'hg38'

    script:
    """
    python /usr/src/app/ref17/Validation_script/merge_vep_amptn.py --vep ${tumor_id}_vs_${normal_id}_filtered_PASS_annot_proccessed.tsv --amp ${tumor_id}_vs_${normal_id}.hg38_multianno.txt.cancervar --output ${tumor_id}_vs_${normal_id}_amp_vep_output.tsv

    """
}

process anno_ampMerge19tn {
    tag "Merge annotation with AMP on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/Annotation", mode: 'copy'
    cpus 4

    input:
    //tuple val(sample_id), path(vcf)
    //val(params.assemblyAMP)
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered_PASS_annot_proccessed.tsv")
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}.hg19_multianno.txt.cancervar")


    output:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_amp_vep_output.tsv"), emit: vep_TSVfhg19

    when:
    params.genome == 'hg19'
    script:
    """
    python /usr/src/app/ref17/Validation_script/merge_vep_amptn.py --vep ${tumor_id}_vs_${normal_id}_filtered_PASS_annot_proccessed.tsv --amp ${tumor_id}_vs_${normal_id}.hg19_multianno.txt.cancervar --output ${tumor_id}_vs_${normal_id}_amp_vep_output.tsv

    """
}

process vcf_annAMPtn {
    tag "VCF information into annotation on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/Annotation", mode: 'copy'
    cpus 4

    input:
    //tuple val(sample_id), path(vcf)
    //val(params.assemblyAMP)
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered_PASS.vcf")
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_amp_vep_output.tsv")


    output:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_amp_vcf_output.tsv"), emit: vcf_vep_TSV

    when:
    params.genome == 'hg38'
    script:
    """
    python /usr/src/app/ref17/Validation_script/vcf_amp_sampleidtn.py --vcf ${tumor_id}_vs_${normal_id}_filtered_PASS.vcf --tsv ${tumor_id}_vs_${normal_id}_amp_vep_output.tsv --output ${tumor_id}_vs_${normal_id}_amp_vcf_output.tsv

    """
}

process vcf_annAMPhg19tn {
    tag "VCF information into annotation on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/Annotation", mode: 'copy'
    cpus 4

    input:
    //tuple val(sample_id), path(vcf)
    //val(params.assemblyAMP)
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered_PASS.vcf")
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_amp_vep_output.tsv")


    output:
    tuple val(sample_id), path("${tumor_id}_vs_${normal_id}_amp_vcf_output.tsv"), emit: vcf_vep_TSV

    when:
    params.genome == 'hg19'
    script:
    """
    python /usr/src/app/ref17/Validation_script/vcf_amp_sampleidtn.py --vcf ${tumor_id}_vs_${normal_id}_filtered_PASS.vcf --tsv ${tumor_id}_vs_${normal_id}_amp_vep_output.tsv --output ${tumor_id}_vs_${normal_id}_amp_vcf_output.tsv

    """
}
