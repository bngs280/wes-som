process ampclasstn {
    tag "AMP class on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/AnnotationAMPtn", mode: 'copy'
    cpus 4
    cache 'false'

    input:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered_PASS.vcf")
    val(params.assemblyAMP)

    output:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}.${params.assemblyAMP}_multianno.txt.cancervar"), path("${tumor_id}_vs_${normal_id}.${params.assemblyAMP}_multianno.txt"), path("${tumor_id}_vs_${normal_id}.${params.assemblyAMP}_multianno.txt.grl_p"), emit: vep_TSVamp

    script:
    """
    python /usr/src/app/ref17/annovarhg38/cancewarhg38/./CancerVar.py  -b ${params.assemblyAMP} -i ${tumor_id}_vs_${normal_id}_filtered_PASS.vcf --input_type=VCF -o ${tumor_id}_vs_${normal_id}
    """
}
