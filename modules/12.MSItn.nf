process msiScoretn {
    tag "MSI Score Calculation on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/MSIScore", mode: 'copy'
    cpus 4

    input:
   tuple val(tumor_id), path(tumor_bam), path(tumor_bai),
          val(normal_id), path(normal_bam), path(normal_bai)
    val(params.ref)
    val(params.mantishg38)

    output:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_kmer_counts.txt"), path("${tumor_id}_vs_${normal_id}_kmer_counts_filtered.txt"), path("${tumor_id}_vs_${normal_id}.txt.status"), path("${tumor_id}_vs_${normal_id}.txt"), emit: msi_score

    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'

    script:
    """
    source activate p3919
    python /usr/src/app/MANTIS/mantis.py --bedfile ${params.mantishg38} --genome ${params.ref} -n ${normal_bam} -t ${tumor_bam} -o ${tumor_id}_vs_${normal_id}.txt
   
    """
}