process SV_somtn {
  
    tag "Structural Variants on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/SV_somatic", mode: 'copy'
    cpus 4

    input:
    tuple val(tumor_id), path(tumor_bam), path(tumor_bai),
          val(normal_id), path(normal_bam), path(normal_bai)
    val(params.ref)

    output:
    path "${tumor_id}_vs_${normal_id}", emit: sv_TUmsomatic

    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'

    script:
    """
    # Activate the environment for Manta
    source activate p2Manta

    # Create a temporary directory in the work directory for Manta to run
    mkdir ${tumor_id}_vs_${normal_id}

    # Run Manta config and workflow
    configManta.py --normalBam ${normal_bam} --tumorBam ${tumor_bam} --referenceFasta ${params.ref} --runDir ${tumor_id}_vs_${normal_id}
    ${tumor_id}_vs_${normal_id}/runWorkflow.py
    """
}