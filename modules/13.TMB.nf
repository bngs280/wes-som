process VCFnormtn {
    tag "Normalize vcf form TMB on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/TMBScore", mode: 'copy'
    cpus 4

    input:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered_PASS.vcf")
    val(params.ref)
    
    output:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}.filtered_norm.vcf.gz"), emit: tmb_score

    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'

    script:
    """
    bcftools norm -m- -f ${params.ref} -o ${tumor_id}_vs_${normal_id}.filtered_norm.vcf.gz ${tumor_id}_vs_${normal_id}_filtered_PASS.vcf
    """
}

process TMBScoretn {

    tag "TMB score calculation on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/TMBScore", mode: 'copy'
    cpus 4

    input:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}.filtered_norm.vcf.gz")
    val(params.bed)
   
    output:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}.TMB_results.log"), emit: tmb_score

    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'

    script:
    """
    source activate pyTMB 
    python /usr/src/app/TMB/bin/pyTMB.py -i ${tumor_id}_vs_${normal_id}.filtered_norm.vcf.gz \
    --dbConfig /usr/src/app/TMB/config/annovar.yml \
	--varConfig /usr/src/app/TMB/config/mutect2.yml \
    --bed ${params.bed}  \
    --sample ${tumor_id} \
	--vaf 0.05 \
	--maf 0.001 \
	--minDepth 50 \
	--minAltDepth 2 \
	--filterLowQual \
	--filterSplice \
	--filterNonCoding \
	--filterSyn \
	--filterPolym \
	--polymDb 1k,gnomad \
	--cancerDb cosmic \
	--export > ${tumor_id}_vs_${normal_id}.TMB_results.log
    """
}