process VCFnorm {
    tag "Normalize vcf form TMB on ${sample_id}"
    publishDir "${params.outdir}/12.TMBScore", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(vcf)
    val(params.ref)
    
    output:
    tuple val(sample_id), path("${sample_id}.filtered_norm.vcf.gz"), emit: tmb_score

    script:
    """
    bcftools norm -m- -f ${params.ref} -o ${sample_id}.filtered_norm.vcf.gz ${vcf}
    """
}

process TMBScore {

    tag "TMB score calculation on ${sample_id}"
    publishDir "${params.outdir}/12.TMBScore", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(vcf)
    val(params.bed)
   
    output:
    tuple val(sample_id), path("${sample_id}.TMB_results.log"), emit: tmb_score

    script:
    """
    source activate pyTMB 
    python /usr/src/app/TMB/bin/pyTMB.py -i ${vcf} \
    --dbConfig /usr/src/app/TMB/config/annovar.yml \
	--varConfig /usr/src/app/TMB/config/mutect2.yml \
    --bed ${params.bed}  \
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
	--export > ${sample_id}.TMB_results.log
    """
}