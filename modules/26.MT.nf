process MTcall {
    tag "MT Variants on ${sample_id}"
    publishDir path: "${params.outdir}/15.MTvariant", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_MT.vcf.gz"), path("${sample_id}_MT.log"), emit: mtCall

    script:
    """
    /usr/src/app/mutserveTool/./mutserve call --reference /usr/src/app/mutserveTool/rCRS.fasta --output ${sample_id}_MT.vcf.gz --threads ${task.cpus} ${bam} > ${sample_id}_MT.log
    """
}