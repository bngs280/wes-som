// Trimming ONT data
process trim_ONT {
    tag "Trimming on ${sample_id} (ONT)"
    publishDir "${params.outdir}/2.Trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(read_pairs_se)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: filt_out_ont

    when:
    params.platform == 'Nanopore' && params.library == 'Single'

    script:
    """
    gunzip -c ${read_pairs_se} | NanoFilt -q 10 -l 500 --headcrop 40 --tailcrop 20 | gzip > ${sample_id}_trimmed.fastq.gz
    """
}

// Validate Trimming ONT data
process validateTrimmedOutONT {
    tag "${sample_id}_TrimmedOutSEvalidation (ONT)"

    input:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz")

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: filt_out_validatedONT

    script:
    """
    python3 ../script/validate_fastp_outSE.py --input ${sample_id}_trimmed.fastq.gz
    """
}

// Trimming of Paired Data
process trim_paired {
    tag "Trimming on ${sample_id} (Paired-End)"
    publishDir "${params.outdir}/2.Trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(read_pairs_pe)

    output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz"), emit: fastp_out_pe, path("${sample_id}_fastp.html")

    when:
    params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI' && params.library == 'Paired' 
    
    script:
    """
    fastp --detect_adapter_for_pe -q 20 \
          -i ${read_pairs_pe[0]} \
          -I ${read_pairs_pe[1]} \
          -o ${sample_id}_R1_trimmed.fastq.gz \
          -O ${sample_id}_R2_trimmed.fastq.gz \
          --html ${sample_id}_fastp.html \
          --report_title "Quality Control for ${sample_id}"
    """
}

// Validate Paired Trimmed data
process validateTrimmedOutPE {
    tag "${sample_id}_TrimmedOutPEvalidation"

    input:
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz")

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: fastp_out_validatedPE

    script:
    """
    python3 ../script/validate_fastp_outPE.py --input1 ${sample_id}_trimmed_R1.fastq.gz --input2 ${sample_id}_trimmed_R2.fastq.gz
    """
}

// Trimming of Single-End data
process trim_single {
    tag "Trimming on ${sample_id} (Single-End)"
    publishDir "${params.outdir}/2.Trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(read_pairs_se)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: fastp_out_se, path("${sample_id}_fastp.html")

    when:
    params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI' || params.platform == 'ThermoFisher' && params.library == 'Single'

    script:
    """
    fastp -q 20 \
          -i ${read_pairs_se} \
          -o ${sample_id}_trimmed.fastq.gz \
          --html ${sample_id}_fastp.html \
          --report_title "Quality Control for ${sample_id}"
    """
}

// Validate Single-End Trimmed data
process validateTrimmedOutSE {
    tag "${sample_id}_TrimmedOutSEvalidation"

    input:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz")

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: fastp_out_validatedSE

    script:
    """
    python3 ../script/validate_fastp_outSE.py --input ${sample_id}_trimmed.fastq.gz
    """
}