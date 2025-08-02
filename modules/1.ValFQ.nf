// validate Paired FASTQ files
process validateFastqPE {
    tag "${sample_id}_InputValidation (Paired-End)"

    input:
    tuple val(sample_id), path(read_pairs_pe)

    output:
    tuple val(sample_id), path(read_pairs_pe), emit: validated_samplesPE

    script:
    """
    python3 ../script/validate_fastqPE.py --input1 ${read_pairs_pe[0]} --input2 ${read_pairs_pe[1]}

    """
}

// validate Single FASTQ files
process validateFastqSE {
    tag "${sample_id}_InputValidation (Single-End)"

    input:
    tuple val(sample_id), path(read_pairs_se)

    output:
    tuple val(sample_id), path(read_pairs_se), emit: validated_samplesSE

    script:
    """
    python3 ../script/validate_fastqSE.py --input ${read_pairs_se} 

    """
}

// validate Nanopore FASTQ files
process validateFastqONT {
    tag "${sample_id}_InputValidation (ONT)"

    input:
    tuple val(sample_id), path(read_pairs_se)

    output:
    tuple val(sample_id), path(read_pairs_se), emit: validated_samplesONT

    script:
    """
    python3 ../script/validate_fastqSE.py --input ${read_pairs_se} 

    """
}