#!/usr/bin/env nextflow
// last update: 2/12/2024
// workd till error intergation
// sudo ~/nextflow run main.nf -params-file config.json
// sudo docker run -v /home/ubuntu/Doc/nf/input:/usr/src/app/input -v /home/ubuntu/Doc/nf/outputSec:/usr/src/app/output -v /home/ubuntu/Doc/nf/nextflow.config:/usr/src/app/nextflow.config -v /home/ubuntu/Doc/nf/config.json:/usr/src/app/config.json rgenxtools:latest
// nextflow run main.nf -with-docker -with-report 1main_1_6_report -with-dag 1main_1_6_DAG.png -with-trace 1main_1_6_trace.txt -with-timeline  1main_1_6_timeline.html
// last update: 21/01/2025 , can generate file with vcf data as well worked for tum_norm mode but tum_norm mode vcf file not consuming by annotation proccess propelry, you have to furtehr wrk on this
// last update: 04022025 worked with tum-norm, TMB, MSI, CNV, SV
nextflow.enable.dsl=2

import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

// Get the current date and time
def now = LocalDateTime.now()
def dateFormat = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss")
def formattedDate = now.format(dateFormat)

// Define ANSI color codes for styling
def cyan = "\033[0;36m"
def green = "\033[0;32m"
def yellow = "\033[1;33m"
def bold = "\033[1m"
def reset = "\033[0m"

// Company name and stylish intro message with current date and time
println """
${bold}${cyan}============================================================
    WELCOME TO THE WES ONCOLOGY PIPELINE!
    Presented by: ${yellow}VGENOMICS INDIA PVT LTD${cyan}
============================================================
   Cutting-edge genomic analysis for precision oncology.
   Powered by ${green}VgenX CLI${cyan} 
============================================================${reset}
"""

// Capture start time for tracking the duration
def startTime = LocalDateTime.now()


//params.outdir = '/home/ubuntu/Doc/nf/illPEhg38' //change output dir
//params.data_dir = '/home/ubuntu/Doc/nf/input' //change input dir
params.genome = params.genome ?: 'hg38'  // Default to hg38 if not specified
params.platform = params.platform ?: 'Illumina'  // Default to Illumina paired-end
params.library = params.library ?: 'Paired' // Default paired end
params.mode = params.mode ?: 'TumorNormal'
params.tumerType = params.tumerType ?: 'All tumor'
params.OKBAPI = '5608b38e-3f27-4b7f-8e4a-13b42a7b7f41'
// Ref Genomes
params.refhg38 = '/usr/src/app/ref17/hg38/hg381_22XYM/Homo_sapiens_assembly38cleaned.fasta' //chr1_22 X_Y_M only
params.refhg37 = '/usr/src/app/ref17/hg19/hg19122XYM/hg19122XYM.fa' //chr1_22 X_Y_M only
params.bedhg38 = '/usr/src/app/ref17/BED/hg38_exome.bed'
params.bedhg37 = '/usr/src/app/ref17/BED/hg37_exome.bed'
params.gnomad38 = '/usr/src/app/ref17/mutect2/hg38/af-only-gnomad.hg38.vcf.gz'
params.gnomad37 = '/usr/src/app/ref17/mutect2/hg19/af-only-gnomad.hg38Tohg19.vcf.gz'
params.db1000g38 = '/usr/src/app/ref17/mutect2/hg38/1000g_pon.hg38.vcf.gz'
params.db1000g37 = '/usr/src/app/ref17/mutect2/hg19/1000g_pon.hg38Tohg19.vcf.gz'

/// MSI score calculation model
params.msimodelhg38 = '/usr/src/app/msisensor2/models_hg38'
params.msimodelhg19 = '/usr/src/app/msisensor2/models_hg19_GRCh37'
params.mantishg38 = '/usr/src/app/ref17/Mantis_MSI/hg38_loci.bed'
params.mantishg19 = '/usr/src/app/ref17/Mantis_MSI/hg38_loci.bed' // need to donload

//// VEP input Databases Hg38
params.dir = '/usr/src/app/vepC/cache'
params.dirPlugin = '/usr/src/app/ensembl-vep/Plugins'
params.dbnsfp38 = '/usr/src/app/vepDB/DBs/hg38/dbNSFP4.7a_grch38.gz'
params.LoFtool = '/usr/src/app/ensembl-vep/Plugins/LoFtool_scores.txt'
params.CADDsnv38 = '/usr/src/app/vepDB/DBs/hg38/whole_genome_SNVs.tsv.gz'
params.CADDindels38 = '/usr/src/app/vepDB/DBs/hg38/gnomad.genomes.r4.0.indel_inclAnno.tsv.gz'
params.dbscSNV38 = '/usr/src/app/vepDB/DBs/hg38/dbscSNV1.1_GRCh38.txt.gz'

//// VEP input Databases Hg37
params.dbnsfp37 = '/usr/src/app/vepDB/DBs/hg19/dbNSFP4.7a_grch37.gz'
params.CADDsnv37 = '/usr/src/app/vepDB/DBs/hg19/whole_genome_SNVs_hg37.tsv.gz'
params.CADDindels37 = '/usr/src/app/vepDB/DBs/hg19/gnomad.genomes-exomes.r4.0.indel_inclAnno_hg37.tsv.gz'
params.dbscSNV37 = '/usr/src/app/vepDB/DBs/hg19/dbscSNV1.1_GRCh37.txt.gz'
params.MaxEntmain = '/usr/src/app/vepDB/DBs/fordownload'

params.maphg38 = '/usr/src/app/ref17/CNV_delly_data/Hg38.map'
params.maphg37 = '/usr/src/app/ref17/CNV_delly_data/hg37/Hg37.map'
params.CNVbaselinecontra = '/usr/src/app/ref17/CNV_delly_data/Contra_CNV_baseline/L1140225.pooled2_TRIM0.2.txt'
//params.CNVbaselinecontraag = '/usr/src/app/ref17/CNV_delly_data/Contra_CNV_baseline/Agilent_SureSelect_All_Exon_50Mb_ICGC.baseline.txt'
params.cosmiccodhg38 = '/usr/src/app/ref17/COSMIC/cosmic_coding_hg38.tsv'
params.cosmiccodhg37 = '/usr/src/app/ref17/COSMIC/cosmic_codingv100_hg37.tsv'
// Genome reference selection
if (params.genome == 'hg38') {
    params.ref = params.refhg38
    params.gnomad = params.gnomad38
    params.db1000g = params.db1000g38
    params.msimodel = params.msimodelhg38
    params.assembly = 'GRCh38'
    params.assemblyAMP = 'hg38'
    params.cachee = params.dir
    params.dirplugin = params.dirPlugin
    params.dbNSFP = params.dbnsfp38
    params.loftool = params.LoFtool
    params.CADDsnv = params.CADDsnv38
    params.CADDindel = params.CADDindels38
    params.dbscSNV = params.dbscSNV38
    params.MaxEnt = params.MaxEntmain
    params.mapCNV = params.maphg38
    params.baselineContra = params.CNVbaselinecontra
    params.cosmiccod = params.cosmiccodhg38
    // Use custom BED file if provided, else use default
    //params.bed = params.bed ?: params.bedhg38
    params.bed = params.bedhg38

} else if (params.genome == 'hg19') {
    params.ref = params.refhg37
    params.gnomad = params.gnomad37
    params.db1000g = params.db1000g37
    params.msimodel = params.msimodelhg19
    params.assembly = 'GRCh37'
    params.assemblyAMP = 'hg19'
    params.cachee = params.dir
    params.dirplugin = params.dirPlugin
    params.dbNSFP = params.dbnsfp37
    params.loftool = params.LoFtool
    params.CADDsnv = params.CADDsnv37
    params.CADDindel = params.CADDindels37
    params.dbscSNV = params.dbscSNV37
    params.MaxEnt = params.MaxEntmain
    params.mapCNV = params.maphg37
    params.baselineContra = params.CNVbaselinecontra
    params.cosmiccod = params.cosmiccodhg37
      // Use custom BED file if provided, else use default
    //params.bed = params.bed ?: params.bedhg37
    params.bed = params.bedhg37

} else {
    error "Unsupported genome reference: ${params.genome}. Please use 'hg38' or 'hg19'."
}

workflow {
 
    // Define channels for input data
    def read_pairs_pe = Channel
        .fromFilePairs("${params.data_dir}/*_R{1,2}.fastq.gz", size: 2, hidden: true)
        .ifEmpty { 
            if ((params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired') {
                error """ 
                No valid FASTQ files found in the input directory (${params.data_dir}).
                Please ensure the following:
                - Paired-end files: *_R1.fastq.gz and *_R2.fastq.gz
                """
            }
        }
    read_pairs_pe.view()

    def read_pairs_se = Channel
        .fromFilePairs("${params.data_dir}/*.fastq.gz", size: 1)
        .ifEmpty { 
            if ((params.platform == 'Illumina' || params.platform == 'Nanopore' || params.platform == 'BGI' || params.platform == 'MGI' || params.platform == 'ThermoFisher') && params.library == 'Single') {
                error """ 
                No valid FASTQ files found in the input directory (${params.data_dir}).
                Please ensure the following:
                - Single-end files: *.fastq.gz (without _R1/_R2 suffix)
                """
            }
        }
    
    // Check for processes to skip from config
    def skipProcesses = params.skip_processes ?: []

    if (!params.library) {
        error """
        Library type not specified.
        Please specified: Single or Paired 
        """
    } 
 
    if ((params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired') {
        if (params.mode == 'TumorOnly') {
            // Paired-end pipeline for Illumina and BGI MGI PE
            // Process steps
            /// Step 1. QC with FastQC
            validated_samplesPE = validateFastqPE(read_pairs_pe)
            quality_check(validated_samplesPE)

            trim_paired(read_pairs_pe)
            fastp_out_pe = trim_paired.out.fastp_out_pe
        
            fastp_out_validatedPE = validateTrimmedOutPE(fastp_out_pe)

            // Step: BAM Analysis (skip if "Alignment" is mentioned in skip_processes)
            if (!skipProcesses.contains('Alignment')) {
                sorted_bams = align_paired(fastp_out_validatedPE, params.ref)
                sorted_bams = align_paired.out.sorted_bams_pe
            } else {
                println "Skipping Alignment analysis as specified in config."
                sorted_bams = null // Set sorted_bams to null to handle dependency
            }
        } else if (params.mode == 'TumorNormal') {
            def sample_manifest = Channel
                .fromPath("${params.data_dir}/manifest.tsv")
                .splitCsv(header: true, sep: '\t')
            sample_manifest.view()
           
            def tumor_normal_pairs = sample_manifest
                .collate(2) // Group rows in pairs of 2
                .map { pair ->
                    def tumor_row = pair.find { it.sample_type == 'tumor' }
                    def normal_row = pair.find { it.sample_type == 'normal' }
                    if (tumor_row && normal_row) {
                        tuple(tumor_row.sample_id, tumor_row.sample_type, normal_row.sample_id, normal_row.sample_type)
                    } else {
                        error "Could not find a tumor-normal pair in: ${pair}"
                    }
                }
            tumor_normal_pairs.view()

            // Validate samples (making sure we get fastq files)
            def validated_samplesTN = tumor_normal_pairs.map { pair ->
                def tumor_fastq1 = file("${params.data_dir}/${pair[0]}_R1.fastq.gz")
                def tumor_fastq2 = file("${params.data_dir}/${pair[0]}_R2.fastq.gz")
                def normal_fastq1 = file("${params.data_dir}/${pair[2]}_R1.fastq.gz")
                def normal_fastq2 = file("${params.data_dir}/${pair[2]}_R2.fastq.gz")
              
                if (!tumor_fastq1.exists() || !tumor_fastq2.exists()) {
                    error "Missing tumor FASTQ files for sample ID: ${pair[0]}"
                }
                if (!normal_fastq1.exists() || !normal_fastq2.exists()) {
                    error "Missing normal FASTQ files for sample ID: ${pair[2]}"
                }
                tuple(pair[0], tumor_fastq1, tumor_fastq2, pair[2], normal_fastq1, normal_fastq2)
            }
            // Add this statement to print the content of the channel
            
            validated_samplesTN.view()
            /// Step 1. QC with FastQC
            //validated_samplesPE = validateFastqPE(validated_samplesTN)
            //quality_check(validated_samplesTN)
            fastqc_outTN = quality_checkTN(validated_samplesTN)
            trimmed_fastqs_tn = fastpTumorNormal(validated_samplesTN)
            //trimmed_fastqs_tn = trim_paired.out.trimmed_fastqs_tn
        
            //fastp_out_validatedPE = validateTrimmedOutPE(fastp_out_pe)
            //aligned_bams_tn = align_tumor_normal(trimmed_fastqs_tn, params.ref)
            sorted_bams = align_tumor_normal(trimmed_fastqs_tn, params.ref)
            sorted_bams = align_tumor_normal.out.sorted_bams_tn
        } else {
            error "Unsupported mode: ${params.mode}. Please choose from 'TumorNormal', 'TumorOnly'."
        }
 
    } else if ((params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI' || params.platform == 'ThermoFisher') && params.library == 'Single') {
        // Single-end pipeline for Illumina SE, Nanopore, and BGI MGI SE
        validated_samplesSE = validateFastqSE(read_pairs_se)
        quality_check(validated_samplesSE)

        trim_single(read_pairs_se)
        fastp_out_se = trim_single.out.fastp_out_se
        
        fastp_out_validatedSE = validateTrimmedOutSE(fastp_out_se)

        // Step: BAM Analysis (skip if "Alignment" is mentioned in skip_processes)
        if (!skipProcesses.contains('Alignment')) {
            sorted_bams = align_single(fastp_out_validatedSE, params.ref)
            sorted_bams = align_single.out.sorted_bams_se
        } else {
            println "Skipping Alignment analysis as specified in config."
            sorted_bams = null // Set sorted_bams to null to handle dependency
        }

    } else if ((params.platform == 'Nanopore') && params.library == 'Single') {
        // Single-end pipeline for Nanopore
        validated_samplesONT = validateFastqONT(read_pairs_se)
        quality_checkONT(validated_samplesONT)

        trim_ONT(read_pairs_se)
        filt_out_ont = trim_ONT.out.filt_out_ont
        
        filt_out_validatedONT = validateTrimmedOutONT(filt_out_ont)

        // Step: BAM Analysis (skip if "Alignment" is mentioned in skip_processes)
        if (!skipProcesses.contains('Alignment')) {
            sorted_bams = align_ONT(filt_out_validatedONT, params.ref)
            sorted_bams = align_ONT.out.sorted_bams_ont
        } else {
            println "Skipping Alignment analysis as specified in config."
            sorted_bams = null // Set sorted_bams to null to handle dependency
        }
    } else {
        error "Unsupported platform: ${params.platform}. Please choose from 'Illumina', 'MGI', 'Nanopore', 'BGI', 'ThermoFisher'."
    }

    // Continue with the rest of your workflow using the `sorted_bams_pe` or `sorted_bams_se`
    if (!skipProcesses.contains('BAMvalidation') && sorted_bams != null) {
        if (params.mode == 'TumorNormal') {
            markdup_bams = markDupTN(sorted_bams)
            valid_markdup_bams = mdBAM_indexTN(markdup_bams)
        } else {
            bam_validation = validateBAM(sorted_bams)
            markdup_bams = markDup(bam_validation)
            valid_markdup_bams = mdvalidate(markdup_bams)
            valid_markdup_bams = mdBAM_index(valid_markdup_bams)
        }
       
    } else {
        println "Skipping BAM file validation analysis as specified in config."
        valid_markdup_bams = null // Set sorted_bams to null to handle dependency
    }

     // Continue variant calling, filtering, etc. with either SE or PE BAMs
    // Step: Variant Calling (skip if "VariantCalling" is mentioned in skip_processes)
    if (!skipProcesses.contains('VariantCalling') && valid_markdup_bams != null) {
        if ((params.platform == 'Nanopore') && params.library == 'Single') {
            // Variant calling for Nanopore platform
            final_vcfs = somVarCall_ont(valid_markdup_bams, params.ref, params.bed)
            //valid_vcfs = validatevcf(raw_vcfs)
            //filtered_vcfs = FilterMT(raw_vcfs, params.ref)
            //final_vcfs = KeepPASS(filtered_vcfs)
            valid_final_vcfs = validateFinalvcf(final_vcfs)
        } else if (params.mode == 'TumorNormal') {
            raw_vcfs = somVarCall_tumor_normal(valid_markdup_bams, params.ref, params.bed, params.gnomad, params.db1000g)
            valid_vcfs = validatevcftn(raw_vcfs)
            filtered_vcfs = FilterMTtn(raw_vcfs, params.ref)
            final_vcfs = KeepPASStn(filtered_vcfs)
            valid_final_vcfs = validateFinalvcftn(final_vcfs)  

        }else {
            raw_vcfs = somVarCall(valid_markdup_bams, params.ref, params.bed, params.gnomad, params.db1000g)
            valid_vcfs = validatevcf(raw_vcfs)
            filtered_vcfs = FilterMT(raw_vcfs, params.ref)
            final_vcfs = KeepPASS(filtered_vcfs)
            valid_final_vcfs = validateFinalvcf(final_vcfs)
        }
        
    } else {
        println "Skipping Variant Calling as specified in config."
        final_vcfs = null // Set final_vcfs to null to handle dependency
    }
    
    // Step: Annotation Analysis (skip if "Annotation" is mentioned in skip_processes)
    if (!skipProcesses.contains('Annotation') && final_vcfs != null) {
        if (params.mode == 'TumorNormal') {
            annotVEP_vcfs = annotationtn(final_vcfs, params.cachee, params.dirplugin, params.ref, params.assembly, params.dbNSFP, params.loftool, params.CADDsnv, params.CADDindel, params.dbscSNV)
            vep_TSV = annot_processingtn(annotVEP_vcfs)
            vep_TSVamp = ampclasstn(final_vcfs, params.assemblyAMP)
            if (params.genome == 'hg38') {
                vep_TSVf = anno_ampMergetn(vep_TSV, vep_TSVamp)
            } else {
                vep_TSVfhg19 = anno_ampMerge19tn(vep_TSV, vep_TSVamp)
            }
            if (params.genome == 'hg38') {
                vcf_vep_TSV = vcf_annAMPtn(final_vcfs, vep_TSVf)
            } else {
                vcf_vep_TSV = vcf_annAMPhg19tn(final_vcfs, vep_TSVfhg19)
            }
        } else {
            annotVEP_vcfs = annotation(final_vcfs, params.cachee, params.dirplugin, params.ref, params.assembly, params.dbNSFP, params.loftool, params.CADDsnv, params.CADDindel, params.dbscSNV)
            vep_TSV = annot_processing(annotVEP_vcfs)
            vep_TSVamp = ampclass(final_vcfs, params.assemblyAMP)
            if (params.genome == 'hg38') {
                vep_TSVf = anno_ampMerge(vep_TSV, vep_TSVamp)
            } else {
                vep_TSVfhg19 = anno_ampMerge19(vep_TSV, vep_TSVamp)
            }
            if (params.genome == 'hg38') {
                vcf_vep_TSV = vcf_annAMP(final_vcfs, vep_TSVf)
            } else {
                vcf_vep_TSV = vcf_annAMPhg19(final_vcfs, vep_TSVfhg19)
            }
            if (params.genome == 'hg38') {
                vcf_vep_TSV1 = cosmic_hg38(params.cosmiccod, vcf_vep_TSV)
            } else {
                vcf_vep_TSV1 = cosmic_hg19(params.cosmiccod, vcf_vep_TSV)
            }
            
            vcf_maf = vcf2maf(final_vcfs, params.cachee, params.ref, params.assembly)
            
            
            oncokbA = oncokb(vcf_maf, params.OKBAPI, params.assembly, params.tumerType)
            
            precision_out = precision_analysis(vcf_vep_TSV1, oncokbA)
           
        }
        
    } else {
        println "Skipping Annotation (annotVEP) because either somVarCall or annotation process is skipped."
    }
   
    // Step: MSI Analysis (skip if "MSI" is mentioned in skip_processes)
    if (!skipProcesses.contains('MSI') && valid_markdup_bams != null) {
        if (params.mode == 'TumorNormal') {
            msi_score = msiScoretn(valid_markdup_bams, params.ref, params.mantishg38)
        } else {
            msi_score = msiScore(valid_markdup_bams, params.msimodel)
        }
        
    } else {
        println "Skipping MSI analysis as specified in config."
    }
    
    // Step: TMB Analysis (skip if "TMB" is mentioned in skip_processes)
    if (!skipProcesses.contains('TMB') && final_vcfs != null) {
        if  (params.mode == 'TumorNormal') {
            tmb_score = VCFnormtn(final_vcfs, params.ref)
            tmb_score = TMBScoretn(tmb_score, params.bed)
        } else {
            tmb_score = VCFnorm(final_vcfs, params.ref)
            tmb_score = TMBScore(tmb_score, params.bed)
        }
    
    } else {
        println "Skipping TMB analysis as specified in config."
    }
    
    // Step: SV Analysis (skip if "SV" is mentioned in skip_processes)
    if (!skipProcesses.contains('SV') && valid_markdup_bams != null) {
        try {
            if (params.mode == 'TumorNormal') {
                sv_TUmsomatic = SV_somtn(valid_markdup_bams, params.ref)
            } else {
                sv_TUmsomatic = SV_som(valid_markdup_bams, params.ref)
            }
            
        } catch (Exception e) {
            println "Error in SV step: ${e.message}"
        }
    } else {
        println "Skipping SV analysis as specified in config."
    }   

    /// Step 13. CNV Somatic Tumor only mode
    if (!skipProcesses.contains('CNV') && valid_markdup_bams != null) {
        if (params.mode == 'TumorNormal') {
            cnvSom = CNV_somtn(valid_markdup_bams, params.ref, params.bed)
        } else {
            cnvSom = CNV_som(valid_markdup_bams, params.ref, params.bed, params.baselineContra)

        }
    } else {
        println "Skipping CNV analysis as specified in config."
    }
    
    // Step 14. MT calling 
    if (!skipProcesses.contains('MT') && valid_markdup_bams != null) {
        mtCall = MTcall(valid_markdup_bams)
    } else {
        println "Skipping MT analysis as specified in config."
    }

}



// Process to validate FASTQ files using validate_fastq.py
process validateFastqPE {
    tag "${sample_id}_InputValidation (Paired-End)"

    input:
    tuple val(sample_id), path(read_pairs_pe)

    output:
    tuple val(sample_id), path(read_pairs_pe), emit: validated_samplesPE

    script:
    """
    python3 /usr/src/app/ref17/Validation_script/validate_fastqPE.py --input1 ${read_pairs_pe[0]} --input2 ${read_pairs_pe[1]}

    """
}
// Process to validate Single FASTQ files using validate_fastq.py
process validateFastqSE {
    tag "${sample_id}_InputValidation (Single-End)"

    input:
    tuple val(sample_id), path(read_pairs_se)

    output:
    tuple val(sample_id), path(read_pairs_se), emit: validated_samplesSE

    script:
    """
    python3 /usr/src/app/ref17/Validation_script/validate_fastqSE.py --input ${read_pairs_se} 

    """
}
// Process to validate Single FASTQ files using validate_fastq.py
process validateFastqONT {
    tag "${sample_id}_InputValidation (ONT)"

    input:
    tuple val(sample_id), path(read_pairs_se)

    output:
    tuple val(sample_id), path(read_pairs_se), emit: validated_samplesONT

    script:
    """
    python3 /usr/src/app/ref17/Validation_script/validate_fastqSE.py --input ${read_pairs_se} 

    """
}
process quality_check {
    tag "Quality Checking on ${sample_id}"
    publishDir "${params.outdir}/1.QC", mode: 'copy'

    input:
    tuple val(sample_id), path(read_files)

    output:
    path "${sample_id}", emit: fastqc_out

    script:
    """
    mkdir ${sample_id}
    fastqc -o ${sample_id} -f fastq -q ${read_files}
    """
}
process quality_checkONT {
    tag "Quality Checking on ${sample_id} (ONT)"
    publishDir "${params.outdir}/1.QC", mode: 'copy'

    input:
    tuple val(sample_id), path(read_files)

    output:
    path "${sample_id}", emit: fastqc_ont_out

    script:
    """
    mkdir ${sample_id}
    NanoPlot --fastq ${read_files} -o ${sample_id}
    """
}

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

// Process to validate FASTQ files using validate_fastq.py
process validateTrimmedOutONT {
    tag "${sample_id}_TrimmedOutSEvalidation (ONT)"

    input:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz")

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: filt_out_validatedONT

    script:
    """
    python3 /usr/src/app/ref17/Validation_script/validate_fastp_outSE.py --input ${sample_id}_trimmed.fastq.gz
    """
}
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
// Process to validate FASTQ files using validate_fastq.py
process validateTrimmedOutPE {
    tag "${sample_id}_TrimmedOutPEvalidation"

    input:
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz")

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: fastp_out_validatedPE

    script:
    """
    python3 /usr/src/app/ref17/Validation_script/validate_fastp_outPE.py --input1 ${sample_id}_trimmed_R1.fastq.gz --input2 ${sample_id}_trimmed_R2.fastq.gz
    """
}

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

// Process to validate FASTQ files using validate_fastq.py
process validateTrimmedOutSE {
    tag "${sample_id}_TrimmedOutSEvalidation"

    input:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz")

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: fastp_out_validatedSE

    script:
    """
    python3 /usr/src/app/ref17/Validation_script/validate_fastp_outSE.py --input ${sample_id}_trimmed.fastq.gz
    """
}
process quality_checkTN {
    tag "QCTN on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/QCTN", mode: 'copy'
    input:
    tuple val(tumor_id), path(tumor_fastq1), path(tumor_fastq2),
          val(normal_id), path(normal_fastq1), path(normal_fastq2)
    
    output:
    path "${tumor_id}_${normal_id}", emit: fastqc_outTN

    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'
    
    script:
    """
    mkdir ${tumor_id}_${normal_id}
    fastqc -o ${tumor_id}_${normal_id} -f fastq -q ${tumor_fastq1} ${tumor_fastq2} ${normal_fastq1} ${normal_fastq2}
    """
}

process fastpTumorNormal {
    tag "FASTP on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/TrimmingTN", mode: 'copy'

    input:
    tuple val(tumor_id), path(tumor_fastq1), path(tumor_fastq2), 
          val(normal_id), path(normal_fastq1), path(normal_fastq2)

    output:
    tuple val(tumor_id), path("${tumor_id}_trimmed_R1.fastq.gz"), path("${tumor_id}_trimmed_R2.fastq.gz"),
          val(normal_id), path("${normal_id}_trimmed_R1.fastq.gz"), path("${normal_id}_trimmed_R2.fastq.gz"), path("${tumor_id}_fastp.html"), path("${normal_id}_fastp.html"), emit: trimmed_fastqs_tn

    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'

    script:
    """
    fastp -q 20 -i ${tumor_fastq1} -I ${tumor_fastq2} -o ${tumor_id}_trimmed_R1.fastq.gz -O ${tumor_id}_trimmed_R2.fastq.gz --html ${tumor_id}_fastp.html --report_title "Quality Control for ${tumor_id}"
    fastp -q 20 -i ${normal_fastq1} -I ${normal_fastq2} -o ${normal_id}_trimmed_R1.fastq.gz -O ${normal_id}_trimmed_R2.fastq.gz --html ${normal_id}_fastp.html --report_title "Quality Control for ${normal_id}"
    """
}

process align_tumor_normal {
    tag "Alignment of ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/AlignmentTN", mode: 'copy'
    cpus 8
    memory '16 GB'

    input:
    tuple val(tumor_id), path(tumor_fastq1), path(tumor_fastq2),
          val(normal_id), path(normal_fastq1), path(normal_fastq2)
    val(params.ref)

    output:
    tuple val(tumor_id), path("${tumor_id}_sorted.bam"),
          val(normal_id), path("${normal_id}_sorted.bam"), emit: sorted_bams_tn
    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'

    script:
    """
    # Align tumor FASTQ
    bwa mem -t ${task.cpus} -M \
            -R "@RG\\tID:${tumor_id}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${tumor_id}" \
            ${params.ref} ${tumor_fastq1} ${tumor_fastq2} | samtools view -@ ${task.cpus} -b - | samtools sort -o ${tumor_id}_sorted.bam

    # Align normal FASTQ
    bwa mem -t ${task.cpus} -M \
            -R "@RG\\tID:${normal_id}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${normal_id}" \
            ${params.ref} ${normal_fastq1} ${normal_fastq2} | samtools view -@ ${task.cpus} -b - | samtools sort -o ${normal_id}_sorted.bam
    """
}

process markDupTN {
    tag "Mark Duplicate on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/MarkDuplicatesTN", mode: 'copy'
    cpus 4

    input:
    tuple val(tumor_id), path("${tumor_id}_sorted.bam"),
          val(normal_id), path("${normal_id}_sorted.bam")

    output:
    tuple val(tumor_id), path("${tumor_id}_sorted_md.bam"), path("${tumor_id}_sorted_md_metrics.txt"),
          val(normal_id), path("${normal_id}_sorted_md.bam"), path("${normal_id}_sorted_md_metrics.txt"), emit: markdup_bams

    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'

    script:
    """
    # MD tumor BAM
    java -jar /usr/src/app/picard.jar MarkDuplicates \
    	I=${tumor_id}_sorted.bam \
    	O=${tumor_id}_sorted_md.bam \
    	M=${tumor_id}_sorted_md_metrics.txt

    # MD normal BAM
    java -jar /usr/src/app/picard.jar MarkDuplicates \
    	I=${normal_id}_sorted.bam \
    	O=${normal_id}_sorted_md.bam \
    	M=${normal_id}_sorted_md_metrics.txt
    """
}

process mdBAM_indexTN {
    tag "Index final md BAM on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/MarkDuplicatesTN", mode: 'copy'
    cpus 4

    input:
    tuple val(tumor_id), path("${tumor_id}_sorted_md.bam"), path("${tumor_id}_sorted_md_metrics.txt"),
          val(normal_id), path("${normal_id}_sorted_md.bam"), path("${normal_id}_sorted_md_metrics.txt")

    output:
    tuple val(tumor_id), path("${tumor_id}_sorted_md.bam"), path("${tumor_id}_sorted_md.bam.bai"),
          val(normal_id), path("${normal_id}_sorted_md.bam"), path("${normal_id}_sorted_md.bam.bai"), emit: valid_markdup_bams

    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'

    script:
    """
    samtools index ${tumor_id}_sorted_md.bam
    samtools index ${normal_id}_sorted_md.bam

    """
}
process somVarCall_tumor_normal {
    tag "Somatic Variant on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/VariantCallingTN", mode: 'copy'
    cpus 4

    input:
    tuple val(tumor_id), path(tumor_bam), path(tumor_bai),
          val(normal_id), path(normal_bam), path(normal_bai)
    val(params.ref)
    val(params.bed)
    val(params.gnomad)
    val(params.db1000g)

    output:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_raw.vcf.gz"), path("${tumor_id}_vs_${normal_id}_raw.vcf.gz.tbi"), path("${tumor_id}_vs_${normal_id}_raw.vcf.gz.stats"), emit: final_vcfs

    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'

    script:
    """
    java -jar /usr/src/app/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar Mutect2 \
        -R ${params.ref} \
        -L ${params.bed} \
        -I ${tumor_bam} \
        -I ${normal_bam} \
        -normal ${normal_id} \
        --germline-resource ${params.gnomad} \
        --panel-of-normals ${params.db1000g} \
        -O ${tumor_id}_vs_${normal_id}_raw.vcf.gz
    """
}

process validatevcftn {
    tag "${tumor_id} vs ${normal_id}_validation"

    input:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_raw.vcf.gz"), path("${tumor_id}_vs_${normal_id}_raw.vcf.gz.tbi")

    output:
    tuple val(tumor_id), val(normal_id),path("${tumor_id}_vs_${normal_id}_raw.vcf.gz"), emit: valid_vcfs

    script:
    """
    python3 /usr/src/app/ref17/Validation_script/vcf_validation.py --vcf ${tumor_id}_vs_${normal_id}_raw.vcf.gz
    """
}


process FilterMTtn {
    tag "Somatic Filter on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/FilterSomtn", mode: 'copy'
    cpus 4

    input:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_raw.vcf.gz"), path("${tumor_id}_vs_${normal_id}_raw.vcf.gz.tbi"), path("${tumor_id}_vs_${normal_id}_raw.vcf.gz.stats")
    val(params.ref)
    
    output:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered.vcf.gz"), emit: filtered_vcfs

    script: 
    """
    java -jar /usr/src/app/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar FilterMutectCalls \
        -R ${params.ref} \
        -V ${tumor_id}_vs_${normal_id}_raw.vcf.gz \
        -O ${tumor_id}_vs_${normal_id}_filtered.vcf.gz
    """
}

process KeepPASStn {
    tag "Extract PASS variants on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/filterMutect", mode: 'copy'
    cpus 4

    input:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered.vcf.gz")//worked
    
    output:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered_PASS.vcf"), emit: final_vcfs

    script: 
    """
    bcftools view -f PASS ${tumor_id}_vs_${normal_id}_filtered.vcf.gz > ${tumor_id}_vs_${normal_id}_filtered_PASS.vcf
    """
}

process validateFinalvcftn {
    tag "${tumor_id} vs ${normal_id}_validation"

    input:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered_PASS.vcf")

    output:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered_PASS.vcf"), emit: valid_final_vcfs

    script:
    """
    python3 /usr/src/app/ref17/Validation_script/vcf_Finalvalidation.py --vcf ${tumor_id}_vs_${normal_id}_filtered_PASS.vcf
    """
}

process annotationtn {
    tag "Annotation on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/Annotation", mode: 'copy'
    cpus 4

    input:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered_PASS.vcf")
    val(params.ref)
    val(params.cachee)
    val(params.dirplugin)
    val(params.dbNSFP)
    val(params.loftool)
    val(params.CADDsnv)
    val(params.CADDindel)
    val(params.dbscSNV)
    val(params.assembly)

    output:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered_PASS_annot.txt"), emit: annotVEP_vcfs

    script:
    """
    /usr/src/app/ensembl-vep/./vep --biotype --buffer_size 500 --offline --cache --check_existing --database \
    --assembly ${params.assembly} \
    --dir ${params.cachee} \
    --dir_plugins ${params.dirplugin} \
    --fasta_dir ${params.ref} --force --fork 4 \
    --input_file ${tumor_id}_vs_${normal_id}_filtered_PASS.vcf --mane \
    --output_file ${tumor_id}_vs_${normal_id}_filtered_PASS_annot.txt --tab \
    --plugin dbNSFP,${params.dbNSFP},aapos,genename,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,Uniprot_acc,Uniprot_entry,HGVSc_ANNOVAR,HGVSp_ANNOVAR,HGVSc_snpEff,HGVSp_snpEff,HGVSc_VEP,HGVSp_VEP,TSL,VEP_canonical,cds_strand,SIFT_score,SIFT_converted_rankscore,SIFT_pred,SIFT4G_score,SIFT4G_converted_rankscore,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,VEST4_score,VEST4_rankscore,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,MetaRNN_score,MetaRNN_rankscore,MetaRNN_pred,M-CAP_score,M-CAP_rankscore,M-CAP_pred,REVEL_score,REVEL_rankscore,MutPred_score,MutPred_rankscore,MutPred_Top5features,MVP_score,MVP_rankscore,gMVP_score,gMVP_rankscore,MPC_score,MPC_rankscore,PrimateAI_score,PrimateAI_rankscore,PrimateAI_pred,DEOGEN2_score,DEOGEN2_rankscore,DEOGEN2_pred,BayesDel_addAF_score,BayesDel_addAF_rankscore,BayesDel_addAF_pred,ClinPred_score,ClinPred_rankscore,ClinPred_pred,LIST-S2_score,LIST-S2_rankscore,LIST-S2_pred,VARITY_R_score,VARITY_R_rankscore,VARITY_ER_score,VARITY_ER_rankscore,VARITY_R_LOO_score,VARITY_R_LOO_rankscore,VARITY_ER_LOO_score,VARITY_ER_LOO_rankscore,ESM1b_score,ESM1b_rankscore,ESM1b_pred,EVE_score,EVE_rankscore,EVE_Class10_pred,EVE_Class20_pred,EVE_Class25_pred,EVE_Class30_pred,EVE_Class40_pred,EVE_Class50_pred,EVE_Class60_pred,EVE_Class70_pred,EVE_Class75_pred,EVE_Class80_pred,EVE_Class90_pred,AlphaMissense_score,AlphaMissense_rankscore,AlphaMissense_pred,Aloft_Fraction_transcripts_affected,Aloft_prob_Tolerant,Aloft_prob_Recessive,Aloft_prob_Dominant,Aloft_pred,Aloft_Confidence,CADD_raw,CADD_raw_rankscore,CADD_phred,CADD_raw_hg19,CADD_raw_rankscore_hg19,CADD_phred_hg19,DANN_score,DANN_rankscore,fathmm-MKL_coding_score,fathmm-MKL_coding_rankscore,fathmm-MKL_coding_pred,fathmm-MKL_coding_group,fathmm-XF_coding_score,fathmm-XF_coding_rankscore,fathmm-XF_coding_pred,Eigen-raw_coding,Eigen-raw_coding_rankscore,Eigen-phred_coding,Eigen-PC-raw_coding,Eigen-PC-raw_coding_rankscore,Eigen-PC-phred_coding,GenoCanyon_score,GenoCanyon_rankscore,integrated_fitCons_score,integrated_fitCons_rankscore,integrated_confidence_value,GM12878_fitCons_score,GM12878_fitCons_rankscore,GM12878_confidence_value,H1-hESC_fitCons_score,H1-hESC_fitCons_rankscore,H1-hESC_confidence_value,HUVEC_fitCons_score,HUVEC_fitCons_rankscore,HUVEC_confidence_value,LINSIGHT,LINSIGHT_rankscore,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP470way_mammalian,phyloP470way_mammalian_rankscore,phyloP17way_primate,phyloP17way_primate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phastCons470way_mammalian,phastCons470way_mammalian_rankscore,phastCons17way_primate,phastCons17way_primate_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore,bStatistic,bStatistic_converted_rankscore,1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,UK10K_AC,UK10K_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,ExAC_nonTCGA_AC,ExAC_nonTCGA_AF,ExAC_nonTCGA_Adj_AC,ExAC_nonTCGA_Adj_AF,ExAC_nonTCGA_AFR_AC,ExAC_nonTCGA_AFR_AF,ExAC_nonTCGA_AMR_AC,ExAC_nonTCGA_AMR_AF,ExAC_nonTCGA_EAS_AC,ExAC_nonTCGA_EAS_AF,ExAC_nonTCGA_FIN_AC,ExAC_nonTCGA_FIN_AF,ExAC_nonTCGA_NFE_AC,ExAC_nonTCGA_NFE_AF,ExAC_nonTCGA_SAS_AC,ExAC_nonTCGA_SAS_AF,ExAC_nonpsych_AC,ExAC_nonpsych_AF,ExAC_nonpsych_Adj_AC,ExAC_nonpsych_Adj_AF,ExAC_nonpsych_AFR_AC,ExAC_nonpsych_AFR_AF,ExAC_nonpsych_AMR_AC,ExAC_nonpsych_AMR_AF,ExAC_nonpsych_EAS_AC,ExAC_nonpsych_EAS_AF,ExAC_nonpsych_FIN_AC,ExAC_nonpsych_FIN_AF,ExAC_nonpsych_NFE_AC,ExAC_nonpsych_NFE_AF,ExAC_nonpsych_SAS_AC,ExAC_nonpsych_SAS_AF,gnomAD_exomes_flag,gnomAD_exomes_AC,gnomAD_exomes_AN,gnomAD_exomes_AF,gnomAD_exomes_POPMAX_AC,gnomAD_exomes_POPMAX_AN,gnomAD_exomes_POPMAX_AF,gnomAD_exomes_AFR_AC,gnomAD_exomes_AFR_AN,gnomAD_exomes_AFR_AF,gnomAD_exomes_AFR_nhomalt,gnomAD_exomes_AMR_AC,gnomAD_exomes_AMR_AN,gnomAD_exomes_AMR_AF,gnomAD_exomes_AMR_nhomalt,gnomAD_exomes_ASJ_AC,gnomAD_exomes_ASJ_AN,gnomAD_exomes_ASJ_AF,gnomAD_exomes_ASJ_nhomalt,gnomAD_exomes_EAS_AC,gnomAD_exomes_EAS_AN,gnomAD_exomes_EAS_AF,gnomAD_exomes_EAS_nhomalt,gnomAD_exomes_FIN_AC,gnomAD_exomes_FIN_AN,gnomAD_exomes_FIN_AF,gnomAD_exomes_FIN_nhomalt,gnomAD_exomes_MID_AC,gnomAD_exomes_MID_AN,gnomAD_exomes_MID_AF,gnomAD_exomes_MID_nhomalt,gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_AN,gnomAD_exomes_NFE_AF,gnomAD_exomes_NFE_nhomalt,gnomAD_exomes_SAS_AC,gnomAD_exomes_SAS_AN,gnomAD_exomes_SAS_AF,gnomAD_exomes_SAS_nhomalt,gnomAD_exomes_non_ukb_AC,gnomAD_exomes_non_ukb_AN,gnomAD_exomes_non_ukb_AF,gnomAD_exomes_non_ukb_nhomalt,gnomAD_exomes_non_ukb_AFR_AC,gnomAD_exomes_non_ukb_AFR_AN,gnomAD_exomes_non_ukb_AFR_AF,gnomAD_exomes_non_ukb_AFR_nhomalt,gnomAD_exomes_non_ukb_AMR_AC,gnomAD_exomes_non_ukb_AMR_AN,gnomAD_exomes_non_ukb_AMR_AF,gnomAD_exomes_non_ukb_AMR_nhomalt,gnomAD_exomes_non_ukb_ASJ_AC,gnomAD_exomes_non_ukb_ASJ_AN,gnomAD_exomes_non_ukb_ASJ_AF,gnomAD_exomes_non_ukb_ASJ_nhomalt,gnomAD_exomes_non_ukb_EAS_AC,gnomAD_exomes_non_ukb_EAS_AN,gnomAD_exomes_non_ukb_EAS_AF,gnomAD_exomes_non_ukb_EAS_nhomalt,gnomAD_exomes_non_ukb_FIN_AC,gnomAD_exomes_non_ukb_FIN_AN,gnomAD_exomes_non_ukb_FIN_AF,gnomAD_exomes_non_ukb_FIN_nhomalt,gnomAD_exomes_non_ukb_MID_AC,gnomAD_exomes_non_ukb_MID_AN,gnomAD_exomes_non_ukb_MID_AF,gnomAD_exomes_non_ukb_MID_nhomalt,gnomAD_exomes_non_ukb_NFE_AC,gnomAD_exomes_non_ukb_NFE_AN,gnomAD_exomes_non_ukb_NFE_AF,gnomAD_exomes_non_ukb_NFE_nhomalt,gnomAD_exomes_non_ukb_SAS_AC,gnomAD_exomes_non_ukb_SAS_AN,gnomAD_exomes_non_ukb_SAS_AF,gnomAD_exomes_non_ukb_SAS_nhomalt,ALFA_European_AC,ALFA_European_AN,ALFA_European_AF,ALFA_African_Others_AC,ALFA_African_Others_AN,ALFA_African_Others_AF,ALFA_East_Asian_AC,ALFA_East_Asian_AN,ALFA_East_Asian_AF,ALFA_African_American_AC,ALFA_African_American_AN,ALFA_African_American_AF,ALFA_Latin_American_1_AC,ALFA_Latin_American_1_AN,ALFA_Latin_American_1_AF,ALFA_Latin_American_2_AC,ALFA_Latin_American_2_AN,ALFA_Latin_American_2_AF,ALFA_Other_Asian_AC,ALFA_Other_Asian_AN,ALFA_Other_Asian_AF,ALFA_South_Asian_AC,ALFA_South_Asian_AN,ALFA_South_Asian_AF,ALFA_Other_AC,ALFA_Other_AN,ALFA_Other_AF,ALFA_African_AC,ALFA_African_AN,ALFA_African_AF,ALFA_Asian_AC,ALFA_Asian_AN,ALFA_Asian_AF,ALFA_Total_AC,ALFA_Total_AN,ALFA_Total_AF,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,Interpro_domain,GTEx_V8_eQTL_gene,GTEx_V8_eQTL_tissue,GTEx_V8_sQTL_gene,GTEx_V8_sQTL_tissue \
    --pubmed \
    --plugin LoFtool,${params.loftool} \
    --plugin CADD,snv=${params.CADDsnv},indels=${params.CADDindel} \
    --plugin dbscSNV,${params.dbscSNV} \
    --plugin LOVD --refseq --quiet --safe --regulatory --species homo_sapiens --tsl \
    --show_ref_allele --stats_text --symbol --transcript_version --sift b --polyphen b --uploaded_allele
    """
}

process annot_processingtn {
    tag "Annot output proccesing on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/AnnotationProctn", mode: 'copy'
    cpus 4

    input:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered_PASS_annot.txt")

    output:
    tuple val(tumor_id), val(normal_id), path("${tumor_id}_vs_${normal_id}_filtered_PASS_annot_proccessed.tsv"), emit: vep_TSV

    script:
    """
    python /usr/src/app/ref17/Validation_script/VEP_postprocessing.py --input ${tumor_id}_vs_${normal_id}_filtered_PASS_annot.txt  --output ${tumor_id}_vs_${normal_id}_filtered_PASS_annot_proccessed.tsv
   
    """
}

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

process CNV_somtn {
    tag "Somatic Copy Number on ${tumor_id} vs ${normal_id}"
    publishDir "${params.outdir}/CNV_somatic", mode: 'copy'
    cpus 4

    input:
    tuple val(tumor_id), path(tumor_bam), path(tumor_bai),
          val(normal_id), path(normal_bam), path(normal_bai)
    val(params.ref)
    val(params.bed)

    output:
    path "${tumor_id}_vs_${normal_id}", emit: cnvSom

    when:
    (params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI') && params.library == 'Paired' && params.mode == 'TumorNormal'
    
    script:
    """
    source activate p2Manta
    python /usr/src/app/CONTRA.v2.0.8/contra.py \
    --target ${params.bed} \
    --test ${tumor_bam} \
    --control ${normal_bam} \
    --fasta ${params.ref} \
    -p --sampleName ${tumor_id}_vs_${normal_id} -o ${tumor_id}_vs_${normal_id}
    """
}
/// end of TN annot

process align_paired {
    tag "Alignment on ${sample_id} (Paired-End)"
    publishDir "${params.outdir}/3.Alignment", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz")
    val(params.ref)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: sorted_bams_pe
    
    when:
    params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI' && params.library == 'Paired'
    
    script:
    """
    bwa mem -t ${task.cpus} -M \
            -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${sample_id}" \
            ${params.ref} ${sample_id}_R1_trimmed.fastq.gz ${sample_id}_R2_trimmed.fastq.gz | samtools view -@ ${task.cpus} -b - | samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam 
    """
}

process align_single {
    tag "Alignment on ${sample_id} (Single-End)"
    publishDir "${params.outdir}/3.Alignment", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz")
    val(params.ref)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: sorted_bams_se

    when:
    params.platform == 'Illumina' || params.platform == 'BGI' || params.platform == 'MGI' || params.platform == 'ThermoFisher' && params.library == 'Single'

    script:
    """
    bwa mem -t ${task.cpus} -M \
            -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${sample_id}" \
            ${params.ref} ${sample_id}_trimmed.fastq.gz | samtools view -@ ${task.cpus} -b - | samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam 
    """
}
process align_ONT {
    tag "Alignment on ${sample_id} (ONT)"
    publishDir "${params.outdir}/3.Alignment", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz")
    val(params.ref)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: sorted_bams_ont

    when:
    params.platform == 'Nanopore' && params.library == 'Single'

    script:
    """
    /usr/src/app/minimap2/./minimap2 -ax map-ont ${params.ref} ${sample_id}_trimmed.fastq.gz | samtools view -@ ${task.cpus} -b - | samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam
    """
}
process validateBAM {
    tag "${sample_id}_BAMvalidation"

    input:
    tuple val(sample_id), path("${sample_id}_sorted.bam")

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: bam_validation

    script:
    """
    python3 /usr/src/app/ref17/Validation_script/bamvalidation.py --bam ${sample_id}_sorted.bam
    """
}

process markDup {
    tag "Mark Duplicate on ${sample_id}"
    publishDir "${params.outdir}/4.MarkDuplicates", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path("${sample_id}_sorted.bam")

    output:
    tuple val(sample_id), path("${sample_id}_sorted_md.bam"), path("${sample_id}_sorted_md_metrics.txt"), emit: markdup_bams

    script:
    """
    java -jar /usr/src/app/picard.jar MarkDuplicates \
    	I=${sample_id}_sorted.bam \
    	O=${sample_id}_sorted_md.bam \
    	M=${sample_id}_sorted_md_metrics.txt
    """
}

process mdvalidate {
    tag "${sample_id}_mdvalidation"

    input:
    tuple val(sample_id), path("${sample_id}_sorted_md.bam")

    output:
    tuple val(sample_id), path("${sample_id}_sorted_md.bam"), emit: valid_markdup_bams

    script:
    """
    python3 /usr/src/app/ref17/Validation_script/picardOutPutvalidation.py --bam ${sample_id}_sorted_md.bam
    """
}

process mdBAM_index {
    tag "Index final md BAM on ${sample_id}"
    publishDir "${params.outdir}/4.MarkDuplicates", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path("${sample_id}_sorted_md.bam")

    output:
    tuple val(sample_id), path("${sample_id}_sorted_md.bam"), path("${sample_id}_sorted_md.bam.bai"), emit: valid_markdup_bams

    script:
    """
    samtools index ${sample_id}_sorted_md.bam

    """
}

process somVarCall {
    tag "Somatic Variant on ${sample_id}"
    publishDir "${params.outdir}/5.Variant", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(bam), path(bai)
    val(params.ref)
    val(params.bed)
    val(params.gnomad)
    val(params.db1000g)

    output:
    tuple val(sample_id), path("${sample_id}_raw.vcf.gz"), path("${sample_id}_raw.vcf.gz.tbi"), path("${sample_id}_raw.vcf.gz.stats"), emit: raw_vcfs

    script:
    """
    java -jar /usr/src/app/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar Mutect2 \
        -R ${params.ref} \
        --native-pair-hmm-threads ${task.cpus} \
        -L ${params.bed} \
        --germline-resource ${params.gnomad} \
        --panel-of-normals ${params.db1000g} \
        -I ${bam} \
        -O ${sample_id}_raw.vcf.gz
    """
}

process validatevcf {
    tag "${sample_id}_validation"

    input:
    tuple val(sample_id), path("${sample_id}_raw.vcf.gz")

    output:
    tuple val(sample_id), path("${sample_id}_raw.vcf.gz"), emit: valid_vcfs

    script:
    """
    python3 /usr/src/app/ref17/Validation_script/vcf_validation.py --vcf ${sample_id}_raw.vcf.gz
    """
}

process FilterMT {
    tag "Somatic Filter on ${sample_id}"
    publishDir "${params.outdir}/6.FilterSom", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(vcf), path(tbi), path(tsv)
    val(params.ref)
    
    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf.gz"), emit: filtered_vcfs

    script: 
    """
    java -jar /usr/src/app/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar FilterMutectCalls \
        -R ${params.ref} \
        -V ${vcf} \
        -O ${sample_id}_filtered.vcf.gz
    """
}

process KeepPASS {
    tag "Extract PASS variants on ${sample_id}"
    publishDir "${params.outdir}/7.FinalFilteredVCF", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(vcf)//worked
    
    output:
    tuple val(sample_id), path("${sample_id}_filtered_PASS.vcf"), emit: final_vcfs

    script: 
    """
    bcftools view -f PASS ${vcf} > ${sample_id}_filtered_PASS.vcf
    """
}
process somVarCall_ont {
    tag "Somatic Variant on ${sample_id} (ONT)"
    publishDir "${params.outdir}/5.Variant", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(bam), path(bai)
    val(params.ref)
    val(params.bed)
 
    output:
    //path "${sample_id}", emit: final_vcfs
    tuple val(sample_id), path("${sample_id}/${sample_id}.vcf"), emit: final_vcfs

    script:
    """
    source activate nanocaller_env 
    mkdir -p ${sample_id}
    /usr/src/app/NanoCaller/./NanoCaller \
        --bam ${bam} \
        --ref ${params.ref} \
        --bed ${params.bed} \
        --cpu 20 \
        --mode all \
        --sequencing ont \
        --output ${sample_id} \
        --prefix ${sample_id}
    """
}
process validateFinalvcf {
    tag "${sample_id}_validation"

    input:
    tuple val(sample_id), path("${sample_id}.vcf")

    output:
    tuple val(sample_id), path("${sample_id}.vcf"), emit: valid_final_vcfs

    script:
    """
    python3 /usr/src/app/ref17/Validation_script/vcf_Finalvalidation.py --vcf ${sample_id}.vcf
    """
}

process annotation {
    tag "Annotation on ${sample_id}"
    publishDir "${params.outdir}/8.Annotation", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(vcf)
    val(params.ref)
    val(params.cachee)
    val(params.dirplugin)
    val(params.dbNSFP)
    val(params.loftool)
    val(params.CADDsnv)
    val(params.CADDindel)
    val(params.dbscSNV)
    val(params.assembly)

    output:
    tuple val(sample_id), path("${sample_id}_filtered_PASS_annot.txt"), emit: annotVEP_vcfs

    script:
    """
    /usr/src/app/ensembl-vep/./vep --biotype --buffer_size 500 --offline --cache --check_existing --database \
    --assembly ${params.assembly} \
    --dir ${params.cachee} \
    --dir_plugins ${params.dirplugin} \
    --fasta_dir ${params.ref} --force --fork 4 \
    --input_file ${vcf} --mane \
    --output_file ${sample_id}_filtered_PASS_annot.txt --tab \
    --plugin dbNSFP,${params.dbNSFP},aapos,genename,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,Uniprot_acc,Uniprot_entry,HGVSc_ANNOVAR,HGVSp_ANNOVAR,HGVSc_snpEff,HGVSp_snpEff,HGVSc_VEP,HGVSp_VEP,TSL,VEP_canonical,cds_strand,SIFT_score,SIFT_converted_rankscore,SIFT_pred,SIFT4G_score,SIFT4G_converted_rankscore,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,VEST4_score,VEST4_rankscore,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,MetaRNN_score,MetaRNN_rankscore,MetaRNN_pred,M-CAP_score,M-CAP_rankscore,M-CAP_pred,REVEL_score,REVEL_rankscore,MutPred_score,MutPred_rankscore,MutPred_Top5features,MVP_score,MVP_rankscore,gMVP_score,gMVP_rankscore,MPC_score,MPC_rankscore,PrimateAI_score,PrimateAI_rankscore,PrimateAI_pred,DEOGEN2_score,DEOGEN2_rankscore,DEOGEN2_pred,BayesDel_addAF_score,BayesDel_addAF_rankscore,BayesDel_addAF_pred,ClinPred_score,ClinPred_rankscore,ClinPred_pred,LIST-S2_score,LIST-S2_rankscore,LIST-S2_pred,VARITY_R_score,VARITY_R_rankscore,VARITY_ER_score,VARITY_ER_rankscore,VARITY_R_LOO_score,VARITY_R_LOO_rankscore,VARITY_ER_LOO_score,VARITY_ER_LOO_rankscore,ESM1b_score,ESM1b_rankscore,ESM1b_pred,EVE_score,EVE_rankscore,EVE_Class10_pred,EVE_Class20_pred,EVE_Class25_pred,EVE_Class30_pred,EVE_Class40_pred,EVE_Class50_pred,EVE_Class60_pred,EVE_Class70_pred,EVE_Class75_pred,EVE_Class80_pred,EVE_Class90_pred,AlphaMissense_score,AlphaMissense_rankscore,AlphaMissense_pred,Aloft_Fraction_transcripts_affected,Aloft_prob_Tolerant,Aloft_prob_Recessive,Aloft_prob_Dominant,Aloft_pred,Aloft_Confidence,CADD_raw,CADD_raw_rankscore,CADD_phred,CADD_raw_hg19,CADD_raw_rankscore_hg19,CADD_phred_hg19,DANN_score,DANN_rankscore,fathmm-MKL_coding_score,fathmm-MKL_coding_rankscore,fathmm-MKL_coding_pred,fathmm-MKL_coding_group,fathmm-XF_coding_score,fathmm-XF_coding_rankscore,fathmm-XF_coding_pred,Eigen-raw_coding,Eigen-raw_coding_rankscore,Eigen-phred_coding,Eigen-PC-raw_coding,Eigen-PC-raw_coding_rankscore,Eigen-PC-phred_coding,GenoCanyon_score,GenoCanyon_rankscore,integrated_fitCons_score,integrated_fitCons_rankscore,integrated_confidence_value,GM12878_fitCons_score,GM12878_fitCons_rankscore,GM12878_confidence_value,H1-hESC_fitCons_score,H1-hESC_fitCons_rankscore,H1-hESC_confidence_value,HUVEC_fitCons_score,HUVEC_fitCons_rankscore,HUVEC_confidence_value,LINSIGHT,LINSIGHT_rankscore,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP470way_mammalian,phyloP470way_mammalian_rankscore,phyloP17way_primate,phyloP17way_primate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phastCons470way_mammalian,phastCons470way_mammalian_rankscore,phastCons17way_primate,phastCons17way_primate_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore,bStatistic,bStatistic_converted_rankscore,1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,UK10K_AC,UK10K_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,ExAC_nonTCGA_AC,ExAC_nonTCGA_AF,ExAC_nonTCGA_Adj_AC,ExAC_nonTCGA_Adj_AF,ExAC_nonTCGA_AFR_AC,ExAC_nonTCGA_AFR_AF,ExAC_nonTCGA_AMR_AC,ExAC_nonTCGA_AMR_AF,ExAC_nonTCGA_EAS_AC,ExAC_nonTCGA_EAS_AF,ExAC_nonTCGA_FIN_AC,ExAC_nonTCGA_FIN_AF,ExAC_nonTCGA_NFE_AC,ExAC_nonTCGA_NFE_AF,ExAC_nonTCGA_SAS_AC,ExAC_nonTCGA_SAS_AF,ExAC_nonpsych_AC,ExAC_nonpsych_AF,ExAC_nonpsych_Adj_AC,ExAC_nonpsych_Adj_AF,ExAC_nonpsych_AFR_AC,ExAC_nonpsych_AFR_AF,ExAC_nonpsych_AMR_AC,ExAC_nonpsych_AMR_AF,ExAC_nonpsych_EAS_AC,ExAC_nonpsych_EAS_AF,ExAC_nonpsych_FIN_AC,ExAC_nonpsych_FIN_AF,ExAC_nonpsych_NFE_AC,ExAC_nonpsych_NFE_AF,ExAC_nonpsych_SAS_AC,ExAC_nonpsych_SAS_AF,gnomAD_exomes_flag,gnomAD_exomes_AC,gnomAD_exomes_AN,gnomAD_exomes_AF,gnomAD_exomes_POPMAX_AC,gnomAD_exomes_POPMAX_AN,gnomAD_exomes_POPMAX_AF,gnomAD_exomes_AFR_AC,gnomAD_exomes_AFR_AN,gnomAD_exomes_AFR_AF,gnomAD_exomes_AFR_nhomalt,gnomAD_exomes_AMR_AC,gnomAD_exomes_AMR_AN,gnomAD_exomes_AMR_AF,gnomAD_exomes_AMR_nhomalt,gnomAD_exomes_ASJ_AC,gnomAD_exomes_ASJ_AN,gnomAD_exomes_ASJ_AF,gnomAD_exomes_ASJ_nhomalt,gnomAD_exomes_EAS_AC,gnomAD_exomes_EAS_AN,gnomAD_exomes_EAS_AF,gnomAD_exomes_EAS_nhomalt,gnomAD_exomes_FIN_AC,gnomAD_exomes_FIN_AN,gnomAD_exomes_FIN_AF,gnomAD_exomes_FIN_nhomalt,gnomAD_exomes_MID_AC,gnomAD_exomes_MID_AN,gnomAD_exomes_MID_AF,gnomAD_exomes_MID_nhomalt,gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_AN,gnomAD_exomes_NFE_AF,gnomAD_exomes_NFE_nhomalt,gnomAD_exomes_SAS_AC,gnomAD_exomes_SAS_AN,gnomAD_exomes_SAS_AF,gnomAD_exomes_SAS_nhomalt,gnomAD_exomes_non_ukb_AC,gnomAD_exomes_non_ukb_AN,gnomAD_exomes_non_ukb_AF,gnomAD_exomes_non_ukb_nhomalt,gnomAD_exomes_non_ukb_AFR_AC,gnomAD_exomes_non_ukb_AFR_AN,gnomAD_exomes_non_ukb_AFR_AF,gnomAD_exomes_non_ukb_AFR_nhomalt,gnomAD_exomes_non_ukb_AMR_AC,gnomAD_exomes_non_ukb_AMR_AN,gnomAD_exomes_non_ukb_AMR_AF,gnomAD_exomes_non_ukb_AMR_nhomalt,gnomAD_exomes_non_ukb_ASJ_AC,gnomAD_exomes_non_ukb_ASJ_AN,gnomAD_exomes_non_ukb_ASJ_AF,gnomAD_exomes_non_ukb_ASJ_nhomalt,gnomAD_exomes_non_ukb_EAS_AC,gnomAD_exomes_non_ukb_EAS_AN,gnomAD_exomes_non_ukb_EAS_AF,gnomAD_exomes_non_ukb_EAS_nhomalt,gnomAD_exomes_non_ukb_FIN_AC,gnomAD_exomes_non_ukb_FIN_AN,gnomAD_exomes_non_ukb_FIN_AF,gnomAD_exomes_non_ukb_FIN_nhomalt,gnomAD_exomes_non_ukb_MID_AC,gnomAD_exomes_non_ukb_MID_AN,gnomAD_exomes_non_ukb_MID_AF,gnomAD_exomes_non_ukb_MID_nhomalt,gnomAD_exomes_non_ukb_NFE_AC,gnomAD_exomes_non_ukb_NFE_AN,gnomAD_exomes_non_ukb_NFE_AF,gnomAD_exomes_non_ukb_NFE_nhomalt,gnomAD_exomes_non_ukb_SAS_AC,gnomAD_exomes_non_ukb_SAS_AN,gnomAD_exomes_non_ukb_SAS_AF,gnomAD_exomes_non_ukb_SAS_nhomalt,ALFA_European_AC,ALFA_European_AN,ALFA_European_AF,ALFA_African_Others_AC,ALFA_African_Others_AN,ALFA_African_Others_AF,ALFA_East_Asian_AC,ALFA_East_Asian_AN,ALFA_East_Asian_AF,ALFA_African_American_AC,ALFA_African_American_AN,ALFA_African_American_AF,ALFA_Latin_American_1_AC,ALFA_Latin_American_1_AN,ALFA_Latin_American_1_AF,ALFA_Latin_American_2_AC,ALFA_Latin_American_2_AN,ALFA_Latin_American_2_AF,ALFA_Other_Asian_AC,ALFA_Other_Asian_AN,ALFA_Other_Asian_AF,ALFA_South_Asian_AC,ALFA_South_Asian_AN,ALFA_South_Asian_AF,ALFA_Other_AC,ALFA_Other_AN,ALFA_Other_AF,ALFA_African_AC,ALFA_African_AN,ALFA_African_AF,ALFA_Asian_AC,ALFA_Asian_AN,ALFA_Asian_AF,ALFA_Total_AC,ALFA_Total_AN,ALFA_Total_AF,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,Interpro_domain,GTEx_V8_eQTL_gene,GTEx_V8_eQTL_tissue,GTEx_V8_sQTL_gene,GTEx_V8_sQTL_tissue \
    --pubmed \
    --plugin LoFtool,${params.loftool} \
    --plugin CADD,snv=${params.CADDsnv},indels=${params.CADDindel} \
    --plugin dbscSNV,${params.dbscSNV} \
    --plugin LOVD --quiet --safe --stats_text --symbol --transcript_version --sift b --polyphen b
    """
}

process annot_processing {
    tag "Annot output proccesing on ${sample_id}"
    publishDir "${params.outdir}/8.Annotation", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path("${sample_id}_filtered_PASS_annot.txt")

    output:
    tuple val(sample_id), path("${sample_id}_filtered_PASS_annot_proccessed.tsv"), emit: vep_TSV

    script:
    """
    python /usr/src/app/ref17/Validation_script/VEP_postprocessing.py --input ${sample_id}_filtered_PASS_annot.txt --output ${sample_id}_filtered_PASS_annot_proccessed.tsv
    """
}

process ampclass {
    tag "AMP class on ${sample_id}"
    publishDir "${params.outdir}/8.Annotation", mode: 'copy'
    cpus 4
    cache 'false'

    input:
    tuple val(sample_id), path(vcf)
    val(params.assemblyAMP)

    output:
    tuple val(sample_id), path("${sample_id}.${params.assemblyAMP}_multianno.txt.cancervar"), path("${sample_id}.${params.assemblyAMP}_multianno.txt"), path("${sample_id}.${params.assemblyAMP}_multianno.txt.grl_p"), emit: vep_TSVamp

    script:
    """
    python /usr/src/app/ref17/annovarhg38/cancewarhg38/./CancerVar.py  -b ${params.assemblyAMP} -i ${vcf} --input_type=VCF -o ${sample_id}
    """
}

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

process cosmic_hg38 {
    tag "COSMIC annotation on ${sample_id}"
    publishDir "${params.outdir}/9.AdvanceAnnot", mode: 'copy'
    cpus 4

    input:
    val(params.cosmiccod)
    tuple val(sample_id), path("${sample_id}_amp_vcf_output.tsv")
 
    output:
    tuple val(sample_id), path("${sample_id}_amp_vcf_output_cosmic.tsv"), emit: vcf_vep_TSV1

    when:
    params.genome == 'hg38'
    script:
    """
    python /usr/src/app/ref17/Validation_script/cosmic_dataTSV_argument.py --cosmic_file ${params.cosmiccod} --input_file ${sample_id}_amp_vcf_output.tsv --output_file ${sample_id}_amp_vcf_output_cosmic.tsv
    """
}

process cosmic_hg19 {
    tag "COSMIC annotation on ${sample_id}"
    publishDir "${params.outdir}/9.AdvanceAnnot", mode: 'copy'
    cpus 4

    input:
    val(params.cosmiccod)
    tuple val(sample_id), path("${sample_id}_amp_vcf_output.tsv")
 
    output:
    tuple val(sample_id), path("${sample_id}_amp_vcf_output_cosmic.tsv"), emit: vcf_vep_TSV1

    when:
    params.genome == 'hg19'
    script:
    """
    python /usr/src/app/ref17/Validation_script/cosmic_dataTSV_argument.py --cosmic_file ${params.cosmiccod} --input_file ${sample_id}_amp_vcf_output.tsv --output_file ${sample_id}_amp_vcf_output_cosmic.tsv
    """
}
process vcf2maf {
    tag "VCF to MAF on ${sample_id}"
    publishDir "${params.outdir}/9.AdvanceAnnot", mode: 'copy'
    cpus 4

    input:
    //tuple val(sample_id), path(vcf)
    //val(params.assemblyAMP)
    tuple val(sample_id), path("${sample_id}_filtered_PASS.vcf")
    val(params.cachee)
    val(params.ref)
    val(params.assembly)

    output:
    tuple val(sample_id), path("${sample_id}_vcf2MAF.maf"), emit: vcf_maf

   
    script:
    """
    perl /usr/src/app/vcf2maf/vcf2maf.pl --vep-path /usr/src/app/ensembl-vep --vep-data ${params.cachee} --ref-fasta ${params.ref} --ncbi-build ${params.assembly} --input-vcf ${sample_id}_filtered_PASS.vcf --output-maf ${sample_id}_vcf2MAF.maf
    """
}

process oncokb {
    tag "Precison Analysis on ${sample_id}"
    publishDir "${params.outdir}/9.AdvanceAnnot", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path("${sample_id}_vcf2MAF.maf")
    val(params.assembly)
    val(params.tumerType)
    val(params.OKBAPI)
    
    output:
    tuple val(sample_id), path("${sample_id}_annotat_oncokb.tsv"), emit: oncokbA

  
    script:
    """
    python /usr/src/app/oncokb-annotator/MafAnnotator.py -i ${sample_id}_vcf2MAF.maf -o ${sample_id}_annotat_oncokb.tsv -b ${params.OKBAPI} -a -d -t ${params.tumerType} -r ${params.assembly}
    """
}

process precision_analysis {
    tag "Precison Analysis on ${sample_id}"
    publishDir "${params.outdir}/10.Precision", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path("${sample_id}_amp_vcf_output_cosmic.tsv")
    tuple val(sample_id), path("${sample_id}_annotat_oncokb.tsv")

    output:
    tuple val(sample_id), path("${sample_id}_precision_final_output.xlsx"), emit: precision_out

    script:
    """
    python /usr/src/app/ref17/Validation_script/cosmic_okb.py --cosmic ${sample_id}_amp_vcf_output_cosmic.tsv --okb ${sample_id}_annotat_oncokb.tsv --output ${sample_id}_precision_final_output.xlsx
    """
}

process msiScore {
    tag "MSI Score Calculation on ${sample_id}"
    publishDir "${params.outdir}/11.MSIScore", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(bam), path(bai) //worked
    val(params.msimodel)

    output:
    path("${sample_id}.MSIscore.log"), emit: msi_score

    script:
    """
    /usr/src/app/msisensor2/msisensor2 msi -b 10 -M ${params.msimodel} -t ${bam} -o ${sample_id}.MSIscore.log
    """
}

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
process SV_som {
  
    tag "Structural Variants on ${sample_id}"
    publishDir "${params.outdir}/13.SV_somatic", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(bam), path(bai)
    val(params.ref)

    output:
    path "${sample_id}", emit: sv_TUmsomatic

    script:
    """
    # Activate the environment for Manta
    source activate p2Manta

    # Create a temporary directory in the work directory for Manta to run
    mkdir ${sample_id}

    # Run Manta config and workflow
    configManta.py --tumorBam ${bam} --referenceFasta ${params.ref} --runDir ${sample_id}
    ${sample_id}/runWorkflow.py
    """
}
process CNV_som {
    tag "Somatic Copy Number on ${sample_id}"
    publishDir "${params.outdir}/14.CNV_somatic", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(bam), path(bai)
    val(params.ref)
    val(params.baselineContra)
    val(params.bed)

    output:
    path "${sample_id}", emit: cnvSom
    
    script:
    """
    source activate p2Manta
    python /usr/src/app/CONTRA.v2.0.8/contra.py \
    --target ${params.bed} \
    --test ${bam} \
    --control ${params.baselineContra} \
    --fasta ${params.ref} \
    --largeDeletion --bed -p --sampleName ${sample_id} -o ${sample_id}
    """
}
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

workflow.onComplete = {

    // any workflow property can be used here
    // Stylish outro message with start and end times
    println """
    ${bold}${green}============================================================
         Thank you for using VgenX CLI - A WES Oncology Pipeline!
    ============================================================
 
    ============================================================
    Your analysis is complete.
    For support, please contact ${yellow}support@vgenomics.co.in${green}.
    ============================================================${reset}
    """

}
/*
// Capture the end time
def endTime = LocalDateTime.now()
def formattedEndTime = endTime.format(dateFormat)

// Stylish outro message with start and end times
println """
${bold}${green}============================================================
         Thank you for using VgenX CLI - A WES Oncology Pipeline!
============================================================
 
  End Time  : ${yellow}${formattedEndTime}${green}
============================================================
  Your analysis is complete.
  For support, please contact ${yellow}support@vgenomics.co.in${green}.
============================================================${reset}
"""
End Time  : ${yellow}${formattedEndTime}${green}


annotVUS_vcfs = annotationVUS(final_vcfs, params.cachee, params.dirplugin, params.ref, params.assembly, params.dbNSFP, params.loftool, params.CADDsnv, params.dbscSNV, params.MaxEnt)
annotVUS_proce = vusprizePre(annotVUS_vcfs)
process annotationVUS {
    tag "VUS Pre-Proccessing on ${sample_id}"
    publishDir "${params.outdir}/16.VUS", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(vcf)
    val(params.ref)
    val(params.cachee)
    val(params.dirplugin)
    val(params.dbNSFP)
    val(params.loftool)
    val(params.CADDsnv)
    val(params.dbscSNV)
    val(params.assembly)
    val(params.MaxEnt)

    output:
    tuple val(sample_id), path("${sample_id}_filtered_PASS_annotVUS.txt"), emit: annotVUS_vcfs

    script:
    """
    /usr/src/app/ensembl-vep/./vep --af --appris --buffer_size 500 --cache --check_existing \
    --distance 5000 --database --offline \
    --assembly ${params.assembly} \
    --fasta ${params.ref} \
    --dir ${params.cachee} \
    --dir_plugins ${params.dirplugin} \
    --mane --pick \
    --fasta_dir ${params.ref} --force --fork 4 \
    --input_file ${vcf} \
    --output_file ${sample_id}_filtered_PASS_annotVUS.txt --tab \
    --plugin  MaxEntScan,${params.MaxEnt} \
    --plugin Blosum62 \
    --plugin dbNSFP,${params.dbNSFP},codon_degeneracy,Eigen-phred_coding,integrated_fitCons_score,GERP++_RS,phyloP100way_vertebrate \
    --plugin dbscSNV,${params.dbscSNV} \
    --plugin LoFtool,${params.loftool} \
    --plugin CADD,snv=${params.CADDsnv} \
    --refseq --quiet --regulatory -sift s --species homo_sapiens --symbol --transcript_version --tsl --safe \
    --show_ref_allele --stats_text --uploaded_allele
    """
}

process vusprizePre {
    tag "VUS Pre-Proccessing on ${sample_id}"
    publishDir "${params.outdir}/16.VUS", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path("${sample_id}_filtered_PASS_annotVUS.txt")

    output:
    tuple val(sample_id), path("${sample_id}_filtered_PASS_annotVUS_proccessed.tsv"), emit: annotVUS_proce

    script:
    """
    python /usr/src/app/ref17/Validation_script/vep_vusprize.py --input ${sample_id}_filtered_PASS_annotVUS.txt --output ${sample_id}_filtered_PASS_annotVUS_proccessed.tsv
   
    """
}
*/