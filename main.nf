#!/usr/bin/env nextflow
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

// Stylish intro message with current date and time
println """
${bold}${cyan}============================================================
    WELCOME TO THE WES ONCOLOGY PIPELINE!
    Powered by: ${yellow}Nextflow + Docker${cyan}
============================================================
    Cutting-edge genomic analysis for precision oncology. 
============================================================${reset}
"""
// Capture start time for tracking the duration
def startTime = LocalDateTime.now()

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
            fastqc_outTN = quality_checkTN(validated_samplesTN)
            trimmed_fastqs_tn = fastpTumorNormal(validated_samplesTN)
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

workflow.onComplete = {
    // Stylish outro message with start and end times
    println """
    ${bold}${green}============================================================
         Thank you for using WES-SOM CLI - A WES Oncology Pipeline!
    ============================================================
 
    ============================================================
    Your analysis is complete.
    For support, please contact ${yellow}d6803148@gmail.com${green}.
    ============================================================${reset}
    """

}