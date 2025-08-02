// Annotation
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