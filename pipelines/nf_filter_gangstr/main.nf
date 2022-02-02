#!/usr/bin/env nextflow

// required parameters
params.sample_ids = null
params.sample_dir = null

// check if all required parameters are set, if not: throw error (should add printing of usage here)
if(params.sample_ids == null || params.sample_dir == null) {
    throw new Exception("Error: one or more required command line arguments were not set")
}

// default output directory
// params.output_dir = "/cfs/earth/scratch/verb/projects/genotype_strs/pipelines/nf_filter_gangstr/nf_out"
params.reference = "/cfs/earth/scratch/verb/projects/genotype_strs/data/1000g_reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.str_panel = "/cfs/earth/scratch/verb/projects/genotype_strs/data/str_loci/tral_panel/all_genes_strs_thresh_PERF.tsv"

// ExpansionHunter denovo default paramters
params.ehdn_min_unit_len = 1
params.ehdn_max_unit_len = 6

// default filtering parameters to use in dumpSTR, can be overwritten at the command line
params.mindepth = 20 
params.maxdepth = 1000
params.minqual = 0.9

Channel
    .fromPath(params.sample_ids)
    .splitCsv(header: true, sep: "\t")
    .map{ row -> tuple(row.sample_id, row.sample_sex, row.alignment) }
    .into{ sample_id_ch1; sample_id_ch2 }

// ----------------------------------------------------
// Process for ExpansionHunter denovo
// ----------------------------------------------------
//  Input:
//      Alignment file, reference genome
//  Output:
//      Annotation of expanded STRs
// ----------------------------------------------------
process ExpansionHunterDN {
    publishDir "${params.sample_dir}/${sample_id}/profiles/EHdn/", mode: "move"
    scratch true
    stageOutMode 'move'

    conda "/cfs/earth/scratch/verb/.conda/envs/EHdn"

    input:
    set sample_id, sample_sex, alignment from sample_id_ch1

    output:
    path "${sample_id}_EHdn*"

    script:
    """
    # Copy data to local scratch on node
    cp ${params.reference}.gz ./reference.fa.gz
    cp ${params.reference}.fai ./reference.fa.fai

    gunzip reference.fa.gz
    
    # cp ${params.sample_dir}/${sample_id}/alignments/${sample_id}.bam* ./    
    cp ${alignment} ./alignment
    cp ${alignment}.crai ./alignment.crai

    # make publishDir, run ExpansionHunter denovo    
    ExpansionHunterDenovo profile \
        --reads alignment \
        --reference reference.fa \
        --output-prefix ${sample_id}_EHdn \
        --min-unit-len ${params.ehdn_min_unit_len} \
        --max-unit-len ${params.ehdn_max_unit_len}
    
    mkdir -p "${params.sample_dir}/${sample_id}/profiles/EHdn/"
    """
}


// ----------------------------------------------------
//  Process for GangSTR
// ----------------------------------------------------
//  Input:
//      Alignment file, reference genome, 
//      STR locus panel
//  Output:
//      Heterozygous genotypes for panel of STR loci
// ----------------------------------------------------
process GangSTR {
    publishDir "${params.sample_dir}/${sample_id}/profiles/GangSTR/raw/", mode: "copy"
    scratch true
    stageOutMode 'move'

    conda "/cfs/earth/scratch/verb/.conda/envs/GangSTR"

    input:
    set sample_id, sample_sex, alignment from sample_id_ch2

    output:
    set sample_id, path("${sample_id}_GangSTR.vcf.gz") into GangSTR_out_ch
    path "${sample_id}_GangSTR*"

    script:
    """
    # Copy data to local scratch on node
    cp ${params.reference}.gz ./reference.fa.gz
    cp ${params.reference}.fai ./reference.fa.fai

    gunzip reference.fa.gz

    cp ${params.str_panel} >(cut -f 1-5 > str_panel.tsv)
    
    # cp ${params.sample_dir}/${sample_id}/alignments/${sample_id}.{cr,b}am* ./    
    cp ${alignment} ./alignment
    cp ${alignment}.crai ./alignment.crai

    # make publishDir, run GangSTR denovo
    GangSTR \
        --bam alignment \
        --bam-samps ${sample_id} \
        --samp-sex ${sample_sex} \
        --ref reference.fa \
        --regions str_panel.tsv \
        --out ${sample_id}_GangSTR \
        --quiet
        
    gzip ${sample_id}_GangSTR.vcf
    mkdir -p "${params.sample_dir}/${sample_id}/profiles/GangSTR/"
    """
}

// ----------------------------------------------------
//  Process for DumpSTR
// ----------------------------------------------------
//  Input:
//      Heterozygous genotypes for panel of STR loci
//      (from GangSTR, in vcf format)
//  Output:
//      Filtered vcf (PASS info field added to VCF 
//      entries that pass all filters)
// ----------------------------------------------------
process dumpSTR {
    publishDir "${params.sample_dir}/${sample_id}/profiles/GangSTR/filt/", mode: "move"
    scratch true
    stageOutMode 'move'

    conda "/cfs/earth/scratch/verb/.conda/envs/trtools"

    input:
    set sample_id, sample_vcf from GangSTR_out_ch
    
    output:
    path "${sample_id}_GangSTR_filt*"

    script:
    """
    mv ${sample_vcf} ./
    gunzip ${sample_id}_GangSTR.vcf.gz

    dumpSTR \
        --vcf ${sample_id}_GangSTR.vcf \
        --out ${sample_id}_GangSTR_filt \
        --gangstr-min-call-DP ${params.mindepth} \
        --gangstr-max-call-DP ${params.maxdepth} \
        --gangstr-min-call-Q ${params.minqual} \
        --gangstr-filter-spanbound-only \
        --gangstr-filter-badCI 

    # compress vcf file, index
    bgzip ${sample_id}_GangSTR_filt.vcf
    tabix -p vcf ${sample_id}_GangSTR_filt.vcf.gz
    """
}
