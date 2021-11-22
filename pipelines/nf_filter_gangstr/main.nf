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

// default filtering parameters to use in dumpSTR, can be overwritten at the command line
params.mindepth = 20 
params.maxdepth = 1000
params.minqual = 0.9

Channel
    .fromPath(params.sample_ids)
    .splitCsv(header: false, sep: "\t")
    .map{ row -> tuple(row[0]) }
    .set{ sample_id_ch }

process dumpSTR {
    publishDir "${params.sample_dir}/${sample_id}/filt/tral_panel/", mode: "move"
    scratch '/data/scratch/$SLURM_JOB_ID'
    stageOutMode 'move'

    conda "/cfs/earth/scratch/verb/.conda/envs/trtools"

    input:
    set sample_id from sample_id_ch
    
    output:
    path "${sample_id}_filt*" into dumpSTR_out_ch

    script:
    """
    # get process launch dir to copy results files into later
    # There has to be some nf variable I could use for this...
    process_launch_dir=\$PWD 

    # go to local scratch dir on node, copy input vcf file
    cd /data/scratch/\$SLURM_JOB_ID

    cp ${params.sample_dir}/${sample_id}/raw/tral_panel/${sample_id}.vcf.gz ./

    gunzip ${sample_id}.vcf.gz

    # make publishDir, run dumpSTR
    mkdir -p "${params.sample_dir}/${sample_id}/filt/tral_panel/"

    dumpSTR \
        --vcf ${sample_id}.vcf \
        --out ${sample_id}_filt \
        --gangstr-min-call-DP ${params.mindepth} \
        --gangstr-max-call-DP ${params.maxdepth} \
        --gangstr-min-call-Q ${params.minqual} \
        --gangstr-filter-spanbound-only \
        --gangstr-filter-badCI \
        --drop-filtered
    
    # generate separate vcf file that only contains non-reference entries
    awk '{if(/#/)print;else exit}' ${sample_id}_filt.vcf > ${sample_id}_filt_nonref.vcf
    grep -v "^#" ${sample_id}_filt.vcf | awk 'BEGIN { OFS="\t" } { if (\$5 != ".") { print \$0 } }' >> ${sample_id}_filt_nonref.vcf

    # compress vcf files, index
    bgzip ${sample_id}_filt.vcf
    tabix -p vcf ${sample_id}_filt.vcf

    bgzip ${sample_id}_filt_nonref.vcf
    tabix -p vcf ${sample_id}_filt_nonref.vcf

    # copy output files to dir that process launched from
    cp ${sample_id}_filt* \$process_launch_dir
    """
}
