#!/usr/bin/bash    
#SBATCH --job-name=gangstr
#SBATCH --partition=single
#SBATCH --qos=single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --export=ALL,LSFM_CLUSTER_KEEP_LOCAL_SCRATCH=true
#SBATCH --chdir=/cfs/earth/scratch/verb/jobs/genotyping/run_dir

set -euo pipefail
# set -euxo pipefail

# earth cluster modules
module load USS/2020 gcc/7.3.0
module load slurm

# GangSTR requirements
module load htslib openblas nlopt

# Check command line argmuments (sort of)
usage="sbatch --export=ALL,LSFM_CLUSTER_KEEP_LOCAL_SCRATCH=true,cl_sample_name=<sample_name>,cl_sex=<sample_sex (M/F)> sbatch_gangstr.sh"

if [[ ! -v cl_sample_name || ! -v cl_sex ]]; then
    echo "Error: need to set variables 'cl_sample_name' and 'cl_sex' at command line"
    echo "Usage: ${usage}"
    exit 1
elif [[ ${cl_sex} != "M" && ${cl_sex} != "F" ]]; then
    echo "Error: specify sample sex as 'M' or 'F'"
    echo "Usage: ${usage}"
    exit 1
fi

# initialize variables
localapps="/cfs/earth/scratch/verb/localapps"
alignment_dir="/cfs/earth/scratch/verb/projects/download_cancer_alignments/data_out/1000g"
reference="/cfs/earth/scratch/verb/projects/genotype_strs/data/1000g_reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"
# tr_regions="/cfs/earth/scratch/verb/projects/genotype_strs/data/str_loci/gangstr_hg38_ver13.bedl"
tr_regions="/cfs/earth/scratch/verb/projects/CRC_STRs/results/db/regions_full_units_gangstr.bedl"
out_dir="/cfs/earth/scratch/verb/projects/genotype_strs/results/raw/gangstr_panel"

# find alignment, does it exist and does it have an index?
align_file=$(find ${alignment_dir} -name "*${cl_sample_name}*.cram" -type f)
if [[ ! -f ${align_file} || ! -f ${align_file}.crai ]]; then
    echo "Error: alignment file '${align_file}' is not a file or has no index"
    exit 1
fi

cd $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH

# copy required data to local scratch on node
cp {${reference},${reference}.fai} ./
cp {${align_file},${align_file}.crai} ./
cp ${tr_regions} ./

# run GangSTR, copy output files to head node
$localapps/bin/GangSTR --bam $(basename ${align_file}) \
--bam-samps ${cl_sample_name} \
--samp-sex ${cl_sex} \
--ref $(basename ${reference}) \
--regions $(basename ${tr_regions}) \
--out ./$cl_sample_name 2> ${cl_sample_name}_gangstr.err

cp ${cl_sample_name}.{insdata.tab,samplestats.tab,vcf} $out_dir
# cp ${cl_sample_name}_gangstr.err $out_dir
