#!/usr/bin/bash    
#SBATCH --job-name=gangstr_dl
#SBATCH --partition=single
#SBATCH --qos=single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --export=ALL,LSFM_CLUSTER_KEEP_LOCAL_SCRATCH=true
#SBATCH --chdir=/cfs/earth/scratch/verb/jobs/genotyping/run_dir

set -euo pipefail

module load USS/2020 gcc/7.3.0
module load slurm
module load htslib openblas nlopt

# check invocation
usage="sbatch --export=ALL,LSFM_CLUSTER_KEEP_LOCAL_SCRATCH=true,cl_ftp_url=<target_ftp_url>,cl_sex=<sample_sex (M/F)> sbatch_download_gangstr.sh"

if [[ ! -v cl_ftp_url || ! -v cl_sex ]]; then
    echo "Error: need to set variables 'cl_sample_name' and 'cl_sex' at command line"
    echo "Usage: ${usage}"
    exit 1
elif [[ ${cl_sex} != "M" && ${cl_sex} != "F" ]]; then
    echo "Error: specify sample sex as 'M' or 'F'"
    echo "Usage: ${usage}"
    exit 1
fi

localapps="/cfs/earth/scratch/verb/localapps"
reference="/cfs/earth/scratch/verb/projects/genotype_strs/data/1000g_reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"
tr_regions="/cfs/earth/scratch/verb/projects/CRC_STRs/results/db/regions_full_units_gangstr.bedl"
out_dir="/cfs/earth/scratch/verb/projects/genotype_strs/results/raw"

cd $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH

cp {$reference,$reference.fai} ./
cp $tr_regions ./

wget -q $cl_ftp_url* 
sample_name=$(basename -s .final.cram $cl_ftp_url)
align_file=$(find ./ -name $sample_name.final.cram -type f)

$localapps/bin/GangSTR --bam $(basename $align_file) \
--bam-samps ${sample_name} \
--samp-sex ${cl_sex} \
--ref $(basename $reference) \
--regions $(basename $tr_regions) \
--out ./$sample_name 2> /dev/null #${sample_name}_gangstr.err

gzip "${sample_name}.vcf"
cp ${sample_name}.{insdata.tab,samplestats.tab,vcf.gz} $out_dir
# cp ${sample_name}_gangstr.err $out_dir
