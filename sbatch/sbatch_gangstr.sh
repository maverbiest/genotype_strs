#!/usr/bin/bash    
#SBATCH --job-name=gangstr
#SBATCH --partition=single
#SBATCH --qos=single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --export=ALL,LSFM_CLUSTER_KEEP_LOCAL_SCRATCH=true
#SBATCH --chdir=/cfs/earth/scratch/verb/jobs/genotyping/run_dir

set -e
set -u
set -o pipefail

module load USS/2020 gcc/7.3.0
module load slurm
module load htslib openblas nlopt

# variable ALIGN_FILE must be specified at the command line
if [[ ! -f $ALIGN_FILE || ! -f $ALIGN_FILE.bai ]]
then
    echo "Specified alignment file '$ALIGN_FILE' is not a file or has no index"
    exit 1
fi

LOCALAPPS="/cfs/earth/scratch/verb/localapps"
REFERENCE="/cfs/earth/scratch/verb/projects/CRC_STRs/data/TCGA_reference_files/GRCh38.d1.vd1.fa"
TR_REGIONS="/cfs/earth/scratch/verb/projects/CRC_STRs/results/db/regions_consensus_lai_sun_thresh.bed"
OUT_DIR="/cfs/earth/scratch/verb/projects/genotype_strs/results/raw/subset_nontarg_new_regions"

cd $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH

cp {$REFERENCE,$REFERENCE.fai} ./
cp {$ALIGN_FILE,$ALIGN_FILE.bai} ./
cp $TR_REGIONS ./

SAMPLE_NAME=$(basename -s _genes.bam $ALIGN_FILE)

$LOCALAPPS/bin/GangSTR --bam $(basename $ALIGN_FILE) \
--ref $(basename $REFERENCE) \
--regions $(basename $TR_REGIONS) \
--out ./$SAMPLE_NAME --coverage 36.5165 2> ${SAMPLE_NAME}_gangstr.err

cp ${SAMPLE_NAME}.{insdata.tab,samplestats.tab,vcf} $OUT_DIR
cp ${SAMPLE_NAME}_gangstr.err $OUT_DIR

# Copy aligned protein coding reads and index from local node scratch to head node
# cp $LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH/*.bam* $OUT_DIR
