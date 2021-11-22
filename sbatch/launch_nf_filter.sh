#!/usr/bin/bash    
#SBATCH --job-name=nf_main
#SBATCH --partition=single
#SBATCH --qos=single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=3G

nextflow run "/cfs/earth/scratch/verb/projects/genotype_strs/pipelines/nf_filter_gangstr/main.nf" \
    --sample_dir "/cfs/earth/scratch/verb/projects/genotype_strs/results/1000g/" \
    --sample_ids "/cfs/earth/scratch/verb/projects/genotype_strs/pipelines/nf_filter_gangstr/info/sample_ids.tsv"
