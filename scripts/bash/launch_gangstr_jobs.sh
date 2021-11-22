#!/usr/bin/bash

set eou -pipefail

# sbatch_script="/cfs/earth/scratch/verb/projects/genotype_strs/sbatch/sbatch_gangstr.sh"
sbatch_script="/cfs/earth/scratch/verb/projects/genotype_strs/sbatch/sbatch_download_gangstr.sh"

# Need to pipe a file into the script where each line contains '<sample_id>\t<sex>'
while IFS="" read -r p || [ -n "$p" ]
do
  line=( ${p} )
  if [[ ${#line[@]} -ne 2 ]]; then
    echo "Error: wrong imput file format detected"
    exit 1
  fi
  # echo ${line[0]}
  # echo ${line[1]}
  echo "Launching GangSTR run for sample: ${line[0]} (${line[1]})"
  sbatch --nodelist=node020 --export=ALL,LSFM_CLUSTER_KEEP_LOCAL_SCRATCH=true,cl_ftp_url=${line[0]},cl_sex=${line[1]} ${sbatch_script}
  
  # sbatch --nodelist=node020 --export=ALL,LSFM_CLUSTER_KEEP_LOCAL_SCRATCH=true,cl_sample_name=${line[0]},cl_sex=${line[1]} ${sbatch_script}
  sleep 60
done
