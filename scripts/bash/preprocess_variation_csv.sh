#!/usr/bin/env bash
set -eou pipefail

CSV_FILE="$1"

if [ -f "${CSV_FILE}" ]; then
    echo "patient,sample_type,repeat_id,chr,start,end,period,ref,alt,tmp_id"
    grep -v "^patient" ${CSV_FILE} | awk 'BEGIN{FS=","; OFS=","} {if($3 == "."){$3="NA"} }  {print $0","$4"_"$5}'
else 
    echo "$CSV_FILE does not exist."
    exit
fi
