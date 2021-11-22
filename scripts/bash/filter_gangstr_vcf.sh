#!/usr/bin/bash
# Quick & dirty script to filter gangstr output vcf file, should incorportate into NF pipeline 
# at some point

# extend command line arguments
set $@ eou -pipefail 

if [ ! $# -eq 3 ] 
then
    echo "Invalid number of command line arguments!"
    echo "Usage: $0 gangstr_output.vcf"
    exit 1
fi

outname=$(basename -s .vcf $1)"_nonref.vcf"

echo "Filtering file ${outname}... "

awk '{if(/#/)print;else exit}' $1 > ${outname}
grep -v "^#" $1  | awk 'BEGIN { OFS="\t" } { if ($5 != ".") { print $0 } }' >> ${outname}

echo "Counting records in file ${outname}..."
grep -v "^#" ${i} | wc -l
