#!/usr/bin/env python3
import argparse
import sys

from cyvcf2 import VCF # conda activate trtools
import numpy as np


def parse_cla():
    parser = argparse.ArgumentParser()

    req = parser.add_argument_group("Required arguments")
    req.add_argument(
        "-v", "--vcf", type=str, required=True, 
        help="Sorted GangSTR output vcf file (1-indexed !!) containing STR locus genotypes of multiple samples merged\
            into one file using mergeSTR (TRTOOLS). Each entry in the VCF will be converted to tsv format and printed to stdout."
    )

    return parser.parse_args()

def vcf_to_tsv(vcf: VCF):
    header = "\t".join([
        "chrom",
        "start",
        "end",
        "motif",
        "motif_len",
        "ref_cn",
        "purity",
        "longest_p_stretch",
        "call_rate",
        "diffs_less",
        "diffs_add",
        "diffs_tot",        
        "avg_diffs"  
    ])
    print(header)
    for variant in vcf:
        ref_cn = variant.INFO.get("REF")
        if not ref_cn:
            # Very strange that this happens sometimes, look into this?
            print(f"No ref_cn found for variant at {variant.CHROM}, {variant.start + 1}. skipping...", file=sys.stderr)    
            continue 
        
        out_line = variant_to_tsv(variant=variant)

        diffs_less = 0  # counter for number of fewer repeat units compared to reference allele, across all samples
        diffs_add = 0   # counter for number of additional repeat units compared to reference allele, across all samples        
        if variant.ALT:
            for i in np.nditer(variant.format("REPCN")):
                if not i >= 0:
                    # values will take -inf if the allele has same cn as ref
                    continue
                # add difference in unit number between reference and current allele to appropriate counter
                if ref_cn > i:
                    diffs_less += ref_cn - i
                elif ref_cn < i:
                    diffs_add += i - ref_cn

        # total number of differences between reference allele and sample alleles, summed over all samples
        diffs = diffs_add + diffs_less
        out_line += f"\t{diffs_less}\t{diffs_add}\t{diffs}\t{round(diffs / (variant.call_rate * len(vcf.samples)), 2)}"
        
        print(out_line)

def variant_to_tsv(variant) -> str:
    """ Take a variant (cyvcf2.Variant) and convert it to a line to put on tsv file.

    Parameters
    variant: cyvcf2.Variant     variant from vcf file to be converted into tsv line

    Returns
    out_line: str               String of tab-separated values representing the variant.
                                Will contain "chrom", "start", "end", "motif", "motif_len", "ref_cn", "call_rate"
    """
    out_line = "\t".join([
            variant.CHROM,
            str(variant.start + 1),
            str(variant.end),
            variant.INFO.get("RU").upper(),
            str(len(variant.INFO.get("RU"))),
            str(variant.INFO.get("REF")),
            str(round(variant.INFO.get("PURITY"), 2)),
            str(variant.INFO.get("LONGEST_P_STRETCH")),
            str(variant.call_rate)
        ])
    return out_line


def main():
    args = parse_cla()

    vcf = args.vcf
    
    # init VCF parser
    vcf_to_tsv(VCF(vcf))



if __name__ == "__main__":
    main()
