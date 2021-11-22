#!/usr/bin/env python3
import argparse
import os

from cyvcf2 import VCF, Writer # conda activate trtools
import numpy as np

def cla_parser():
    parser = argparse.ArgumentParser()

    req = parser.add_argument_group("Required arguments")
    req.add_argument(
        "-v", "--vcf", type=str, required=True, 
        help="Sorted GangSTR output vcf file (1-indexed !!) containing STR locus genotypes. Only those that are located within\
            regions specified in the region file will be written to output"
    )
    req.add_argument(
        "-f", "--fold", type=float, required=True,
        help="Minimum fold difference of STR locus length required for locus to be retained.\
            e.g. if set to 2.0, only STRs where a genotype twice or more as long/short \
                as the reference will be retained (regardless of zygosity of call).\
                    Setting fold to 1.0 keep reference genotypes as well."
    )

    nonreq = parser.add_argument_group("Optional arguments")
    nonreq.add_argument(
        "-r", "--regions", type=str, required=False,
        help="Sorted bed file (0-indexed !!) containing \
            genomic regions of interest. All variant calls for loci not in or overlapping these \
                regions will be discarded"
    )    
    nonreq.add_argument(
        "-o", "--output", type=str, required=False,
        help="Filename to use for filtered output file. Will default to <input_file>_filt.<extension>\
            if nothing is specified"
    )

    return parser.parse_args()

def get_output_name(input_path: str, suffix: str=None) -> str:
    in_file = os.path.basename(input_path)
    if suffix:
        out_file = os.path.splitext(in_file)[0] + suffix + os.path.splitext(in_file)[1]
    else:
        out_file = os.path.splitext(in_file)[0] + "_filt" + os.path.splitext(in_file)[1]

    return os.path.join(os.path.dirname(input_path), out_file)


def main():
    args = cla_parser()

    if not args.output:
        output_file = get_output_name(args.vcf)
    else:
        output_file = args.output

    if not args.regions:
        vcf = VCF(args.vcf)
        w = Writer(output_file, vcf)    
        for variant in vcf:
            # FORMAT field REF, copy number in reference (e.g. 4)
            # INFO field REPCN, diploid genotype in sample (e.g. 4,5)
            ref_cn = variant.INFO.get("REF")
            for i in np.nditer(variant.format("REPCN")):
                if i >= (ref_cn * args.fold) or i <= (ref_cn / args.fold):
                    # print(variant, end="")
                    w.write_record(variant)
                    break
        vcf.close()
        w.close()
    else:
        vcf = VCF(args.vcf)
        w = Writer(output_file, vcf)

        variant = next(vcf)
        with open(args.regions, "r") as f:
            for line in f:
                line_split = line.strip().split("\t")
                region = {
                    "chrom": line_split[0],
                    "begin": int(line_split[1]),
                    "end": int(line_split[2])
                }
                try:
                    while variant.CHROM != region["chrom"]:
                        # we are on the wrong chromosome, keep getting neext variant
                        # until we are on the right one
                        variant = next(vcf)
                    
                    while variant.end < (region["begin"] + 1):
                        # we are upstream of region, keep getting neext variant until
                        # we are in the region
                        variant = next(vcf)

                    while variant.start <= (region["end"] + 1):
                        # we are in the region, keep writing variants until we
                        # go out of range
                        ref_cn = variant.INFO.get("REF")
                        for i in np.nditer(variant.format("REPCN")):
                            if i >= (ref_cn * args.fold) or i <= (ref_cn / args.fold):
                                w.write_record(variant)
                                break
                        variant = next(vcf)
                        
                except StopIteration:
                    # we are out of variants, break
                    break
        vcf.close()
        w.close()


if __name__ == "__main__":
    main()
