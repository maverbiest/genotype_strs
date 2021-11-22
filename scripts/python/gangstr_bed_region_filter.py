#!/usr/bin/env python3

import argparse 

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-r", "--regions", type=str, required=True, 
        help="Sorted bed file (0-indexed !!) containing \
            genomic regions for which STR loci will be retained in the \
            GangSTR input file"
    )
    parser.add_argument(
        "-l", "--loci", type=str, required=True, 
        help="Sorted bed-like file (1-indexed !!) containing STR loci. Only those that are located within\
            regions specified in the region file will be written to output"
    )

    return parser.parse_args()


def get_locus(file_generator) -> dict:
    locus = next(file_generator)
    locus_split = locus.strip().split("\t")
    locus_d = {
        "full_line": locus,
        "chrom": locus_split[0],
        "begin": int(locus_split[1]) - 1,
        "end": int(locus_split[2]) - 1
    }
    return locus_d


def main():
    args = cla_parser()

    loci = open(args.loci, "r")

    # get the first locus
    locus = get_locus(loci)

    with open(args.regions, "r") as f:
        for line in f:
            line_split = line.strip().split("\t")
            region = {
                "chrom": line_split[0],
                "begin": int(line_split[1]),
                "end": int(line_split[2])
            }
            try:
                while locus["chrom"] != region["chrom"]:
                    # we are on the wrong chromosome, keep getting new locus
                    # until we are on the right one
                    locus = get_locus(loci)
                
                while locus["end"] < region["begin"]:
                    # we are upstream of region, keep getting new locus until
                    # we are in the region
                    locus = get_locus(loci)

                while locus["begin"] <= region["end"]:
                    # we are in the region, keep printing loci until we
                    # go out of range
                    print(locus["full_line"], end="")
                    locus = get_locus(loci)
                    
            except StopIteration:
                # we are out of STR loci, break
                break

    loci.close()

if __name__ == "__main__":
    main()
