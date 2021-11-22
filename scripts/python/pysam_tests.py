#!/usr/bin/env python3

import pysam


def main():
    reference = "/cfs/earth/scratch/verb/projects/genotype_strs/data/1000g_reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    cram_file = "/cfs/earth/scratch/verb/projects/genotype_strs/results/testing/HG00118.final.cram"

    with pysam.AlignmentFile(cram_file, 'rc', reference_filename=reference) as cram_reader:
        if not cram_reader.check_index():
            raise FileNotFoundError("No index was found for specified alignment file")

        for read in cram_reader.fetch("chr1", 39477030, 39477091):
            print(read)
        region_offset = 1000
        region_length = 5000

        chromosome = "chr1"
        begin = 64096 + region_offset
        end = 64102 + region_offset + region_length
        for align_col in cram_reader.pileup(chromosome, begin, end):
            if not begin <= align_col.pos <= end:
                continue
            print(f"\nCoverage at base {align_col.pos} = {align_col.n}")
            for align_read in align_col.pileups:
                if not align_read.is_del and not align_read.is_refskip:
                    print(f"\tbase in read {align_read.alignment.query_name} = {align_read.alignment.query_sequence[align_read.query_position]}")
                    break

if __name__ == "__main__":
    main()
