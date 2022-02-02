import argparse 
import os
import sys

from cyvcf2 import VCF, Writer # conda activate trtools

def parse_cla():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-v", "--vcf", type=str, required=True, help="GangSTR VCF output file. DB ID and purity of repeat wil be added for each entry here (can be gzipped)"
    )
    parser.add_argument(
        "-b", "--bedl", type=str, required=False, help="Bed-like file of STR loci in GangSTR input format (1-based coordinates). Needs an additional column \
            containing a <db_id>:<repeat_units>:<str_purity>"
    )

    return parser.parse_args()

def process_gangstr_vcf(vcf: str, bedl: str=None):
    vcf = VCF(vcf)
    vcf.add_info_to_header({
        'ID': 'PURITY', 
        'Description': 'Fraction of nucleotides in STR locus that match the consensus unit', 
        'Type': 'Float', 
        'Number': '1'
    })
    vcf.add_info_to_header({
        'ID': 'LONGEST_P_STRETCH', 
        'Description': 'Length of the longest stretch of consecutive consensus unit for this TR locus', 
        'Type': 'Integer', 
        'Number': '1'
    })
    w = Writer("-", vcf)    # we will be writing to stdout
    
    try:            
        with open(bedl, "r") as f:
            variant = get_next_variant(vcf) # collect the first variant    
            for line in f:
                line_split = line.strip().split("\t")    
                if bedl_matches_vcf(variant, line_split):      
                    str_metainf = line_split[6].split(":")
                    variant.ID = str_metainf[0]
                    variant.INFO['PURITY'] = str_metainf[2]

                    longest_p_stretch = get_longest_p_stretch(line_split[4], str_metainf[1])
                    variant.INFO['LONGEST_P_STRETCH'] = longest_p_stretch

                    w.write_record(variant)
                    variant = get_next_variant(vcf)
    except StopIteration:
        # we are out of variants
        pass

    vcf.close()

def get_next_variant(vcf):
    variant = next(vcf)
    while not variant.INFO.get('RU'):
        print(f"WARNING: no RU found for variant at {variant.CHROM}, {variant.start + 1}. Skipping...", file=sys.stderr)
        variant = next(vcf)
    return variant

def bedl_matches_vcf(variant, line_split: list) -> bool:   
    """ Consecutive checks of whether the chromosome, start position and repeat unit
    in the supplied variant and (split + stripped) bedl line match. If all matches, return True,
    otherwise return False
    """    
    if not line_split[0] == variant.CHROM:
        return False
    if not int(line_split[1]) == variant.start + 1:
        return False    
    try: 
        if not variant.INFO.get('RU').upper() == line_split[4]:
            return False
    except AttributeError:
        return False
    return True

def get_longest_p_stretch(cons_unit: str, msa_string: str) -> int:
    cur_p_stretch = 0
    max_p_stretch = 0
    for unit in msa_string.split(","):
        if unit == cons_unit:
            cur_p_stretch += 1
            if cur_p_stretch > max_p_stretch:
                max_p_stretch = cur_p_stretch
        else:
            cur_p_stretch = 0
    return max_p_stretch

def main():
    args = parse_cla()
    if not os.path.isfile(args.vcf) or not os.path.isfile(args.bedl):
        raise ValueError(f"One or both of the specified input files do not exist")

    process_gangstr_vcf(args.vcf, args.bedl)



if __name__ == "__main__":
    main()
