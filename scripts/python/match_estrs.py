#!/usr/bin/env python3

def get_consensus_unit(unit: str) -> str:
    unit2 = unit + unit
    options = []

    for i in range(0, len(unit)):
         options.append( unit2[i : i+len(unit)] )

    return sorted(options)[0]


def main():
    fotsing = "/Users/maxverbiest/PhD/projects/genotype_strs/data/fotsing_etal_data_sheet_hg38.tsv"
    panel = "/Users/maxverbiest/PhD/data/str_panels/tral_and_perf_panel.tsv"

    with open(panel, 'r') as p:        
        with open(fotsing, 'r') as f:
            try:
                pSTR = next(p).strip().split("\t")
                p_start, p_end = int(pSTR[1]), int(pSTR[2])
                for eSTR in f:
                    eSTR = eSTR.strip().split("\t")
                    e_start, e_end = int(eSTR[1]), int(eSTR[2])
                    while pSTR[0] != eSTR[0]:
                        pSTR = next(p).strip().split("\t")
                        p_start, p_end = int(pSTR[1]), int(pSTR[2])
                    while p_end < e_start:
                        while pSTR[0] != eSTR[0]:
                            eSTR = next(f).strip().split("\t")
                            e_start, e_end = int(eSTR[1]), int(eSTR[2])
                        pSTR = next(p).strip().split("\t")
                        p_start, p_end = int(pSTR[1]), int(pSTR[2])
                        # while p_end < e_start:
                        #     pSTR = next(p).strip().split("\t")
                        #     p_start, p_end = int(pSTR[1]), int(pSTR[2])
                    if not p_start < e_end:
                        continue
                    p_cons = get_consensus_unit(pSTR[4])
                    print(("\t").join(eSTR))
                    
            except StopIteration:
                exit()

if __name__ == "__main__":
    main()
