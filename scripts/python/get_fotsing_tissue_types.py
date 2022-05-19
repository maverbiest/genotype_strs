#!/usr/bin/env python3

def main():
    all_tissues = set()
    with open("/Users/maxverbiest/PhD/projects/genotype_strs/data/fotsing_etal_data_sheet.csv", 'r') as f:
        next(f) # skip header
        for line in f:
            line = line.strip().split(',')
            tissue_info = line[7].split(";")            
            tissues = {i.split("_")[0] for i in tissue_info}
            all_tissues = all_tissues.union(tissues)
    for i in all_tissues:
        print(i)


if __name__ == "__main__":
    main()
