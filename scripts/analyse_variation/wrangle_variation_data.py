import argparse

import numpy as np
import pandas as pd

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-e", "--healthy", type=str, required=True, help="Path to csv file containing genotype calls for healthy samples"
    )
    parser.add_argument(
        "-t", "--tumor", type=str, required=True, help="Path to csv file containing genotype calls for tumor samples"
    )
    parser.add_argument(
        "-o", "--output_file", type=str, required=True, help="File that the output data frame will be written to"
    )

    return parser.parse_args()

def to_wide_df(df: pd.DataFrame, sample_type: str) -> pd.DataFrame:
    if not sample_type in {"healthy", "tumor"}:
        raise ValueError("Sample type must be either 'healthy' or 'tumor'")
    df_to_dupe = (
        df[["patient", "tmp_id", "alt"]]
            .groupby(["patient", "tmp_id"])
            .filter(lambda x: len(x) == 1)
    )
    df_wide = (
        pd.concat([df_to_dupe, df[["patient", "tmp_id", "alt"]]])
            .sort_values(by=["patient", "tmp_id"])
            .reset_index(drop=True)
    )
    df_wide["allele"] = [f'allele_a_{sample_type}', f'allele_b_{sample_type}'] * int((df_wide.shape[0] / 2))
    df_wide = (
        df_wide
            .pivot(index=["patient", "tmp_id"], columns="allele", values="alt")
            .reset_index()
            .merge(df.drop(["patient", "alt"], axis=1).drop_duplicates(), how="left", on="tmp_id")
    )
    
    return df_wide
    

def main():
    args = cla_parser()
    
    df_healthy = pd.read_csv(args.healthy, sep=",").drop("sample_type", axis=1).drop_duplicates()
    patients = df_healthy.patient.unique()
    
    df_tumor = pd.read_csv(args.tumor, sep=",").drop("sample_type", axis=1).drop_duplicates()
    df_tumor = df_tumor[df_tumor["patient"].isin(patients)]

    df_healthy_wide = to_wide_df(df_healthy, sample_type="healthy")
    df_tumor_wide = to_wide_df(df_tumor, sample_type="tumor")

    df_merged_wide = df_healthy_wide.merge(
        df_tumor_wide[["patient", "tmp_id", "allele_a_tumor", "allele_b_tumor"]],
        how="left",
        on=["patient", "tmp_id"]
    ).dropna(subset=["allele_a_healthy", "allele_b_healthy", "allele_a_tumor", "allele_b_tumor"], axis=0)

    # reorder columns
    df_merged_wide = df_merged_wide[['tmp_id', 'repeat_id', 'chr', 'start', 'end', 'period', 'ref', 'patient', 'allele_a_healthy', 'allele_b_healthy', 'allele_a_tumor','allele_b_tumor']]

    # make sure relevant columns are int. repeat_id column should also be int, but it contains 'NA' values so conversion does not work
    int_cols = ["start", "end", "period", "ref", "allele_a_healthy", "allele_b_healthy", "allele_a_tumor", "allele_b_tumor"]
    df_merged_wide[int_cols] = df_merged_wide[int_cols].astype(int)

    df_merged_wide.to_csv(args.output_file, header=True, index=False)



if __name__ == "__main__":
    main()
