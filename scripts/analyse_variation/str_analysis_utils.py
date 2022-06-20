#!/usr/bin/env python3
import pandas as pd

def get_difference_options(row) -> int:
    a_to_a = abs(row["allele_a_healthy"] - row["allele_a_tumor"]) + abs(row["allele_b_healthy"] -  row["allele_b_tumor"])
    a_to_b = abs(row["allele_b_healthy"] - row["allele_b_tumor"]) + abs(row["allele_b_healthy"] -  row["allele_a_tumor"])

    return a_to_a, a_to_b

def is_frameshift(allele_pair: tuple, period: int) -> bool:
    if (abs(allele_pair[0] - allele_pair[1]) * period) % 3 == 0:
        return False
    return True

def count_frameshifts(row) -> int:
    a_to_a, a_to_b = get_difference_options(row)

    if a_to_a <= a_to_b:
        pair1 = (row["allele_A_healthy"], row["allele_A_tumor"])
        pair2 = (row["allele_B_healthy"], row["allele_B_tumor"])
    else:
        pair1 = (row["allele_A_healthy"], row["allele_B_tumor"])
        pair2 = (row["allele_B_healthy"], row["allele_A_tumor"])
    
    return is_frameshift(pair1, row["period"]) + is_frameshift(pair2, row["period"])


def filter_variation_df(df_variation: pd.DataFrame, df_str_info: pd.DataFrame, min_period=0) -> pd.DataFrame:
    df_variation = (
        df_variation
            .merge(df_str_info, how="left", on="tmp_id", suffixes=(None, "_tmp"))
            .query(f"not in_segdup and neighbour_type == 'no_neighbour' and period >= {min_period}")
            .loc[:, df_variation.columns] 
    )
    
    return df_variation
