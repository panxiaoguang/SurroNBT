import os
import pandas as pd
import numpy as np


def enrich_depletion(case_stat_file, control_stat_file, outfile):
    if not os.path.exists(case_stat_file) or not os.path.exists(control_stat_file):
        raise FileNotFoundError(f"{case_stat_file} dosen't exists!")
    df = pd.read_csv(case_stat_file, sep="\t")
    df_2 = pd.read_csv(control_stat_file, sep="\t")
    case_total_cleans = df.clean.sum()
    control_total_cleans = df_2.clean.sum()
    df.depth = df.clean / case_total_cleans * 1000000
    df_2.depth = df_2.clean / control_total_cleans * 1000000
    result = pd.merge(df, df_2, on="library", how="outer")
    result.fillna(0, inplace=True)
    result["filter"] = (
        (result["depth_x"] == 0)
        | (result["depth_y"] == 0)
        | (np.abs(np.log2(result["depth_x"] / result["depth_y"])) > 2)
    )
    result.to_csv(outfile, sep="\t", header=True, index=False)


if __name__ == "__main__":
    case_stat_file = snakemake.input["Case"] + "/statistic.txt"
    control_stat_file = snakemake.input["Control"] + "/statistic.txt"
    enrich_depletion(case_stat_file, control_stat_file, snakemake.output[0])
