"""
Generates a tab separated file of the desired functional information
"""
import argparse
import gzip
import json
import logging
import os
import pandas as pd

gzopen = lambda f: gzip.open(f, "rt") if f.endswith(".gz") else open(f)
TAX_LEVELS = {
    "superkingdom": 1,
    "kingdom": 1,
    "phylum": 2,
    "class": 3,
    "order": 4,
    "family": 5,
    "genus": 6,
    "species": 7,
}
FUNCTION_COLS = ["ko", "ec", "product"]


# def parse_kegg_json(json_file):
#     kegg_dict = dict()
#     with gzopen(json_file) as filehandle:
#         json_class = json.load(filehandle)
#         for root in json_class["children"]:
#             for level_1 in root["children"]:
#                 for level_2 in level_1["children"]:
#                     try:
#                         for level_3 in level_2["children"]:
#                             ko = "ko:" + level_3["name"].partition(" ")[0]
#                             if ko in kegg_dict:
#                                 kegg_dict[ko]["level_1"] += ";" + root["name"]
#                                 kegg_dict[ko]["level_2"] += ";" + level_1["name"]
#                                 kegg_dict[ko]["level_3"] += ";" + level_2["name"]
#                             else:
#                                 kegg_dict[ko] = {
#                                     "level_1": root["name"],
#                                     "level_2": level_1["name"],
#                                     "level_3": level_2["name"],
#                                 }
#                     except:
#                         ko = "ko:K" + level_2["name"].partition(" ")[0]
#                         if ko in kegg_dict:
#                             kegg_dict[ko]["level_1"] += ";" + root["name"]
#                             kegg_dict[ko]["level_2"] += ";" + level_1["name"]
#                             kegg_dict[ko]["level_3"] += ";" + level_2["name"]
#                         else:
#                             kegg_dict[ko] = {
#                                 "level_1": root["name"],
#                                 "level_2": level_1["name"],
#                                 "level_3": level_2["name"],
#                             }
#     df = pd.DataFrame.from_dict(kegg_dict, orient="index").reset_index()
#     df = df.rename(columns={"index": "ko"})
#     return df[["ko", "level_1", "level_2", "level_3"]]


def df_from_classifications(tbl_path, group_on, split_idx, min_perc_id, min_len):
    #metric_cols = ["aa_percent_id", "aa_alignment_length", "kaiju_alignment_length"]
    #sample = os.path.basename(tbl_path).partition("_classifications.txt")[0]
    # grab columns of interest
    logging.debug(f"Parsing {tbl_path}")
    df = pd.read_table(tbl_path, usecols=group_on + metric_cols)
    logging.debug(f"Initial table size: {len(df)}")
    df = df.query("kaiju_length != 0")
    logging.debug(f"After filtering table size: {len(df)}")
    # functional assignment filters
    if any(FUNCTION_COLS) in group_on:
        df = df[
            (df["aa_percent_id"] > min_perc_id) & (df["aa_alignment_length"] > min_len)
        ]
        logging.debug(f"After functional filters: {len(df)}")
    # taxonomy assignment filters
    if "kaiju_classification" in group_on:
        df = df[(df["kaiju_alignment_length"] > min_len)]
        logging.debug(f"After taxonomy filters: {len(df)}")
        # with a split index of 2 (phylum), converts:
        # Bacteria; Proteobacteria; NA; NA; NA;  --->  Bacteria; Proteobacteria
        df["kaiju_classification"] = df["kaiju_classification"].apply(
            lambda x: "; ".join(x.split("; ", split_idx)[0:split_idx])
        )
    df = df.drop(metric_cols, axis=1)
    df = df.groupby(group_on).size().reset_index(name=sample)
    logging.debug(f"Final table length: {len(df)}")
    return df


def main(json_file, output, tables, tax_level, min_perc_id, min_len, group_on):
    hmm_libs=["tigrfams","hamap","dbcan","kaiju"]
    hmms = [df.filter(like=hmm_lib).columns.values.tolist() for hmm_lib in hmm_libs]
    samples = df.filter(regex="[0-9]+").columns.values.tolist()
    tax_split_level = TAX_LEVELS[tax_level]
    sample_df = None
    logging.debug(f"Preparing to parse {len(tables)} tables")
    for tbl in tables:
        logging.info(f"Parsing {tbl}")
        df = df_from_classifications(
            tbl, group_on, tax_split_level, min_perc_id, min_len
        )
        if sample_df is None:
            sample_df = df.copy()
            continue
        sample_df = sample_df.merge(df, on=group_on, how="outer").fillna("0")
    sample_df = sample_df.rename(
        columns={"kaiju_classification": f"taxonomy_{tax_level}"}
    )
    if "ko" in group_on:
        logging.info(f"Adding KEGG hierarchy from {json_file}")
        kegg_pd = parse_kegg_json(json_file)
        logging.debug(
            f"Table shape prior to merging with KEGG hierarchy: {sample_df.shape}"
        )
        sample_df = sample_df.merge(kegg_pd, on="ko", how="inner")
        logging.debug(f"After performing inner join: {sample_df.shape}")
    sample_df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("json", help="KEGG hierarchy")
    parser.add_argument("output", help="Output table")
    parser.add_argument("tables", nargs="+", help="classifications.txt")
    parser.add_argument(
        "--tax-level",
        choices=[
            "superkingdom",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ],
        default="phylum",
        help="tax level to show",
    )
    parser.add_argument(
        "--min-evalue",
        type=int,
        default=0.001,
        help="lowest evalue to retain per sequence per sample",
    )
    parser.add_argument(
        "--min-score",
        type=int,
        default=40,
        help="lowest hmmsearch score to retain per sequence per sample",
    )
    parser.add_argument(
        "--group-on",
        # choices=["kaiju_classification", "ko", "ec", "product"],
        nargs="+",
        default="kaiju_taxonomy",
        help="cluster tables on taxonomy and/or function; headers must match in tables",
    )
    parser.add_argument(
        "--verbose", action="store_true", help="increase output verbosity"
    )
    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    main(
        args.json,
        args.output,
        args.tables,
        args.tax_level,
        args.min_id,
        args.min_len,
        args.group_on,
    )
