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
FUNCTION_COLS = {
    "tigrfams": [],
    "hamap":[],
    "dbcan":[]
}




# def df_from_classifications(tbl_path, group_on, split_idx, min_perc_id, min_len):
#     #metric_cols = ["aa_percent_id", "aa_alignment_length", "kaiju_alignment_length"]
#     #sample = os.path.basename(tbl_path).partition("_classifications.txt")[0]
#     # grab columns of interest
#     logging.debug(f"Parsing {tbl_path}")
#     df = pd.read_table(tbl_path, usecols=group_on + metric_cols)
#     logging.debug(f"Initial table size: {len(df)}")
#     df = df.query("kaiju_length != 0")
#     logging.debug(f"After filtering table size: {len(df)}")
#     # functional assignment filters
#     if any(FUNCTION_COLS) in group_on:
#         df = df[
#             (df["aa_percent_id"] > min_perc_id) & (df["aa_alignment_length"] > min_len)
#         ]
#         logging.debug(f"After functional filters: {len(df)}")
#     # taxonomy assignment filters
#     if "kaiju_taxonomy" in group_on:
#         df = df[(df["kaiju_length"] > min_len)]
#         logging.debug(f"After taxonomy filters: {len(df)}")
#         # with a split index of 2 (phylum), converts:
#         # Bacteria; Proteobacteria; NA; NA; NA;  --->  Bacteria; Proteobacteria
#         df["kaiju_taxonomy"] = df["kaiju_taxonomy"].apply(
#             lambda x: "; ".join(x.split("; ", split_idx)[0:split_idx])
#         )
#     df = df.drop(metric_cols, axis=1)
#     df = df.groupby(group_on).size().reset_index(name=sample)
#     logging.debug(f"Final table length: {len(df)}")
#     return df

def single_output(df,tax_level,samples,min_evlaue,min_score,min_len,group_on,output):
    split_idx = TAX_LEVELS[tax_level]
    # group_one=["tigrfams","hamap","dbcan","kaiju"]
    # for item in group_one:
    hmm_cols = df.filter(like=group_on[0]).columns.values.tolist()
    #TODO might not need
    test_df = df[hmm_cols+samples].copy()
    if "kaiju" in item:
        test_df = test_df[(df["kaiju_length"] > min_len)]
        test_df["kaiju_taxonomy"] = df["kaiju_taxonomy"].apply(lambda x: "; ".join(x.split("; ", split_idx)[0:split_idx]))
        test_df = test_df.drop(columns=["kaiju_length"]).groupby(["kaiju_taxonomy"]).sum(axis=1).reset_index()
        test_df = test_df.rename(
            columns={"kaiju_taxonomy": f"taxonomy_{tax_level}"}
        )
    else:
        metric_cols = test_df.filter(items=[group_on[0]+"_score",group_on[0]+"_evalue"]).columns.values.tolist()
        test_df = test_df[(test_df[group_on[0].lower()+"_score"] > min_score) & (test_df[group_on[0].lower()+"_evalue"] < min_evalue)]
        test_df = test_df.drop(columns=metric_cols).groupby(hmm_cols[0:3]).sum(axis=1).reset_index().head()
    test_df.to_csv(item+output, sep="\t", index=False)


def combined_output(df,tax_level,samples,min_evalue,min_score,min_len,group_on,output):
    split_idx = TAX_LEVELS[tax_level]
    # group_one=["tigrfams","hamap","dbcan"]
    kaiju_cols = df.filter(like="kaiju").columns.values.tolist()
    # for item in group_one:
    hmm_cols = df.filter(like=group_on[0].lower()).columns.values.tolist()
    #TODO might not need
    test_df = df[kaiju_cols+hmm_cols+samples].copy()
    test_df = test_df[(df["kaiju_length"] > min_len) & (test_df[group_on[0]+"_score"] > min_score) & (test_df[group_on[0]+"_evalue"] < min_evalue)]
    test_df["kaiju_taxonomy"] = df["kaiju_taxonomy"].apply(lambda x: "; ".join(x.split("; ", split_idx)[0:split_idx]))
    metric_cols = test_df.filter(items=[group_on[0].lower()+"_score",group_on[0].lower()+"_evalue"]).columns.values.tolist()
    metric_cols.append("kaiju_length")
    test_df = test_df.drop(columns=metric_cols).groupby(["kaiju_taxonomy"]+hmm_cols[0:3]).sum(axis=1).reset_index()
    test_df = test_df.rename(
        columns={"kaiju_taxonomy": f"taxonomy_{tax_level}"}
    )
    test_df.to_csv(output, sep="\t", index=False)


def main(output, table, tax_level, min_evalue, min_score, min_len,group_on):
    df = pd.read_table(table)
    # remove the rows w/ no information
    df = df.dropna(
    subset=[
        "kaiju_taxonomy", "tigrfams_ec", "tigrfams_gene", "tigrfams_product",
        "tigrfams_score", "tigrfams_evalue", "hamap_ec", "hamap_gene",
        "hamap_product", "hamap_score", "hamap_evalue", "dbcan_ec",
        "dbcan_enzyme_class", "dbcan_enzyme_class_subfamily"
    ],
    thresh=1)
    samples = df.filter(regex="[0-9]+").columns.values.tolist()
    if len(group_on) > 1:
        combined_output(df,tax_level,min_evalue, min_score, min_len, group_on,output)
    else:
        single_output(df,tax_level,samples,min_evlaue,min_score,min_len,group_on,output)
    # split_idx = TAX_LEVELS[tax_level]
    # group_one=["tigrfams","hamap","dbcan","kaiju"]
    # for item in group_one:
    #     print(item)
    #     use_these = df.filter(like=item).columns.values.tolist()
    #     print(use_these)
    #     #TODO might not need
    #     test_df = df[use_these+samples].copy()
    #     if "kaiju" in item:
    #         test_df = test_df[(df["kaiju_length"] > min_len)]
    #         test_df["kaiju_taxonomy"] = df["kaiju_taxonomy"].apply(lambda x: "; ".join(x.split("; ", split_idx)[0:split_idx]))
    #         test_df = test_df.drop(columns=["kaiju_length"]).groupby(["kaiju_taxonomy"]).sum(axis=1).reset_index()
    #         test_df = test_df.rename(
    #             columns={"kaiju_taxonomy": f"taxonomy_{tax_level}"}
    #         )
    #     else:
    #         metric_cols = test_df.filter(items=[item+"_score",item+"_evalue"]).columns.values.tolist()
    #         test_df = test_df[(test_df[item+"_score"] > min_score) & (test_df[item+"_evalue"] < min_evalue)]
    #         test_df = test_df.drop(columns=metric_cols).groupby(use_these[0:3]).sum(axis=1).reset_index().head()
    #     test_df.to_csv(item+output, sep="\t", index=False)
    # tax_split_level = TAX_LEVELS[tax_level]
    # sample_df = None
    # logging.debug(f"Preparing to parse {len(tables)} tables")
    # for hmm in tables:
    #     logging.info(f"Parsing {tbl}")
    #     df = df_from_classifications(
    #         tbl, group_on, tax_split_level, min_perc_id, min_len
    #     )
    #     if sample_df is None:
    #         sample_df = df.copy()
    #         continue
    #     sample_df = sample_df.merge(df, on=group_on, how="outer").fillna("0")
    # sample_df = sample_df.rename(
    #     columns={"kaiju_classification": f"taxonomy_{tax_level}"}
    # )
    # if "ko" in group_on:
    #     logging.info(f"Adding KEGG hierarchy from {json_file}")
    #     kegg_pd = parse_kegg_json(json_file)
    #     logging.debug(
    #         f"Table shape prior to merging with KEGG hierarchy: {sample_df.shape}"
    #     )
    #     sample_df = sample_df.merge(kegg_pd, on="ko", how="inner")
    #     logging.debug(f"After performing inner join: {sample_df.shape}")
    # sample_df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", help="Output table")
    parser.add_argument("table",help="annotations.txt")
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
        "--min-len",
        type=int,
        default=0.001,
        help="minimum kaiju alignment length to retain",
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
        args.output,
        args.table,
        args.tax_level,
        args.min_evalue,
        args.min_score,
        args.min_len,
        args.group_on,
    )
