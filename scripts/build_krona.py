#!/usr/bin/env python
import argparse
import re
import sys
from collections import defaultdict

import pandas as pd


def build_tax(tax_file, output):
    """
    Takes in the file to the order level in classification_txt rule
    and parses those to put them into a form then calls Krona Plot in order to build
    a hierarchical graph of the taxa present in the sample
    """
    # read in tax file to order level from build_tax_table rule
    test_tax = pd.read_table(tax_file)
    col_names = list(test_tax)[1:]
    taxonomies = ["Kingdom", "Phylum", "Class", "Order"]
    for i, item in enumerate(taxonomies, 0):
        test_tax[item] = test_tax.taxonomy_order.str.split("; ").str.get(i)
    # proposed cleanup for above
    # I don't think this proposal works so I tried a for loop above
    # test_tax[["Kingdom", "Phylum", "Class", "Order"]] = test_tax.taxonomy_order.str.split("; ").str[0:4]
    # generate individual files per sample
    for item in col_names:
        df_name = f"test_tax_{item}"
        df_name = test_tax[[f"{item}", "Kingdom", "Phylum", "Class", "Order"]].dropna(
            thresh=2
        )
        df_name.to_csv(
            output + "/" + f"{item}" + "_tax.txt", sep="\t", index=False, header=False
        )


def build_ec_dict(ec_file, dat_file):
    """
    Parses enzyme commission numbers in their hierarchy
    """

    with open(ec_file) as file, open(dat_file) as dat_file:
        d = defaultdict(list)
        # build hierarchy of ec numbers
        for item in file:
            toks = item.strip().partition(";")
            path = toks[0].replace(" ", "").strip()
            pieces = path.split("-")[0].split(".")[:-1]
            desc = toks[2][:-1]

            if len(pieces) == 1:
                d[path].append(desc)
                first_level = desc

            elif len(pieces) == 2:
                d[path].append(desc)
                second_level = desc
                d[path].insert(0, first_level)
            else:
                d[path].append(desc)
                d[path].insert(0, first_level)
                d[path].insert(1, second_level)

        # need to grab the lowest level of ec numbers
        for line in dat_file:
            if line.startswith("ID"):
                line = line.strip().split("  ")[1].strip()
                desc = next(dat_file).split("  ")[1].strip()
                first_levels = d[line.rpartition(".")[0] + ".-"]
                for level in first_levels:
                    d[line].append(level)
                d[line].append(desc)
    return d


def parse_ec_file(ec_file_from_summaries):
    """
    Split sequences that aligned to multiple ec numbers into last common ec number
    shared amongst matches
    """
    with open(ec_file_from_summaries) as ec_file_sum:

        next(ec_file_sum)
        new_ec_dict = {"ec": []}
        # split ec numbers that mapped to multiple
        for item in ec_file_sum:
            toks = item.partition("\t")[0]
            ec = toks.split(";")
            # print(ec)
            first_item = ec[0][:5]
            if len(ec) == 1:
                new_ec_dict["ec"].append(ec[0])
            else:
                n = 0
                for item in ec[1:]:
                    if item[:5] == first_item:
                        new_ec_dict["ec"].append(item[:5] + ".-")
                        # print('same',item,first_item)
                    elif item[:3] == first_item[:3]:
                        new_ec_dict["ec"].append(item[:3] + ".-.-")
                        # print('same to high level',item,first_item)
                    elif item[:1] == first_item[:1]:
                        new_ec_dict["ec"].append(item[:1] + ".-.-.-")
    return new_ec_dict


def build_ec_output(new_ec_dict, ec_file_from_summaries, parsed_dict, output):
    """
    Splits the ec number counts and hierarchy by sample
    """
    ec_replace = pd.DataFrame.from_dict(new_ec_dict)
    ec_table = pd.read_table(ec_file_from_summaries)
    col_names = list(ec_table)[1:]
    ec_table["ec"] = ec_replace["ec"]
    ec_hier_dict = pd.DataFrame.from_dict(parsed_dict, orient="index").reset_index()
    ec_hier_dict = ec_hier_dict.rename(columns={"index": "ec"})
    # print(ec_hier_dict.head(5))
    grouped_ec_tbl = ec_table.merge(ec_hier_dict, on="ec", how="inner")
    # print(grouped_ec_tbl.head(5))
    grouped_ec_tbl = grouped_ec_tbl.rename(
        columns={0: "level_1", 1: "level_2", 2: "level_3", 3: "level_4"}
    )
    # split the samples into separate files
    for item in col_names:
        df_name = f"test_tax_{item}"
        df_name = grouped_ec_tbl[
            [f"{item}", "level_1", "level_2", "level_3", "level_4"]
        ]

        df_name.to_csv(
            output + "/" + f"{item}" + "_ec.txt", sep="\t", index=False, header=False
        )


def main(tax_file, output, ec_file, ec_file_from_summaries, dat_file):
    if len(sys.argv) == 2:
        build_tax(tax_file, output)
    else:
        new_ec_dict = parse_ec_file(ec_file_from_summaries)
        parsed_dict = build_ec_dict(ec_file, dat_file)
        build_ec_output(new_ec_dict, ec_file_from_summaries, parsed_dict, output)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--tax-file",
        type=str,
        help="file that contains the order information from Kaiju alignments",
    )
    p.add_argument("--output", type=str, help="output directory name")
    p.add_argument(
        "--ec-file", type=str, help="converter file that contains full ec hierarchy"
    )
    p.add_argument(
        "--ec-file-from-summaries",
        type=str,
        help="file that contains the ec numbers from the diamond alignment",
    )
    p.add_argument(
        "--dat-file",
        type=str,
        help="enzyme.dat file that has the lowest level of ec classification",
    )

    args = p.parse_args()
    # logging.basicConfig(
    #     level=logging.INFO, datefmt="%Y-%m-%d %H:%M", format="[%(asctime)s] %(message)s"
    # )
    main(
        args.tax_file,
        args.output,
        args.ec_file,
        args.ec_file_from_summaries,
        args.dat_file,
    )
