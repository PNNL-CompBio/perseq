#!/usr/bin/env python
import argparse
import gzip
from collections import defaultdict

import pandas as pd

gzopen = lambda f: gzip.open(f, "rt") if f.endswith(".gz") else open(f)


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
    # split taxonomy column into hierarchy
    test_tax[
        ["Kingdom", "Phylum", "Class", "Order"]
    ] = test_tax.taxonomy_order.str.split("; ", expand=True)
    # generate individual files per sample
    for item in col_names:
        df_name = f"test_tax_{item}"
        df_name = test_tax[[f"{item}", "Kingdom", "Phylum", "Class", "Order"]].dropna(
            thresh=2
        )
        df_name.to_csv(
            output + "/" + f"{item}" + "_tax.txt", sep="\t", index=False, header=False
        )


def build_ec_dict(ec_filename, dat_filename):
    """
    Parses enzyme commission numbers in their hierarchy

    From the config file:

    # curl ftp://ftp.expasy.org/databases/enzyme/enzyme.dat | gzip > enzyme.dat.gz
    enzyme_dat_file: resources/enzyme.dat.gz
    # curl ftp://ftp.expasy.org/databases/enzyme/enzclass.txt | gzip > enzclass.txt.gz
    ec_converter: resources/enzclass.txt.gz
    """
    with gzopen(ec_filename) as ec_file, gzopen(dat_filename) as dat_file:
        d = defaultdict(list)
        # build hierarchy of ec numbers
        for line in ec_file:
            if not line[0].isdigit():
                continue

            # all ECs are 9 characters wide and include spaces as padding
            # and the lowest level is never defined
            # 1. 1.99.-    With other acceptors." -> 1. 1.99.-
            path = line[0:9].replace(" ", "")
            # 1. 1.99.- -> ['1', ' 1', '99']
            pieces = path.split("-")[0].split(".")[:-1]
            # 1. 1.99.-    With other acceptors." -> With other acceptors
            desc = line[9:].strip(" \r\n.")
            # depends on the input always being sorted by increasing complexity
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
                # line = line.strip().split("  ")[1].strip()
                path = line.strip().partition("  ")[-1].strip()
                # advance to the next line
                desc_line = next(dat_file)
                assert desc_line.startswith("DE")
                desc = desc_line.split("  ", 2)[1].strip()
                d[path] = [i for i in d[path.rpartition(".")[0] + ".-"]] + [desc]
    return d


def parse_ec_file(ec_file_from_summaries):
    """
    Split sequences that aligned to multiple ec numbers into last common ec number
    shared amongst matches
    """
    with gzopen(ec_file_from_summaries) as ec_file_sum:

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
    if tax_file is not None:
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
    p.add_argument("output", type=str, help="output directory name")
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
