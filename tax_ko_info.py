#!/usr/bin/env python

import os
import pandas as pd
import json
import argparse

# @click.command()
# @click.argument('json_file')


def parse_kegg_json(json_file_path):
    with open(json_file_path) as json_file:
        json_class = json.load(json_file)
        kegg_dict = dict()
        for root in json_class["children"]:
            #         print("1", children)
            #         print("1 name", test_list["name"])
            for level_1 in root["children"]:
                # print("sg1", subgroup_1)
                for level_2 in level_1["children"]:
                    # print(level_2)
                    try:
                        for level_3 in level_2["children"]:
                            # print("sg3", subgroup_3["name"], subgroup_2["name"], subgroup_1["name"], children["name"])
                            # final print
                            ko = "ko:" + level_3["name"].partition(" ")[0]
                            if ko in kegg_dict:
                                # if kegg_dict[ko]['level_1'].find(root['name']) == -1:
                                kegg_dict[ko]["level_1"] += ";" + root["name"]
                                # kegg_dict[ko]['level_2']+='/'+level_1['name']

                                # if kegg_dict[ko]['level_2'].find(level_1['name']) == -1:
                                kegg_dict[ko]["level_2"] += ";" + level_1["name"]
                                # if kegg_dict[ko]['level_3'].find(level_2['name']) == -1:
                                kegg_dict[ko]["level_3"] += ";" + level_2["name"]

                            else:

                                # func=level_3['name'].partition(' ')[2]
                                kegg_dict[ko] = {
                                    "level_1": root["name"],
                                    "level_2": level_1["name"],
                                    "level_3": level_2["name"],
                                }

                    except:

                        ko = "ko:K" + level_2["name"].partition(" ")[0]
                        if ko in kegg_dict:
                            # if ko in kegg_dict:
                            # if kegg_dict[ko]['level_1'].find(root['name']) == -1:
                            kegg_dict[ko]["level_1"] += ";" + root["name"]
                            # kegg_dict[ko]['level_2']+='/'+level_1['name']

                            # if kegg_dict[ko]['level_2'].find(level_1['name']) == -1:
                            kegg_dict[ko]["level_2"] += ";" + level_1["name"]
                            # if kegg_dict[ko]['level_3'].find(level_2['name']) == -1:
                            kegg_dict[ko]["level_3"] += ";" + level_2["name"]
                        else:
                            # print(ko,level_2['name'])
                            kegg_dict[ko] = {
                                "level_1": root["name"],
                                "level_2": level_1["name"],
                                "level_3": level_2["name"],
                            }
        kegg_pd = pd.DataFrame.from_dict(kegg_dict, orient="index").reset_index()
        kegg_pd = kegg_pd.rename(columns={"index": "KO"})
        return kegg_pd


# @click.option('--taxa',default='Phylum',help='level of taxonomy to show in final table default: Phylum')
# @click.option('--min_perc_id',default=40,help='minimum cut off of aa percent id to retain read default: 40' )
# @click.option('--min_len',default=50, help='minimum aa alignment length to retain read default:50')
# @click.argument('file_path')
# @click.option('--function',default='ko',help='options include: ec, ko or product default: ko')


# files=open('/Users/zavo603/Documents/Nicki_files/perseq/tables/t140m_classifications.txt')
# file_original='/Users/zavo603/Documents/Nicki_files/perseq/tables/t140m_classifications.txt'
def main(
    json_path,
    classification_file_path,
    tax_level,
    min_perc_id,
    min_len,
    group_on,
    output,
):

    # function='ko'
    # tax_level='Phylum'
    # perc_id=50
    # align_len=40
    # json_path='/Users/zavo603/Documents/Nicki_files/perseq/ko00001.json'
    # path='/Users/zavo603/Documents/Nicki_files/perseq/tables/'
    if type(group_on) == str:
        group_on = [group_on]
    group_with = list(group_on)
    group_on = list(group_on)
    group_on.extend(["aa_alignment_length", "aa_percent_id"])
    tax_levels = {
        "superkingdom": 1,
        "kingdom": 1,
        "phylum": 2,
        "class": 3,
        "order": 4,
        "family": 5,
        "genus": 6,
        "species": 7,
    }
    grouped_sample_tbl = None
    samples = []
    ## THIS WILL NEED TO BE REWORKED TO ACCOMODATE SNAKEMAKE
    kegg_pd = parse_kegg_json(json_path)
    for f in os.listdir(classification_file_path):
        f = os.path.join(classification_file_path, f)
        sample = os.path.basename(f).partition("_classifications.txt")[0]
        samples.append(sample)
        with open(f) as file:
            print("working on", f)
            tax_ko = {"tax": []}
            next(file)
            for line in file:
                toks = line.strip().split("\t")

                try:
                    tex = toks[7].split(";")
                    tax = ";".join(tex[: tax_levels[tax_level]])
                    tax_ko["tax"].append(tax)
                except IndexError:
                    tax_ko["tax"].append("NA")
            # print(group_on)
            # print(group_with)
            if grouped_sample_tbl is None:
                tax_ko_tbl = pd.DataFrame(tax_ko)
                full_tbl = pd.read_table(
                    f, usecols=group_on
                )  # grab only the two columns of interest(can be changed later)

                full_tbl[sample] = "NA"  # add an empty column
                if "tax_classification" in group_with:
                    full_tbl["tax_classification"] = tax_ko_tbl[
                        "tax"
                    ]  ##this will assign the tax classification to the edited taxonomy from above
                else:
                    pass
                final_tbl = full_tbl.fillna("NA")

                final_tbl = final_tbl[final_tbl["aa_percent_id"] > min_perc_id]
                final_tbl = final_tbl[final_tbl["aa_alignment_length"] > min_len]
                final_tbl = final_tbl.drop(
                    ["aa_percent_id", "aa_alignment_length"], axis=1
                )
                # print(group_with)
                grouped_sample_tbl = final_tbl.groupby(group_with).count().reset_index()
                continue

            # print(sample)
            tax_ko_tbl = pd.DataFrame(tax_ko)
            full_tbl = pd.read_table(
                f, usecols=group_on
            )  # grab only the two columns of interest(can be changed later)
            full_tbl[sample] = "NA"  # add an empty column
            if "tax_classification" in group_with:
                full_tbl["tax_classification"] = tax_ko_tbl[
                    "tax"
                ]  ##this will assign the tax classification to the edited taxonomy from above
            else:
                pass
            final_tbl = full_tbl.fillna("NA")
            final_tbl = final_tbl[final_tbl["aa_percent_id"] > min_perc_id]
            final_tbl = final_tbl[final_tbl["aa_alignment_length"] > min_len]
            final_tbl = final_tbl.drop(["aa_percent_id", "aa_alignment_length"], axis=1)
            grouped_tbl = final_tbl.groupby(group_with).count().reset_index()
            grouped_sample_tbl = grouped_sample_tbl.merge(
                grouped_tbl, on=group_with, how="outer"
            ).fillna("0")
    grouped_sample_tbl = grouped_sample_tbl.rename(
        columns={"tax_classification": "taxonomy_%s" % tax_level, "ko": "KO"}
    )

    try:
        # print(grouped_sample_tbl.shape)
        print(kegg_pd.shape)
        grouped_sample_tbl = grouped_sample_tbl.merge(kegg_pd, on="KO")
    except:
        # print('didnt work')
        pass

    grouped_sample_tbl.to_csv(output, sep="\t", index=False)


# def main():
#      test = build_func_tax_tbl(json_path,tax_level,min_perc_id,min_len,file_path,group_with,group_on)
#      test.to_csv(out,sep='\t')
#      print('done')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generates a tab seperated file of the desired functional information"
    )
    parser.add_argument("-jfp", "--json-file-path", help="path to json file of KOs")
    parser.add_argument(
        "-cfp",
        "--classification-file-path",
        help="path to the _classification.txt files",
    )
    parser.add_argument(
        "-t",
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
        help="taxa level to show default: phylum",
    )
    parser.add_argument(
        "-mp",
        "--min-perc-id",
        type=int,
        default=50,
        help="lowest percent ID to retain default: 50",
    )
    parser.add_argument(
        "-ml",
        "--min-len",
        type=int,
        default=40,
        help="lowest alignment length to retain default: 40",
    )
    # parser.add_argument('--function',
    #                    help='functional aspect of species to reatin(ko,ec,or product) default: ko')
    parser.add_argument(
        "-g",
        "--group-on",
        choices=["tax_classification", "ko", "ec", "product"],
        nargs="+",
        default="tax_classification",
        help="Information to be retained(taxonomy, function, etc.). Must match the file headers default: tax_classifcation",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="tax_ko_out.txt",
        type=str,
        help="Output table name default: tax_ko_out.txt",
    )
    args = parser.parse_args()
    main(
        args.json_file_path,
        args.classification_file_path,
        args.tax_level,
        args.min_perc_id,
        args.min_len,
        args.group_on,
        args.output,
    )
