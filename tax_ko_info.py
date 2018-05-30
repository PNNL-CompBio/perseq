##this works
# user provides the given taxa(class, order, etc..),aa_percent_id filter, aa_alignment_length filter,function to retain(if any), and path to files to be analyzed
import os
import pandas as pd
import json

# test_list["children"]
# def recursive_func(dictionary, min_depth=4):
def parse_kegg_json(json_file_path):
    with open(json_file_path) as json_file:
        json_class = json.load(json_file)
        kegg_dict = dict()
        for root in json_class["children"]:

            for level_1 in root["children"]:

                for level_2 in level_1["children"]:

                    try:
                        for level_3 in level_2["children"]:
                            # print("sg3", subgroup_3["name"], subgroup_2["name"], subgroup_1["name"], children["name"])
                            # final print
                            ko = "ko:" + level_3["name"].partition(" ")[0]
                            func = level_3["name"].partition(" ")[2]
                            kegg_dict[ko] = {
                                "level_1": level_1["name"],
                                "level_2": level_2["name"],
                                "level_3": func,
                            }
                    except KeyError:
                        ko = "ko:KO" + level_2["name"].partition(" ")[0]
                        kegg_dict[ko] = {
                            "level_1": level_1["name"],
                            "level_2": level_2["name"],
                            "level_3": "NA",
                        }
        kegg_pd = pd.DataFrame.from_dict(kegg_dict, orient="index").reset_index()
        kegg_pd = kegg_pd.rename(columns={"index": "KO"})
        return kegg_pd


# files=open('/Users/zavo603/Documents/Nicki_files/perseq/tables/t140m_classifications.txt')
# file_original='/Users/zavo603/Documents/Nicki_files/perseq/tables/t140m_classifications.txt'
def build_func_tax_tbl(json_path, function, tax_level, min_perc_id, min_len, file_path):

    # function='ko'
    # tax_level='Phylum'
    # perc_id=50
    # align_len=40
    # json_path='/Users/zavo603/Documents/Nicki_files/perseq/ko00001.json'
    # path='/Users/zavo603/Documents/Nicki_files/perseq/tables/'
    tax_levels = {
        "Super Kingdom": 1,
        "Kingdom": 1,
        "Phylum": 2,
        "Class": 3,
        "Order": 4,
        "Family": 5,
        "Genus": 6,
        "Species": 7,
    }
    grouped_sample_tbl = None
    samples = []
    for f in os.listdir(path):
        f = os.path.join(path, f)
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
            if grouped_sample_tbl is None:
                tax_ko_tbl = pd.DataFrame(tax_ko)
                full_tbl = pd.read_table(
                    f,
                    usecols=[
                        function,
                        "tax_classification",
                        "aa_alignment_length",
                        "aa_percent_id",
                    ],
                )  # grab only the two columns of interest(can be changed later)

                full_tbl["Count_%s" % sample] = "NA"  # add an empty column
                full_tbl["tax_classification"] = tax_ko_tbl[
                    "tax"
                ]  ##this will assign the tax classification to the edited taxonomy from above
                final_tbl = full_tbl.fillna("NA")

                final_tbl = final_tbl[final_tbl["aa_percent_id"] > min_perc_id]
                final_tbl = final_tbl[final_tbl["aa_alignment_length"] > min_len]
                final_tbl = final_tbl.drop(
                    ["aa_percent_id", "aa_alignment_length"], axis=1
                )
                grouped_sample_tbl = (
                    final_tbl.groupby(["tax_classification", function])
                    .count()
                    .reset_index()
                )
                continue

            print(sample)
            tax_ko_tbl = pd.DataFrame(tax_ko)
            full_tbl = pd.read_table(
                f,
                usecols=[
                    function,
                    "tax_classification",
                    "aa_percent_id",
                    "aa_alignment_length",
                ],
            )  # grab only the two columns of interest(can be changed later)
            full_tbl["Count_%s" % sample] = "NA"  # add an empty column
            full_tbl["tax_classification"] = tax_ko_tbl[
                "tax"
            ]  ##this will assign the tax classification to the edited taxonomy from above
            final_tbl = full_tbl.fillna("NA")
            final_tbl = final_tbl[final_tbl["aa_percent_id"] > min_perc_id]
            final_tbl = final_tbl[final_tbl["aa_alignment_length"] > min_len]
            final_tbl = final_tbl.drop(["aa_percent_id", "aa_alignment_length"], axis=1)
            grouped_tbl = (
                final_tbl.groupby(["tax_classification", function])
                .count()
                .reset_index()
            )
            grouped_sample_tbl = grouped_sample_tbl.merge(
                grouped_tbl, on=["tax_classification", "ko"], how="outer"
            ).fillna("0")
    grouped_sample_tbl = grouped_sample_tbl.rename(
        columns={
            "tax_classification": "Tax. Classification(%s)" % tax_level,
            "ko": "KO",
        }
    )
    kegg_pd = parse_kegg_json(json_path)
    final_tbl = grouped_sample_tbl.merge(kegg_pd, on="KO")
    return final_tbl
