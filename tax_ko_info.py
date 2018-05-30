import os
import pandas as pd
import json
import argparse

# @click.command()
# @click.argument('json_file')
parser = argparse.ArgumentParser()
parser.add_argument('--json_file_path',
                    help='path to json file of KOs')
parser.add_argument('--file_path',
                    help='path to the _classification.txt files')
parser.add_argument('--tax_level', default='Phylum',
                    help='taxa level to show')
parser.add_argument('--min_perc_id',type=int, default=50,
                    help='lowest percent ID to retain default: 50')
parser.add_argument('--min_len', type=int, default=40,
                    help='lowest alingment length to reatin default: 40')
#parser.add_argument('--function',
#                    help='functional aspect of species to reatin(ko,ec,or product) default: ko')
parser.add_argument('--group_on', nargs='+',default='tax_classification',help='Information to be retained(taxonomy, function, etc.). Must match the file headers default: tax_classifcation')
parser.add_argument('--output',default='tax_ko_out.txt',type=str,help='Output table name')

args = parser.parse_args()
json_path=args.json_file_path
file_path=args.file_path
tax_level=args.tax_level
min_perc_id=args.min_perc_id
min_len=args.min_len
out=args.output
group_on=args.group_on

group_with=list(group_on)

group_on.extend(["aa_alignment_length","aa_percent_id"])

def parse_kegg_json(json_file):
    with open(json_file) as json_file:
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
# @click.option('--taxa',default='Phylum',help='level of taxonomy to show in final table default: Phylum')
# @click.option('--min_perc_id',default=40,help='minimum cut off of aa percent id to retain read default: 40' )
# @click.option('--min_len',default=50, help='minimum aa alignment length to retain read default:50')
# @click.argument('file_path')
# @click.option('--function',default='ko',help='options include: ec, ko or product default: ko')



# files=open('/Users/zavo603/Documents/Nicki_files/perseq/tables/t140m_classifications.txt')
# file_original='/Users/zavo603/Documents/Nicki_files/perseq/tables/t140m_classifications.txt'
def build_func_tax_tbl(json_path,taxa,min_perc_id,min_len,file_path,group_with,group_on):

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
    ## THIS WILL NEED TO BE REWORKED TO ACCOMODATE SNAKEMAKE
    for f in os.listdir(file_path):
        f = os.path.join(file_path, f)
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
                    usecols=group_on,
                )  # grab only the two columns of interest(can be changed later)

                full_tbl["Count_%s" % sample] = "NA"  # add an empty column
                try:
                    full_tbl["tax_classification"] = tax_ko_tbl[
                    "tax"
                    ]  ##this will assign the tax classification to the edited taxonomy from above
                except:
                    pass
                final_tbl = full_tbl.fillna("NA")

                final_tbl = final_tbl[final_tbl["aa_percent_id"] > min_perc_id]
                final_tbl = final_tbl[final_tbl["aa_alignment_length"] > min_len]
                final_tbl = final_tbl.drop(
                    ["aa_percent_id", "aa_alignment_length"], axis=1
                )
                print(group_with)
                grouped_sample_tbl = (
                    final_tbl.groupby(group_with)
                    .count()
                    .reset_index()
                )
                continue

            print(sample)
            tax_ko_tbl = pd.DataFrame(tax_ko)
            full_tbl = pd.read_table(
                f,
                usecols=group_on,
            )  # grab only the two columns of interest(can be changed later)
            full_tbl["Count_%s" % sample] = "NA"  # add an empty column
            try:
                full_tbl["tax_classification"] = tax_ko_tbl[
                "tax"
            ]  ##this will assign the tax classification to the edited taxonomy from above
            except:
                pass
            final_tbl = full_tbl.fillna("NA")
            final_tbl = final_tbl[final_tbl["aa_percent_id"] > min_perc_id]
            final_tbl = final_tbl[final_tbl["aa_alignment_length"] > min_len]
            final_tbl = final_tbl.drop(["aa_percent_id", "aa_alignment_length"], axis=1)
            grouped_tbl = (
                final_tbl.groupby(group_with)
                .count()
                .reset_index()
            )
            grouped_sample_tbl = grouped_sample_tbl.merge(
                grouped_tbl, on=group_with, how="outer"
            ).fillna("0")
    grouped_sample_tbl = grouped_sample_tbl.rename(
        columns={
            "tax_classification": "Tax. Classification(%s)" % tax_level,
            "ko": "KO",
        }
    )
    kegg_pd = parse_kegg_json(json_path)
    try:
        grouped_sample_tbl = grouped_sample_tbl.merge(kegg_pd, on="KO")
    except:
        pass

    return grouped_sample_tbl

def main():
     test = build_func_tax_tbl(json_path,tax_level,min_perc_id,min_len,file_path,group_with,group_on)
     test.to_csv(out,sep='\t')
     print('done')


if __name__ == '__main__':
    main()
