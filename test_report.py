import argparse
import os
import pandas as pd
import plotly.graph_objs as go
import plotly
from plotly.offline import offline, iplot

from snakemake.utils import report

PLOTLY_PARAMS = dict(
    include_plotlyjs=False, show_link=False, output_type="div", image_height=700
)


def parse_log_files(path, expect_one_value=True):

    key_val = {
        "merge_sequences.log": "Input:",
        "decontamination.log": "Reads Used:",
        "deduplicate_reads.log": "Reads Processed:",
    }
    df_table = {}
    for f in os.listdir(path):
        f = os.path.join(path, f)

        # print(f)
        if f.endswith("_decontamination_by_reference.log"):
            continue
        if f.endswith("_deduplicate_reads.log"):
            continue
        sample = os.path.basename(f).partition("_")[0]
        ### NEED TO FIND A WAY TO DO THIS THAT IS NOT HARDCODED
        if not sample in df_table:
            df_table[sample] = [0, 1, 2, 3, 4]
            if sample == "t140m":
                df_table[sample][0] = "2707038"
            elif sample == "t164m":
                df_table[sample][0] = "8276745"
            elif sample == "t174m":
                df_table[sample][0] = "4238294"
        else:
            pass
        file_type = os.path.basename(f).partition("_")[2]
        keyword = key_val[file_type]

        content = open(f).read()
        pos = content.find(keyword)
        if pos == -1:
            # raise Exception("Didn't find {} in file:\n\n{}".format(keyword, path))
            continue
        # else:
        elif keyword == "Reads In:":  # the deduplication file
            if expect_one_value:
                out = content[pos:].split()[2][:-1] + "000"
                df_table[sample][1] = out

                # out=out +'000'
                # print(out)

            else:
                out = content[pos:].split()[2:][:-1] + "000"
                df_table[sample][1] = out

        elif keyword == "Input:":  # the merged file
            if expect_one_value:
                out = content[pos:].split()[1]
                df_table[sample][2] = out
                keyword = "Joined:"

                pos = content.find(keyword)
                out = content[pos:].split()[1]
                df_table[sample][3] = out
                out = content[pos:].split()[2]
                df_table[sample][4] = out

            else:
                out = content[pos:].split()[1:]
                df_table[sample][2] = out
                keyword = "Joined:"

                pos = content.find(keyword)
                out = content[pos:].split()[1]
                df_table[sample][3] = out
                out = content[pos:].split()[2]
                df_table[sample][4] = out
        elif keyword == "Reads Used:":  # the deduplication file(ignore?)
            if expect_one_value:
                out = content[pos:].split()[2]
                df_table[sample][1] = out

            else:
                out = content[pos:].split()[2:]
                df_table[sample][1] = out
        # print(df_table)
    final = pd.DataFrame.from_dict(df_table, orient="index")
    final.columns = [
        "Initial \nRead Count",
        "Decontamination",
        "Deduplication",
        "Merged",
        "Pairs Joined",
    ]
    html_tbl = final.to_html().replace("\n", "\n" + 10 * " ")
    return html_tbl


# def parse_log_files(path, expect_one_value=True):
#
#     key_val={'merge_sequences.log': 'Input:','decontamination.log':'Reads Used:', 'deduplicate_reads.log':'Reads Processed:'}
#     df_table={}
#     for f in os.listdir(path):
#         f = os.path.join(path, f)
#
#         # print(f)
#         if f.endswith("_decontamination_by_reference.log") or (): continue
#         sample = os.path.basename(f).partition("_")[0]
#         if not sample in df_table:
#             df_table[sample]=[0,1,2,3]
#             if sample == 't140m':
#                 df_table[sample][0]='2707038'
#             if sample == 't164m':
#                 df_table[sample][0]='8276745'
#             else:
#                 df_table[sample][0]='4238294'
#         else:
#             pass
#         file_type=os.path.basename(f).partition("_")[2]
#         keyword=key_val[file_type]
#         #print(keyword)
#         content = open(f).read()
#         pos = content.find(keyword)
#         if pos == -1:
#             #raise Exception("Didn't find {} in file:\n\n{}".format(keyword, path))
#             continue
#         #else:
#         elif keyword == 'Reads Processed:':
#             if expect_one_value:
#                 out=content[pos:].split()[2][:-1]+'000'
#                 df_table[sample][2]=out
#                 #out=out +'000'
#                 #print(out)
#
#             else:
#                 out=content[pos:].split()[2:][:-1]+'000'
#                 df_table[sample][2]=out
#         elif keyword == 'Input:':
#             if expect_one_value:
#                 out=content[pos:].split()[1]
#                 df_table[sample][3]=out
#
#             else:
#                 out=content[pos:].split()[1:]
#                 df_table[sample][3]=out
#         elif keyword =='Reads Used:':
#             if expect_one_value:
#                 out=content[pos:].split()[2]
#                 df_table[sample][1]=out
#
#
#             else:
#                 out=content[pos:].split()[2:]
#                 df_table[sample][1]=out
#         #print(df_table)
#     final=pd.DataFrame.from_dict(df_table, orient='index')
#     final.columns=['Decontamination','Dedup','Merged']
#     html_tbl=final.to_html().replace('\n', '\n' + 10 * ' ')
#
#     return html_tbl

# parses log_file for extra info

# def parse_log_file(log_file, keyword, expect_one_value=True):
#     content = open(log_file).read()
#     pos = content.find(keyword)
#     if pos == -1:
#         raise Exception("Didn't find {} in file:\n\n{}".format(keyword, log_file))

#     else:
#         if expect_one_value:
#             return content[pos:].split()[2]

#         else:
#             return content[pos:].split()[2:]

# dependent functions for the others
# function to parse kaiju output

# parse the kaiju output
def parse_align(alignment_file):
    total_lines = 0
    order_dict = {}
    class_dict = {}
    phylum_dict = {}
    total_align = []
    # for alignment_file in file_list:
    with open(alignment_file) as file:
        classified = 0
        unc = 0
        for line in file:

            total_lines += 1
            toks = line.strip().split("\t")
            # print(toks)
            classification = toks[0]

            if classification == "C":
                classified += 1
                bits = toks[7].split(";")
                order = bits[3].strip()
                clas = bits[2]
                phylum = bits[1]
                if order in order_dict:
                    order_dict[order] += 1
                else:
                    order_dict[order] = 1
                if clas in class_dict:
                    class_dict[clas] += 1
                else:
                    class_dict[clas] = 1
                if phylum in phylum_dict:
                    phylum_dict[phylum] += 1
                else:
                    phylum_dict[phylum] = 1
            else:
                unc += 1
    total = classified + unc
    # or should this be a list. prob list
    # print('family',family)
    return order_dict, class_dict, phylum_dict, total


# calculates the order for diversity
def gini(x, corr=False):
    """Calculates Gini coefficient, the measure of inequality among values
    of a frequency distribution. A Gini coefficient of zero expresses
    perfect equality.

    Port from ineq package in R::

            > library(ineq)
            > t <- c(1,2,6,7,8)
            > Gini(t)
            [1] 0.3166667
            > Gini(t, corr=TRUE)
            [1] 0.3958333

    Args:
        x (list): list or array of numbers
        corr (Optional[bool]): finite sample correction

    Returns:
        float

    >>> import numpy as np
    >>> t = [1,2,6,7,8]
    >>> gini(t) # doctest: +ELLIPSIS
    0.3166...
    >>> gini(t, corr=True) # doctest: +ELLIPSIS
    0.3958...
    >>> gini([]) # doctest: +ELLIPSIS
    Traceback (most recent call last):
     ...
    AssertionError: x is empty
    >>> t = [1,2,6,7,"A"]
    >>> gini(t) # doctest: +ELLIPSIS
    Traceback (most recent call last):
     ...
    ValueError: could not convert...
    """
    x = np.array(x, dtype=float)
    # filter out nan values as list is coming from merged dataframe
    x = x[~np.isnan(x)]
    n = len(x)
    assert n > 0, "x is empty"
    x.sort(kind="mergesort")
    G = sum(np.arange(1, n + 1) * x)
    G = 2 * G / sum(x) - (n + 1)
    if corr:
        return G / (n - 1)
    else:
        return G / n


# step 2
# so this will make the counts information that we need. but this should be put into a function because it needs to be done
# for all three of the taxa
import numpy as np


def process_reads(tax_level, pandas_df, total_counts):
    """This function will take the merged reads from the file and put them into a graphable
    form for both counts and percents. It will return both. """
    # cols = [i for i in m.columns.tolist() if i.startswith("Percent_")]
    header = pandas_df.columns.tolist()
    header.pop(0)
    full = pandas_df.sort_values(by=header, ascending=False)  # organize these by count

    sample_order = []
    for read_col in pandas_df.columns.tolist():  # this will make a list of the gini coefficient of where these should be located by diversity of the sample
        if read_col.startswith("Reads_"):
            sample_order.append([read_col.strip("Reads_"), gini(pandas_df[read_col])])
    sample_order = sorted(sample_order, key=lambda x: x[1])

    sub = full[
        [tax_level]
        + ["Reads_%s" % sample_order[i][0] for i in range(len(sample_order))]
    ]
    sub = sub.head(100)
    cols = sub[tax_level].tolist()

    # find the percent
    sub_perc = sub[header].div(total_counts)
    sub_perc_t = sub_perc.transpose()
    sub_perc_t.columns = cols  # reset the names

    # these are just the counts
    sub_t = sub.transpose()  # needs to be in this form or it wont work
    sub_t.columns = sub_t.loc[tax_level]  # set the column names as the taxoomy levels
    sub_t.drop([tax_level], inplace=True)  # this just takes out the duplicate header

    return sub_t, sub_perc_t


# step 1
def df_by_merge(path, tax_level):
    """
    Reads in multiple sample alignments from diamond in a given directory and merges them into
    a single pandas.DataFrame. It returns a pandas dataframe for each of th
    phylum that is ready to plug into the processing function. Also returns total counts
    which is necessary to calculate the percentage of total that is being represented
    """
    merged_phy = None
    samples = []
    total_counts = []
    for f in os.listdir(path):
        f = os.path.join(path, f)

        # print(f)
        # if not f.endswith("%s_summary.txt" % tax_level): continue
        sample = os.path.basename(f).partition("_aln_names.txt")[0]
        samples.append(sample)  # making a list of the sample headers
        header = [tax_level, "Reads_%s" % sample]
        if merged_phy is None:
            # print('yes')
            # merged = pd.read_table(f, header=0, names=header, comment="-")
            phy, ordr, clas, tc = parse_align(f)
            # print(phy)
            total_counts.append(tc)
            merged_phy = pd.DataFrame(data=phy, index=[0]).transpose()
            merged_phy.reset_index(inplace=True)
            merged_phy.columns = header
            merged_clas = pd.DataFrame(data=clas, index=[0]).transpose()
            merged_clas.reset_index(inplace=True)
            merged_clas.columns = header
            merged_ord = pd.DataFrame(data=ordr, index=[0]).transpose()
            merged_ord.reset_index(inplace=True)
            merged_ord.columns = header
            # cols = merged.columns.tolist() #make a list of the column headers from the merged file
            # cols.remove(tax_level) #remove the order tag
            # cols.insert(0, tax_level) #insert it at the front.
            # merged = merged[cols]
            continue
        phy, ordr, clas, tc = parse_align(f)
        total_counts.append(tc)
        df1 = pd.DataFrame(data=phy, index=[0]).transpose()
        df1.reset_index(inplace=True)
        df1.columns = header
        df2 = pd.DataFrame(data=ordr, index=[0]).transpose()
        df2.reset_index(inplace=True)
        df2.columns = header
        df3 = pd.DataFrame(data=clas, index=[0]).transpose()
        df3.reset_index(inplace=True)
        df3.columns = header
        merged_phy_final = merged_phy.merge(df1, on=tax_level, how="outer")
        merged_ord_final = merged_ord.merge(df2, on=tax_level, how="outer")
        merged_class_final = merged_clas.merge(df3, on=tax_level, how="outer")
        # print(total_counts)

    # return merged_phy_final,merged_class_final,merged_order_final
    p_done, p_perc_done = process_reads(tax_level, merged_phy_final, total_counts)
    c_done, c_perc_done = process_reads(tax_level, merged_class_final, total_counts)
    o_done, o_perc_done = process_reads(tax_level, merged_ord_final, total_counts)
    return p_done, p_perc_done, c_done, c_perc_done, o_done, o_perc_done


# def complete_processing(merged_phy_final,merged_class_final,merged_order_final):
#     p_done,p_perc_done=process_reads(tax_level,merged_phy_final)
#     c_done,c_perc_done=process_reads(tax_level,merged_class_final)
#     o_done,o_perc_done=process_reads(tax_level,merged_order_final)
#
#     return p_done,p_perc_done,c_done, c_perc_done,o_done,o_perc_done


def make_plots(p_done, c_done, o_done, y_axis_title, variable):
    ##code in this block works

    # data traces are taxonomies across samples
    data = [
        go.Bar(x=p_done.index, y=p_done[tax], name=tax, visible=True)
        for tax in p_done.columns.tolist()
    ] + [
        go.Bar(x=c_done.index, y=c_done[tax], name=tax, visible=False)
        for tax in c_done.columns.tolist()
    ] + [
        go.Bar(x=o_done.index, y=o_done[tax], name=tax, visible=False)
        for tax in o_done.columns.tolist()
    ]

    # the number of taxa
    trace_length = len(o_done.columns)

    # plot buttons
    updatemenus = list(
        [
            dict(
                type="buttons",
                active=0,
                buttons=list(
                    [
                        dict(
                            label="Phylum",
                            method="update",
                            args=[
                                {
                                    "visible": [True]
                                    * trace_length
                                    + [False]
                                    * trace_length
                                    + [False]
                                    * trace_length
                                },
                                {"yaxis": {"title": y_axis_title}},
                            ],
                        ),
                        dict(
                            label="Class",
                            method="update",
                            args=[
                                {
                                    "visible": [False]
                                    * trace_length
                                    + [True]
                                    * trace_length
                                    + [False]
                                    * trace_length
                                },
                                {"yaxis": {"title": y_axis_title}},
                            ],
                        ),
                        dict(
                            label="Order",
                            method="update",
                            args=[
                                {
                                    "visible": [False]
                                    * trace_length
                                    + [False]
                                    * trace_length
                                    + [True]
                                    * trace_length
                                },
                                {"yaxis": {"title": y_axis_title}},
                            ],
                        ),
                    ]
                ),
                direction="left",
                pad={"r": 0, "t": 0},
                showactive=True,
                x=0,
                xanchor="left",
                y=1.2,
                yanchor="top",
            )
        ]
    )

    # initial layout
    layout = dict(
        title="Assignments per Sample By " + variable + "(Top 20)",
        updatemenus=updatemenus,
        barmode="stack",
        margin={"b": 200},
        xaxis={"tickangle": -60},
        yaxis=dict(title=y_axis_title),
        showlegend=False,
        hovermode="closest",
    )

    fig = go.Figure(data=data, layout=layout)

    return fig


def main(path_to_files, report_out, path_to_log_files):
    tax_level = "taxon"
    p_done, p_perc_done, c_done, c_perc_done, o_done, o_perc_done = df_by_merge(
        path_to_files, tax_level
    )
    div = {}  # this this appears to be what the actual html is going to call on(because it contains the plot info)
    labels = {
        "Percentage": "Percent of Total Reads Aligned", "Counts": "Count of Alignments"
    }
    for variable in ["Percentage", "Counts"]:
        if variable == "Counts":
            y_axis_label = labels[variable]
            div[variable] = offline.plot(
                make_plots(p_done, c_done, o_done, y_axis_label, variable),
                include_plotlyjs=False,
                show_link=False,
                output_type="div",
                image_height=700,
            )
        else:
            y_axis_label = labels[variable]
            div[variable] = offline.plot(
                make_plots(
                    p_perc_done, c_perc_done, o_perc_done, y_axis_label, variable
                ),
                include_plotlyjs=False,
                show_link=False,
                output_type="div",
                image_height=700,
            )
    html_tbl = parse_log_files(path_to_log_files)
    report_str = """

.. raw:: html

    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>


=============================================================
PerSeq_ - Something
=============================================================

.. _PerSeq: https://github.com/pnnl


.. contents::
    :backlinks: none


Summary
-------

Counts by Taxonomy
******************

.. raw:: html

    {div[Counts]}


Percentage by Taxonomy
**********************

.. raw:: html

    {div[Percentage]}


Reads Surviving
****************

.. raw:: html

    {html_tbl}


Downloads
---------


"""

    report(
        report_str,
        report_out,
        stylesheet="/Users/zavo603/Documents/Nicki_files/perseq/perseq/resources/style.css",
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--path_to_alignment_files")
    p.add_argument("--path_to_log_files")
    p.add_argument("--report-out")

    args = p.parse_args()
    main(args.path_to_alignment_files, args.report_out, args.path_to_log_files)
