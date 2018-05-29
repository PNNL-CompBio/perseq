import argparse
import os
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go
from plotly.offline import iplot, offline
from snakemake.utils import logger, report

PLOTLY_PARAMS = dict(
    include_plotlyjs=False, show_link=False, output_type="div", image_height=700
)
STYLE = """
    <style type="text/css">
    body{font-family:Helvetica,arial,sans-serif;font-size:14px;line-height:1.6;padding-bottom:10px;background-color:#fff;color:#333;margin:0}body>div .section::before{content:"";display:block;height:80px;margin:-80px 0 0}#summary::before{margin:0}.topic-title{font-size:18pt}body>div>.section{margin-left:22%;margin-bottom:3em}div.section{margin-right:20px}#contents>p{display:none}button,li p.first{display:inline-block}#contents{margin-top:80px;padding-left:0;width:20%;background-color:#f1f1f1;height:100%;position:fixed;overflow:auto}#contents ul{list-style-type:none}#contents ul>li{font-size:14pt}#contents ul>li a:hover{color:#151d26}button,h1.title{color:#fff;background-color:#151d26}#contents ul>li>ul>li{font-size:12pt}h1.title{margin-top:0;position:fixed;z-index:10;padding:20px;width:100%}code,table tr:nth-child(2n),tt{background-color:#f8f8f8}.one-col{min-width:310px;height:500px;margin:0 auto}.two-col-left{height:300px;width:49%;float:left}.two-col-right{height:300px;width:49%;float:right}button{margin:0 5px 0 0;padding:5px 25px;font-size:18px;line-height:1.8;appearance:none;box-shadow:none;border-radius:3px;border:none}button:focus{outline:0}button:hover{background-color:#4183C4}button:active{background-color:#27496d}.legend-rect{width:20px;height:20px;margin-right:8px;margin-left:20px;float:left;-webkit-border-radius:2px;border-radius:2px}a{color:#4183C4;text-decoration:none}a.absent{color:#c00}a.anchor{padding-left:30px;margin-left:-30px;cursor:pointer;position:absolute;top:0;left:0;bottom:0}dl,dl dt,dl dt:first-child,hr,table,table tr{padding:0}table tr td,table tr th{border:1px solid #ccc;text-align:left;padding:6px 13px}h1,h2,h3,h4,h5,h6{margin:20px 0 10px;padding:0;font-weight:700;-webkit-font-smoothing:antialiased;cursor:text;position:relative}h1:hover a.anchor,h2:hover a.anchor,h3:hover a.anchor,h4:hover a.anchor,h5:hover a.anchor,h6:hover a.anchor{text-decoration:none}h1 code,h1 tt,h2 code,h2 tt,h3 code,h3 tt,h4 code,h4 tt,h5 code,h5 tt,h6 code,h6 tt{font-size:inherit}h1{font-size:28px;color:#151d26;border-bottom:1px solid #ccc}h2{font-size:24px;color:#000}h3{font-size:18px}h4{font-size:16px}dl dt,h5,h6{font-size:14px}h6{color:#777}blockquote,dl,li,ol,p,pre,table,ul{margin:15px 0}hr{background:url(http://tinyurl.com/bq5kskr) repeat-x;border:0;color:#ccc;height:4px}a:first-child h1,a:first-child h2,a:first-child h3,a:first-child h4,a:first-child h5,a:first-child h6{margin-top:0;padding-top:0}h1 p,h2 p,h3 p,h4 p,h5 p,h6 p{margin-top:0}dl dt{font-weight:700;font-style:italic;margin:15px 0 5px}blockquote>:first-child,dl dd>:first-child,dl dt>:first-child,table tr td :first-child,table tr th :first-child{margin-top:0}blockquote>:last-child,dl dd>:last-child,dl dt>:last-child{margin-bottom:0}dl dd{margin:0 0 15px;padding:0 15px}blockquote{border-left:4px solid #ddd;padding:0 15px;color:#777}table{border-spacing:0;border-collapse:collapse}table tr{border-top:1px solid #ccc;background-color:#fff;margin:0}table tr th{font-weight:700;margin:0}table tr td{margin:0}table tr td :last-child,table tr th :last-child{margin-bottom:0}img{max-width:100%}span.frame{display:block;overflow:hidden}span.frame>span{border:1px solid #ddd;display:block;float:left;overflow:hidden;margin:13px 0 0;padding:7px;width:auto}span.frame span img{display:block;float:left}span.frame span span{clear:both;color:#333;display:block;padding:5px 0 0}span.align-center{display:block;overflow:hidden;clear:both}span.align-center>span{display:block;overflow:hidden;margin:13px auto 0;text-align:center}span.align-center span img{margin:0 auto;text-align:center}span.align-right{display:block;overflow:hidden;clear:both}span.align-right>span{display:block;overflow:hidden;margin:13px 0 0;text-align:right}span.align-right span img{margin:0;text-align:right}span.float-left{display:block;margin-right:13px;overflow:hidden;float:left}span.float-left span{margin:13px 0 0}span.float-right{display:block;margin-left:13px;overflow:hidden;float:right}span.float-right>span{display:block;overflow:hidden;margin:13px auto 0;text-align:right}code,tt{margin:0 2px;padding:0 5px;white-space:nowrap;border:1px solid #eaeaea;border-radius:3px}pre code{margin:0;padding:0;white-space:pre;background:0 0}.highlight pre,pre{background-color:#f8f8f8;border:1px solid #ccc;font-size:13px;line-height:19px;overflow:auto;padding:6px 10px;border-radius:3px}pre code,pre tt{background-color:transparent;border:none}div#metadata{text-align:right}
    </style>
"""


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
    x = x[~ np.isnan(x)]
    n = len(x)
    assert n > 0, "x is empty"
    x.sort(kind="mergesort")
    G = sum(np.arange(1, n + 1) * x)
    G = 2 * G / sum(x) - (n + 1)
    if corr:
        return G / (n - 1)

    else:
        return G / n


def parse_classifications_for_taxonomy(path):
    """
    SN1035:381:h3cv7bcx2:2:2202:3165:73750	-1	-1
    SN1035:381:h3cv7bcx2:1:1105:7280:88523	67.5	40				111	Archaea; Thaumarchaeota; Nitrososphaerales; Nitrososphaeria; Nitrososphaeraceae; Candidatus Nitrosocosmicus; Candidatus Nitrocosmicus oleophilus;
    SN1035:381:h3cv7bcx2:1:1204:20140:83036	78.4	37	ko:K00820	glmS, GFPT; glucosamine---fructose-6-phosphate aminotransferase (isomerizing)	2.6.1.16	229	Bacteria; Actinobacteria; NA; NA; NA; NA; Actinobacteria bacterium 13_1_40CM_66_12;
    """
    logger.info("Parsing {}".format(path))
    # hardcoded tax levels :(
    taxonomy_counter = {"order": Counter(), "class": Counter(), "phylum": Counter()}
    summary_counter = Counter()
    with open(path) as fh:
        # skip the header
        next(fh)
        for i, line in enumerate(fh, start=1):
            toks = line.strip("\r\n").split("\t")
            if toks[3]:
                summary_counter.update(["Assigned Function"])
                if toks[7]:
                    summary_counter.update(["Assigned Both"])
            if toks[7]:
                taxonomy = [j.strip() for j in toks[7].split(";")]
                taxonomy_counter["order"].update([taxonomy[3]])
                taxonomy_counter["class"].update([taxonomy[2]])
                taxonomy_counter["phylum"].update([taxonomy[1]])
                summary_counter.update(["Assigned Taxonomy"])
    return taxonomy_counter, i, summary_counter


def get_df_at_tax_level(count_obj, sample):
    df = pd.DataFrame(data=count_obj, index=[0]).transpose()
    df.reset_index(inplace=True)
    return df


def process_reads(tax_level, pandas_df, total_counts):
    """
    This function will take the merged reads from the file and put them into a graphable
    form for both counts and percents. It will return both.
    """
    header = pandas_df.columns.tolist()
    header.pop(0)
    # organize these by count
    full = pandas_df.sort_values(by=header, ascending=False)
    sample_order = []
    # this will make a list of the gini coefficient of where these should be
    # located by diversity of the sample
    for read_col in pandas_df.columns.tolist():
        if not read_col == tax_level:
            sample_order.append([read_col, gini(pandas_df[read_col])])
    sample_order = sorted(sample_order, key= lambda x: x[1])
    sub = full[[tax_level] + [sample_order[i][0] for i in range(len(sample_order))]]
    cols = sub[tax_level].tolist()
    # find the percent
    sub_perc = sub[header].div(total_counts)
    sub_perc_t = sub_perc.transpose()
    # reset the names
    sub_perc_t.columns = cols
    # these are just the counts
    # needs to be in this form or it wont work
    sub_t = sub.transpose()
    # set the column names as the taxonomy levels
    sub_t.columns = sub_t.loc[tax_level]
    # this just takes out the duplicate header
    sub_t.drop([tax_level], inplace=True)
    return sub_t, sub_perc_t


def compile_summary_df(classification_tables, tax_levels= ["phylum", "class", "order"]):
    """
    Reads in multiple sample alignments from diamond in a given directory and merges them into
    a single pandas.DataFrame. It returns a pandas dataframe for each of th
    phylum that is ready to plug into the processing function. Also returns total counts
    which is necessary to calculate the percentage of total that is being represented
    """
    samples = []
    total_counts = []
    dfs = {}
    classifications_per_sample = {}
    # summary_counts = Counter()
    for classification_table in classification_tables:
        sample = get_sample(classification_table, "_classifications.txt")
        samples.append(sample)
        counts, observed_seqs, assignment_summary = parse_classifications_for_taxonomy(classification_table)
        # assigned #'s in summary table
        classifications_per_sample[sample] = assignment_summary
        total_counts.append(observed_seqs)
        if len(dfs) == 0:
            for tax_level in tax_levels:
                dfs[tax_level] = get_df_at_tax_level(counts[tax_level], sample)
                dfs[tax_level].columns = [tax_level, sample]
            continue

        for tax_level in tax_levels:
            df = get_df_at_tax_level(counts[tax_level], sample)
            df.columns = [tax_level, sample]
            dfs[tax_level] = dfs[tax_level].merge(df, on=tax_level, how="outer")
    observations_at_levels = {"Counts": dict(), "Percentage": dict()}
    for tax_level in tax_levels:
        c, p = process_reads(tax_level, dfs[tax_level], total_counts)
        observations_at_levels["Counts"][tax_level] = c
        observations_at_levels["Percentage"][tax_level] = p
    return observations_at_levels, pd.DataFrame.from_dict(classifications_per_sample, orient="index")


def make_plots(observations, summary_type):
    # data traces are taxonomies across samples
    labels = {"Percentage": "Percent of Total Reads", "Counts": "Count of Alignments"}
    # tax levels are hardcoded at this point
    data = [
        go.Bar(
            x=observations[summary_type]["phylum"].index,
            y=observations[summary_type]["phylum"][tax],
            name=tax,
            visible=True,
        )
        for tax in observations[summary_type]["phylum"].columns.tolist()
    ] + [
        go.Bar(
            x=observations[summary_type]["class"].index,
            y=observations[summary_type]["class"][tax],
            name=tax,
            visible=False,
        )
        for tax in observations[summary_type]["class"].columns.tolist()
    ] + [
        go.Bar(
            x=observations[summary_type]["order"].index,
            y=observations[summary_type]["order"][tax],
            name=tax,
            visible=False,
        )
        for tax in observations[summary_type]["order"].columns.tolist()
    ]
    # the number of taxa
    trace_length_phy = len(observations[summary_type]["phylum"].columns)
    trace_length_cla = len(observations[summary_type]["class"].columns)
    trace_length_ord = len(observations[summary_type]["order"].columns)
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
                                    "visible": [True] *
                                    trace_length_phy +
                                    [False] *
                                    trace_length_cla +
                                    [False] *
                                    trace_length_ord
                                },
                                {"yaxis": {"title": labels[summary_type]}},
                            ],
                        ),
                        dict(
                            label="Class",
                            method="update",
                            args=[
                                {
                                    "visible": [False] *
                                    trace_length_phy +
                                    [True] *
                                    trace_length_cla +
                                    [False] *
                                    trace_length_ord
                                },
                                {"yaxis": {"title": labels[summary_type]}},
                            ],
                        ),
                        dict(
                            label="Order",
                            method="update",
                            args=[
                                {
                                    "visible": [False] *
                                    trace_length_phy +
                                    [False] *
                                    trace_length_cla +
                                    [True] *
                                    trace_length_ord
                                },
                                {"yaxis": {"title": labels[summary_type]}},
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
        title="Assignments per Sample By {}".format(summary_type),
        updatemenus=updatemenus,
        barmode="stack",
        margin={"b": "auto", "r": "auto"},
        yaxis=dict(title=labels[summary_type]),
        showlegend=False,
        hovermode="closest",
    )
    fig = go.Figure(data=data, layout=layout)
    return fig


def get_sample(path, key):
    return os.path.basename(path).partition(key)[0]


def parse_merge_file(path):
    counts = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("Input:"):
                toks = line.strip().split("\t")
                counts.append(int(toks[1].partition(" ")[0]))
            if line.startswith("Joined:"):
                toks = line.strip().split("\t")
                # joined
                counts.append(int(toks[1]))
                # join rate
                counts.append(toks[2])
            if line.startswith("Avg Insert:"):
                counts.append(float(line.strip().split("\t")[1]))
    return counts


def parse_classifications_tables(tables):
    logger.info("Parsing summary tables")
    final = {}
    for table_file in tables:
        sample = os.path.basename(table_file).partition("_classifications.txt")[0]
        with open(table_file) as fh:
            counter = Counter()
            header = next(fh)
            for line in fh:
                toks = line.strip("\r\n").split("\t")
                if toks[3]:
                    counter.update(["Assigned Function"])
                    if toks[7]:
                        counter.update(["Assigned Both"])
                if toks[7]:
                    counter.update(["Assigned Taxonomy"])
        final[sample] = counter
    df = pd.DataFrame.from_dict(final, orient="index")
    return df


def parse_log_files(
    deduplication_logs, decontamination_logs, merge_logs, classifications_per_sample
):
    logger.info("Parsing log files for summary table")
    header = [
        "Initial Read Count",
        "Deduplication",
        "Decontamination",
        "Pairs Joined",
        "Join Rate",
        "Average Insert Size",
    ]
    count_table = defaultdict(list)
    # initial read count from first step's log
    for deduplication_log in deduplication_logs:
        sample = get_sample(deduplication_log, "_deduplicate_reads.log")
        with open(deduplication_log) as fh:
            for line in fh:
                if line.startswith("Reads In:"):
                    count_table[sample].append(int(line.strip("\r\n").split(" ")[-1]))
    # deduplication count from decontamination input
    for decontamination_log in decontamination_logs:
        sample = get_sample(decontamination_log, "_decontamination.log")
        with open(decontamination_log) as fh:
            for line in fh:
                if line.startswith("Reads Used:"):
                    count_table[sample].append(int(line.strip("\r\n").split("\t")[1]))
    # decontamination count, join count, join rate, insert size from merge step
    for merge_log in merge_logs:
        sample = get_sample(merge_log, "_merge_sequences.log")
        count_table[sample].extend(parse_merge_file(merge_log))
    log_df = pd.DataFrame.from_dict(count_table, orient="index")
    log_df.columns = header
    # parse the summary tables for assignments
    # tables_df = parse_classifications_tables(summary_tables)
    log_df = log_df.merge(classifications_per_sample, left_index=True, right_index=True)
    header.extend(["Assigned Function", "Assigned Taxonomy", "Assigned Both"])
    return log_df[header].to_html().replace("\n", "\n" + 10 * " ")


def main(
    decontamination_logs,
    deduplication_logs,
    merge_logs,
    summary_tables,
    html,
):
    observations_at_levels, classifications_per_sample = compile_summary_df(summary_tables)
    div = {}
    for v in ["Percentage", "Counts"]:
        div[v] = offline.plot(make_plots(observations_at_levels, v), **PLOTLY_PARAMS)
    html_tbl = parse_log_files(
        deduplication_logs, decontamination_logs, merge_logs, classifications_per_sample
    )
    report_str = """

.. raw:: html

    {STYLE}
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>


=============================================================
PerSeq_ - Per sequence functional and taxonomic assignments
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
    report(report_str, html, stylesheet="")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--decontamination-logs", nargs="+")
    p.add_argument("--deduplication-logs", nargs="+")
    p.add_argument("--merge-logs", nargs="+")
    p.add_argument("--summary-tables", nargs="+")
    p.add_argument("--html")
    args = p.parse_args()
    main(
        args.decontamination_logs,
        args.deduplication_logs,
        args.merge_logs,
        args.summary_tables,
        args.html,
    )
