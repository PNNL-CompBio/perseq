import argparse
import csv
import os
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import plotly
import plotly.figure_factory as ff
import plotly.graph_objs as go
from plotly.offline import iplot, offline
from snakemake.utils import logger, report

PLOTLY_PARAMS = dict(
    include_plotlyjs=False, show_link=False, output_type="div", image_height=700
)

# data tables and then overall page style
STYLE = """
    <link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css"/>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.18/css/dataTables.bootstrap.min.css"/>
    <style type="text/css">
    body{font-family:Helvetica,arial,sans-serif;font-size:14px;line-height:1.6;padding-bottom:10px;background-color:#fff;color:#333;margin:0}body>div .section::before{content:"";display:block;height:80px;margin:-80px 0 0}#summary::before{margin:0}.topic-title{font-size:18pt}body>div>.section{margin-left:22%;margin-bottom:3em}div.section{margin-right:20px}#contents>p{display:none}button,li p.first{display:inline-block}#contents{margin-top:80px;padding-left:0;width:20%;background-color:#f1f1f1;height:100%;position:fixed;overflow:auto}#contents ul{list-style-type:none}#contents ul>li{font-size:14pt}#contents ul>li a:hover{color:#151d26}button,h1.title{color:#fff;background-color:#151d26}#contents ul>li>ul>li{font-size:12pt}h1.title{margin-top:0;position:fixed;z-index:10;padding:20px;width:100%}code,table tr:nth-child(2n),tt{background-color:#f8f8f8}.one-col{min-width:310px;height:500px;margin:0 auto}.two-col-left{height:300px;width:49%;float:left}.two-col-right{height:300px;width:49%;float:right}button{margin:0 5px 0 0;padding:5px 25px;font-size:18px;line-height:1.8;appearance:none;box-shadow:none;border-radius:3px;border:none}button:focus{outline:0}button:hover{background-color:#4183C4}button:active{background-color:#27496d}.legend-rect{width:20px;height:20px;margin-right:8px;margin-left:20px;float:left;-webkit-border-radius:2px;border-radius:2px}a{color:#4183C4;text-decoration:none}a.absent{color:#c00}a.anchor{padding-left:30px;margin-left:-30px;cursor:pointer;position:absolute;top:0;left:0;bottom:0}dl,dl dt,dl dt:first-child,hr,table,table tr{padding:0}table tr td,table tr th{border:1px solid #ccc;text-align:left;padding:6px 13px}h1,h2,h3,h4,h5,h6{margin:20px 0 10px;padding:0;font-weight:700;-webkit-font-smoothing:antialiased;cursor:text;position:relative}h1:hover a.anchor,h2:hover a.anchor,h3:hover a.anchor,h4:hover a.anchor,h5:hover a.anchor,h6:hover a.anchor{text-decoration:none}h1 code,h1 tt,h2 code,h2 tt,h3 code,h3 tt,h4 code,h4 tt,h5 code,h5 tt,h6 code,h6 tt{font-size:inherit}h1{font-size:28px;color:#151d26;border-bottom:1px solid #ccc}h2{font-size:24px;color:#000}h3{font-size:18px}h4{font-size:16px}dl dt,h5,h6{font-size:14px}h6{color:#777}blockquote,dl,li,ol,p,pre,table,ul{margin:15px 0}hr{background:url(http://tinyurl.com/bq5kskr) repeat-x;border:0;color:#ccc;height:4px}a:first-child h1,a:first-child h2,a:first-child h3,a:first-child h4,a:first-child h5,a:first-child h6{margin-top:0;padding-top:0}h1 p,h2 p,h3 p,h4 p,h5 p,h6 p{margin-top:0}dl dt{font-weight:700;font-style:italic;margin:15px 0 5px}blockquote>:first-child,dl dd>:first-child,dl dt>:first-child,table tr td :first-child,table tr th :first-child{margin-top:0}blockquote>:last-child,dl dd>:last-child,dl dt>:last-child{margin-bottom:0}dl dd{margin:0 0 15px;padding:0 15px}blockquote{border-left:4px solid #ddd;padding:0 15px;color:#777}table{border-spacing:0;border-collapse:collapse}table tr{border-top:1px solid #ccc;background-color:#fff;margin:0}table tr th{font-weight:700;margin:0}table tr td{margin:0}table tr td :last-child,table tr th :last-child{margin-bottom:0}img{max-width:100%}span.frame{display:block;overflow:hidden}span.frame>span{border:1px solid #ddd;display:block;float:left;overflow:hidden;margin:13px 0 0;padding:7px;width:auto}span.frame span img{display:block;float:left}span.frame span span{clear:both;color:#333;display:block;padding:5px 0 0}span.align-center{display:block;overflow:hidden;clear:both}span.align-center>span{display:block;overflow:hidden;margin:13px auto 0;text-align:center}span.align-center span img{margin:0 auto;text-align:center}span.align-right{display:block;overflow:hidden;clear:both}span.align-right>span{display:block;overflow:hidden;margin:13px 0 0;text-align:right}span.align-right span img{margin:0;text-align:right}span.float-left{display:block;margin-right:13px;overflow:hidden;float:left}span.float-left span{margin:13px 0 0}span.float-right{display:block;margin-left:13px;overflow:hidden;float:right}span.float-right>span{display:block;overflow:hidden;margin:13px auto 0;text-align:right}code,tt{margin:0 2px;padding:0 5px;white-space:nowrap;border:1px solid #eaeaea;border-radius:3px}pre code{margin:0;padding:0;white-space:pre;background:0 0}.highlight pre,pre{background-color:#f8f8f8;border:1px solid #ccc;font-size:13px;line-height:19px;overflow:auto;padding:6px 10px;border-radius:3px}pre code,pre tt{background-color:transparent;border:none}div#metadata{text-align:right}h1{line-height:1.6}.simple{padding-left:20px}.docutils.container{width:100%}
    </style>
"""
SCRIPT = """
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <script src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
    <script src="https://cdn.datatables.net/1.10.18/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/1.10.18/js/dataTables.bootstrap4.min.js"></script>
    <script>
    $(document).ready( function () {
        $('#summaryTable').DataTable();
    } );
    </script>
"""


def shannon_index(arr):
    """
    >>> counter = Counter({'28DPL6Y5ZLJQ': 20.0,
             '5JJ3I433OD4D': 3.0,
             '66ZS8L2FIBHE': 189.0,
             '7W0EEVTOR5VV': 1.0,
             '8NZCUXPEUIKE': 2.0,
             '8XTZF810ENMI': 486.0,
             'C8PP69BX09ZZ': 4.0,
             'CMULPOCAFROS': 2.0,
             'E3BOBZO5MIX3': 2.0,
             'E517B87AQ3XM': 2.0,
             'EUDVJH6XLSY3': 3.0,
             'F7O31BBBUDCE': 2.0,
             'JWWBB7Y2PTQ4': 2.0,
             'JYKOAK2U228V': 3.0,
             'NAMWMPHCUBJC': 2.0,
             'O5N7ADO5K9I4': 3.0,
             'P34M5UNS6P9B': 6.0,
             'RLU5W9AMJHFT': 123.0,
             'RMPWYKVRSXNL': 1.0,
             'T13Z885J5HJ4': 2.0,
             'TUK6CZH506C1': 2.0,
             'ZX7LUKOMKLA4': 3.0})
    >>> shannon_index(np.array(list(counter.values())))
    1.32095...
    """
    total = sum(arr)
    return -1 * sum([(x / total) * np.log(x / total) for x in arr if x > 0])


def parse_classifications_for_taxonomy(path):
    """
    SN1035:381:h3cv7bcx2:2:2202:3165:73750	-1	-1
    SN1035:381:h3cv7bcx2:1:1105:7280:88523	67.5	40				111	Archaea; Thaumarchaeota; Nitrososphaerales; Nitrososphaeria; Nitrososphaeraceae; Candidatus Nitrosocosmicus; Candidatus Nitrocosmicus oleophilus;
    SN1035:381:h3cv7bcx2:1:1204:20140:83036	78.4	37	ko:K00820	glmS, GFPT; glucosamine---fructose-6-phosphate aminotransferase (isomerizing)	2.6.1.16	229	Bacteria; Actinobacteria; NA; NA; NA; NA; Actinobacteria bacterium 13_1_40CM_66_12;
    """
    logger.info("Parsing {}".format(path))
    # hardcoded tax levels :(
    taxonomy_level_counter = {
        "order": Counter(), "class": Counter(), "phylum": Counter()
    }
    taxonomy_counter = Counter()
    summary_counter = Counter()
    with open(path) as fh:
        # skip the header
        next(fh)
        for i, line in enumerate(fh, start=1):
            toks = line.strip("\r\n").split("\t")
            if toks[3]:
                summary_counter.update(["Assigned\nFunction"])
                if toks[7]:
                    summary_counter.update(["Assigned\nBoth"])
            if toks[7]:
                taxonomy_counter.update([toks[7]])
                taxonomy = [j.strip() for j in toks[7].split(";")]
                taxonomy_level_counter["order"].update([taxonomy[3]])
                taxonomy_level_counter["class"].update([taxonomy[2]])
                taxonomy_level_counter["phylum"].update([taxonomy[1]])
                summary_counter.update(["Assigned\nTaxonomy"])
    idx = shannon_index(np.array(list(taxonomy_counter.values())))
    return dict(
        taxonomy_level_counter=taxonomy_level_counter,
        summary_counter=summary_counter,
        shannon=idx,
    )


def get_df_at_tax_level(count_obj, sample, tax_level):
    df = pd.DataFrame(data=count_obj, index=[0]).transpose()
    df.reset_index(inplace=True)
    df.columns = [tax_level, sample]
    return df


def get_sample_order(lst):
    """
    >>> lst = [["j", 2],["o", 1],["e", 3]]
    ['e', 'j', 'o']
    """
    return [i[0] for i in sorted(lst, key= lambda x: x[1], reverse=True)]


def process_reads(df, tax_level, sample_order):
    """
    This function will take the merged reads from the file and put them into a graphable
    form for both counts and percents. It will return both.
    """
    header = df.columns.tolist()
    df["summed_col"] = df[sample_order].sum(axis=1)
    # sort df by taxon abundance which organizes the traces in the plots
    df.sort_values(by="summed_col", ascending=False, inplace=True)
    # drop the extra column
    df = df[header]
    # sort by diversity
    sub = df[[tax_level] + sample_order].copy()
    cols = sub[tax_level].tolist()
    # find the relative percent
    sub_perc_t = sub[sample_order].div(sub[sample_order].sum())
    sub_perc_t = sub_perc_t * 100
    sub_perc_t = sub_perc_t.transpose()
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
    dfs = {}
    classifications_per_sample = {}
    for classification_table in classification_tables:
        sample = get_sample(classification_table, "_classifications.txt")
        parsed_taxonomy = parse_classifications_for_taxonomy(classification_table)
        samples.append([sample, parsed_taxonomy["shannon"]])
        # assigned #'s in summary table
        classifications_per_sample[sample] = parsed_taxonomy["summary_counter"]
        if len(dfs) == 0:
            for tax_level in tax_levels:
                dfs[tax_level] = get_df_at_tax_level(
                    parsed_taxonomy["taxonomy_level_counter"][tax_level],
                    sample,
                    tax_level,
                )
            continue

        for tax_level in tax_levels:
            df = get_df_at_tax_level(
                parsed_taxonomy["taxonomy_level_counter"][tax_level], sample, tax_level
            )
            dfs[tax_level] = dfs[tax_level].merge(df, on=tax_level, how="outer")
    # most diverse to least
    sample_order = get_sample_order(samples)
    observations_at_levels = {"Counts": dict(), "Percentage": dict()}
    for tax_level in tax_levels:
        c, p = process_reads(dfs[tax_level], tax_level, sample_order)
        observations_at_levels["Counts"][tax_level] = c
        observations_at_levels["Percentage"][tax_level] = p
    return observations_at_levels, pd.DataFrame.from_dict(
        classifications_per_sample, orient="index"
    )


def make_plots(observations, summary_type):
    # data traces are taxonomies across samples
    labels = {"Percentage": "Relative Abundance", "Counts": "Counts"}
    # tax levels are hardcoded at this point
    data = [
        go.Bar(
            x=observations[summary_type]["phylum"].index,
            y=observations[summary_type]["phylum"][tax],
            name=tax,
            text=tax,
            hoverinfo="text+y",
            visible=True,
        )
        for tax in observations[summary_type]["phylum"].columns.tolist()
    ] + [
        go.Bar(
            x=observations[summary_type]["class"].index,
            y=observations[summary_type]["class"][tax],
            name=tax,
            text=tax,
            hoverinfo="text+y",
            visible=False,
        )
        for tax in observations[summary_type]["class"].columns.tolist()
    ] + [
        go.Bar(
            x=observations[summary_type]["order"].index,
            y=observations[summary_type]["order"][tax],
            name=tax,
            text=tax,
            hoverinfo="text+y",
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
        height=700,
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


def parse_log_files(
    deduplication_logs, decontamination_logs, merge_logs, classifications_per_sample
):
    logger.info("Parsing log files for summary table")
    header = [
        "Sequences", "Unique", "Clean", "Pairs\nJoined", "Join\nRate", "Average\nInsert"
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
    log_df = log_df.merge(classifications_per_sample, left_index=True, right_index=True)
    log_df.reset_index(inplace=True)
    header.insert(0, "Sample")
    header.extend(["Assigned\nFunction", "Assigned\nTaxonomy", "Assigned\nBoth"])
    log_df.columns = header
    sample_summary_table = log_df[header].to_html(
        index=False, bold_rows=False, classes=["table", "table-bordered"], table_id="summaryTable"
    )
    sample_summary_table = sample_summary_table.replace("\n", "\n" + 10 * " ")
    return sample_summary_table


def build_quality_plot(r1_quality_files):
    logger.info("Building Sequence Quality plots")
    raw_qual_stats = defaultdict( lambda: defaultdict(list))
    for r1_ee_file in r1_quality_files:
        sample_name = get_sample(r1_ee_file, "_R1_eestats.txt")
        r2_ee_file = "_R2".join(r1_ee_file.rsplit("_R1", 1))
        with open(r1_ee_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                raw_qual_stats["R1"][sample_name].append(float(row["Mean_Q"]))
        with open(r2_ee_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                q = float(row["Mean_Q"])
                raw_qual_stats["R2"][sample_name].append(q)
    data = []
    for read_index, sample_data in raw_qual_stats.items():
        color = "rgba(31,119,180,0.3)" if read_index == "R1" else "rgba(214,39,40,0.3)"
        for sample_name, sample_stats in sample_data.items():
            data.append(
                go.Scatter(
                    x=list(range(1, len(sample_stats))),
                    y=sample_stats,
                    name=sample_name,
                    text=sample_name,
                    hoverinfo="text+x+y",
                    legendgroup=read_index,
                    mode="lines",
                    line=dict(color=color, dash="solid"),
                )
            )
    layout = go.Layout(
        title='Mean Quality Scores for R1 and R2',
        margin={"b": "auto", "r": "auto"},
        xaxis={"title": "Position"},
        yaxis={"title": "Quality (Phred score)"},
        hovermode="closest",
        showlegend=False,
        autosize=True,
        annotations=[
            dict(
                x=0,
                y=1.1,
                xref="paper",
                yref="paper",
                text="Forward",
                showarrow=False,
                font=dict(size=16, color="#ffffff"),
                align="left",
                borderpad=4,
                bgcolor="#1f77b4",
            ),
            dict(
                x=0.15,
                y=1.1,
                xref="paper",
                yref="paper",
                text="Reverse",
                showarrow=False,
                font=dict(size=16, color='#ffffff'),
                align='left',
                borderpad=4,
                bgcolor="#d62728",
            ),
        ],
    )
    fig = go.Figure(data=data, layout=layout)
    quality_plot = offline.plot(fig, **PLOTLY_PARAMS)
    return quality_plot


def get_conda_env_str(conda_env_file):
    # conda environment
    conda_env = ""
    with open(conda_env_file) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                conda_env += line
            else:
                conda_env += "    " + line
    return conda_env


def main(
    decontamination_logs,
    deduplication_logs,
    merge_logs,
    summary_tables,
    r1_quality_files,
    html,
    conda_env,
    function_table,
    taxonomy_table,
    taxonomy_function_table,
):
    observations_at_levels, classifications_per_sample = compile_summary_df(
        summary_tables
    )
    div = {}
    for v in ["Percentage", "Counts"]:
        div[v] = offline.plot(make_plots(observations_at_levels, v), **PLOTLY_PARAMS)
    html_tbl = parse_log_files(
        deduplication_logs, decontamination_logs, merge_logs, classifications_per_sample
    )
    quality_plot = build_quality_plot(r1_quality_files)
    conda_env = get_conda_env_str(conda_env)
    report_str = """

.. raw:: html

    {STYLE}
    {SCRIPT}

=============================================================
PerSeq_ - Per sequence functional and taxonomic assignments
=============================================================

.. _PerSeq: https://github.com/pnnl


.. contents::
    :backlinks: none


Summary
-------

Sequence Counts
***************

.. raw:: html

    <div style="overflow-x:auto;">
    {html_tbl}
    </div>

Sequence Quality
****************

.. raw:: html

    {quality_plot}

Taxonomy by Count
*****************

Samples are sorted based on their Shannon index calculated from taxonomically
annotated sequences. The order is most to least diverse.

.. raw:: html

    {div[Counts]}

Taxonomy by Percent
*******************

.. raw:: html

    {div[Percentage]}

Methods
-------

Paired-end sequences were evaluated for quality using VSEARCH [1]. Sequence
reads are deduplicated using the clumpify tool [2] then filtered of phiX and
rRNA using bbsplit [2]. Passing sequences are quality trimmed after successful
merging using bbmerge [2]. Sequences are allowed to be extended up 300 bp
during the merging process to account for non-overlapping R1 and R2 sequences
(``k=60 extend2=60 iterations=5 qtrim2=t``). Functional annotation and taxonomic
classification were performed on the merged sequences.

Functional Annotation
*********************

The blastx algorithm of DIAMOND [3] was used to align nucleotide sequences to
the KEGG protein reference database [4] consisting of non-redundant, family
level fungal eukaryotes and genus level prokaryotes (``--strand=both``). The
highest scoring alignment per sequence was used for functional annotation.

Taxonomic Annotation
********************

Kmer-based taxonomic classification was performed on the merged reads using
Kaiju [5] in greedy mode (``-a greedy -E 0.05``). NCBI's nr database [6]
containing reference sequences for archaea, bacteria, viruses, fungi, and
microbial eukaryotes was used as the reference index for Kaiju.

References
**********

1. Rognes T, Flouri T, Nichols B, Quince C, Mahé F. VSEARCH: a versatile open source tool for metagenomics. PeerJ. PeerJ Inc; 2016;4:e2584.
2. Bushnell B. BBTools [Internet]. Available from: https://sourceforge.net/projects/bbmap/
3. Buchfink B, Xie C, Huson DH. Fast and sensitive protein alignment using DIAMOND. Nat. Methods. Nature Publishing Group; 2015;12:59–60.
4. Kanehisa M, Sato Y, Kawashima M, Furumichi M, Tanabe M. KEGG as a reference resource for gene and protein annotation. Nucleic Acids Res. 2016;44:D457–62.
5. Menzel P, Ng KL, Krogh A. Fast and sensitive taxonomic classification for metagenomics with Kaiju. Nat Commun. Nature Publishing Group; 2016;7:11257.
6. NCBI Resource Coordinators. Database resources of the National Center for Biotechnology Information. Nucleic Acids Res. 2018;46:D8–D13.


Execution Environment
---------------------

::

    {conda_env}

Output
------

Per sample classifications in tables/ contain:

.. table::
    :name: classificationtable

    ====================  ==========================================================================================================================================
    Header ID             Definition
    ====================  ==========================================================================================================================================
    aa_alignment_length   The length of the DIAMOND blastx hit
    aa_percent_id         The percent ID of the DIAMOND blastx hit; could be used to increase post-processing stringency
    ec                    Enzyme Commission number from KEGG; semicolon delimited where multiple
    ko                    KEGG entry ID
    product               KEGG gene ID <semicolon> KEGG product
    read_id               The sequence identifier (unique)
    tax_alignment_length  The length of the Kaiju hit
    tax_classification    The Kaiju classification in order of superkingdom, phylum, order, class, family, genus, species; "NA" for each taxonomic level not defined
    ====================  ==========================================================================================================================================

Per taxonomy assignments in tables named taxonomy_<level>.txt contain:

.. table::
    :name: taxonomydeftable

    ====================  ======================================================
    Header ID             Definition
    ====================  ======================================================
    taxonomy_<level>      taxonomic level into which counts have been summed
    samples names         Non-normalized, per sample sum at this taxonomic level
    ====================  ======================================================

Per function assignments in tables named function_<type>.txt contain:

.. table::
    :name: functiondeftable

    ====================  ===================================================================
    Header ID             Definition
    ====================  ===================================================================
    <type>                either KO or EC into which counts have been summed
    samples names         Non-normalized, per sample sum for this particular functional group
    level_1               KEGG hierarchy [level 1] for KO defined in first column
    level_2               KEGG hierarchy [level 2] for KO defined in first column
    level_3               KEGG hierarchy [level 3] for KO defined in first column
    ====================  ===================================================================


Downloads
---------

"""
    report(
        report_str,
        html,
        file1=function_table,
        file2=taxonomy_table,
        file3=taxonomy_function_table,
        stylesheet="",
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--decontamination-logs", nargs="+")
    p.add_argument("--deduplication-logs", nargs="+")
    p.add_argument("--merge-logs", nargs="+")
    p.add_argument("--summary-tables", nargs="+")
    p.add_argument("--r1-quality-files", nargs="+")
    p.add_argument("--html")
    p.add_argument("conda_env")
    p.add_argument("function_table")
    p.add_argument("taxonomy_table")
    p.add_argument("taxonomy_function_table")
    args = p.parse_args()
    main(
        args.decontamination_logs,
        args.deduplication_logs,
        args.merge_logs,
        args.summary_tables,
        args.r1_quality_files,
        args.html,
        args.conda_env,
        args.function_table,
        args.taxonomy_table,
        args.taxonomy_function_table,
    )
