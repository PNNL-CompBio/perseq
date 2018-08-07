import argparse
import csv
import os
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import plotly
import relatively
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



def get_sample_order(lst):
    """
    >>> lst = [["j", 2],["o", 1],["e", 3]]
    ['e', 'j', 'o']
    """
    return [i[0] for i in sorted(lst, key=lambda x: x[1], reverse=True)]


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
#
#
# def compile_summary_df(classification_tables, tax_levels=["phylum", "class", "order"]):
#     """
#     Reads in multiple sample alignments from diamond in a given directory and merges them into
#     a single pandas.DataFrame. It returns a pandas dataframe for each of th
#     phylum that is ready to plug into the processing function. Also returns total counts
#     which is necessary to calculate the percentage of total that is being represented
#     """
#     samples = []
#     dfs = {}
#     classifications_per_sample = {}
#     for classification_table in classification_tables:
#         sample = get_sample(classification_table, "_classifications.txt")
#         parsed_taxonomy = parse_classifications_for_taxonomy(classification_table)
#         samples.append([sample, parsed_taxonomy["shannon"]])
#         # assigned #'s in summary table
#         classifications_per_sample[sample] = parsed_taxonomy["summary_counter"]
#         if len(dfs) == 0:
#             for tax_level in tax_levels:
#                 dfs[tax_level] = get_df_at_tax_level(
#                     parsed_taxonomy["taxonomy_level_counter"][tax_level],
#                     sample,
#                     tax_level,
#                 )
#             continue
#
#         for tax_level in tax_levels:
#             df = get_df_at_tax_level(
#                 parsed_taxonomy["taxonomy_level_counter"][tax_level], sample, tax_level
#             )
#             dfs[tax_level] = dfs[tax_level].merge(df, on=tax_level, how="outer")
#     return dfs
    # # most diverse to least
    # sample_order = get_sample_order(samples)
    # observations_at_levels = {"Counts": dict(), "Percentage": dict()}
    # for tax_level in tax_levels:
    #     c, p = process_reads(dfs[tax_level], tax_level, sample_order)
    #     observations_at_levels["Counts"][tax_level] = c
    #     observations_at_levels["Percentage"][tax_level] = p
    # return (
    #     observations_at_levels,
    #     pd.DataFrame.from_dict(classifications_per_sample, orient="index"),
    # )

def build_taxonomy_plot(txt, value_cols, height=900):
    df = pd.read_table(txt)
    levels = ["kingdom", "phylum", "class", "order"]
    df[levels] = df["taxonomy_order"].str.split(";", expand=True)
    df.head()
    hierarchy = ["phylum", "class", "order"]
    df[hierarchy] = df[hierarchy].fillna("NA")
    df[value_cols] = df[value_cols].fillna(0)
    df = df[hierarchy + value_cols]
    dfs = relatively.get_dfs_across_hierarchy(df, hierarchy, value_cols, reorder='shannon')
    fig = relatively.get_abundance_figure_from_dfs(
        dfs, hierarchy, "Assigned Taxonomy Per Sample", height=height
)
    return fig

##NEED to change back to work with snakemake
def get_sample(path, key):
    return [os.path.basename(item).partition(key)[0]for item in path]
#
# def get_sample(path, key):
#     return [os.path.basename(item).partition(key)[0] for item in os.listdir(path)]


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


def parse_log_files(merge_logs, unique_logs, clean_logs, classifications_per_sample):
    logger.info("Parsing log files for summary table")
    header = [
        "Sequences",
        "Pairs\nJoined",
        "Join\nRate",
        "Average\nInsert",
        "Unique",
        "Clean",
    ]
    count_table = defaultdict(list)
    # initial read count, join count, join rate, insert size from merge step
    for merge_log in merge_logs:
        sample = get_sample(merge_log, "_merge_sequences.log")
        count_table[sample].extend(parse_merge_file(merge_log))

    # unique count after merging and deduplication
    for unique_log in unique_logs:
        sample = get_sample(unique_log, "_02_unique_readlengths.txt")
        with open(unique_log) as fh:
            for line in fh:
                if line.startswith("#Reads:"):
                    count_table[sample].append(int(line.strip("\r\n").split("\t")[-1]))

    # clean count after merging, deduplication, and decontamination
    for clean_log in clean_logs:
        sample = get_sample(clean_log, "_03_clean_readlengths.txt")
        with open(clean_log) as fh:
            for line in fh:
                if line.startswith("#Reads:"):
                    count_table[sample].append(int(line.strip("\r\n").split("\t")[-1]))

    log_df = pd.DataFrame.from_dict(count_table, orient="index")
    log_df.columns = header
    log_df = log_df.merge(classifications_per_sample, left_index=True, right_index=True)
    log_df.reset_index(inplace=True)
    header.insert(0, "Sample")
    header.extend(["Assigned\nFunction", "Assigned\nTaxonomy", "Assigned\nBoth"])
    log_df.columns = header
    sample_summary_table = log_df[header].to_html(
        index=False,
        bold_rows=False,
        classes=["table", "table-bordered"],
        table_id="summaryTable",
    )
    sample_summary_table = sample_summary_table.replace("\n", "\n" + 10 * " ")
    return sample_summary_table


def build_quality_plot(r1_quality_files):
    logger.info("Building Sequence Quality plots")
    raw_qual_stats = defaultdict(lambda: defaultdict(list))
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
        title="Mean Quality Scores for R1 and R2",
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
                font=dict(size=16, color="#ffffff"),
                align="left",
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
    summary_tables,
    html,
    function_table,
    taxonomy_table,
    taxonomy_function_table,
    krona_tax,
    krona_ec
):
    value_cols = get_sample(summary_tables,'_classifications.txt')
    fig = build_taxonomy_plot(taxonomy_table,value_cols)
    plots = offline.plot(fig,**PLOTLY_PARAMS)
    # div = {}
    # for v in ["Percentage", "Counts"]:
    #     div[v] = offline.plot(make_plots(observations_at_levels, v), **PLOTLY_PARAMS)
    html_tbl = parse_log_files(
        merge_logs, unique_logs, clean_logs, classifications_per_sample
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
    :depth: 2


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


Taxonomy Assignment Summary
***************************

Samples are sorted based on their Shannon index calculated from taxonomically
annotated sequences. The order is most to least diverse.

.. raw:: html

    {plots}

Methods
-------

Paired-end sequences were evaluated for quality using VSEARCH [1]. Sequence
reads are quality trimmed after successful merging using bbmerge [2].
Sequences are allowed to be extended up 300 bp
during the merging process to account for non-overlapping R1 and R2 sequences
(``k=60 extend2=60 iterations=5 qtrim2=t``). Merged sequences are deduplicated
using the clumpify tool [2] then, by default, filtered of PhiX and
rRNA using bbsplit [2]. An arbitrary number of Name:FASTA pairs may be
specified during the decontamination process. Functional annotation and
taxonomic classification were performed following the decontamination step.

Functional Annotation
*********************

The blastx algorithm of DIAMOND [3] was used to align nucleotide sequences to
the KEGG protein reference database [4] consisting of non-redundant, family
level fungal eukaryotes and genus level prokaryotes
(``--strand=both --evalue 0.00001``). The highest scoring alignment per
sequence was used for functional annotation.

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

Classification Tables
*********************

Per sample classifications in tables/ contain:

.. table::
    :name: classificationtable

    =========================  ==========================================================================================================================================
    Header ID                  Definition
    =========================  ==========================================================================================================================================
    aa_alignment_length        The length of the DIAMOND blastx hit
    aa_percent_id              The percent ID of the DIAMOND blastx hit; could be used to increase post-processing stringency
    ec                         Enzyme Commission number from KEGG; semicolon delimited where multiple
    ko                         KEGG entry ID
    product                    KEGG gene ID <semicolon> KEGG product
    read_id                    The sequence identifier (unique)
    kaiju_alignment_length     The length of the Kaiju hit
    kaiju_classification       The Kaiju classification in order of superkingdom, phylum, class, order, family, genus, species; "NA" for each taxonomic level not defined
    blastx_lca_classification  The LCA result from the blastx HSPs
    =========================  ==========================================================================================================================================

Summary Tables
**************

Taxonomy
````````

Per taxonomy assignments in tables named **summaries/taxonomy/<level>.txt**
contain:

.. table::
    :name: taxonomydeftable

    ====================  ======================================================
    Header ID             Definition
    ====================  ======================================================
    taxonomy_<level>      taxonomic level into which counts have been summed
    samples names         non-normalized, per sample sum at this taxonomic level
    ====================  ======================================================

Function
````````

Per function assignments in tables named **summaries/function/<type>.txt**
contain:

.. table::
    :name: functiondeftable

    ====================  ===================================================================
    Header ID             Definition
    ====================  ===================================================================
    <type>                either KO, EC, or product into which counts have been summed
    samples names         non-normalized, per sample sum for this particular functional group
    level_1               KEGG hierarchy [level 1] if KO defined in first column
    level_2               KEGG hierarchy [level 2] if KO defined in first column
    level_3               KEGG hierarchy [level 3] if KO defined in first column
    ====================  ===================================================================

Combined
````````

Per taxonomy+function assignments in tables named
**summaries/combined/<type>_<level>.txt** contain:

.. table::
    :name: combineddeftable

    ====================  ====================================================================
    Header ID             Definition
    ====================  ====================================================================
    <type>                either KO, EC, or product; counts are summed using <type>+<taxonomy>
    taxonomy_<level>      taxonomic level; counts are summed using <type>+<taxonomy>
    sample names          non-normalized, per sample sum for this particular functional group
    level_1               KEGG hierarchy [level 1] if KO defined in first column
    level_2               KEGG hierarchy [level 2] if KO defined in first column
    level_3               KEGG hierarchy [level 3] if KO defined in first column
    ====================  ====================================================================

Downloads
---------

"""
    report(
        report_str,
        html,
        file1=function_table,
        file2=taxonomy_table,
        file3=taxonomy_function_table,
        kronaplot_tax=krona_tax,
        kronaplot_ec=krona_ec,
        stylesheet="",
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--summary-tables", nargs="+")
    p.add_argument("--html")
    p.add_argument("function_table")
    p.add_argument("taxonomy_table")
    p.add_argument("taxonomy_function_table")
    p.add_argument("krona_tax")
    p.add_argument("krona_ec")
    args = p.parse_args()
    main(
        args.summary_tables,
        args.html,
        args.function_table,
        args.taxonomy_table,
        args.taxonomy_function_table,
        args.krona_tax,
        args.krona_ec,
    )
