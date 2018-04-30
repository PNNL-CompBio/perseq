#!/usr/bin/env python
# coding=utf-8
"""
Important types:
    click.File(mode='r', encoding=None, errors='strict', lazy=None, atomic=False)
    click.Path(exists=True, writable=True)|click.Choice(["1", "2"])
Choices:
    @click.option('--item', nargs=2, type=click.Tuple([unicode, int]))
Specify multiple:
    @click.option('--message', '-m', multiple=True)
Boolean:
    @click.option('--shout', is_flag=True)

genes_ko.list.gz
hsa:578	ko:K14021
hsa:9373	ko:K14018

20180419_ko_list.txt
ko:K00001	E1.1.1.1, adh; alcohol dehydrogenase [EC:1.1.1.1]


read_id, aa_percent_id, aa_alignment_length, ko, product, ec, tax_alignment_length, tax_classification


"""

import click
import gzip
from itertools import groupby


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.argument("genes_ko")
@click.argument("ko_list")
@click.argument("kaiju_names")
@click.argument("diamond_hits")
@click.argument("outtable")
def main(genes_ko, ko_list, kaiju_names, diamond_hits, outtable):
    gene_map = {}
    with gzip.open(genes_ko, "rt") as fh:
        for line in fh:
            toks = line.strip().split("\t")
            gene_map[toks[0]] = toks[1]
    function_map = {}
    with open(ko_list) as fh:
        for line in fh:
            toks = line.strip().split("\t")
            if "[EC:" in toks[1]:
                product = toks[1].partition(" [EC:")[0]
                ec = toks[1].partition("[EC:")[-1].strip("]").replace(" ", ";")
            else:
                product = toks[1]
                ec = ""
            function_map[toks[0]] = (ec, product)
    tax_classifications = {}
    with open(kaiju_names) as fh:
        for line in fh:
            toks = line.strip().split("\t")
            if toks[0] == "U":
                tax_classifications[toks[1]] = ["", ""]
            else:
                tax_classifications[toks[1]] = [toks[3], toks[7]]
    with open(diamond_hits) as ifh, open(outtable, "w") as ofh:
        print("read_id", "aa_percent_id", "aa_alignment_length", "ko",
            "product", "ec", "tax_alignment_length", "tax_classification",
            sep="\t", file=ofh)
        for seqid, seqgroup in groupby(ifh, key=lambda i: i.partition("\t")[0]):
            for line in seqgroup:
                toks = line.strip().split("\t")
                read_id = toks[0]
                # no functional hit
                if toks[1] == "*":
                    aa_percent_id = -1
                    aa_alignment_length = -1
                    ko = ""
                    product = ""
                    ec = ""
                else:
                    try:
                        ko = gene_map[toks[1]]
                    except KeyError:
                        ko = ""
                    if ko:
                        ec, product = function_map[ko]
                    else:
                        ec = ""
                        product = ""
                    aa_percent_id = toks[2]
                    aa_alignment_length = toks[3]
                # get taxonomy; this should never cause a KeyError
                tax_alignment_length, tax_classification = tax_classifications[read_id]
                print(read_id, aa_percent_id, aa_alignment_length, ko, product,
                    ec, tax_alignment_length, tax_classification, sep="\t", file=ofh)
                # just print best HSP from DIAMOND
                break


if __name__ == "__main__":
    main()
