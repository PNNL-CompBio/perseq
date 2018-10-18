"""
per sequence classifications
"""
import argparse
import gzip
import logging
import os
import sys
from itertools import groupby

gzopen = lambda f: gzip.open(f, "rt") if f.endswith(".gz") else open(f)



def parse_kaiju_output(kaiju):
    kaiju_classifications = dict()
    with gzopen(kaiju) as fh:
        for line in fh:
            toks = line.strip("\r\n").split("\t")
            if toks[0] == "U":
                kaiju_classifications[toks[1]] = ["", ""]
            else:
                kaiju_classifications[toks[1]] = [toks[3], toks[7]]
    return kaiju_classifications


def main(
    kaiju, output, lca_threshold
):
    hamap_hits = parse_hamap(...)
    dbcan_hits = parse_dbcan(...)
    tax_classifications = parse_kaiju_output(kaiju)

    with gzopen() as ifh, open(output, "w") as ofh:
        print(
            "seq_id",
            "kaiju_length",
            "kaiju_score",
            "kaiju_classification",
            "hamap...stuff"
            "dbcan...stuff"
            sep="\t",
            file=ofh,
        )
        for seqid, seqgroup in groupby(ifh, key=lambda i: i.partition("\t")[0]):


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("kaiju", help="kaiju addTaxonNames output")
    p.add_argument("blastx", help="BLAST outfmt=6")
    p.add_argument(
        "gene2ko", help="convert blastx references to KO (hsa:9373<TAB>ko:K14018)"
    )
    p.add_argument("code2id", help="http://rest.kegg.jp/list/genome")
    p.add_argument(
        "kolist",
        help="convert KO to definition (ko:K00204<TAB>fwdH; 4Fe-4S ferredoxin)",
    )
    p.add_argument("namesdmp")
    p.add_argument("nodesdmp")
    p.add_argument("output", help="classification table")
    p.add_argument(
        "--lca-threshold",
        type=float,
        default=1,
        help="allow the lca to be determined at a node with this fraction of the hits",
    )
    args = p.parse_args()
    logging.basicConfig(
        level=logging.INFO, datefmt="%Y-%m-%d %H:%M", format="[%(asctime)s] %(message)s"
    )
    main(
        args.kaiju,
        args.blastx,
        args.gene2ko,
        args.code2id,
        args.kolist,
        args.namesdmp,
        args.nodesdmp,
        args.output,
        args.lca_threshold,
    )
