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


def parse_hmm_hits(filename):
    # contains list of split ACC [4], full seq score [1], and e-value [1]
    annotations = dict()
    with gzopen(filename) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for hit_id, hits in groupby(reader, key=lambda i: i[0]):
            for hit in hits:
                annotations[hit_id] = row[4].split("~~~").extend(row[6:8])
                break
    return annotations


def main(kaiju, hamap, dbcan, tigrfams, output):
    # row[4].split("~~~") -> ec, gene, product.replace("^", " "), HMM ID
    hamap_hits = parse_hmm_hits(hamap)
    # row[4].split("~~~") -> ec, enzyme class, enzyme class subfamily, HMM ID
    # ECs may be a comma delimited list
    dbcan_hits = parse_hmm_hits(dbcan)
    # row[4].split("~~~") -> ec, gene, product.replace("^", " "), HMM ID
    # ECs may be a comma delimited list
    tigrfams_hits = parse_hmm_hits(tigrfams)

    with gzopen(kaiju) as ifh, open(output, "w") as ofh:
        print(
            "seq_id",
            "kaiju_length",
            "kaiju_score",
            "kaiju_classification",
            "hamap_ec",
            "hamap_gene",
            "hamap_product",
            "hamap_score",
            "hamap_evalue",
            "dbcan...",
            "tigrfams...",
            sep="\t",
            file=ofh,
        )
        for seqid, seqgroup in groupby(ifh, key=lambda i: i.partition("\t")[0]):


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("kaiju", help="kaiju addTaxonNames output")
    p.add_argument("hamap")
    p.add_argument("dbcan")
    p.add_argument("tigrfams")
    p.add_argument("output", help="classification table")
    args = p.parse_args()
    logging.basicConfig(
        level=logging.INFO, datefmt="%Y-%m-%d %H:%M", format="[%(asctime)s] %(message)s"
    )
    main(
        args.kaiju,
        args.hamap,
        args.dbcan,
        args.tigrfams,
        args.output,
    )
