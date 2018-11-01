"""
per sequence classifications
"""
import argparse
import csv
import gzip
import logging
import os
import sys
from collections import Counter, defaultdict
from glob import glob
from itertools import groupby

gzopen = lambda f: gzip.open(f, "rt") if f.endswith(".gz") else open(f)


def parse_hmm_hits(filename):
    # contains list of split ACC [4], full seq score [1], and e-value [1]
    annotations = dict()
    with gzopen(filename) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for hit_id, hits in groupby(reader, key=lambda i: i[0]):
            for row in hits:
                print(row)
                # drops the HMM ID
                annotation = row[4].split("~~~")[0:3]
                print(annotation)
                # product strings within the ACC can't have spaces
                annotation[2] = annotation[2].replace("^", " ")
                annotations[hit_id] = annotation.extend(row[6:8])
                break
    return annotations


def parse_diamond_outputs(filenames):
    sequence_counts = defaultdict(Counter)
    for filename in filenames:
        sample = os.path.basename(filename).rstrip(".tsv")
        with gzopen(filename) as fh:
            reader = csv.reader(fh, delimiter="\t")
            for row in reader:
                sequence_counts[row[1]].update([sample])
    return sequence_counts


def main(kaiju, hamap, dbcan, tigrfams, diamond, output):
    diamond = glob(diamond)
    samples = [os.path.basename(i).rstrip(".tsv") for i in diamond]
    sequence_counts = parse_diamond_outputs(diamond)

    # ECs may be a comma delimited list
    # sequence_id -> [ec, gene, product, hmm_id, score, evalue]
    hamap_hits = parse_hmm_hits(hamap)
    # sequence_id -> [ec, enzyme class, enzyme class subfamily, hmm_id, score, evalue]
    dbcan_hits = parse_hmm_hits(dbcan)
    # sequence_id -> [ec, gene, product, hmm_id, score, evalue]
    tigrfams_hits = parse_hmm_hits(tigrfams)

    output_header = [
        "seq",
        "kaiju_length",
        "kaiju_taxonomy",
        "tigrfams_ec",
        "tigrfams_gene",
        "tigrfams_product",
        "tigrfams_score",
        "tigrfams_evalue",
        "hamap_ec",
        "hamap_gene",
        "hamap_product",
        "hamap_score",
        "hamap_evalue",
        "dbcan_ec",
        "dbcan_enzyme_class",
        "dbcan_enzyme_class_subfamily",
        "dbcan_score",
        "dbcan_evalue",
    ].extend(samples)
    # iterate over the kaiju hits and print table
    with gzopen(kaiju) as kaiju_hits, open(output, "w") as out_fh:
        print(*output_header, sep="\t", file=out_fh)
        reader = csv.reader(kaiju_hits, delimiter="\t")
        for hit in reader:
            seq = hit[1]
            if hit[0] == "U":
                kaiju_length = 0
                kaiju_taxonomy = ""
            else:
                kaiju_length = hit[3]
                kaiju_taxonomy = hit[7]
            try:
                tigrfams_ec, tigrfams_gene, tigrfams_product, tigrfams_score, tigrfams_evalue = tigrfams_hits[seq]
            except KeyError:
                tigrfams_ec, tigrfams_gene, tigrfams_product, tigrfams_score, tigrfams_evalue = [""] * 5
            try:
                hamap_ec, hamap_gene, hamap_product, hamap_score, hamap_evalue = hamap_hits[seq]
            except KeyError:
                hamap_ec, hamap_gene, hamap_product, hamap_score, hamap_evalue = [""] * 5
            try:
                dbcan_ec, dbcan_enzyme_class, dbcan_enzyme_class_subfamily, dbcan_score, dbcan_evalue = dbcan_hits[seq]
            except KeyError:
                dbcan_ec, dbcan_enzyme_class, dbcan_enzyme_class_subfamily, dbcan_score, dbcan_evalue = [""] * 5
            hit_counts = [sequence_counts[seq][sample] for sample in samples]
            print(
                seq,
                kaiju_length,
                kaiju_taxonomy,
                tigrfams_ec,
                tigrfams_gene,
                tigrfams_product,
                tigrfams_evalue,
                hamap_ec,
                hamap_gene,
                hamap_product,
                hamap_score,
                hamap_evalue,
                dbcan_ec,
                dbcan_enzyme_class,
                dbcan_enzyme_class_subfamily,
                dbcan_score,
                dbcan_evalue,
                *hit_counts,
                sep="\t",
                file=out_fh,
            )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("kaiju", help="kaiju addTaxonNames output")
    p.add_argument("hamap")
    p.add_argument("dbcan")
    p.add_argument("tigrfams")
    p.add_argument("diamond")
    p.add_argument("--output", help="annotation table with sample counts")
    args = p.parse_args()
    logging.basicConfig(
        level=logging.INFO, datefmt="%Y-%m-%d %H:%M", format="[%(asctime)s] %(message)s"
    )
    main(args.kaiju, args.hamap, args.dbcan, args.tigrfams, args.diamond, args.output)
