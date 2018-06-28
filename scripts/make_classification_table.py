"""
per sequence classifications
"""
import argparse
import errno
import gzip
import logging
import os
import pickle
import sys
from collections import Counter, OrderedDict, defaultdict, deque
from itertools import groupby

gzopen = lambda f: gzip.open(f, "rt") if f.endswith(".gz") else open(f)


class Node(object):
    def __init__(self, taxonomy, node_id, parent_id, tax_level):
        """Represents a node within a tree.
        Args:
            taxonomy (str): taxonomy name or ID
            node_id (str): taxonomy ID
            parent_id (str): taxonomy ID of parent
            tax_level (str): the taxonomic level for this node_id
        """
        # the current node's string ID
        self.taxonomy = taxonomy
        # the current node's digit ID
        self.node_id = node_id
        self.parent_id = parent_id
        self.tax_level = tax_level


class Tree(object):
    def __init__(self, namesdmp, nodesdmp):
        """Builds reference dictionary of Taxonomy Name, Taxonomy ID, and Parent Taxonomy ID."""
        self.tax_levels = [
            "superkingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ]
        self.tax_levels_abb = ["k", "p", "c", "o", "f", "g", "s"]
        self.tree = defaultdict(dict)
        tax_to_name = {}
        logging.info("Reading in %s" % namesdmp)
        with open(namesdmp) as dmp:
            for tax_id, group in groupby(
                dmp, key=lambda x: [i.strip() for i in x.strip().split("|")][0]
            ):
                scientific_name = ""
                synonym = ""
                for line in group:
                    toks = [
                        x.strip().replace("'", "").replace('"', "")
                        for x in line.strip().split("|")
                    ]
                    if toks[-2] == "scientific name":
                        tax_to_name[toks[0]] = toks[1]
                        # we really just want the scientific name
                        break
                    elif toks[-2] == "misspelling":
                        continue
                    # prioritize synonyms over authority, type material, etc.
                    elif toks[-2] == "synonym":
                        synonym = toks[1]
                        tax_to_name[toks[0]] = synonym
                    else:
                        if synonym:
                            tax_to_name[toks[0]] = synonym
                        else:
                            tax_to_name[toks[0]] = toks[1]
        logging.info("Reading in %s" % nodesdmp)
        with open(nodesdmp) as dmp:
            for line in dmp:
                toks = [x.strip() for x in line.strip().split("|")]
                if not toks[0] == "1" and not toks[1] == "1":
                    assert not toks[0] == toks[1]
                self.add_node(tax_to_name[toks[0]], toks[0], toks[1], toks[2])

    def add_node(self, taxonomy, node_id, parent_id, tax_level):
        """
        Adds node to tree dictionary.

        Args:
            taxonomy (str): the taxonomy name
            node_id (str): the taxonomy id
            parent_id (str): the parent's taxonomy id
            tax_level (str): the taxonomic level for this node_id
        """
        # taxonomy id to node mapping; ensures unique nodes despite non-unique names
        self.tree[node_id] = Node(taxonomy, node_id, parent_id, tax_level)

    def lca(self, taxonomies, threshold=1.):
        """
        Returns the taxonomy of the LCA and optionally only use the top
        fraction of hits.

        Args:
            taxonomies (list): list of node id of corresponding taxonomic
                assignment; when using threshold < 1 they should be ordered
                by decreasing bitscore
            threshold (Optional[float]): 0-1; threshold fraction of hits to be
                factored into lca

        Returns:
            str: node id of LCA
        """
        if threshold > 1:
            threshold = 1
        elif threshold < 0.01:
            # 10% as the minimum
            threshold = 0.1
        count_target = len(taxonomies) * threshold
        count_taxonomies = Counter()
        for taxonomy in taxonomies:
            try:
                current_taxonomy = self.tree[taxonomy].node_id
            except AttributeError:
                # dict when key not present
                # taxonomy represented in the reference database, but is not present in the tree
                continue

            while not current_taxonomy == "1":
                count_taxonomies.update([current_taxonomy])
                if count_taxonomies[current_taxonomy] >= count_target:
                    return self.tree[current_taxonomy].node_id
                # traverse up tree
                current_taxonomy = self.tree[current_taxonomy].parent_id
        return "1"

    def taxonomic_lineage(self, taxonomy):
        if taxonomy == "1":
            return [taxonomy]

        lineage = [taxonomy]
        while not taxonomy == "1":
            taxonomy = self.tree[taxonomy].parent_id
            # prepend
            lineage.insert(0, taxonomy)
        return lineage

    def get_lineage_str(self, taxonomy, sep=";"):
        """
        Args:
            tax_id (str): Lowest common ancestor. Can be determined by the tree.lca function
        Returns:
            str: Taxonomic classification in this format k__;p__;c__;o__;__;g__;s__
        """
        lineage = {}
        for item in self.taxonomic_lineage(taxonomy):
            node = self.tree[item]
            # no rank class subclass infraclass superclass cohort tribe
            # varietas forma subtribe family subfamily superfamily genus subgenus
            # order infraorder parvorder superorder suborder phylum subphylum
            # superphylum species species group species subgroup subspecies
            # kingdom superkingdom subkingdom
            if node.tax_level in self.tax_levels:
                # does not account for "no rank" and some other cases of "unclassified"
                lineage[
                    "k" if node.tax_level == "superkingdom" else node.tax_level[0]
                ] = node.taxonomy
        valid_lineage = []
        for idx in self.tax_levels_abb:
            valid_lineage.append(
                "%s__%s" % (idx, lineage.get(idx, "?").replace(",", ""))
            )
        return sep.join(valid_lineage)


def parse_gene2ko(gene2ko):
    """
    hsa:9373  ko:K14018
    hsa:9344  ko:K04429
    hsa:5894  ko:K04366
    """
    logging.info("Parsing %s" % gene2ko)
    gene_map = dict()
    with gzopen(gene2ko) as fh:
        for line in fh:
            toks = line.strip("\r\n").split("\t")
            gene_map[toks[0]] = toks[1]
    return gene_map


def parse_kolist(kolist):
    """
    ko:K00001   E1.1.1.1, adh; alcohol dehydrogenase [EC:1.1.1.1]
    ko:K00002   AKR1A1, adh; alcohol dehydrogenase (NADP+) [EC:1.1.1.2]
    """
    logging.info("Parsing %s" % kolist)
    ko_map = dict()
    with gzopen(kolist) as fh:
        for line in fh:
            toks = line.strip("\r\n").split("\t")
            if "[EC:" in toks[1]:
                # iolG; myo-inositol 2-dehydrogenase / D-chiro-inositol 1-dehydrogenase [EC:1.1.1.18 1.1.1.369]
                product = toks[1].partition(" [EC:")[0]
                # [EC:1.1.1.18 1.1.1.369]
                ec = toks[1].partition("[EC:")[-1].strip("]").replace(" ", ";")
            else:
                product = toks[1]
                ec = ""
            ko_map[toks[0]] = [ec, product]
    return ko_map


def parse_kegg_code_to_ncbi_tax(species2ncbi):
    """
    Takes a converter file that builds a dictionary where the key is the
    KEGG ID and the value is a tuple of the NCBI ID and the associated node
    input example string:

    'gn:T00001    hin, HAEIN, 71421; Haemophilus influenzae Rd KW20 (serotype d)'

    Can be downloaded from http://rest.kegg.jp/list/genome

    Returns:
        dict of kegg species to ncbi tax id
    """
    logging.info("Parsing %s" % species2ncbi)
    namemap = dict()
    with gzopen(species2ncbi) as fh:
        for line in fh:
            if not ";" in line:
                continue
            toks = line.partition("; ")[0].partition("\t")[-1]
            toks = [i.strip() for i in toks.split(",")]
            namemap[toks[0]] = toks[-1]
    return namemap


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
    kaiju, blastx, gene2ko, code2id, kolist, namesdmp, nodesdmp, output, lca_threshold
):
    gene_map = parse_gene2ko(gene2ko)
    function_map = parse_kolist(kolist)
    kegg_to_ncbi = parse_kegg_code_to_ncbi_tax(code2id)
    tax_classifications = parse_kaiju_output(kaiju)
    tree = Tree(namesdmp, nodesdmp)

    with gzopen(blastx) as ifh, open(output, "w") as ofh:
        print(
            "read_id",
            "aa_percent_id",
            "aa_alignment_length",
            "ko",
            "product",
            "ec",
            "kaiju_alignment_length",
            "kaiju_classification",
            "blastx_lca_classification",
            sep="\t",
            file=ofh,
        )
        for seqid, seqgroup in groupby(ifh, key=lambda i: i.partition("\t")[0]):
            tax_ids = []
            for line in seqgroup:
                toks = line.strip("\r\n").split("\t")
                read_id = toks[0]
                # no functional hit
                if toks[1] == "*":
                    aa_percent_id = -1
                    aa_alignment_length = -1
                    ko = ""
                    product = ""
                    ec = ""
                    lca_tax = ""
                else:
                    try:
                        ko = gene_map[toks[1]]
                        ec, product = function_map[ko]
                    except KeyError:
                        ko = ""
                        ec = ""
                        product = ""
                    try:
                        # 2058097,498374,1968910,1968909,1643678,651137,*1783275*,2157,131567,1
                        tax_id = kegg_to_ncbi[toks[1].partition(":")[0]]
                        tax_ids.append(tax_id)
                    except KeyError:
                        # species code was not able to be converted
                        pass
                    aa_percent_id = toks[2]
                    aa_alignment_length = toks[3]
            if tax_ids:
                lca = tree.lca(tax_ids, threshold=lca_threshold)
                lineage = tree.get_lineage_str(lca)
            else:
                lineage = ""
            tax_alignment_length, tax_classification = tax_classifications[read_id]
            print(
                read_id,
                aa_percent_id,
                aa_alignment_length,
                ko,
                product,
                ec,
                tax_alignment_length,
                tax_classification,
                lineage,
                sep="\t",
                file=ofh,
            )


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
