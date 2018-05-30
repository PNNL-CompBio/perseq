"""
per sequence classifications
"""

import os
import sys
from snakemake.logging import logger
from snakemake.utils import report
from collections import Counter, defaultdict, deque, OrderedDict
##tree class
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

    def __init__(self, tree_file):
        """Builds reference dictionary of Taxonomy Name, Taxonomy ID, and Parent Taxonomy ID."""
        self.tree = defaultdict(dict)
        with open(tree_file) as tree_fh:
            for line in tree_fh:
                toks = line.strip().split("\t")
                if not toks[0] == '1' and not toks[2] == '1':
                    assert not toks[0] == toks[2]
                if not len(toks) == 4:
                    logging.warning("Line [%s] does not have ID, NAME, PARENTID, TAX LEVEL" % line.strip())
                    continue
                self.add_node(toks[1], toks[0], toks[2], toks[3])




    def add_node(self, taxonomy, node_id, parent_id, tax_level):
        """Adds node to tree dictionary.
        Args:
            taxonomy (str): the taxonomy name
            node_id (str): the taxonomy id
            parent_id (str): the parent's taxonomy id
            tax_level (str): the taxonomic level for this node_id
        """

        # taxonomy id to node mapping; ensures unique nodes despite non-unique names
        self.tree[node_id] = Node(taxonomy, node_id, parent_id, tax_level)


    def lca(self, taxonomies, threshold=1.):

        if threshold > 1:
            threshold = 1
        elif threshold < 0.01:
            # 10% as the minimum
            threshold = 0.1

        count_target = len(taxonomies) * threshold
        count_taxonomies = defaultdict(int)

        for taxonomy in taxonomies:
            #print('tax',taxonomy)

            try:
                current_taxonomy = self.tree[taxonomy].node_id
                #print(current_taxonomy)
            except AttributeError:
                # dict when key not present
                # taxonomy represented in the reference database, but is not present in the tree
                continue

            while not current_taxonomy == "1":
                #print('curtx',current_taxonomy)
                count_taxonomies[current_taxonomy] += 1
                if count_taxonomies[current_taxonomy] >= count_target:
                    #print('target found')
                    return self.tree[current_taxonomy].node_id

                # traverse up tree
                current_taxonomy = self.tree[current_taxonomy].parent_id
        return "1"

    def taxonomic_lineage(self, taxonomy):

        # a large portion of ORFs
        if taxonomy == "1":
            return [taxonomy]

        lineage = [taxonomy]
        while not taxonomy == "1":
            taxonomy = self.tree[taxonomy].parent_id
            # prepend
            lineage.insert(0, taxonomy)
        return lineage




TAX_LEVELS=["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
def validate_lineage(lineage, sep=";"):

    levels = ["k" if tax_level == "superkingdom" else tax_level[0] for tax_level in TAX_LEVELS]
    #print(levels)
    valid_lineage = []
    for idx in levels:
        # removes commas in tax names
        #print(idx)
        valid_lineage.append("%s__%s" % (idx, lineage.get(idx, "?").replace(",", "")))
    return sep.join(valid_lineage)

#calls validate lineage
def lineage_form(lca):

    lineage = {}
    for item in test.taxonomic_lineage(lca):
        node = test.tree[item]
        #print(node)
        if node.tax_level in TAX_LEVELS:
            #print(node.tax_level)
            # does not account for "no rank" and some other cases of "unclassified"
            lineage["k" if node.tax_level == "superkingdom" else node.tax_level[0]] = node.taxonomy
            print(node.tax_level)
            #print(lineage)
    print(lineage)
    lineage = validate_lineage(lineage)
    return lineage




#this converts the kegg ID into the NCBI ID
def convert_tax(converter_file):

    with open(converter_file,'r') as tax_id:
        converter={}
        for line in tax_id:

            toks=line.strip().split('\t')[1]

            bits=toks.split(',')

            toks=line.strip().split('\t')[1]

            bits=toks.split(';')


            lil_bits=bits[0].split(',')
            #print(lil_bits)
            if len(lil_bits) ==3:
                kegg_id=lil_bits[0]
                node_id=lil_bits[2]
                NCBI_id=bits[1].strip()
            elif len(lil_bits)<3:
                lil_bits.insert(1,'X')
                kegg_id=lil_bits[0]
                node_id=lil_bits[2]
                NCBI_id=bits[1].strip()
            else:
                break

            converter[kegg_id]=(NCBI_id,node_id)
        return converter

def grab_node(aligned):

        toks=aligned.strip().split('\t')
        kegg_id=str(toks[1].split(':')[0])
        #print(kegg_id)
        try:

            if kegg_id in converter:

                return converter[kegg_id](1)

        except KeyError:

            pass


def get_kaiju_db_dir(config):
    db_dir = ""
    if "kaijudb" in config:
        db_dir = config["kaijudb"]
    elif os.environ.get("KAIJUDB"):
        db_dir = os.environ.get("KAIJUDB")
    else:
        logger.error("Either set KAIJUDB or pass --config kaijudb=/path")
        sys.exit(1)
    if not os.path.exists(os.path.join(db_dir, "kaiju_db.fmi")):
        logger.error("kaiju_db.fmi does not exist in your database directory")
        sys.exit(1)
    if not os.path.exists(os.path.join(db_dir, "names.dmp")):
        logger.error("names.dmp does not exist in your database directory")
        sys.exit(1)
    if not os.path.exists(os.path.join(db_dir, "nodes.dmp")):
        logger.error("nodes.dmp does not exist in your database directory")
        sys.exit(1)
    return db_dir


KAIJUDB = get_kaiju_db_dir(config)


rule all:
    input:
        expand("tables/{sample}_classifications.txt", sample=config["samples"].keys())


rule deduplicate_reads:
        input:
            r1 = lambda wildcards: config["samples"][wildcards.sample][0],
            r2 = lambda wildcards: config["samples"][wildcards.sample][1]
        output:
            r1 = "quality_control/{sample}_00_deduplicate_R1.fastq.gz",
            r2 = "quality_control/{sample}_00_deduplicate_R2.fastq.gz"
        log:
            "logs/{sample}_deduplicate_reads.log"
        # conda:
        #     "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        resources:
            java_mem = config.get("java_mem", 60)
        shell:
            """
            clumpify.sh in={input.r1} in2={input.r2} \
                out={output.r1} out2={output.r2} \
                dedupe=t optical=t threads={threads} \
                -Xmx{resources.java_mem}G 2> {log}
            """


# rule apply_quality_filter:
#     input:
#         r1 = "quality_control/{sample}_00_deduplicate_R1.fastq.gz",
#         r2 = "quality_control/{sample}_00_deduplicate_R2.fastq.gz"
#     output:
#         r1 = "quality_control/{sample}_01_trimmed_R1.fastq.gz",
#         r2 = "quality_control/{sample}_01_trimmed_R2.fastq.gz"
#     params:
#         trimq = config.get("minimum_base_quality", 10),
#         qtrim = config.get("qtrim", "rl"),
#         minlength = config.get("minimum_passing_read_length", 51),
#         minbasefrequency = config.get("minimum_base_frequency", 0.05)
#     log:
#         "logs/{sample}_apply_quality_filter.log"
#     # conda:
#     #     "%s/required_packages.yaml" % CONDAENV
#     threads:
#         config.get("threads", 1)
#     resources:
#         java_mem = config.get("java_mem", 60)
#     shell:
#         """
#         bbduk.sh in={input.r1} in2={input.r2} \
#             out={output.r1} out2={output.r2} \
#             qout=33 trd=t trimq={params.trimq} qtrim={params.qtrim} \
#             threads={threads} minlength={params.minlength} \
#             minbasefrequency={params.minbasefrequency} \
#             -Xmx{resources.java_mem}G 2> {log}
#         """


rule build_decontamination_db:
    output:
        "ref/genome/1/summary.txt"
    # conda:
    #     "%s/required_packages.yaml" % CONDAENV
    params:
        k = config.get("contaminant_kmer_length", 13),
        refs_in = " ".join("ref_%s=%s" % (n, fa) for n, fa in [["PhiX", config["phix"]], ["rRNA", config["rRNA"]]])
    resources:
        java_mem = config.get("java_mem", 60)
    threads:
        config.get("threads", 1)
    shell:
        """
        bbsplit.sh -Xmx{resources.java_mem}G {params.refs_in} threads={threads} k={params.k} local=t
        """


rule run_decontamination:
    input:
        r1 = "quality_control/{sample}_00_deduplicate_R1.fastq.gz",
        r2 = "quality_control/{sample}_00_deduplicate_R2.fastq.gz",
        db = "ref/genome/1/summary.txt"
    output:
        r1 = "quality_control/{sample}_02_decontamination_R1.fastq.gz",
        r2 = "quality_control/{sample}_02_decontamination_R2.fastq.gz",
        stats = "logs/{sample}_decontamination_by_reference.log"
    params:
        maxindel = config.get("contaminant_max_indel", 5),
        minratio = config.get("contaminant_min_ratio", 0.80),
        minhits = config.get("contaminant_minimum_hits", 3),
        ambiguous = config.get("contaminant_ambiguous", "best"),
        k = config.get("contaminant_kmer_length", 15),
    log:
        "logs/{sample}_decontamination.log"
    # conda:
    #     "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        java_mem = config.get("java_mem", 60),
    shell:
        """
        bbsplit.sh in1={input.r1} in2={input.r2} \
            outu1={output.r1} outu2={output.r2} refstats={output.stats} \
            maxindel={params.maxindel} minratio={params.minratio} \
            minhits={params.minhits} ambiguous={params.ambiguous} \
            threads={threads} k={params.k} \
            local=t -Xmx{resources.java_mem}G 2> {log}
        """


rule merge_sequences:
    input:
        r1 = "quality_control/{sample}_02_decontamination_R1.fastq.gz",
        r2 = "quality_control/{sample}_02_decontamination_R2.fastq.gz"
    output:
        merged = "quality_control/{sample}_03_merged.fastq.gz",
        r1 = "quality_control/{sample}_03_unmerged_R1.fastq.gz",
        r2 = "quality_control/{sample}_03_unmerged_R2.fastq.gz"
    params:
        adapters = "" if not config.get("adapters") else "adapter=%s" % config.get("adapters")
    log:
        "logs/{sample}_merge_sequences.log"
    threads:
        config.get("threads", 1)
    resources:
        java_mem = config.get("java_mem", 60),
    shell:
        """
        bbmerge.sh threads={threads} k=60 extend2=60 iterations=5 \
            reassemble=t shave rinse prealloc=t -Xmx{resources.java_mem}G \
            loose=t qtrim2=t in={input.r1} in2={input.r2} \
            {params.adapters} out={output.merged} \
            outu={output.r1} outu2={output.r2} 2> {log}
        """


rule run_taxonomic_classification:
    input:
        "quality_control/{sample}_03_merged.fastq.gz"
    output:
        temp("kaiju/{sample}_aln.txt")
    params:
        evalue = config.get("kaiju_evalue", 0.05)
    threads:
        config.get("threads", 1)
    shell:
        """
        kaiju -t {KAIJUDB}/nodes.dmp \
            -f {KAIJUDB}/kaiju_db.fmi \
            -i {input} -o {output} -z {threads} \
            -a greedy -x -v -E {params.evalue}
        """


rule add_full_taxonomy:
    input:
        "kaiju/{sample}_aln.txt"
    output:
        "kaiju/{sample}_aln_names.txt"
    shell:
        """
        addTaxonNames -t {KAIJUDB}/nodes.dmp -n {KAIJUDB}/names.dmp \
            -i {input} -o {output} \
            -r superkingdom,phylum,order,class,family,genus,species
        """


rule build_diamond_index:
    input:
        config["diamonddb"]
    output:
        config["diamonddb"] + ".dmnd"
    threads:
        config.get("threads", 1)
    shell:
        """
        diamond makedb --in {input} --threads {threads} --db {output}
        """


rule run_functional_classification:
    input:
        fq = "quality_control/{sample}_03_merged.fastq.gz",
        db = config["diamonddb"] + ".dmnd"
    output:
        "diamond/{sample}_aln.txt"
    threads:
        config.get("threads", 1)
    shell:
        """
        diamond blastx --threads {threads} --db {input.db} \
            --out {output} --outfmt 6 --query {input.fq} --strand both \
            --unal 1 --top 3 --block-size 4 --index-chunks 1
        """


rule combine_sample_output:
    input:
        tax = "kaiju/{sample}_aln_names.txt",
        fun = "diamond/{sample}_aln.txt"
    output:
        "tables/{sample}_classifications.txt"
    params:
        gene2ko = config["gene2ko"],
        kodefs = config["kodefs"]
    run:
        import gzip
        from itertools import groupby

        genes_ko = params["gene2ko"]
        ko_list = params["kodefs"]
        kaiju_names = input["tax"]
        diamond_hits = input["fun"]
        outtable = output[0]

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
        # kegg_map[three_letter_code] = ncbi_tax_id
        kegg_map = convert_tax(converter_file) #this will return the whole dict
        with open(diamond_hits) as ifh, open(outtable, "w") as ofh:
            print("read_id", "aa_percent_id", "aa_alignment_length", "ko",
                "product", "ec", "kaiju_alignment_length",
                "kaiju_classification", "lca_classification",
                sep="\t", file=ofh)
            for seqid, seqgroup in groupby(ifh, key=lambda i: i.partition("\t")[0]):
                taxonomies = []
                for line in seqgroup:
                    toks = line.strip().split("\t")
                    read_id = toks[0]
                    # qseqid
                    sseqid= toks[1]
                    try:

                    tax_id = kegg_map[toks["sseqid"].partition(":")[0]](1)
                except KeyError:
                    pass

                    # no functional hit
                    if toks[1] == "*":
                        aa_percent_id = -1
                        aa_alignment_length = -1
                        ko = ""
                        product = ""
                        ec = ""
                        # TODO
                        lca_classification = "k__?;p__?;c__?;o__?;f__?;g__?;s__?"
                        print(read_id, aa_percent_id, aa_alignment_length, ko, product,
                            ec, tax_alignment_length, tax_classification, lca_classification, sep="\t", file=ofh)
                        # just print the best HSP
                        break
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
                        # TODO
                        # taxonomies = lca_finder(line)
                        taxonomies.append(tax_id)

                # get taxonomy; this should never cause a KeyError
                lca=test.lca(taxonomies)
                lca_classification=lineage_form(lca)
                tax_alignment_length, tax_classification = tax_classifications[read_id]
                print(read_id, aa_percent_id, aa_alignment_length, ko, product,
                    ec, tax_alignment_length, tax_classification, lca_classification, sep="\t", file=ofh)
                # just print the best HSP
                break
rule build_tax_ko:
    input:
        json='/Users/zavo603/Documents/Nicki_files/perseq/ko00001.json',
        path_to_class_tables = "tables/{sample}_classifications.txt"
    output:
        'tax_ko_tbl.txt'
    threads:
        config.get("threads", 1)
    shell:
        """
        tax_ko_info.py --json_file_path {input.json} \
        --file_path {input.path_to_class_tables} \
        --tax_level --min_perc_id --min_len --group_on --output {output}

        """
