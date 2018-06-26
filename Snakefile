"""
per sequence classifications
"""
import os
import sys

from collections import Counter, defaultdict, deque, OrderedDict

from snakemake.logging import logger
from snakemake.utils import report


TAX_LEVELS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

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
        """
        taxonomies (list): FIXME
        """
        if threshold > 1:
            threshold = 1
        elif threshold < 0.01:
            # 10% as the minimum
            threshold = 0.1
        count_target = len(taxonomies) * threshold
        # FIXME this should be a Counter
        count_taxonomies = defaultdict(int)
        for taxonomy in taxonomies:
            try:
                current_taxonomy = self.tree[taxonomy].node_id
            except AttributeError:
                # dict when key not present
                # taxonomy represented in the reference database, but is not present in the tree
                continue

            while not current_taxonomy == "1":
                count_taxonomies[current_taxonomy] += 1
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


def validate_lineage(lineage, sep=";"):
    levels = ["k" if tax_level == "superkingdom" else tax_level[0] for tax_level in TAX_LEVELS]
    valid_lineage = []
    for idx in levels:
        valid_lineage.append("%s__%s" % (idx, lineage.get(idx, "?").replace(",", "")))
    return sep.join(valid_lineage)


def lineage_form(lca):
    """
    calls validate lineage <-- technically, yes, but something about what
    it accepts as a parameter and maybe something about what this method
    returns would be more appropriate
    """
    lineage = {}
    # FIXME bug
    for item in test.taxonomic_lineage(lca):
        # FIXME
        node = test.tree[item]
        if node.tax_level in TAX_LEVELS:
            # does not account for "no rank" and some other cases of "unclassified"
            lineage["k" if node.tax_level == "superkingdom" else node.tax_level[0]] = node.taxonomy
    lineage = validate_lineage(lineage)
    return lineage


def convert_tax(converter_file):
    """converts the kegg ID into the NCBI ID"""

    with open(converter_file) as tax_id:
        converter = {}
        for line in tax_id:

            # FIXME something doesn't seem right here
            toks=line.strip().split('\t')[1]
            bits=toks.split(',')
            toks=line.strip().split('\t')[1]
            bits=toks.split(';')
            lil_bits=bits[0].split(',')

            if len(lil_bits) == 3:
                kegg_id = lil_bits[0]
                node_id = lil_bits[2]
                NCBI_id = bits[1].strip()
            elif len(lil_bits) < 3:
                # FIXME why?
                lil_bits.insert(1, 'X')
                kegg_id = lil_bits[0]
                node_id = lil_bits[2]
                NCBI_id = bits[1].strip()
            else:
                break

            converter[kegg_id] = (NCBI_id, node_id)
        return converter


def grab_node(aligned):
    # FIXME bug
    toks = aligned.strip().split('\t')
    kegg_id = str(toks[1].split(':')[0])
    # FIXME why a try/except if we're also checking for the key's existence?
    try:
        # FIXME another bug
        if kegg_id in converter:
            return converter[kegg_id][1]
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


def get_samples_from_dir(config):
    """
    Key assumptions:
        + names still end in .fastq and they are gzip compressed (.fastq.gz)
        + files are paired-end
        + names end with _R?_001.fastq.gz or _R?.fastq.gz
        + names start with SAMPLE_ and everything before the first underscore maintains a unique set of sample names
    """
    # groups = {key: set(value) for key, value in groupby(sorted(mylist, key = lambda e: os.path.splitext(e)[0]), key = lambda e: os.path.splitext(e)[0])}
    fastq_dir = config.get("data")
    if not fastq_dir:
        logger.error("'data' dir with FASTQs has not been set; pass --config data=/path")
        sys.exit(1)

    logger.info("Finding samples in %s" % fastq_dir)
    samples = dict()
    seen = set()
    for fname in os.listdir(fastq_dir):
        if not ".fastq" in fname and not ".fq" in fname: continue
        if not "_R1" in fname and not "_r1" in fname: continue
        fq_path = os.path.join(fastq_dir, fname)
        sample_id = fname.partition(".fastq")[0]
        if ".fq" in sample_id:
            sample_id = fname.partition(".fq")[0]
        sample_id = sample_id.replace("_R1", "").replace("_r1", "")
        sample_id = sample_id.replace(".", "_").replace(" ", "_").replace("-", "_")

        # sample ID after rename
        if sample_id in seen:
            # but FASTQ has yet to be added
            # if one sample has a dash and another an underscore, this
            # is a case where we should warn the user that this file
            # is being skipped
            if not fq_path in seen:
                logger.warning("Duplicate sample %s was found after renaming; skipping..." % sample_id)
            continue
        # simple replace of right-most read index designator
        if fname.find("_R1") > fname.find("_r1"):
            r2 = os.path.join(fastq_dir, "_R2".join(fname.rsplit("_R1", 1)))
        else:
            r2 = os.path.join(fastq_dir, "_r2".join(fname.rsplit("_r1", 1)))
        # not paired-end?
        if not os.path.exists(r2):
            logger.error("File [%s] for %s was not found. Exiting." % (r2, sample_id))
            sys.exit(1)
        seen.add(fq_path)
        seen.add(sample_id)
        samples[sample_id] = {"R1": fq_path, "R2": r2}

    if len(samples) == 0:
        logger.error("No samples were found for processing.")
        sys.exit(1)
    logger.info("Found %d samples for processing:\n" % len(samples))
    samples_str = ""
    for k, v in samples.items():
        samples_str += "%s: %s; %s\n" % (k, v["R1"], v["R2"])
    logger.info(samples_str)
    # add sample into config
    config["samples"] = samples


def get_summaries():
    function = ["ec", "ko", "product"]
    taxonomy = ["phylum", "class", "order"]
    # code to generate the possible files
    file_paths = expand("summaries/combined/{function}_{taxonomy}.txt",
        function=function, taxonomy=taxonomy)
    file_paths.extend(expand("summaries/function/{function}.txt",
        function=function))
    file_paths.extend(expand("summaries/taxonomy/{taxonomy}.txt",
        taxonomy=taxonomy))
    return file_paths


get_samples_from_dir(config)
KAIJUDB = get_kaiju_db_dir(config)
CONDAENV = "envs/environment.yml"


rule all:
    input:
        expand("logs/{sample}_{idx}_eestats.txt",
            sample=config["samples"].keys(), idx=["R1", "R2"]),
        expand("quality_control/{sample}_03_{db}.fasta.gz",
            sample=config["samples"].keys(),
            db=config["contaminant_references"].keys()),
        get_summaries(),
        "summary.html"


rule get_raw_fastq_qualities:
    input:
        lambda wildcards: config["samples"][wildcards.sample][wildcards.idx]
    output:
        "logs/{sample}_{idx}_eestats.txt"
    conda:
        CONDAENV
    threads:
        1
    group:
        "sample_group"
    shell:
        """
        vsearch --threads 1 --fastq_eestats {input} --output {output}
        """


rule merge_sequences:
    input:
        unpack(lambda wildcards: config["samples"][wildcards.sample])
    output:
        merged = "quality_control/{sample}_01_merged.fastq.gz",
        R1 = "quality_control/{sample}_01_unmerged_R1.fastq.gz",
        R2 = "quality_control/{sample}_01_unmerged_R2.fastq.gz",
        log = "logs/{sample}_merge_sequences.log"
    params:
        adapters = "" if not config.get("adapters") else "adapter=%s" % config.get("adapters")
    threads:
        config.get("threads", 1)
    resources:
        java_mem = config.get("java_mem", 60)
    conda:
        CONDAENV
    group:
        "sample_group"
    shell:
        """
        bbmerge.sh threads={threads} k=60 extend2=60 iterations=5 \
            ecctadpole=t reassemble=t shave rinse prealloc=t \
            prefilter=10 -Xmx{resources.java_mem}G \
            loose=t qtrim2=t in={input.R1} in2={input.R2} \
            {params.adapters} out={output.merged} \
            outu={output.R1} outu2={output.R2} 2> {output.log}
        """


rule deduplicate_reads:
    input:
        "quality_control/{sample}_01_merged.fastq.gz"
    output:
        fa = "quality_control/{sample}_02_unique.fasta.gz",
        log = "logs/{sample}_unique_reads.log"
    threads:
        config.get("threads", 1)
    resources:
        java_mem = config.get("java_mem", 60)
    conda:
        CONDAENV
    group:
        "sample_group"
    shell:
        """
        clumpify.sh in={input} out={output.fa} dedupe=t threads={threads} \
            -Xmx{resources.java_mem}G 2> {output.log}
        """


rule build_decontamination_db:
    output:
        os.path.join(os.path.dirname(config.get("diamonddb")),
            "ref", "genome", "1", "summary.txt")
    params:
        k = config.get("contaminant_kmer_length", 13),
        refs_in = " ".join(["ref_%s=%s" % (n, fa) for n, fa in config["contaminant_references"].items()]),
        path = os.path.dirname(config.get("diamonddb"))
    resources:
        java_mem = config.get("java_mem", 60)
    threads:
        config.get("threads", 1)
    conda:
        CONDAENV
    shell:
        """
        bbsplit.sh -Xmx{resources.java_mem}G {params.refs_in} \
            threads={threads} k={params.k} local=t path={params.path}
        """


rule run_decontamination:
    input:
        fa = "quality_control/{sample}_02_unique.fasta.gz",
        db = rules.build_decontamination_db.output
    output:
        fa = "quality_control/{sample}_03_clean.fasta.gz",
        contaminants = expand("quality_control/{{sample}}_03_{db}.fasta.gz",
            db=list(config["contaminant_references"].keys())),
        stats = "logs/{sample}_decontamination_by_reference.log",
        log = "logs/{sample}_decontamination.log"
    params:
        path = os.path.join(os.path.dirname(config.get("diamonddb"))),
        prefix = lambda wc, output: "".join(output.contaminants[0].rpartition("_03_")[0:2]),
        maxindel = config.get("contaminant_max_indel", 5),
        minratio = config.get("contaminant_min_ratio", 0.80),
        minhits = config.get("contaminant_minimum_hits", 3),
        ambiguous = config.get("contaminant_ambiguous", "best"),
        k = config.get("contaminant_kmer_length", 12),
    threads:
        config.get("threads", 1)
    resources:
        java_mem = config.get("java_mem", 60)
    conda:
        CONDAENV
    group:
        "sample_group"
    shell:
        """
        bbsplit.sh in={input.fa} outu={output.fa} fastareadlen=300 \
            refstats={output.stats} basename={params.prefix}%.fasta.gz \
            maxindel={params.maxindel} minratio={params.minratio} \
            minhits={params.minhits} ambiguous={params.ambiguous} \
            threads={threads} k={params.k} path={params.path} \
            local=t -Xmx{resources.java_mem}G 2>&1 | tee {output.log}
        """


rule calculate_read_lengths:
    input:
        "quality_control/{file}.fasta.gz"
    output:
        "logs/{file}_readlengths.txt"
    conda:
        CONDAENV
    shell:
        """
        readlength.sh in={input} out={output}
        """


rule run_taxonomic_classification:
    input:
        "quality_control/{sample}_03_clean.fasta.gz"
    output:
        temp("kaiju/{sample}_aln.txt")
    params:
        evalue = config.get("kaiju_evalue", 0.05)
    threads:
        config.get("threads", 1)
    conda:
        CONDAENV
    group:
        "sample_group"
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
    conda:
        CONDAENV
    group:
        "sample_group"
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
    conda:
        CONDAENV
    shell:
        """
        diamond makedb --in {input} --threads {threads} --db {output}
        """


rule run_functional_classification:
    input:
        fq = "quality_control/{sample}_03_clean.fasta.gz",
        db = config["diamonddb"] + ".dmnd"
    output:
        "diamond/{sample}_aln.txt"
    params:
        evalue = config.get("evalue", 0.00001)
    threads:
        config.get("threads", 1)
    conda:
        CONDAENV
    group:
        "sample_group"
    shell:
        """
        diamond blastx --threads {threads} --db {input.db} \
            --evalue {params.evalue} --out {output} --outfmt 6 \
            --query {input.fq} --strand both --unal 1 --top 3 \
            --block-size 4 --index-chunks 1
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
                    # just print the best HSP
                    break


rule download_kegg_hierarchy:
    output:
        os.path.join(os.path.dirname(config.get("diamonddb")), "kegg_hierarchy.json")
    message:
        "Attempting to download KEGG's hierarchy which will write to %s" % os.path.join(os.path.dirname(config.get("diamonddb")), "kegg_hierarchy.json")
    shell:
        """
        curl 'http://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=json' \
            --location --output '{output}'
        """


rule build_functional_table:
    input:
        tables = expand('tables/{sample}_classifications.txt', sample=config["samples"].keys()),
        json = rules.download_kegg_hierarchy.output
    output:
        "summaries/function/{function}.txt"
    params:
        min_id = config.get("min_percent_id", 50),
        min_len = config.get("min_alignment_length", 40)
    threads:
        1
    conda:
        CONDAENV
    shell:
        """
        python scripts/summarize_classifications.py \
            --group-on {wildcards.function} --min-id {params.min_id} \
            --min-len {params.min_len} {input.json} {output} \
            {input.tables}
        """


rule build_tax_table:
    input:
        tables = expand('tables/{sample}_classifications.txt', sample=config["samples"].keys()),
        json = rules.download_kegg_hierarchy.output
    output:
        "summaries/taxonomy/{tax_classification}.txt"
    params:
        min_id = config.get("min_percent_id", 50),
        min_len = config.get("min_alignment_length", 40)
    threads:
        1
    conda:
        CONDAENV
    shell:
        """
        python scripts/summarize_classifications.py \
            --group-on tax_classification --min-len {params.min_len} \
            --min-id {params.min_id} --tax-level {wildcards.tax_classification} \
            {input.json} {output} {input.tables}
        """


rule build_functional_and_tax_table:
    input:
        tables = expand('tables/{sample}_classifications.txt', sample=config["samples"].keys()),
        json = rules.download_kegg_hierarchy.output
    output:
        "summaries/combined/{function}_{tax_classification}.txt"
    params:
        min_id = config.get("min_percent_id", 50),
        min_len = config.get("min_alignment_length", 40)
    threads:
        1
    conda:
        CONDAENV
    shell:
        """
        python scripts/summarize_classifications.py \
            --group-on {wildcards.function} tax_classification \
            --tax-level {wildcards.tax_classification} --min-id {params.min_id} \
            --min-len {params.min_len} {input.json} {output} {input.tables}
        """


rule build_report:
    input:
        classifications = expand("tables/{sample}_classifications.txt", sample=config["samples"].keys()),
        ee_stats = expand("logs/{sample}_R1_eestats.txt", sample=config["samples"].keys()),
        clean_length_logs = expand("logs/{sample}_03_clean_readlengths.txt", sample=config["samples"].keys()),
        unique_length_logs = expand("logs/{sample}_02_unique_readlengths.txt", sample=config["samples"].keys()),
        clean_logs = expand("logs/{sample}_decontamination.log", sample=config["samples"].keys()),
        merge_logs = expand("logs/{sample}_merge_sequences.log", sample=config["samples"].keys()),
        function = "summaries/function/ko.txt",
        taxonomy = "summaries/taxonomy/phylum.txt",
        combined = "summaries/combined/ko_phylum.txt"
    output:
        "summary.html"
    shell:
        """
        python scripts/build_report.py --merge-logs {input.merge_logs} \
            --unique-logs {input.unique_length_logs} \
            --clean-logs {input.clean_length_logs} \
            --summary-tables {input.classifications} \
            --r1-quality-files {input.ee_stats} \
            --html {output} \
            {CONDAENV} {input.function} {input.taxonomy} {input.combined}
        """
