"""
per sequence classifications
"""

import os
import sys
from snakemake.logging import logger
from snakemake.utils import report


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
CONDAENV = "environment.yml"


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
    conda:
        CONDAENV
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
    conda:
        CONDAENV
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
        java_mem = config.get("java_mem", 60)
    conda:
        CONDAENV
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
        java_mem = config.get("java_mem", 60)
    conda:
        CONDAENV
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
    conda:
        CONDAENV
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
        fq = "quality_control/{sample}_03_merged.fastq.gz",
        db = config["diamonddb"] + ".dmnd"
    output:
        "diamond/{sample}_aln.txt"
    threads:
        config.get("threads", 1)
    conda:
        CONDAENV
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
        with open(diamond_hits) as ifh, open(outtable, "w") as ofh:
            print("read_id", "aa_percent_id", "aa_alignment_length", "ko",
                "product", "ec", "kaiju_alignment_length",
                "kaiju_classification", "lca_classification",
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
                        # TODO
                        lca_classification = "k__?;p__?;c__?;o__?;f__?;g__?;s__?"
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
                        lca_classification = your_lca_method()
                    # get taxonomy; this should never cause a KeyError
                    tax_alignment_length, tax_classification = tax_classifications[read_id]
                    print(read_id, aa_percent_id, aa_alignment_length, ko, product,
                        ec, tax_alignment_length, tax_classification, lca_classification, sep="\t", file=ofh)
                    # just print the best HSP
                    break
