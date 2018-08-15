"""
Per sequence characterization
"""
import os
import sys

from snakemake.logging import logger


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
    # groups = {key: set(value) for key, value in groupby(sorted(mylist, key = lambda e: os.path.splitext(e)[0]), key = lambda e: os.path.splitext(e)[0])}
    fastq_dir = config.get("data")
    if not fastq_dir:
        logger.error("'data' dir with FASTQs has not been set; pass --config data=/path")
        sys.exit(1)

    raw_fastq_files = list()
    for fq_dir in fastq_dir.split(","):
        if "*" in fastq_dir:
            from glob import glob
            logger.info("Finding samples matching %s" % fq_dir)
            raw_fastq_files.extend(glob(fq_dir))
        else:
            logger.info("Finding samples in %s" % fq_dir)
            raw_fastq_files.extend([os.path.join(fq_dir, i) for i in os.listdir(fq_dir)])

    samples = dict()
    seen = set()
    for fq_path in raw_fastq_files:
        fname = os.path.basename(fq_path)
        fastq_dir = os.path.dirname(fq_path)
        if not ".fastq" in fname and not ".fq" in fname:
            continue
        if not "_R1" in fname and not "_r1" in fname:
            continue

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


def subsample(wildcards):
    # FIXME what if the user puts in "No"
    sample = config.get("subsample", -1)
    if sample < 0:
        return config["samples"][wildcards.sample]
    else:
        return {"R1": "subsampled/{wildcards.sample}_R1.fastq".format(wildcards=wildcards),
            "R2":"subsampled/{wildcards.sample}_R2.fastq".format(wildcards=wildcards)}


get_samples_from_dir(config)
KAIJUDB = get_kaiju_db_dir(config)
CONDAENV = "envs/environment.yml"


rule all:
    input:
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


# TODO write this rule so that it operates on a single index
rule subsample:
    input:
        unpack(lambda wildcards: config["samples"][wildcards.sample])
    output:
        R1 = "subsampled/{sample}_R1.fastq.gz",
        R2 = "subsampled/{sample}_R2.fastq.gz"
    params:
        subsample = config.get("subsample", 60000)
    threads:
        1
    conda:
        CONDAENV
    group:
        "sample_group"
    shell:
        """
          seqtk sample -s100 {input.R1} {params.subsample} | gzip > {output.R1}
          seqtk sample -s100 {input.R2} {params.subsample} | gzip > {output.R2}
        """


rule merge_sequences:
    input:
        unpack(subsample)
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
    group:
        "sample_group"
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
            -r superkingdom,phylum,class,order,family,genus,species
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
        kaiju = "kaiju/{sample}_aln_names.txt",
        blastx = "diamond/{sample}_aln.txt",
        gene2ko = config["gene2ko"],
        kolist = config["ko_list"],
        code2id = config["genome_list"],
        names = os.path.join(config["kaijudb"], "names.dmp"),
        nodes = os.path.join(config["kaijudb"], "nodes.dmp")
    output:
        "tables/{sample}_classifications.txt"
    params:
        lca_threshold = config.get("lca_threshold", 1.0)
    shell:
        """
        python scripts/make_classification_table.py \
            --lca-threshold {params.lca_threshold} {input.kaiju} {input.blastx} \
            {input.gene2ko} {input.code2id} {input.kolist} {input.names} \
            {input.nodes} {output}
        """


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
        tables = expand("tables/{sample}_classifications.txt", sample=config["samples"].keys()),
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
        tables = expand("tables/{sample}_classifications.txt", sample=config["samples"].keys()),
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
            --group-on kaiju_classification --min-len {params.min_len} \
            --min-id {params.min_id} --tax-level {wildcards.tax_classification} \
            {input.json} {output} {input.tables}
        """


rule build_functional_and_tax_table:
    input:
        tables = expand("tables/{sample}_classifications.txt", sample=config["samples"].keys()),
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
            --group-on {wildcards.function} kaiju_classification \
            --tax-level {wildcards.tax_classification} --min-id {params.min_id} \
            --min-len {params.min_len} {input.json} {output} {input.tables}
        """


rule build_krona_ec_input:
    input:
        ec_file = "summaries/function/ec.txt",
        ec_converter = config["ec_converter"],
        ec_dat_file = config["enzyme_dat_file"]
    output:
        expand("krona_plots/{sample}_ec.txt", sample=config["samples"].keys())
    threads:
        1
    conda:
        CONDAENV
    shell:
        """
        python scripts/krona_plots.py --ec-file {input.ec_converter} \
            --dat-file {input.ec_dat_file} --ec-file-from-summaries {input.ec_file} \
            krona_plots
        """


rule build_krona_taxonomy_input:
    input:
        tax_file="summaries/taxonomy/order.txt"
    output:
        expand("krona_plots/{sample}_tax.txt", sample=config["samples"].keys())
    threads:
        1
    conda:
        CONDAENV
    shell:
        """
        python scripts/krona_plots.py {input.tax_file} krona_plots
        """


rule build_krona_plots:
    input:
        tax = expand("krona_plots/{sample}_tax.txt", sample=config["samples"].keys()),
        ec = expand("krona_plots/{sample}_ec.txt", sample=config["samples"].keys())
    output:
        tax = "krona_plots/tax.krona.html",
        ec = "krona_plots/ec.krona.html"
    threads:
        1
    conda:
        CONDAENV
    shell:
        """
        ktImportText {input.tax} -o {output.tax}
        ktImportText {input.ec} -o {output.ec}
        """


# FIXME I think these need to be passed as patterns like 'logs/*_R1_eestats.txt'
# and resolved in build_report.py using `glob()` [from glob import glob]
rule build_report:
    input:
        classifications = expand("tables/{sample}_classifications.txt", sample=config["samples"].keys()),
        ee_stats = expand("logs/{sample}_{idx}_eestats.txt", sample=config["samples"].keys(), idx=["R1", "R2"]),
        clean_length_logs = expand("logs/{sample}_03_clean_readlengths.txt", sample=config["samples"].keys()),
        unique_length_logs = expand("logs/{sample}_02_unique_readlengths.txt", sample=config["samples"].keys()),
        clean_logs = expand("logs/{sample}_decontamination.log", sample=config["samples"].keys()),
        merge_logs = expand("logs/{sample}_merge_sequences.log", sample=config["samples"].keys()),
        function = "summaries/function/ko.txt",
        taxonomy = "summaries/taxonomy/phylum.txt",
        combined = "summaries/combined/ko_phylum.txt",
        krona_tax = "krona_plots/tax.krona.html",
        krona_ec = "krona_plots/ec.krona.html"
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
            {CONDAENV} {input.function} {input.taxonomy} {input.combined} \
            {input.krona_tax} {input.krona_ec}
        """
