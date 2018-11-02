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


def readfx(fp):
    """
    Generator for FAST{A,Q}. See https://github.com/lh3/readfq.
    """
    last = None
    while True:
        if not last:
            for l in fp:
                if l[0] in ">@":
                    last = l[:-1]
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:
            if l[0] in "@+>":
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != "+":
            yield name, "".join(seqs), None
            if not last:
                break
        else:
            seq, leng, seqs = "".join(seqs), 0, []
            for l in fp:
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):
                    last = None
                    yield name, seq, "".join(seqs)
                    break
            if last:
                yield name, seq, None
                break


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
    logger.info("Found %d samples for processing" % len(samples))
    # samples_str = ""
    # for k, v in samples.items():
    #     samples_str += "%s: %s; %s\n" % (k, v["R1"], v["R2"])
    # logger.info(samples_str)
    # add sample into config
    config["samples"] = samples


def get_summaries():
    function = ["ec", "product"]
    taxonomy = ["phylum", "class", "order"]
    # code to generate the possible files
    file_paths = expand("summaries/combined/{function}_{taxonomy}.txt",
        function=function, taxonomy=taxonomy)
    file_paths.extend(expand("summaries/function/{function}.txt",
        function=function))
    file_paths.extend(expand("summaries/taxonomy/{taxonomy}.txt",
        taxonomy=taxonomy))
    return file_paths


def get_merge_input(wildcards):
    subsample = config.get("subsample", -1)
    if isinstance(subsample, int):
        if subsample < 1:
            # logger.info("No subsampling performed.")
            files = config["samples"][wildcards.sample]
        else:
            files = {
                "R1": "subsampled/{wc.sample}_R1.fastq.gz".format(wc=wildcards),
                "R2": "subsampled/{wc.sample}_R2.fastq.gz".format(wc=wildcards),
            }
    else:
        logger.error(f"Invalid argument provided to subsample: {subsample}")
        sys.exit(1)
    return files


def get_hmm(wildcards):
    if wildcards.hmm == "HAMAP":
        hmm = config["hamap_hmm"]
    elif wildcards.hmm == "dbCAN":
        hmm = config["dbcan_hmm"]
    elif wildcards.hmm == "TIGRFAMs":
        hmm = config["tigrfams_hmm"]
    else:
        logger.error("Unsure which HMM is currently selected.")
        sys.exit(status=1)
    return dict(
                hmm=hmm,
                h3f="%s.h3f" % hmm,
                h3i="%s.h3i" % hmm,
                h3m="%s.h3m" % hmm,
                h3p="%s.h3p" % hmm
    )


get_samples_from_dir(config)
KAIJUDB = get_kaiju_db_dir(config)
CONDAENV = "envs/environment.yml"


localrules: all
rule all:
    input:
        "SAMPLES.txt",
        expand("quality_control/{sample}_03_{db}.fasta.gz",
            sample=config["samples"].keys(),
            db=config["contaminant_references"].keys()),
        get_summaries(),
        # expand("translated/{sample}_kaiju.txt",sample=config["samples"].keys()),
        # dynamic("fasta_chunks/{sample}_03_clean_{chunk}.fasta"),
        # expand("full_hmmscan/{sample}_hmmscan.txt",sample=config["samples"].keys()),
        # expand("full_hmmscan/{sample}_hmmscan.txt",sample=config["samples"].keys()),

        "summary.html"


localrules: print_samples
rule print_samples:
    output:
        "SAMPLES.txt"
    params:
        samples = config["samples"]
    run:
        with open(output[0], "w") as fh:
            for k, v in params.samples.items():
                print("%s: %s; %s" % (k, v["R1"], v["R2"]), file=fh)


rule get_raw_fastq_qualities:
    input:
        lambda wildcards: config["samples"][wildcards.sample][wildcards.idx]
    output:
        "logs/{sample}_{idx}_eestats.txt"
    conda:
        CONDAENV
    threads:
        1
    shell:
        # TODO seqtk fqchk
        """
        vsearch --threads 1 --fastq_eestats {input} --output {output}
        """


rule subsample_sequences:
    input:
        lambda wildcards: config["samples"][wildcards.sample][wildcards.idx]
    output:
        "subsampled/{sample}_{idx}.fastq.gz"
    params:
        subsample = config.get("subsample", 60000)
    threads:
        1
    conda:
        CONDAENV
    shell:
        """
        seqtk sample -s100 {input} {params.subsample} | gzip > {output}
        """


rule merge_sequences:
    input:
        unpack(get_merge_input)
    output:
        merged = temp("quality_control/{sample}_01_merged.fastq.gz"),
        R1 = temp("quality_control/{sample}_01_unmerged_R1.fastq.gz"),
        R2 = temp("quality_control/{sample}_01_unmerged_R2.fastq.gz"),
        log = "logs/{sample}_merge_sequences.log"
    params:
        adapters = "" if not config.get("adapters") else "adapter=%s" % config.get("adapters")
    threads:
        config.get("threads", 1)
    resources:
        java_mem = config.get("java_mem", 60)
    conda:
        CONDAENV
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
        fa = temp("quality_control/{sample}_02_unique.fasta.gz"),
        log = "logs/{sample}_unique_reads.log"
    threads:
        config.get("threads", 1)
    resources:
        java_mem = config.get("java_mem", 60)
    conda:
        CONDAENV
    shell:
        """
        clumpify.sh in={input} out={output.fa} dedupe=t threads={threads} \
            -Xmx{resources.java_mem}G 2> {output.log}
        """


rule build_decontamination_db:
    output:
        os.path.join(os.path.dirname(config["hamap_hmm"]),
            "ref", "genome", "1", "summary.txt")
    params:
        k = config.get("contaminant_kmer_length", 13),
        refs_in = " ".join(["ref_%s=%s" % (n, fa) for n, fa in config["contaminant_references"].items()]),
        path = os.path.dirname(config["hamap_hmm"])
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
        path = os.path.join(os.path.dirname(config["hamap_hmm"])),
        prefix = lambda wc, output: "".join(output.contaminants[0].rpartition("_03_")[0:2]),
        maxindel = config.get("contaminant_max_indel", 5),
        minratio = config.get("contaminant_min_ratio", 0.95),
        minhits = config.get("contaminant_minimum_hits", 3),
        ambiguous = config.get("contaminant_ambiguous", "best"),
        k = config.get("contaminant_kmer_length", 12),
    threads:
        config.get("threads", 1)
    resources:
        java_mem = config.get("java_mem", 60)
    conda:
        CONDAENV
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


rule run_prodigal:
    input:
        "quality_control/{sample}_03_clean.fasta.gz"
    output:
        "gene_catalog/prodigal/{sample}.faa"
    params:
        null = os.devnull
    conda:
        CONDAENV
    shell:
        """
        gunzip -c {input} | prodigal -q -p meta -a {output} -o {params.null}
        """


localrules: aggregate_all_genes
rule aggregate_all_genes:
    input:
        faa = expand(
            "gene_catalog/prodigal/{sample}.faa",
            sample=config["samples"].keys()
        )
    output:
        faa = "gene_catalog/all_genes.faa"
    run:
        name_index = 1
        # renames to short, numeric ID
        with open(output.faa, "w") as out_faa:
            for f in input.faa:
                with open(f) as fh:
                    for name, seq, _ in readfx(fh):
                        # per sequence, choose first only
                        if name.endswith("_2"):
                            continue
                        print(">%d" % name_index, seq.replace("*", ""), sep="\n", file=out_faa)
                        name_index += 1


rule build_gene_db:
    input:
        faa = "gene_catalog/all_genes.faa"
    output:
        db = temp("gene_catalog/preclustered_genes_db"),
        extras = temp([
            "gene_catalog/preclustered_genes_db.dbtype",
            "gene_catalog/preclustered_genes_db.index",
            "gene_catalog/preclustered_genes_db.lookup",
            "gene_catalog/preclustered_genes_db_h",
            "gene_catalog/preclustered_genes_db_h.index"
        ]),
        clustered_db = temp("gene_catalog/clustered_genes_db"),
        cluster_extras = temp("gene_catalog/clustered_genes_db.index"),
        representative_seqs = temp("gene_catalog/representative_seqs"),
        representative_seqs_extras = temp([
            "gene_catalog/representative_seqs.dbtype",
            "gene_catalog/representative_seqs.index"
        ]),
        fasta = "gene_catalog/clustered_genes.faa"
    params:
        cluster_id = config.get("clustering_threshold", 0.90),
        tmpdir = lambda wildcards, input: os.path.join(os.path.dirname(input.faa), "TMP")
    conda:
        CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """
        mmseqs createdb {input.faa} {output.db}
        mmseqs linclust --threads {threads} -v 1 --min-seq-id 0.90 \
            {output.db} {output.clustered_db} {params.tmpdir}
        mmseqs result2repseq {output.db} {output.clustered_db} \
            {output.representative_seqs}
        mmseqs result2flat {output.db} {output.db} \
            {output.representative_seqs} {output.fasta} --use-fasta-header
        rm -r {params.tmpdir}
        """


rule index_representative_sequences:
    input:
        "gene_catalog/clustered_genes.faa"
    output:
        "gene_catalog/clustered_genes.dmnd"
    conda:
        CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """
        diamond makedb --threads {threads} --db {output} --in {input}
        """


rule index_hmm_libraries:
    input:
        hamap = config["hamap_hmm"],
        tigrfams = config["tigrfams_hmm"],
        dbcan = config["dbcan_hmm"]
    output:
        expand(
            "{hmm}.{exts}",
            hmm=[config[i] for i in ["hamap_hmm", "tigrfams_hmm", "dbcan_hmm"]],
            exts=["h3f", "h3i", "h3m", "h3p"]
        )
    conda:
        CONDAENV
    shell:
        """
        hmmpress -f {input.hamap}
        hmmpress -f {input.tigrfams}
        hmmpress -f {input.dbcan}
        """


localrules: split_fasta
rule split_fasta:
    input:
        fasta = "gene_catalog/clustered_genes.faa"
    output:
        temp((dynamic("gene_catalog/tmp/clustered_genes_{chunk}.faa")))
    params:
        # consider a smaller chunk size
        chunk_size = config.get("chunk_size", 1000000)
    run:
        fasta_chunk = None
        with open(input.fasta) as fasta_fh:
            for lineno, (name, seq, _) in enumerate(readfx(fasta_fh)):
                if lineno % params.chunk_size == 0:
                    if fasta_chunk:
                        fasta_chunk.close()
                    fasta_chunk_filename = f"gene_catalog/tmp/clustered_genes_{lineno + params.chunk_size}.faa"
                    fasta_chunk = open(fasta_chunk_filename, "w")
                fasta_chunk.write(f">{name}\n")
                fasta_chunk.write(f"{seq}\n")
            if fasta_chunk:
                fasta_chunk.close()


rule run_hmmsearch:
    # output is sorted by the target HMM library
    input:
        unpack(get_hmm),
        faa = "gene_catalog/tmp/clustered_genes_{chunk}.faa"
    output:
        hits = temp("gene_catalog/{hmm}/unsorted_alignments_{chunk}.txt")
    params:
        evalue = config.get("evalue", 0.05),
        null = os.devnull
    conda:
        CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """
        hmmsearch --noali --notextw --acc --cpu {threads} -E {params.evalue} \
            --domtblout {output.hits} -o {params.null} {input.hmm} {input.faa}
        """


rule sort_hmm_hits:
    # Remove the header, remove spacing, replace spaces with tabs, sort by
    # query then score. Best hit will be first of group.
    input:
        hits = dynamic("gene_catalog/{hmm}/unsorted_alignments_{chunk}.txt")
    output:
        # column[4] contains annotation data
        hits = "gene_catalog/{hmm}/alignments.tsv"
    conda:
        CONDAENV
    shell:
        """
        cat {input.hits} | \
            grep -v '^#' | \
            tr -s ' ' | \
            tr ' ' '\t' | \
            sort --buffer-size=50% -k1,1 -k8,8nr > {output.hits}
        """


rule align_sequences_to_clusters:
    # reports best hit only per sequence
    input:
        dmnd = "gene_catalog/clustered_genes.dmnd",
        faa = "gene_catalog/prodigal/{sample}.faa"
    output:
        "gene_catalog/diamond/{sample}.tsv"
    params:
        sequence_id = config.get("sequence_threshold", 0.75),
        block_size = config.get("block_size", 8)
    conda:
        CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """
        diamond blastp --threads {threads} --id {params.sequence_id} \
            --max-target-seqs 1 --db {input.dmnd} --out {output} --outfmt 6 \
            --query {input.faa} --block-size {params.block_size}
        """


rule run_taxonomic_classification:
    input:
        faa = "gene_catalog/clustered_genes.faa",
        fmi = f"{KAIJUDB}/kaiju_db.fmi",
        nodes = f"{KAIJUDB}/nodes.dmp"
    output:
        temp("gene_catalog/kaiju/alignments_no_names.txt")
    params:
        evalue = config.get("kaiju_evalue", 0.05)
    threads:
        config.get("threads", 1)
    conda:
        CONDAENV
    shell:
        """
        kaiju -t {input.nodes} -f {input.fmi} -p -i {input} -o {output} \
            -z {threads} -a greedy -x -v -E {params.evalue}
        """


rule add_full_taxonomy:
    input:
        alignments = "gene_catalog/kaiju/alignments_no_names.txt",
        nodes = f"{KAIJUDB}/nodes.dmp",
        names = f"{KAIJUDB}/names.dmp"
    output:
        "gene_catalog/kaiju/alignments.tsv"
    conda:
        CONDAENV
    shell:
        """
        addTaxonNames -t {input.nodes} -n {input.names} -i {input} \
            -o {output} -r superkingdom,phylum,class,order,family,genus,species
        """


rule combine_sample_output:
    input:
        kaiju = "gene_catalog/kaiju/alignments.tsv",
        # row[4].split("~~~") -> ec, gene, product.replace("^", " "), HMM ID
        hamap = "gene_catalog/HAMAP/alignments.tsv",
        # row[4].split("~~~") -> ec, enzyme class, enzyme class subfamily, HMM ID
        dbcan = "gene_catalog/dbCAN/alignments.tsv",
        # row[4].split("~~~") -> ec, gene, product.replace("^", " "), HMM ID
        tigrfams = "gene_catalog/TIGRFAMs/alignments.tsv",
        hsps = expand("gene_catalog/diamond/{sample}.tsv", sample=config["samples"].keys())
    output:
        "gene_catalog/annotations.txt"
    conda:
        CONDAENV
    shell:
        """
        python scripts/make_classification_table.py --output {output} \
            {input.kaiju} {input.hamap} {input.dbcan} {input.tigrfams} \
            'gene_catalog/diamond/*.tsv'
        """


rule build_functional_table:
    input:
        "gene_catalog/annotations.txt"
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
            --min-len {params.min_len} {output} {input}
        """


rule build_tax_table:
    input:
        "gene_catalog/annotations.txt"
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
            {output} {input}
        """


rule build_functional_and_tax_table:
    input:
        "gene_catalog/annotations.txt"
    output:
        "summaries/combined/{function}_{tax_classification}.txt"
    params:
        min_id = config.get("min_percent_id", 50),
        min_len = config.get("min_alignment_length", 40)
    conda:
        CONDAENV
    shell:
        """
        python scripts/summarize_classifications.py \
            --group-on {wildcards.function} kaiju_classification \
            --tax-level {wildcards.tax_classification} --min-id {params.min_id} \
            --min-len {params.min_len} {output} {input}
        """


rule build_krona_ec_input:
    input:
        ec_file = "summaries/function/ec.txt",
        ec_converter = config["enzyme_classes"],
        ec_dat_file = config["enzyme_nomenclature"]
    output:
        expand("krona_plots/{sample}_ec.txt", sample=config["samples"].keys())
    conda:
        CONDAENV
    shell:
        """
        python scripts/build_krona.py --ec-file {input.ec_converter} \
            --dat-file {input.ec_dat_file} --ec-file-from-summaries \
            {input.ec_file} krona_plots
        """


rule build_krona_taxonomy_input:
    input:
        "summaries/taxonomy/order.txt"
    output:
        expand("krona_plots/{sample}_tax.txt", sample=config["samples"].keys())
    conda:
        CONDAENV
    shell:
        """
        python scripts/build_krona.py --tax-file {input} krona_plots
        """


rule build_krona_plots:
    input:
        tax = expand("krona_plots/{sample}_tax.txt", sample=config["samples"].keys()),
        ec = expand("krona_plots/{sample}_ec.txt", sample=config["samples"].keys())
    output:
        tax = "krona_plots/tax.krona.html",
        ec = "krona_plots/ec.krona.html"
    conda:
        CONDAENV
    shell:
        """
        ktImportText {input.tax} -o {output.tax}
        ktImportText {input.ec} -o {output.ec}
        """


rule zip_attachments:
    input:
        function = "summaries/function/ko.txt",
        taxonomy = "summaries/taxonomy/order.txt",
        combined = "summaries/combined/ko_phylum.txt",
        krona_tax = "krona_plots/tax.krona.html",
        krona_ec = "krona_plots/ec.krona.html"
    output:
        temp("perseq_downloads.zip")
    shell:
        """
        zip {output} {input.function} {input.taxonomy} {input.combined}
        """


rule build_report:
    input:
        annotations = "gene_catalog/annotations.txt",
        ee_stats = expand("logs/{sample}_{idx}_eestats.txt", sample=config["samples"].keys(), idx=["R1", "R2"]),
        clean_length_logs = expand("logs/{sample}_03_clean_readlengths.txt", sample=config["samples"].keys()),
        unique_length_logs = expand("logs/{sample}_02_unique_readlengths.txt", sample=config["samples"].keys()),
        clean_logs = expand("logs/{sample}_decontamination.log", sample=config["samples"].keys()),
        merge_logs = expand("logs/{sample}_merge_sequences.log", sample=config["samples"].keys()),
        taxonomy = "summaries/taxonomy/order.txt",
        zipped_files = "perseq_downloads.zip"
    output:
        "summary.html"
    conda:
        CONDAENV
    shell:
        """
        python scripts/build_report.py \
            --clean-logs 'logs/*_03_clean_readlengths.txt' \
            --unique-logs 'logs/*_02_unique_readlengths.txt' \
            --merge-logs 'logs/*_merge_sequences.log' \
            --summary-tables {input.annotations} \
            --r1-quality-files 'logs/*_R1_eestats.txt' \
            --html {output} \
            {CONDAENV} {input.taxonomy} {input.zipped_files}
        """
