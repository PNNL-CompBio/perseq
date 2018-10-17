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


get_samples_from_dir(config)
KAIJUDB = get_kaiju_db_dir(config)
CONDAENV = "envs/environment.yml"


localrules: all
rule all:
    input:
        expand("quality_control/{sample}_03_{db}.fasta.gz",
            sample=config["samples"].keys(),
            db=config["contaminant_references"].keys()),
        get_summaries(),
        expand("translated/{sample}_kaiju.txt",sample=config["samples"].keys()),
        # dynamic("fasta_chunks/{sample}_03_clean_{chunk}.fasta"),
        # expand("full_hmmscan/{sample}_hmmscan.txt",sample=config["samples"].keys()),
        expand("full_hmmscan/{sample}_hmmscan.txt",sample=config["samples"].keys()),
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
    group:
        "sample_group"
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
        fa = temp("quality_control/{sample}_02_unique.fasta.gz"),
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


localrules: parse_kaiju_for_prot
rule parse_kaiju_for_prot:
    input:
        "kaiju/{sample}_aln_names.txt"
    output:
        "kaiju/{sample}.faa"
    run:
        with open(input, "r") as in_file, open(output, "w") as out_file:
            for line in in_file:
                if line.startswith("C"):
                    toks = line.split("\t")
                    name = toks[1]
                    # multiple sequences, same ID, best hit is chosen later
                    for seq in toks[6].split(sep=","):
                        if seq:
                            print(">%s" % name, file=out_file)
                            print(seq, file=out_file)


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


rule split_fasta:
    input:
        #"quality_control/{sample}_03_clean.fasta.gz"
        "translated/{sample}_kaiju.txt"
    output:
        #dynamic("fasta_chunks/{sample}_03_clean_{chunk}.fasta")
        dynamic("fasta_chunk/{sample}_kaiju_{chunk}.fasta")
    # conda:
    #     CONDAENV
    params:
        # 20 MB chunks
        chunk_size=config.get("chunk_size",1048576),
        sample='{sample}'
    run:

        with open(input, 'r') as file:
            n=0
            #file_out = "sample_03_clean_{INDEX}.fasta"
            while True:
                #file_out={params.sample}
                # continue to write to one file
                #with open("fasta_chunks/"+{params.sample}+"_03_clean_"+str(n)+".fasta", 'w') as out:
                with open("fasta_chunks/"+{params.sample}+"_kaiju"+str(n)+".fasta", 'w') as out:
                    current_bytes = 0
                    for_loop = False
                    for x, (name, seq, other) in enumerate(readfx(file)):
                        if current_bytes > {params.chunk_size}:
                            for_loop = True
                            n+=1
                            break
                        current_bytes += len(name) + len(seq)
                        out.write(">" + name + '\n')
                        out.write(seq + '\n')
                    if not for_loop:
                        break
#
#
# rule translate_nuc_prot:
#     input:
#         #dynamic(expand("fasta_chunks/{sample}_03_clean_{{chunk}}.fasta", sample=config["samples"].keys()))
#         "fasta_chunks/{sample}_03_clean_{chunk}.fasta"
#     output:
#         tranlated = temp("translated/{sample}_03_clean_{chunk}.faa"),
#         log = "logs/{sample}_prodigal_{chunk}.log"
#     conda:
#         CONDAENV
#     shell:
#         """
#         prodigal -i {input} -o {output.log} -a {output.translated} -q
#         """
#
#
rule hmmscan:
    input:
        #"translated/{sample}_03_clean_{chunk}.faa"
        "fasta_chunk/{sample}_kaiju_{chunk}.fasta"
    output:
        hmm_out =temp("hmm_scan/{sample}_hmmscan_{chunk}.txt"),
        log = "logs/{sample}_hmmscan_{chunk}_out.log"
    threads:
        config.get("threads", 1)
    params:
        hmm_db = config.get("hmm_db"),
        evalue = config.get("evalue",10)
    conda:
        CONDAENV
    shell:
        """
        hmmscan --noali --notextw -o {output.log} --acc -E {params.evalue} --cpu {threads} --domtblout {output.hmm_out}
        """
#
rule merge_chunks:
    input:
        dynamic(expand("hmm_scan/{sample}_hmmscan_{{chunk}}.txt",sample=config["samples"].keys()))
    output:
        "full_hmmscan/{sample}_hmmscan.txt"
    shell:
        "cat {input} > {output}"
#
#
# rule parse_hmmscan_output:
#     input:
#         "full_hmmscan/{sample}_hmmscan.txt"
#     output:
#         "hmm/{sample}_hmmscan.txt"
#     run:
#         with open(input,"r") as file, open(output,"w") as out:
#             print("seq_name","evalue","ec","gene","desc",sep='\t',file=out)
#             for line in file:
#                 if line.startswith('#'):
#                     continue
#                 toks= line.strip().split()
#                 seq_name=toks[3]
#                 try:
#                     evalue = toks[6]
#                 except:
#                     print(toks)
#                     break
#                 info = toks[22].split('~~~')
#                 ec = info[0]
#                 gene = info[1]
#                 desc = info[2]
#                 print(seq_name,evalue,ec,gene,desc,sep='\t',file=out)


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
    conda:
        CONDAENV
    shell:
        """
        python scripts/make_classification_table.py \
            --lca-threshold {params.lca_threshold} {input.kaiju} {input.blastx} \
            {input.gene2ko} {input.code2id} {input.kolist} {input.names} \
            {input.nodes} {output}
        """


rule build_functional_table:
    input:
        tables = expand("tables/{sample}_classifications.txt", sample=config["samples"].keys()),
        json = config["ko_hierarchy"]
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
        json = config["ko_hierarchy"]
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
        json = config["ko_hierarchy"]
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
        python scripts/build_krona.py --ec-file {input.ec_converter} \
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
        python scripts/build_krona.py --tax-file {input.tax_file} krona_plots
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


rule build_report:
    input:
        classifications = expand("tables/{sample}_classifications.txt", sample=config["samples"].keys()),
        ee_stats = expand("logs/{sample}_{idx}_eestats.txt", sample=config["samples"].keys(), idx=["R1", "R2"]),
        clean_length_logs = expand("logs/{sample}_03_clean_readlengths.txt", sample=config["samples"].keys()),
        unique_length_logs = expand("logs/{sample}_02_unique_readlengths.txt", sample=config["samples"].keys()),
        clean_logs = expand("logs/{sample}_decontamination.log", sample=config["samples"].keys()),
        merge_logs = expand("logs/{sample}_merge_sequences.log", sample=config["samples"].keys()),
        function = "summaries/function/ko.txt",
        taxonomy = "summaries/taxonomy/order.txt",
        combined = "summaries/combined/ko_phylum.txt",
        krona_tax = "krona_plots/tax.krona.html",
        krona_ec = "krona_plots/ec.krona.html"
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
            --summary-tables 'tables/*_classifications.txt' \
            --r1-quality-files 'logs/*_R1_eestats.txt' \
            --html {output} \
            {CONDAENV} {input.function} {input.taxonomy} {input.combined} \
            {input.krona_tax} {input.krona_ec}
        """
