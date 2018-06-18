# PerSeq: Per sequence functional and taxonomic assignments

PerSeq is an annotation workflow implemented in [Snakemake](https://snakemake.readthedocs.io/en/stable/) and is designed to be copied to your analysis directory. Dependencies are defined in
``envs/required.yaml`` and are installed at runtime via [Bioconda](https://bioconda.github.io/)
(``snakemake --use-conda``) or [Biocontainers](https://biocontainers.pro/) (``snakemake --use-singularity``).

To report a bug or suggest changes, please visit the [GitHub repository](https://github.com/brwnj/perseq).

An introduction
===============

+ Outline the methods included and what is performed at each step
+ Workflow and replicability
+ Version control and sharing experimental archive

+ deduplication
+ decontamination
+ merge seqs
+ kaiju
+ diamond
+ report

# Installation

```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install python3 snakemake
```

+ Install prerequisites
+ download conda environments or containers
+ download reference data
+ basic command to start the protocol

Using the workflow
==================

- htseq=0.9.1

```
mkdir kaiju_db
cd kaiju_db
makeDB.sh -e -t {threads}
```

In the config, set `kaiju_db` as the absolute path to this directory:

```
kaijudb: /full/path/kaiju_db
```

Expected files after building the database are:

+ names.dmp
+ nodes.dmp
+ kaiju_db.fmi


+ Cloning the repo
+ Placing data in the right spot and expectations of format
+ Executing snakemake with basic commands; link to snakemake docs
+ Expected output
