# PerSeq: Per sequence functional and taxonomic assignments

PerSeq is an annotation workflow implemented in [Snakemake](https://snakemake.readthedocs.io/en/stable/) and is designed to be copied to your analysis directory. Dependencies are defined in
``envs/required.yaml`` and are installed at runtime via [Bioconda](https://bioconda.github.io/)
(``snakemake --use-conda``) or [Biocontainers](https://biocontainers.pro/) (``snakemake --use-singularity``).

To report a bug or suggest changes, please visit the [GitHub repository](https://github.com/brwnj/perseq).

# Quality Control

Sequences are merged (step 01) using `bbmerge.sh` of BBTools. Merging allows
for read extensions up to 300 bp and will quality trim the sequences after
successfully merging R1 and R2. Merged sequences are then deduplicated
(step 02) using `clumpify.sh` of BBTools. Unique, merged sequences are then
mapped against contaminant libraries (step 03) in the form of key:value pairs
in the configuration. To filter against rRNA and PhiX, the specification looks
like:

```
contaminant_references:
    PhiX: resources/phix.fa
    rRNA: resources/rrna.fa.gz
```

Additional references are added to the list like:

```
contaminant_references:
    PhiX: resources/phix.fa
    rRNA: resources/rrna.fa.gz
    Ecoli: /local/path/ecoli.fa
```

Hits to contaminant references are detailed in
`quality_control/<sample>_03_<contaminant key>.fasta.gz`.

# Taxonomic Annotation

Kmer-based taxonomic classification is performed on the merged reads using
Kaiju [5] in greedy mode (``-a greedy -E 0.05``) with the user supplied
reference index. To create a Kaiju reference of NCBI's nr database
containing reference sequences for archaea, bacteria, viruses, fungi, and
microbial eukaryotes execute `makeDB.sh -e` of Kaiju's executable library.

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

# Functional Annotation

The blastx algorithm of DIAMOND is used to align nucleotide sequences to
a KEGG protein reference database consisting of non-redundant, family
level fungal eukaryotes and genus level prokaryotes
(``--strand=both --evalue 0.00001``). Due to KEGG licensing we cannot
distribute this reference database. The highest scoring alignment per
sequence is used for functional annotation.

# Installation

```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install python3 snakemake
```

# Updating References

To update the KO, pathway, and KEGG hierarchy reference data, execute:

```
./scripts/update_public_kegg_references.sh
```
