Per sequence functional and taxonomic assignments
=================================================

.. raw:: html

    <div style="clear: both"></div>
    <div class="container-fluid">
    <div class="row">


.. raw:: html

    <div class="col-sm-3">
    <h3>Documentation</h3>


.. toctree::
   :maxdepth: 1

   introduction
   installing
   usage
   example


.. raw:: html

    </div>
    <div class="col-md-9">
    <br>


PerSeq is an annotation workflow implemented in Snakemake_ and is designed to
be copied to your analysis directory. Dependencies are defined in
``envs/required.yaml`` and are installed at runtime via Bioconda_
(``snakemake --use-conda``) or Biocontainers_ (``snakemake --use-singularity``).

To report a bug or suggest changes, please visit the `GitHub repository`_.


.. raw:: html

    </div>
    </div>
    </div>


.. _Snakemake: https://snakemake.readthedocs.io/en/stable/
.. _Bioconda: https://bioconda.github.io/
.. _Biocontainers: https://biocontainers.pro/
.. _GitHub repository: https://github.com/brwnj/perseq
