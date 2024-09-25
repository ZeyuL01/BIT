Installation Guide
==================

BIT is built by Rcpp to achieve fast computation. The following packages are required for the usage of BIT, please check the installation before running BIT.

.. note::

   Rcpp, RcppArmadillo, RcppDist, RcppProgress, rtracklayer, data.table, R.utils, ggsci, basicPlotteR


You can directly install BIT from github by typing following command in the R terminal:

.. code-block::

   import(devtools)
   devtools::install_github("ZeyuL01/BIT")
   Downloading GitHub repo ZeyuL01/BIT@HEAD
   ── R CMD build ──────────────────────────────────────────────────────────────────────────────────────────
   ✔  checking for file ‘/private/var/folders/cm/9xs2krg94ygcyk5m9xsr31qc0000gn/T/Rtmp0KrI6F/remotes3c8c5b5ed5f3/ZeyuL01-BIT-2e57689/DESCRIPTION’ ...
   ─  preparing ‘BIT’:
   ✔  checking DESCRIPTION meta-information ...
   ...
   ...
   ...
   ** testing if installed package can be loaded from temporary location
   ** checking absolute paths in shared objects and dynamic libraries
   ** testing if installed package can be loaded from final location
   ** testing if installed package keeps a record of temporary installation path
   * DONE (BIT)

Next, run the `load_chip_data()` function to connect BIT to the reference database. You can download the data from the `Zenodo<https://zenodo.org/records/13732877>`_. Assume the local directory for the folder containing the reference data for genome hg38 is `/../hg38/`. Enter the following code in the R terminal:

.. code-block::

   load_chip_data("/../hg38/", bin_width = 1000, genome="hg38")
   [1] "ChIP-seq data successfully loaded, please run BIT with input to check!"
   NULL


You can reconnect to the updated or modified database at any time.

If you encounter compilation problems, please refer to the following resources:

- **Windows:** `CRAN R Development Guide <https://cran.r-project.org/bin/windows/base/howto-R-devel.html>`_
- **macOS:** `R Compiler Tools for Rcpp on macOS <https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/>`_

Alternatively, you can `submit an issue <https://github.com/ZeyuL01/BIT/issues>`_ on GitHub.

