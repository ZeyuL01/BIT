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

Next you have to run load_chip_data() function to connect BIT to the reference database. You can download the data from `Box folder <https://smu.box.com/s/dswrvsz4chh7ygkjpwdq3lex2gvrz2gi>`_, assume the location of the reference 1000 bps database folder is /../chip_seq_data_1000/. Just type the following code in the R terminal:

.. code-block::

   load_chip_data("/../chip_seq_data_1000/", bin_width = 1000)
   [1] "ChIP-seq data successfully loaded, please run BIT with input to check!"
   NULL


You can reconnect the updated or modified database anytime later.

If you have compilation problem, please refer to the following posts.

windows: https://cran.r-project.org/bin/windows/base/howto-R-devel.html

macos: https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/

Or you can `submit <https://github.com/ZeyuL01/BIT/issues>`_ the issue on GitHub.