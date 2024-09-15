API References
==============

.. function:: BIT (file, output_path, show=TRUE, plot.bar=TRUE, format=NULL, N = 5000 ,bin_width = 1000, burnin=NULL, genome= c("hg38","mm10"))

	*Main interface to run BIT method, please set the input file path, input file format, number of iterations and bin width*.

	 **file**: file path to the user-input file.

	 **output_path**: absoluate or relative directory to save the output.

   **show**: `TRUE` / `FALSE`. Whether to display the ranking table. Default: `TRUE`.

   **plot.bar**: `TRUE` / `FALSE`. Whether to plot the top 10 TRs BIT scores in a horizontal bar plot. Default: `TRUE`.

   **format**: One of `"bed"`, `"narrowPeak"`, `"broadPeak"`, `"bigNarrowPeak"`, `"csv"`, or `NULL`. Specifies the format of the input file. Default: `NULL`. If set to `NULL`, **BIT** will automatically determine the file format based on its extension.

   **N**: Integer. The number of iterations in the Gibbs sampler. Default: `5000`.

   **bin_width**: Integer. The width of the bin used to divide the chromatin into non-overlapping bins. Default: `1000`. Only change this if you compile a different reference database.

   **burnin**: `NULL` or an integer. Specifies the burn-in period when deriving the Bayesian inference for TR-level parameters. Default: `NULL`, which will use `N/2`.

   **genome**: `hg38` or `mm10` for TR ChIP-seq data collected from different genome.



.. function:: BIT_compare (file1, file2, output_path, show=TRUE, plot_scatter=TRUE, format=c(NULL,NULL), N = 5000, bin_width = 1000, burnin=NULL, genome=c("hg38","mm10"))

	 *compare BIT identifid TRs for two user input epigenomic region sets*.

	 **file1**: file path to the user-input file 1.

	 **file2**: file path to the user-input file 2.

	 **output_path**: absoluate or relative directory to store the Gibbs sampler data.

   **format**: One of `"bed"`, `"narrowPeak"`, `"broadPeak"`, `"bigNarrowPeak"`, `"csv"`, or `NULL`. Specifies the format of the input file. Default: `NULL`. If set to `NULL`, **BIT** will automatically determine the file format based on its extension.

   **N**: Integer. The number of iterations in the Gibbs sampler. Default: `5000`.

   **bin_width**: Integer. The width of the bin used to divide the chromatin into non-overlapping bins. Default: `1000`. Only change this if you compile a different reference database.

   **burnin**: `NULL` or an integer. Specifies the burn-in period when deriving the Bayesian inference for TR-level parameters. Default: `NULL`, which will use `N/2`.

   **genome**: `hg38` or `mm10` for TR ChIP-seq data collected from different genome.



.. function:: display_tables (file_path, output_path, burnin=NULL)

	*To show the ranking table by inspecting the results of Gibbs sampler*.

	**file_path**: path to the saved BIT Gibbs sampling results.

	**output_path**: path to save the rank table.

	**burnin**: number of samples used for burn-in. If not specify, BIT will use the half of the iterations as burn in.



.. function:: rank_plot (file_path=NULL, output_path, burnin=NULL, n=10, colors="NPG", main=NULL, xlab="BIT score", ylab="TR symbols")

	*To draw a barplot for the top n TRs*.

    **file_path**: path to the saved BIT Gibbs sampling results.

    **output_path**: path to save the barplot.

    **burnin**: number of samples used for burn-in. If not specify, BIT will use the half of the iterations as burn in.

    **n**: top n TRs will show in the barplot, default: 10.

    **colors**: colors for each bar, default "NPG" for n<=10, has to be manually specified if n>10.

    **main**: main title for the barplot, default: NULL.

    **xlab**: x axis label, default: BIT score.

    **ylab**: y axis label, default: TR symbols.



.. function:: compare_scatter_plot (file1_path, file2_path, output_path, burnin=NULL)

	*To draw a scatterplot for the comparison bewteen two input region sets*.

    **file1_path**: path to the saved BIT Gibbs sampling results of input 1.

    **file2_path**: path to the saved BIT Gibbs sampling results of input 2.

    **output_path**: path to save the barplot.

    **burnin**: number of samples used for burn-in. If not specify, BIT will use the half of the iterations as burn in.


.. function:: load_chip_data (data_path, bin_width, genome = c("hg38", "mm10"))

	*load the pre-compiled chip-seq data*.

    **data_path**: path to the ChIP-seq data folder, can be absolute or relative path.

    **bin_width** width of bin, which should be in 100/500/1000 and map with your ChIP-seq data.

    **genome**: `hg38` or `mm10` for TR ChIP-seq data collected from different genome.


.. function:: import_input_regions (file, format = NULL, bin_width = 1000, genome=c("hg38", "mm10"))

  *Transform the input regions to binary vector*.

    **file**: file path to the user-input.

    **bin_width**: desired width of bin, default: 1000.

    **genome**: the genome of TR ChIP-seq data, either as "hg38" or "mm10".


.. function:: alignment_wrapper (input_vec, bin_width, genome=c("hg38", "mm10"))

  *Count the 'good' and 'informative' cases by comparing input with the reference database.*

    **input_vec**: A input vector contains index of transformed regions by applying import_input_regions.

    **bin_width**: desired width of bin, default: 1000.

    **genome**: the genome of TR ChIP-seq data, either as "hg38" or "mm10".
