API References
==============

.. function:: BIT (file, output_path, show=TRUE, plot.bar=TRUE, format=NULL, N = 5000 ,bin_width = 1000, option="ALL",burnin=NULL)

	*Main interface to run BIT method, please set the input file path, input file format, number of iterations and bin width*.
	
	**file**: file path to the user-input.

	**output_path**: absoluate or relative directory to store the Gibbs sampler data.

	**format**: if specify as NULL, BIT will automatically read and judge the file type based on extension.

	**N**: number of iterations for gibbs sampler, recommended for >= 5000, default: 5000.

	**bin_width**: desired width of bin, should be in 100/500/1000, default: 1000.

	**option**: option to filter peaks with candidate cis-regulatory elements from ENCODE, default use "ALL", can be "PLS" for all promoters, "ELS" for all enhancers.



.. function:: BIT_compare (file1, file2, output_path, show=TRUE, plot.scatter=TRUE, format=c(NULL,NULL), N = 5000, bin_width = 1000, option="ALL", burnin=NULL)

	*compare BIT identifid TRs for two user input epigenomic region sets*.

	**file1**: file path to the user-input file 1.

	**file2**: file path to the user-input file 2.

	**output_path**: absoluate or relative directory to store the Gibbs sampler data.

	**format**: if specify as NULL, BIT will automatically read and judge the file type based on extension, default: NULL.

	**N**: number of iterations for gibbs sampler, recommended for >= 5000, default: 5000.

	**bin_width**: desired width of bin, should be in 100/500/1000, default: 1000.

	**option**: option to filter peaks with candidate cis-regulatory elements from ENCODE, default use "ALL", can be "PLS" for all promoters, "ELS" for all enhancers.



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


.. function:: load_chip_data (data_path, bin_width)

	*load the pre-compiled chip-seq data*.

    **data_path**: path to the ChIP-seq data folder, can be absolute or relative path.

    **bin_width** width of bin, which should be in 100/500/1000 and map with your ChIP-seq data.





