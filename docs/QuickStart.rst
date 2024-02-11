Quick Start
===========

After successfully installing BIT and configuring the reference ChIP-seq database, you're ready to use BIT for your analyses. The process is straightforward if you have any of the following file types indicating the chromosomal start and end positions of your epigenomic regions: .bed, .narrowPeak, .broadPeak, .bigNarrowPeak (.bb), or .csv. 

Assuming your input file is located at "file_path/file_name (acceptable formats: .bed, .bb, etc.)", and you aim to store the output in "output_path", the command below will first start a comprehensive comparison of your input regions against all available ChIP-seq reference datasets. Next it initiates the BIT model using the Gibbs sampler to give the Bayesian inference of TR level BIT score, the output data from each iteration will be saved in a R data object named "file_name.rds".

.. code-block:: R

	> BIT("file_path/file_name(.bed / .bb, etc..)","output_path")
	[1] "Load and map peaks to bins..."
	[1] "Done."
	[1] "Compare the input regions with the pre-compiled reference ChIP-seq data, bin width used: 1000 bps"
	==================================================
	[1] "Done."
	[1] "Start BIT Gibbs sampler for file 1, iterations: 5000"
	0%   10   20   30   40   50   60   70   80   90   100%
	[----|----|----|----|----|----|----|----|----|----|
	**************************************************|
	[1] "Done."
	[1] "Output data saved as output_path/file_name.rds"

There are several default parameters you can change:

.. note::

	show = TRUE / FALSE, whether to show the ranking table, default: TRUE.

	plot.bar = TRUE / FALSE, whether to plot the top 10 TRs BIT score in a horizontal barplot, default TRUE. 

	format = c("bed", "narrowPeak", "broadPeak", "bigNarrowPeak", "csv", NULL), to specify the format of input file, default: NULL. When setting as NULL, BIT will automatically judge the file format by its extension.

	N = 5000 , number of iterations in Gibbs sampler.

	bin_width = c(100, 500, 1000), width of bin used to separate the chromatin into non-overlapping bins, default: 1000.

	option = c("ALL", "PLS", "ELS"), Whether to filter the input regions with ENCODE defined cis-regulatory elements, default: "ALL" to use all input regions.

	burnin = NULL / integer, burn in when derive the Bayesian inference for TR level parameters, default: NULL will use N/2.


If you set show = FALSE and later want to generate the rank table, you can use the following command:

.. code-block:: R

	Results_Table<-display_tables("file_path/file_name.rds","output_path")
	> Results_Table
	         TR   Theta_i     lower     upper  BIT_score BIT_score_lower BIT_score_upper Rank
	1      CTCF -2.010571 -2.011593 -2.009676 0.11809745      0.11799114      0.11819079    1
	2     RAD21 -2.028610 -2.031747 -2.025619 0.11623164      0.11590978      0.11653925    2
	3      SMC3 -2.110100 -2.120542 -2.100907 0.10811898      0.10711622      0.10900866    3
	4     SMC1A -2.181305 -2.193447 -2.170246 0.10144192      0.10034048      0.10245443    4
	5      PHF2 -2.222166 -2.633869 -1.990377 0.09777755      0.06699025      0.12021697    5
	6    GABPB1 -2.301673 -2.689461 -2.062208 0.09098450      0.06359813      0.11282466    6
	7      CHD2 -2.317842 -2.370804 -2.270597 0.08965600      0.08542628      0.09358759    7
	8      NFYA -2.353494 -2.396030 -2.313238 0.08678843      0.08347594      0.09003250    8
	9    NELFCD -2.371794 -2.476460 -2.276018 0.08534898      0.07752502      0.09312872    9
	10     NFYC -2.373795 -2.767204 -2.122314 0.08519290      0.05912233      0.10694685   10