# Bayesian Set Benchmark Dose (BS-BMD) estimation 
This repository holds the code and data to reproduce the results of our manuscript, "Bayesian Gene Set Benchmark Dose Estimation for "omic'' responses."  A production version of the software suitable for non-experts is being developed and will be released as either a standalone package or as a feature of an existing package.  This repo will eventually redirect to the production repo of https://github.com/NIEHS/transcryptR, which is currently internal to the NIEHS.

The code has not be prepared for novice users but some basic instructions are provided below for reference.

1. Open the BSBMD-manuscript.Rproj project file to automatically have the correct working directory and local references.
2. Run the MCMC script, specifying which file to run by index.  The default index is 18, referring to Fenofibrate Liver.  The number of iterations is set to 2000 and the BMD samples for each hallmark set are saved.  The sampled parameters are not saved by default due to potential storage issues, but there is a boolean that can be set to True and the parameters will also be saved.  Note that for 5000 iterations, the parameter file will take about 3GB for each input file.
4. The BMD_plots.R script can reproduce most of the plots of the paper.  Starting at around line 750, the output from the MCMC can be aggregated into a single large data frame, "hallmark_bigset_df," and used to reproduce the posterior ridge plots.
5.  Other plots, such as the gene contribution for the test statistic, can be run after the MCMC has run a least a few iterations and uses the current parameter values in memory to create the plot.



