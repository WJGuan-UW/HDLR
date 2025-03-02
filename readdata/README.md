Unzip the file debinf_0208.zip, which contains all the latest simulation results.

Once unzipped, go to the directory, copy the folders 'debias_res_xx', 'debiascf_res_xx', 'SIHR_res_xx' and 'refit_res_xx' to the same directory as the readdata file.
After you read in the results, you can generate all these plots.
Please manually save all the plots in .pdf file.
For QQ plots, save in 6in * 6in, for the remaining plots, save in 10in * 3.5in.
After generating the QQ plots, don't forget to clear the plot screen. Or you can add a line par(mfrow=c(1,1)) to change back to the original.
Otherwise the shape for the remaining plots will be wrong and too small.
