This repository contains code and the data for the manuscript "Spatial modeling and future projection of extreme
precipitation extents" 

To reproduce the results shown in the manuscript, one needs to
* First, specify the working directory, which should be *~/ExtremePrecip/* in our case
* Then, one need to install the packages mvPotST as well as other packages listed in the R script files.
* Run the marginal fit and the dependence model fit via the following code, *Rscript code/bootstrap.R "idx.region=1;bootstrap.ind=0;computer=\"local\""*. Here we fit the model for the region 1 (Danube) by setting *idx.region=1*. 
* Generate the plots using the code in *code/plots.R*, in some cases, one can directly load the plots into memory by load the R Data file data/plot_\<name of the plot\>.RData
* Other R scripts are used to pre-process the raw data (provided upon request). 