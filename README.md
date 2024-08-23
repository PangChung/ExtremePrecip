This repository contains code and the data for the manuscript "Spatial modeling and future projection of extreme
precipitation extents". 

The dataset we studied in the project consists of observed daily precipitation data in millimeters (mm) from 125 monitoring stations in the Danube river basin (Europe) and from 2229 monitoring stations in the Mississippi river basin (North America) over the period from 1965 to 2020. We divided the Mississippi river basin into 7 subregions according to Watershed Boundary Dataset, see the figure below. The temperature covariate (Celsius degree) used to fit the model in this project was derived from the ERA5-Land reanalysis data for the corresponding region we selected, which is a global land-surface dataset with a high spatial resolution of 9km. The projected temperature covariate is derived from climate models outputs of the sixth Coupled Model Intercomparison Project (CMIP6), namely AWI, MIROC, and MPI.

We have saved the data used in our marginal modeling and dependence modeling in the format of `R` data file (.RData), and we host those data on Github[https://github.com/PangChung/ExtremePrecip] under the folder `data`, which are publicly accessible. Though the raw data of the ERA5-Land Reanalysis data from which we derived the temperature covariate in Celsius degree can be downloaded from the website https://www.copernicus.eu/en/access-data, this raw data will be provided upon email request as the size of it is hundreds of gigabytes and it will take great efforts for someone to download directly from the Copernicus website. The email address is peng[dot]zhong[at]unsw[dot]edu[dot]au. The descriptions of each `R` data file are following:

* `data/precip.RData`: The precipitation data for the eight subregions, which contain `R` objects:
   + `precip`: Lists of lists, raw data of the precipitation in millimeters (mm) in 8 subregions. This is the response variable used in the marginal fit and dependence fit.
   + `region.name`: a vector, names of the 8 subregions
   + `region.id`: a vector, regional ID number that corresponding to the regional number in the shapefiles.
   + `station`: a data frame, contains geographical information about the monitoring stations, where the precipitation data were recorded. The ``X'' column is latitude degrees, the ``Y'' column is the longitude in degrees, the ``start'' column is the start measuring year, the ``end'' is the last measuring year, the ``elev'' is the elevation of the monitoring stations in meters, and the ``group.id'' corresponds to the region.id for identifying each of the 8 subregions.
   + `START.date and END.date`: dates, between which the data were used in our analysis. 
* `data/temperature.RData`: The derived temperature covariate over the period 1965--2020, which contains `R` objects:
    + `date.df`: data frame, contains the date between 1965--2020, and its corresponding season of the year. 
    + `temperature.covariate`: list of vectors, the derived temperature covariate (30 day moving averages) used in marginal modeling and dependence modeling.
    + `temperature`: list of vectors, daily spatial temperature avarages in Celsius degrees, used only during data preprocessing. 
    + `loc_df`: geographical information about the locations of the ERA5-Land temperature data, used only during data preprocessing. 
* `data/temperature_pred.RData`: The derived temperature covariate from the climate models over the period 2015--2100 under different shared socioeconomic pathways (SSP 2-4.5 or SSP 5-8.5), which contains `R` objects:
    + `date.245 and date.585`: vector of dates, corresponding to the temperature covariate from the climate models. 
    + `idx.models`: 5 climate models, we used the first, the 3rd, and the 4th for future projections, which are denoted by AWI, MIROC, and MPI. 
    + `temperature.245.avg and temperature.585.avg`: lists of list, each list contains the derived temperature covariate in Celsius degrees from one climate model for future projections before realign with the temperature covariate derived from the ERA5-Land data under SSP 2-4.5 (or SSP 5-8.5). SSP 2-4.5 is denoted by 245, and the same goes for SSP 5-8.5.    
* `data/marginal_fit_quantiles.RData`: Transformed margins based on the marginal fit,  which contains `R` objects:
    + `theretical.quantiles`: list of matrices, each matrix corresponding to the fitted marginal data based on the marginal model fit for each subregion, this is only used to generate the QQplots in the manuscript. 
* `data/dep.fit.boot.results3.RData`: Fitted results from the bootstrap scheme for the dependence model, which contains `R` objects:
    + `boot.result.df`: data frame, contains the 95\% confidence interval (e.g, the column low.shape, high.shape) based on the bootstrap for each season, each region, and each risk functional.
    + `boot.result.list`:  lists of lists for each risk functional (first dim), each season (second dim) and each subregion (third dim). In each list, it contains the standard error of the three parameter estimates, the estimates using the whole data set, and the 300 nonparametric bootstrap estimates, named by jack.
* `data/era5_geoinfo.RData`: Shapefiles for the 8 subregions, which contains `R` objects:
    + `shape1`: shapefile of the US river-basins.
    + `shape3`: shapefile of Danube river-basin.
* `data/transformed_coordinates.RData`: Transformed coordinates for the eight subregions, which transformed the latitude/longitude coordinate to the Euclidean coordinate and contains `R` objects:  
    + `loc.trans.list`: list of matrices: transformed coordinates for each of the 8 subregions in meters. 

![Subregions of the two river basins.](/figures/precip_loc.png)

The objective of this project is to modeling the precipitation data in Mississippi river basin and Danube river basins and use the temperature covariate derived as described above to study the dynamic changes in the marginal and dependence structure of the precipitation data, i.e., the marginal return level and the spatial extent of the extreme precipitation events.  

To reproduce the results shown in the manuscript, please follow the instructions in the R markdown file `ResultsReproduction.pdf`.