# softwarePackageSimplify2021

Data repository and scripts associated with "A Software Package to Simplify Tikhonov Regularization with Examples for Matrix-Based Multi-Charge Inversion of SMPS and HTDMA Data", Markus D. Petters, 2021

## Data

### Raw Data 
```data/raw/bbmlsmps/```: Preprocessed response functions from SMPS systems deployed at Bodega Bay Marine Laboratory. The response functions comprise the measured CPC concentration as a function of diameter on a 120 bin diameter grid. The CSV format is self-explanatory. Diameter bin edges and midpoints are given in the first three rows. RH denotes relative humidity during the scan. The data in the file were derived from preprocessed data archive on Zenodo

Petters, Markus D., Rothfuss, Nicholas E., Taylor, Hans, Kreidenweis, Sonia M., DeMott, Paul J., and Atwood, Samuel A. (2019). Size-resolved cloud condensation nuclei data collected during the CalWater 2015 field campaign (Version v1.0) [Data set]. Zenodo. [http://doi.org/10.5281/zenodo.2605668](http://doi.org/10.5281/zenodo.2605668).


```data/sgpaoshtdma/```: Data files in netCDF formar downloaded from the DOE ARM repository and pruned to the days analyzed in the manuscript.

Atmospheric Radiation Measurement (ARM) user facility. 2017, updated hourly. Humidified Tandem Differential Mobility Analyzer (AOSHTDMA). 2020-01-01 to 2020-02-22, Southern Great Plains (SGP) Lamont, OK (Extended and Co-located with C1) (E13). Compiled by J. Uin, C. Salwen and G. Senum. ARM Data Center. Data set accessed 2020-09-25 at [http://dx.doi.org/10.5439/1095581](http://dx.doi.org/10.5439/1095581).


```data/sgpaossmps/```: Data files in netCDF format downloaded from the DOE ARM repository and pruned to the days analyzed in the manuscript.

Atmospheric Radiation Measurement (ARM) user facility. 2016, updated hourly. Scanning mobility particle sizer (AOSSMPS). 2020-01-01 to 2020-09-27, Southern Great Plains (SGP) Lamont, OK (Extended and Co-located with C1) (E13). Compiled by C. Kuang, C. Salwen, M. Boyer and A. Singh. ARM Data Center. Data set accessed 2020-09-29 at [http://dx.doi.org/10.5439/1095583](http://dx.doi.org/10.5439/1095583).


### Processed Data

Processed data are generated from the raw data files using the scripts in ```src/processing/```

```processed/sgpaossmps.csv```: Data from netCDF archive converted to common CSV format for easier handling.

```processed/methodsummary.csv```: Output from the 60000 Monte-Carlo simulations for plotting

```processed/bbmlsmpsL0x0B.csv```: Inverted SMPS data from Bodega Bay Marine Laboratory using the L0x0B method, stored in common CSV format. 

```processed/bbmlsmpsL2B.csv```: Inverted SMPS data from Bodega Bay Marine Laboratory using the L2B method, stored in common CSV format. 

```processed/sgpaoshtdmainverted.csv```: Summary of inverted HTDMA data 

## Scripts

### Processing Scripts
```src/processing/```: Self contained set of scripts for data processing. These scripts should be executed from within the folder. The Project.toml and Manifest.toml files summarize the environment (package versions) that was active when running the scripts. 

### Figure Scripts
```src/figures/```: Self contained set of scripts to generate the figures in the manuscript. These scripts should be executed from within the folder. The Project.toml and Manifest.toml files summarize the environment (package versions) that was active when running the scripts.
