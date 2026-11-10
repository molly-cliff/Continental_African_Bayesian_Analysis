# Spatiotemporal Modelling of Bacterial Meningitis Outbreaks using Climatic Drivers in Africa

[![DOI](https://zenodo.org/badge/1354748639.svg)](https://doi.org/10.5281/zenodo.23803897)


This repository contains the code used to develop a district-level Bayesian spatiotemporal model to estimate monthly meningitis risk across 18 African countries using climatic and demographic variables. Building on our previous analysis, the current study investigates how the relationship between climatic factors and meningitis has evolved spatially and temporally across the belt. Within this study we analysed district-month (ADMN2) epidemic occurrences from 2003 to 2022, assessing the effects of temperature, specific humidity, wind speed and direction, dust, rainfall, and land cover. We also examined the influence of population density and the occurrence of MenAfriVac vaccination campaigns. The study uses the Integrated Nested Laplace Approximation (INLA) approach via the R-INLA package to account for spatial and temporal heterogeneity in meningitis outbreak occurrence.

# Data availability

The district-month outbreak data used in this study were provided by WHO-AFRO under a data-sharing agreement and cannot be shared publicly. The data provided to us were anonymised weekly counts per district, which we cleaned in Stata into a binary outbreak indicator (1/0) per district-month.  Additionally, as the analysis relies on the use of district-level identifying geometry, which would reveal exact outbreak locations, the code cannot be run from this repository alone. To support transparency and reuse, we have provided an anonymised version of the final dataset at xxxxx with the same structure and variables as the original analysis data frame.

As such this repository serves as an overview of the code used within each step of the analysis, as detailed within the respository contents below. Access to the original data can be requested via Clément Lingani who can be contacted at linganic@who.int for further advice. Additionally, as the underlying anonymised dataset is very large we reccomend running any analyses using it on an HPC.

Spatial boundaries (ADMN2) are publicly available from [GADM v 4.2 / (https://gadm.org/data.html) ; the script that builds the spatial adjacency matrix from them is included (folder 2).
Public input data for the climate and covariate variables are not stored in this repository. Sources:

| Descriptor | Variable summary | INLA data source |
|---|---|---|
| Specific humidity | Monthly mean specific humidity (g/kg). Specific humidity refers to the mass of water vapour per unit mass of moist air. | ERA5 reanalysis – ECMWF [1] |
| Rainfall | Monthly mean atmospheric precipitation (mm). | CHIRPS [2] |
| Aerosol optical depth (dust) | Monthly mean measure of the light extinction effect of atmospheric aerosols. | NASA MODIS [3] |
| Land-use type | Yearly ADM2 district-level proportions of land classified as barren land, forest, or cropland. | LUH2 [4] |
| Wind speed | Monthly mean wind speed (m/s) at 10 m above ground level. | ERA5 reanalysis – ECMWF [1] |
| Meridional wind component | Monthly mean meridional wind component (m/s) at 10 m above ground level. Meridional refers to north–south wind flow. | ERA5 reanalysis – ECMWF [1] |
| Surface-level temperature | Monthly mean surface air temperature. | NASA AIRS [5] |
| Population density | Population density at five-year intervals (2005, 2010, 2015, and 2020). Calculated by dividing the population count grid for each target year by the corresponding gridded land area. | SEDAC [6] |

1. Hersbach H, Bell B, Berrisford P, Biavati G, Horányi A, Muñoz Sabater J, et al. ERA5 monthly averaged data on single levels from 1940 to present; 2023 [cited 2026 Aug 27]. Database: Copernicus Climate Change Service (C3S) Climate Data Store (CDS) [Internet]. doi: 10.24381/cds.f17050d7
2. Funk C, Peterson P, Landsfeld M, Pedreros D, Verdin J, Shukla S, et al. The climate hazards infrared precipitation with stations-a new environmental record for monitoring extremes. Sci Data. 2015 Dec 8:2:150066. doi: 10.1038/sdata.2015.66
3. Levy R, Hsu C, et al. MODIS Atmosphere L2 Aerosol Product [dataset]. NASA MODIS Adaptive Processing System, Goddard Space Flight Center, USA; 2015. doi:10.5067/MODIS/MOD04_L2.006. Accessed 27 Aug 2026.
4. Hurtt GC, Chini L, Sahajpal R, Frolking S, Bodirsky BL, Calvin K, et al. Harmonization of global land use change and management for the period 850-2100 (LUH2) for CMIP6. Geoscientific Model Development. 2020 Nov 10;13(11):5425-64. doi: 10.5194/gmd-13-5425-2020
5. AIRS Science Team, Texeira J. Aqua AIRS Level 3 Standard Monthly Product using AIRS IR-only V6; 2013 [cited 2026 Aug 27]. Database: NASA Goddard Earth Sciences Data and Information Services Center (GES DISC) [Internet]. doi: 10.5067/Aqua/AIRS/DATA321
6. Center for International Earth Science Information Network (CIESIN), Columbia University. Gridded Population of the World, Version 4 (GPWv4): Population Density, Revision 11 (Version 4.11) [dataset]. Palisades (NY): Socioeconomic Data and Applications Center (SEDAC); 2017. doi:10.7927/H49C6VHW. Accessed 27 Aug 2026.
Gridded values were aggregated to ADMN2-month mean values using terra::extract. Continuous covariates (temperature, humidity, rainfall, meridional and zonal wind, wind speed, AOD and population density) were standardised to z scores using the scale  function in base R

# System Requirements

Within this analysis we used R version 4.6, and R-INLA version 26.06.08  (not on CRAN; version can be installed via the following command: remotes::install_version("INLA", version = "26.06.08", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/testing"), dep = TRUE). Additional key packages are as follows:  dplyr 1.2.1 , sf 1.1-1,sp 2.2-1 ,spdep 1.4-2 ,ggplot24.0.3 ,pROC1.19.0.1, terra 1.9-34, Matrix 1.7-5.

# Repository contents

1. [Data_processing](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/1.Data_processing): First, climatic variables were extracted on a district (ADM2) month level following the code within the [extraction folder](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/1.Data_processing/a-Extraction) for each variable. Here, the file shape2 refers to the African continent-wide shapefile built using country-level shapefiles from GADM. The full outbreak matched file represents district-month-level outbreak occurrence from Africa between 2003 and 2022. Alongside this within the landcover file LUH2 state data representing the fractions of each grid-cell occupied by various land uses in a given year were aggregated into categories representing land-use states associated with meningitis outbreak occurrence/non-occurence within the literature (barren, forested, cropland).

Following this, the joining [sub-folder](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/1.Data_processing/b-Joining) was used to combine the extracted climatic variables from the different data sources into a single outbreak–climate dataset. This was first done through running [Climate_data_joining](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/1.Data_processing/b-Joining/Climate_data_joining.R) Following on from this we created an additional variable (MenAfriVac campaign status), providing the district-month locations and timing of MenAfriVac campaigns, using campaign data provided by WHO AFRO. Due to the sensitive nature of this data, the associated R file cannot be shared. Finally, the [pop_density_joining script] ((https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/1.Data_processing/b-Joining/Pop_density_joining.R) was run to extract and merge population density data with the existing dataset, producing the final analytical dataset used in the subsequent analyses.

edit source of pop density data.

2. Space_time_analysis: This folder contains the examination and selection of spatial and temporal random-effects structures for the climatic INLA model. Temporal structures (AR1, RW1 and RW2) were compared, considering month as both a continuous and cyclic variable, alongside four spatial structures (Besag, Besag proper, BYM and BYM2). The temporal and spatial components were assessed separately using WAIC/DIC, with the temporal analysis conducted first, followed by the spatial analysis. The results of these comparisons are provided in the minimal dataset.

Based on these comparisons, an RW2 temporal effect, with month treated as continuous, and a BYM2 spatial effect were selected for subsequent analyses, which were restricted to meningitis belt countries. The spatial_random_effects script also constructs the spatial adjacency matrix using GADM ADM2 boundaries and queen contiguity (spdep::poly2nb()), before converting it into an INLA-compatible graph using nb2INLA(). This graph was subsequently used throughout the remaining climatic analyses.

3. Outbreak_cluster_analysis: This folder contains the initial outbreak analysis informing  Figures 1a and 1b of the manuscript, examining the number, timing, and geographic distribution of meningitis outbreak events across Africa. Here we pull in the outbreak dataset, examine the number of district-month affected by an outbreak and use this to create Figure 1a.  Following on from this we created Figure 1b, which defined the number of seperate outbreak events on a ADMN district level across Africa from 2003-2022. To achieve this we built an  adjacency matrix to define neighbouring districts affected by the same outbreak event and defined seperate outbreak events through creating a lag effect defining the time gapp between outbreaks within the same district.

4. Univariate analysis: Here we carry out a univariate analysis of climatic and socioeconomic variable relationships with district-month level bacterial meningitis outbreaks. Using the results from the space-time analysis here we have included a BYM2 spatial effect, a continuous RW2 temporal effect for month, and inverse probability weighting due to the rareity of outbreak events. Continuous covariates (temperature, humidity, rainfall, meridional and zonal wind, wind speed, dust and population density have been standardised to z scores). The univariate analysis provides both the odds and log odds of variable association with district-month outbreak occurrence. Prior to the univariate analysis, this R script also assesses for variable collinearity through the calculation of the Pearson correlation coefficients (Supplementary Information) and VIF analysis for collinearity (threshold VIF > 5]). This analysis led to the inclusion of temperature as the first variable within our stepwise analysis.

5. Multivariate analysis: Here we follow a stepwise analysis adding in climatic/socioeconomic variables linearly based on changes in WAIC and DIC, statistical significance as well as epidemiological relevance. The minimal dataset supporting this analysis provies the WAIC and DIC values for all variables at each step. Models are run in order a-d, with zonal wind, dust and humidity being added to the model in each sequential step.


6. We first ran the rw2_exploration file to evaluate non-linear relationships using RW2 effects for the included explanatory variables, alongside an assessment of the interaction between dust and zonal wind. In File A, all explanatory variables were tested for potential non-linear RW2 effects, with temperature also supported as a non-linear effect in the minimal dataset. The same approach was repeated in File B; however, the inclusion of RW2 effects did not produce any noticeable improvement in model fit, with no appreciable reduction in WAIC or DIC.

Following this, we ran the zonal_wind_dust interaction analysis to investigate the relationship between dust and zonal wind. This involved examining the interaction in relation to zonal wind strength and direction, with alternative model specifications compared using WAIC and DIC to assess whether including the interaction improved model fit.

7. Cross-validation: Here we conducted stratified five-fold cross-validation of the final model, with predictive performance evaluated using the Brier score, log loss, and area under the receiver operating characteristic curve (ROC AUC) (Figure 5).  Before splitting the data, consecutive outbreak months within each district are grouped into outbreak episodes. All months belonging to the same episode are assigned to the same fold, preventing a single outbreak episode. For each fold outcomes in the test folds are maskeed and INLA model is fitted using the remaining observations, Predictions are generated for the held-out observations and evaluated using brier score, log loss, roc cruve, precision and recall. ROC curves from all five folds are combined into a single plot, and the fold-level metrics are combined into an overall cross-validation results file. 

8. Plots: Final model and plots for random effects (temperature, space, and time), as well as the interaction effect between dust and wind (Figures 2-4 and Supplementary). A GIF of monthly predicted probability across Africa is also included here (Supplementary Information) [or: hosted in the supplementary repository, DOI].



9. Sensitivity analysis: Prior sensitivity analyses assessing the robustness of model estimates to alternative prior specifications, including their effects on fixed-effect estimates, WAIC, and DIC. This includes the fixed effect plot and WAIC/DIC table for different prior setups (Supplementary Information).

# LICENSE 

# MIT License

Copyright (c) [2026] [Molly Cliff]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

