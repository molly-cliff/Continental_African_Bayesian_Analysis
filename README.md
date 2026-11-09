# Spatiotemporal Modelling of Bacterial Meningitis Outbreaks using Climatic Drivers in Africa

[![DOI](https://zenodo.org/badge/1354748639.svg)](https://doi.org/10.5281/zenodo.23803897)


This repository contains the code used to develop a district-level Bayesian spatiotemporal model to estimate monthly meningitis risk across Africa using climatic and demographic variables. Building on our previous analysis, the current study investigates how the relationship between climatic factors and meningitis has evolved spatially and temporally across the continent.
This repository contains all author-generated code underpinning the findings in the manuscript, shared without restriction. 
Within this study we analysed district-month (ADMN2) epidemic occurrences from 2003 to 2022, assessing the effects of specific humidity, wind speed and direction, dust, rainfall, and land cover. We also examined the influence of population density and the occurrence of MenAfriVac vaccination campaigns. The study uses the Integrated Nested Laplace Approximation (INLA) approach to account for spatial and temporal heterogeneity in meningitis outbreak occurrence.
The analysis requires confidential outbreak data, so the code cannot be run end to end from this repository alone. The data provided to us were anonymised weekly counts per district, which we cleaned in Stata into a binary outbreak indicator (1/0) per district-month. An anonymised version of the final dataset is available at [link].
Data availability
The district-month outbreak data used in this study were provided by WHO-AFRO under a data-sharing agreement and cannot be shared publicly. To support transparency and reuse, we provide:
Anonymised example dataset ([data/example_data.csv]) with the same structure, column names and variable types as the real analysis dataframe. It shows the exact format the code expects. It is not intended to reproduce the analysis, but is provided in full so that every step of the analysis can be inspected. We would recommend running any analyses with the anonymised example dataset on an HPC as the dataset is very large.
The analysis dataframe has one row per ADMN2 district per month across Africa (2003-2022). It includes a 1/0 indicator for whether the district experienced an outbreak, and a 1/0 indicator for whether the district received a MenAfriVac campaign within the given month.
Access to the original data can be requested via Clément Lingani who can be contacted at linganic@who.int for further advice. 
Code for this manuscript is designed to be run in the order of folders.
Final model outputs ([outputs/]) underlying the manuscript tables and figures (fixed-effect summaries, WAIC/DIC tables, cross-validation metrics), so that reported estimates can be checked without access to the confidential data.


Spatial boundaries (ADMN2) are publicly available from [GADM v 4.2 / (https://gadm.org/data.html) ; the script that builds the spatial adjacency matrix from them is included (folder 2).
Public input data for the climate and covariate variables are not stored in this repository. Sources:

| Descriptor | Variable summary | INLA data source |
|---|---|---|
| Specific humidity | Monthly mean specific humidity (g/kg). Specific humidity refers to the mass of water vapour in a unit mass of moist air. | ERA5 reanalysis-ECMWF [1] |
| Rainfall | Monthly mean atmospheric precipitation (mm). | CHIRPS [2] |
| Aerosol optical depth (dust) | Monthly mean measure of the extinction effect of atmospheric aerosols. | NASA MODIS [3] |
| Land-use type | Yearly ADMN2 district-level proportions of land classified as barren land, forest, or cropland. | LUH2 [4] |
| Wind speed | Monthly mean wind speed (m/s) at 10 m above ground level. | ERA5 reanalysis-ECMWF [1] |
| Meridional wind component | Monthly mean meridional component of wind (m/s) at 10 m above ground level. Meridional refers to north/south (longitudinal flow). | ERA5 reanalysis-ECMWF [1] |


1. Hersbach H, Bell B, Berrisford P, Biavati G, Horányi A, Muñoz Sabater J, et al. ERA5 monthly averaged data on single levels from 1940 to present; 2023 [cited 2026 Aug 27]. Database: Copernicus Climate Change Service (C3S) Climate Data Store (CDS) [Internet]. doi: 10.24381/cds.f17050d7
2. Funk C, Peterson P, Landsfeld M, Pedreros D, Verdin J, Shukla S, et al. The climate hazards infrared precipitation with stations-a new environmental record for monitoring extremes. Sci Data. 2015 Dec 8:2:150066. doi: 10.1038/sdata.2015.66
3. Levy R, Hsu C, et al. MODIS Atmosphere L2 Aerosol Product [dataset]. NASA MODIS Adaptive Processing System, Goddard Space Flight Center, USA; 2015. doi:10.5067/MODIS/MOD04_L2.006. Accessed 27 Aug 2026.
4. Hurtt GC, Chini L, Sahajpal R, Frolking S, Bodirsky BL, Calvin K, et al. Harmonization of global land use change and management for the period 850-2100 (LUH2) for CMIP6. Geoscientific Model Development. 2020 Nov 10;13(11):5425-64. doi: 10.5194/gmd-13-5425-2020
5. AIRS Science Team, Texeira J. Aqua AIRS Level 3 Standard Monthly Product using AIRS IR-only V6; 2013 [cited 2026 Aug 27]. Database: NASA Goddard Earth Sciences Data and Information Services Center (GES DISC) [Internet]. doi: 10.5067/Aqua/AIRS/DATA321
6. Center for International Earth Science Information Network (CIESIN), Columbia University. Gridded Population of the World, Version 4 (GPWv4): Population Density, Revision 11 (Version 4.11) [dataset]. Palisades (NY): Socioeconomic Data and Applications Center (SEDAC); 2017. doi:10.7927/H49C6VHW. Accessed 27 Aug 2026.
Gridded values were aggregated to ADMN2-month mean values using terra::extract. Continuous covariates (temperature, humidity, rainfall, meridional and zonal wind, wind speed, AOD and population density) were standardised to z scores using the scale  function in base R


# System Requirements

1.R version 4.6 

2.R-INLA version 26.06.08  (not on CRAN; version can be installed via the following command: remotes::install_version("INLA", version = "26.06.08", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/testing"), dep = TRUE) 

3. Additional key packages: dplyr 1.2.1   , sf 1.1-1, sp 2.2-1 , spdep 1.4-2 , ggplot24.0.3 , pROC1.19.0.1, terra1.9-34, Matrix 1.7-5.

# Repository contents

1. Data_processing: Code used to process and extract district-month level data on specific humidity, wind speed and direction, dust, rainfall, and land cover for inclusion in the spatiotemporal model. Data are first extracted for each climatic variable using the code in the Extraction folder, then joined using the Climate data joining and Population data joining files. Data are joined onto a larger outbreak dataframe, described under Data availability. Because this step uses confidential data, it cannot be run with the example dataset.
   
2. Space_time_analysis: Examination of spatial and temporal random-effects structures to determine which components should be included in the climatic model. Candidate structures compared: [e.g. BYM2, iid, RW1/RW2, seasonal], selected using [WAIC/DIC]. Spatial adjacency matrix. The adjacency structure was built from ADMN2 boundaries ([GADM v x.x / geoBoundaries]) using queen contiguity (spdep::poly2nb) and converted to an INLA graph with nb2INLA. 

3. Outbreak cluster analysis: Analysis informing Figures 1a and 1b of the manuscript, examining the number, timing, and geographic distribution of meningitis outbreak events across Africa. Here we pull in the outbreak dataset, examine the number of district-month affected by an outbreak for Figure 1a. Alongside this we also create an  adjacency matrix, clustering outbreaks by districts adjoining and concurrent month to create a map of number of districts affected by an outbreak
4. Univariate analysis: Assessment of the univariate relationships between climatic and socioeconomic variables and district-month-level bacterial meningitis outbreak occurrence. This also includes the Pearson correlation coefficients (Supplementary Information) and VIF analysis for collinearity ([threshold, e.g. VIF > 5]). Alongside this both we provide code to calculate both the odds and log odds of variable association with district-month outbreak occurrence
5. Multivariate analysis: Stepwise multivariate analysis examining the relationships between climatic and socioeconomic variables and district-month-level bacterial meningitis outbreak occurrence. Variables were retained if [criterion, e.g. improved WAIC].
6. Non-linearity and interactions: Evaluation of non-linear relationships using RW2 effects for included explanatory variables, as well as assessment of the interaction between dust and zonal wind. These are run in order with the multivariate step_1 file preceding the step 2 file.
7. Cross-validation: Stratified five-fold cross-validation of the final model, with predictive performance evaluated using the Brier score, log loss, and area under the receiver operating characteristic curve (ROC AUC) (Figure 5). Folds were stratified by [outcome / district / time].
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






