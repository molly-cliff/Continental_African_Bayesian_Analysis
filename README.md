# Spatiotemporal Modelling of Bacterial Meningitis Outbreaks using Climatic Drivers in Africa

[![DOI](https://zenodo.org/badge/1354748639.svg)](https://doi.org/10.5281/zenodo.23803897)


This repository contains the code used to develop a district-level Bayesian spatiotemporal model to estimate monthly meningitis risk across 18 African countries using climatic and demographic variables. Building on our previous analysis, the current study investigates how the relationship between climatic factors and meningitis has evolved spatially and temporally across the belt. Within this study we analysed district-month (ADMN2) epidemic occurrences from 2003 to 2022, assessing the effects of temperature, specific humidity, wind speed and direction, dust, rainfall, and land cover. We also examined the influence of population density and the occurrence of MenAfriVac vaccination campaigns. The study uses the Integrated Nested Laplace Approximation (INLA) approach via the R-INLA package to account for spatial and temporal heterogeneity in meningitis outbreak occurrence.

# Data availability

The district-month outbreak data used in this study were provided by WHO-AFRO under a data-sharing agreement and cannot be shared publicly. The data provided to us were anonymised weekly counts per district, which we cleaned in Stata into a binary outbreak indicator (1/0) per district-month.  Additionally, as the analysis relies on the use of district-level identifying geometry, which would reveal exact outbreak locations, the code cannot be run from this repository alone. To support transparency and reuse, we have provided an anonymised version of the final dataset at xxxxx with the same structure and variables as the original analysis data frame.

As such this repository serves as an overview of the code used within each step of the analysis, as detailed within the respository contents below. Access to the original data can be requested via Clément Lingani who can be contacted at linganic@who.int for further advice. Additionally, as the underlying anonymised dataset is very large we reccomend running any analyses using it on an HPC.

Spatial boundaries (ADMN2) are publicly available from [GADM v 4.2 / (https://gadm.org/data.html) ; the script that builds the spatial adjacency matrix from them is included (folder 2).

# System Requirements

Within this analysis we used R version 4.6, and R-INLA version 26.06.08  (not on CRAN; version can be installed via the following command: remotes::install_version("INLA", version = "26.06.08", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/testing"), dep = TRUE). Additional key packages are as follows:  dplyr 1.2.1 , sf 1.1-1,sp 2.2-1 ,spdep 1.4-2 ,ggplot24.0.3 ,pROC1.19.0.1, terra 1.9-34, Matrix 1.7-5.

# Repository contents

1. [Data_processing](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/1.Data_processing): First, climatic variables were extracted on a district (ADM2) month level following the code within the [extraction folder](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/1.Data_processing/a-Extraction) for each variable. Here, the file shape2 refers to the African continent-wide shapefile built using country-level shapefiles from GADM. The full outbreak matched file represents district-month-level outbreak occurrence from Africa between 2003 and 2022. Alongside this within the landcover file LUH2 state data representing the fractions of each grid-cell occupied by various land uses in a given year were aggregated into categories representing land-use states associated with meningitis outbreak occurrence/non-occurence within the literature (barren, forested, cropland). Following this, the joining [sub-folder](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/1.Data_processing/b-Joining) was used to combine the extracted climatic variables from the different data sources into a single outbreak–climate dataset. This was first done through running [Climate_data_joining](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/1.Data_processing/b-Joining/Climate_data_joining.R) Following on from this we created an additional variable (MenAfriVac campaign status), providing the district-month locations and timing of MenAfriVac campaigns, using campaign data provided by WHO AFRO. Due to the sensitive nature of this data, the associated R file cannot be shared. Finally, the [pop_density_joining script] ((https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/1.Data_processing/b-Joining/Pop_density_joining.R) was run to extract and merge population density data with the existing dataset, producing the final analytical dataset used in the subsequent analyses.


2. [Space_time_analysis](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/2.Space_time_analysis): This folder contains the examination and selection of spatial and temporal random-effects structures for the climatic INLA model. Temporal structures (AR1, RW1 and RW2) were compared, considering month as both a continuous and cyclic variable, alongside four spatial structures (Besag, Besag proper, BYM and BYM2). The temporal and spatial components were assessed separately using WAIC/DIC, with the [temporal analysis conducted first](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/2.Space_time_analysis/a- Temporal_random_effects.R), followed by the [spatial analysis](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/2.Space_time_analysis/b-Spatial_random_effects.R). The results of these comparisons are provided in the minimal dataset. The spatial_random_effects script also constructs the spatial adjacency matrix using GADM ADM2 boundaries and queen contiguity (spdep::poly2nb()), before converting it into an INLA-compatible graph using nb2INLA(). This graph was subsequently used throughout the remaining climatic analyses.

3. [Outbreak_cluster_analysis] (https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/3.Outbreak_cluster_analysis): This folder contains the initial outbreak analysis informing  Figures 1a and 1b of the manuscript, examining the  geographic distribution and timing of meningitis outbreak events across Africa. Here we pull in the outbreak dataset, examine the number of district-month affected by an outbreak and use this to create Figure 1a.  Following on from this we created Figure 1b, which defined the number of seperate outbreak events on a ADMN district level across Africa from 2003-2022. 

4. [Univariate_analysis] (https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/4. Univariate_analysis): Here we carry out a univariate analysis of climatic and socioeconomic variable relationships with district-month level bacterial meningitis outbreaks.  Continuous covariates (temperature, humidity, rainfall, meridional and zonal wind, wind speed, dust and population density have been standardised to z scores). The univariate analysis provides both the odds and log odds of variable association with district-month outbreak occurrence. This R script also assesses for variable collinearity through the calculation of the Pearson correlation coefficients and VIF (threshold VIF > 5]). This analysis led to the inclusion of temperature as the first variable within our stepwise analysis.

5. [Multivariate_analysis] (https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/5. Multivariate_analysis): Here we follow a stepwise analysis adding in variables linearly based on changes in WAIC and DIC, statistical significance as well as epidemiological relevance. The minimal dataset supporting this analysis provies the WAIC and DIC values for all variables at each step. Models are run in order a-d, with zonal wind, dust and humidity being added to the model in each sequential step.


6. [Non_linearity] (https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/6.Non_linearity): We first evaluated non-linear relationships using [RW2 effects for the included explanatory variables](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/6.Non_linearity/a- RW2_exploration.R). First all linear model variables were examined for potential non-linear RW2 effects, with this analysis being carried out in a stepwise fashion. This led to the inclusion of an RW2 effect for temperature, but in the second step of the analysis the inclusion of any additional RW2 effects did not improve model fit. Following this, we carried out a [secondary analysis](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/6.Non_linearity/b- Zonal_wind_dust_interaction.R) examining  potential interactions between zonal wind and dust, through the creation of additional zonal wind strength and direction variables. 

7. [Cross-validation](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/7.Cross_validation): Here we conducted stratified five-fold cross-validation of the final model, with predictive performance evaluated using the Brier score, log loss, recall and the ROC curve (Figure 5).  Before splitting the data, consecutive outbreak months within each district were grouped into outbreak episodes, with all months belonging to the same episode assigned to the same fold. In turn, outcomes within the test fold is masked the and the INLA model is fitted onto the remaining train folds, generating predictions for the held-out observations.

8.[Plots:](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/8.Plots) The code here was used to generate Figures 2-4, and represents the extracted non-linear effects from our INLA model, through holding other variables as constant. We first extract [our modelled RW2 temperature] (https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/8.Plots/a-non_linear_temperature_plots.R) results, generating the log odds of outbreak occurence within a given range of temperature, with the variable converted back to degrees celsius. Alongside this we have also generated a histogram with the distribuion of district-months affected by outbreaks across a range of temperatures. Following on from this we have extracted out and generated plots for the [RW2 temporal effect, BYM2 spatial effect and dust and directional zonal wind interaction] (https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/8.Plots/b-non_linear_effect_space_time_interaction.R). Finally the [(]predicted probability R file] (https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/8.Plots/c-predicted_probability_map.R) uses a shapefile to alongside the predicted probability of outbreak occurence within a given district-month to create an animated map of weighted outbreak risk over time.

9. [Sensitivity analysis:](https://github.com/molly-cliff/Continental_African_Bayesian_Analysis/tree/main/9.Sensitivity_analysis) Finally we carried out a prior sensitivity analyses assessing the robustness of model estimates to alternative prior specifications, including their effects on fixed-effect estimates, WAIC, and DIC.

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

