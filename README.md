# Spatiotemporal Modelling of Bacterial Meningitis Outbreaks using Climatic Drivers in Africa

This repository contains the code used to develop a district-level Bayesian spatiotemporal model to estimate monthly meningitis risk across Africa using climatic and demographic variables. Building on our previous analysis, the current study investigates how the relationship between climatic factors and meningitis has evolved spatially and temporally across the continent.

We analysed district-month (ADMN2) epidemic occurrences from 2003 to 2022, assessing the effects of specific humidity, wind speed and direction, dust, rainfall, and land cover. We also examined the influence of population density and the occurrence of MenAfriVac vaccination campaigns. The study uses the Integrated Nested Laplace Approximation (INLA) approach to account for spatial and temporal heterogeneity in meningitis outbreak occurrence.

# Repository contents

1. Data preparation: Code used to process and extract district-month level data on specific humidity, wind speed and direction, dust, rainfall, and land cover for inclusion in the spatiotemporal model. 
2. Spatiotemporal analysis: Examination of spatial and temporal random-effects structures to determine which components should be included in the climatic model.
3. Outbreak cluster analysis: Analysis informing Figures 1a and 1b of the manuscript, examining the number, timing, and geographic distribution of meningitis outbreak events across Africa.
4. Univariate analysis: Assessment of the univariate relationships between climatic and socioeconomic variables and district-month-level bacterial meningitis outbreak occurrence. This also includes the Pearsons correlation coeffient (Supplementary Information) and VIF analysis for collinearity. 
5. Multivariate analysis: Stepwise multivariate analysis examining the relationships between climatic and socioeconomic variables and district-month-level bacterial meningitis outbreak occurrence.
6. Non-linearity and interactions: Evaluation of non-linear relationships using RW2 effects for included explanatory variables, as well as assessment of the interaction between dust and zonal wind.
7. Cross-validation: Stratified five-fold cross-validation of the final model, with predictive performance evaluated using the Brier score, log loss, and area under the receiver operating characteristic curve (ROC AUC).
8. Sensitivity analysis: Prior sensitivity analyses assessing the robustness of model estimates to alternative prior specifications, including their effects on fixed-effect estimates, WAIC, and DIC.
R/Stata version
9. Additional plots

# R/Stata Version

This project was carried out using R version 4.4.2 and Stata SE 18

