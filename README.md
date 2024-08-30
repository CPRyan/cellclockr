# cellclockr

## Description

This package is designed to streamline the analysis of associations between epigenetic clocks and immune cell composition across various exposures and contexts. It provides a comprehensive and reproducible pipeline, allowing users to harmonize data products across studies without sharing raw data. The package supports data preparation, model fitting, and output generation, including stratified analyses and visualizations. It is built in R and relies on a minimal set of well-established packages to ensure broad compatibility. Users are guided through setup and execution, with detailed instructions for input data and required configurations. For those preferring not to run the analysis themselves, the package developers offer to perform the calculations upon request.

## Details

This package currently requires the following: 

A dataset containing columns with the following:
 * Epigenetic clock measures (e.g. Horvath, PCPhenoAge, DunedinPACE, etc.)
 * Estimated cell composition measures (e.g. Bcell, CD4Tnaive, etc.)
 * Chronological age (optional, but used in most functions)
 * (Optional) Groups of interest with one or more categories (e.g. sex, self-reported ethnicity, disease status, etc.)
 * (Optional) Continuous variables of interest (e.g. BMI, SES, etc.)
 * (Optional) Covariates to control for unwanted technical variation (e.g. Batch, Array, Plate, Row, Column etc.)
