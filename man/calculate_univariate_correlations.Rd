% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_univariate_correlations.R
\name{calculate_univariate_correlations}
\alias{calculate_univariate_correlations}
\title{Calculate Univariate Correlations between Cell Types and Clock Measures}
\usage{
calculate_univariate_correlations(
  data,
  study,
  cell_types,
  all_clocks,
  age_column,
  stratify_by = NULL,
  control_covariates = NULL,
  output_dir = NULL,
  save_results = FALSE
)
}
\arguments{
\item{data}{A data frame containing cell types, clock measures, age, and optional
control covariates and stratification variables.}

\item{study}{A character string specifying the study name, used for naming output files.}

\item{cell_types}{A character vector of column names in \code{data} that contain cell type measurements.}

\item{all_clocks}{A character vector of column names in \code{data} that contain clock measurements.}

\item{age_column}{A character string specifying the name of the age column in \code{data}.}

\item{stratify_by}{A character vector of column names in \code{data} to use for stratification.
Default is NULL (no stratification).}

\item{control_covariates}{A character vector of column names in \code{data} to use as additional control
covariates (Batch, Array, Plate, Row, Column etc.) in the regression models. Default is \code{NULL} (no additional control covariates).}

\item{output_dir}{A character string specifying the output directory. If NULL, the current
working directory is used. Default is NULL.}

\item{save_results}{Logical; if TRUE, results are saved as a CSV file. Default is FALSE.}
}
\value{
A data frame with columns:
\item{Clock_Variable}{Name of the clock measure}
\item{Cell_Variable}{Name of the cell type}
\item{Beta}{Regression coefficient}
\item{Lower_CI}{Lower bound of the 95\% confidence interval for Beta}
\item{Upper_CI}{Upper bound of the 95\% confidence interval for Beta}
\item{Categorical_Variable}{Name of the stratification variable (if applicable)}
\item{Categorical_Level}{Level of the stratification variable (if applicable)}
}
\description{
This function performs regressions with each clock variable against each cell variable,
including specified control covariates and optionally stratified by the levels of specified
categorical variables. It calculates the beta coefficients and confidence intervals
for each correlation.
}
\examples{
\dontrun{
result <- calculate_univariate_correlations(
  data = my_data,
  study = "my_study",
  cell_types = c("Neutrophils", "Monocytes"),
  all_clocks = c("GrimAge", "PhenoAge"),
  age_column = "Age",
  stratify_by = c("Sex", "SmokingStatus"),
  control_covariates = c("BMI", "Batch"),
  save_results = TRUE
)
}

}
