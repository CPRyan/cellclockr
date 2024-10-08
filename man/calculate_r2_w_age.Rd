% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_r2_w_age.R
\name{calculate_r2_w_age}
\alias{calculate_r2_w_age}
\title{Calculate R2 for Age, Accounting for Control Covariates}
\usage{
calculate_r2_w_age(
  data,
  study,
  cell_types,
  all_clocks,
  age_column,
  control_covariates,
  output_dir = NULL,
  save_results = FALSE
)
}
\arguments{
\item{data}{A data frame containing cell types, clock measures, age, and control covariates.}

\item{study}{A string specifying the study name, used for naming output files.}

\item{cell_types}{A character vector of column names in \code{data} that contain cell type estimates.}

\item{all_clocks}{A character vector of column names in \code{data} that contain clock estimates.}

\item{age_column}{A character string specifying the name of the age column in \code{data}.}

\item{control_covariates}{A character vector of column names in \code{data} to use as
additional control covariates (Batch, Array, Plate, Row, Column etc.) in the regression models.}

\item{output_dir}{A character string specifying the output directory. If \code{NULL}, the current
working directory is used. Default is \code{NULL.}}

\item{save_results}{Logical; if TRUE, results are saved as a CSV file. Default is FALSE.}
}
\value{
A data frame with columns:
\item{Cell_or_Clock_Outcome}{Name of the cell type or clock measure}
\item{R2_with_age}{R-squared value for the model with age only}
\item{R2_w_age_and_control_covariates}{R-squared value for the model with age and control covariates}
\item{Delta_R2}{Difference between the two R-squared values}
\item{n_obs}{Number of complete observations used in the model}
}
\description{
This function calculates R-squared values for the relationship between age and
various cell type and clock measures, both with and without control covariates.
It produces a table with the R2 values for a basic model (age only), a full model
(age plus control covariates), and the difference between them (delta-R2).
}
\examples{
\dontrun{
result <- calculate_r2_w_age(
  data = data,
  study = "my_study",
  cell_types = c("Neutrophils", "Monocytes"),
  all_clocks = c("PCGrimAge", "PCPhenoAge"),
  age_column = "Age",
  control_covariates = c("Sex", "BMI"),
  save_results = TRUE
)
}

}
