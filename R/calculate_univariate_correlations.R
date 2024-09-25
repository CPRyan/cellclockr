#' Calculate Univariate Correlations between Cell Types and Clock Measures
#'
#' This function performs regressions with each clock variable against each cell variable,
#' including specified control covariates and optionally stratified by the levels of specified
#' categorical variables. It calculates the beta coefficients and confidence intervals
#' for each correlation.
#'
#' @param data A data frame containing cell types, clock measures, age, and optional
#'   control covariates and stratification variables.
#' @param study A character string specifying the study name, used for naming output files.
#' @param cell_types A character vector of column names in `data` that contain cell type measurements.
#' @param all_clocks A character vector of column names in `data` that contain clock measurements.
#' @param age_column A character string specifying the name of the age column in `data`.
#' @param stratify_by A character vector of column names in `data` to use for stratification.
#'   Default is NULL (no stratification).
#' @param control_covariates A character vector of column names in `data` to use as additional control
#'   covariates (Batch, Array, Plate, Row, Column etc.) in the regression models. Default is `NULL` (no additional control covariates).
#' @param output_dir A character string specifying the output directory. If NULL, the current
#'   working directory is used. Default is NULL.
#' @param save_results Logical; if TRUE, results are saved as a CSV file. Default is FALSE.
#'
#' @return A data frame with columns:
#'   \item{Clock_Variable}{Name of the clock measure}
#'   \item{Cell_Variable}{Name of the cell type}
#'   \item{Beta}{Regression coefficient}
#'   \item{Lower_CI}{Lower bound of the 95% confidence interval for Beta}
#'   \item{Upper_CI}{Upper bound of the 95% confidence interval for Beta}
#'   \item{Categorical_Variable}{Name of the stratification variable (if applicable)}
#'   \item{Categorical_Level}{Level of the stratification variable (if applicable)}
#'
#' @import dplyr tidyr
#' @importFrom broom tidy
#' @importFrom rlang .data
#' @importFrom xfun dir_create
#'
#' @examples
#' \dontrun{
#' result <- calculate_univariate_correlations(
#'   data = my_data,
#'   study = "my_study",
#'   cell_types = c("Neutrophils", "Monocytes"),
#'   all_clocks = c("GrimAge", "PhenoAge"),
#'   age_column = "Age",
#'   stratify_by = c("Sex", "SmokingStatus"),
#'   control_covariates = c("BMI", "Batch"),
#'   save_results = TRUE
#' )
#' }
#'
#' @export

calculate_univariate_correlations <- function(data, study, cell_types, all_clocks, age_column, stratify_by = NULL, control_covariates = NULL, output_dir=NULL, save_results=FALSE) {
  
  # Error-Checking
  raise_errors_calculate_univariate_correlations(data, study, cell_types, all_clocks, age_column, stratify_by, control_covariates, output_dir, save_results)
  
  # Create output directory
  if (save_results) {
    output_dirs <- create_output_directories_calculate_univariate_correlations(output_dir=output_dir, save_results=save_results)
  } else {
    output_dirs <- NULL
  }
  
  # Define all columns for residuals
  all_columns <- c(all_clocks, cell_types)
  
  # Create an empty data frame to store the results
  results_df <- data.frame(Clock_Variable = character(),
                           Categorical_Variable = character(),
                           Categorical_Level = character(),
                           Cell_Variable = character(),
                           Beta = numeric(),
                           Lower_CI = numeric(),
                           Upper_CI = numeric(),
                           stringsAsFactors = FALSE)
  
  # Standardize residuals for each column
  standardized_residuals_list <- lapply(all_columns, function(column) {
    regression_and_standardization_calculate_univariate_correlations(column, data, age_column, all_columns, control_covariates)
  })  
  
  # Name the columns of residualized outcomes
  names(standardized_residuals_list) <- paste0(all_columns, "_resids")
  
  cell_clock_df <-bind_cols(data, as_tibble(standardized_residuals_list))
  
  # If no control covariates are provided, set =1 to fit the regressions
  if (is.null(control_covariates)){
    control_covariates = 1
  }
  
  # Calculate coefficient estimates and Wald interval across all cell-clock pairs
  # And across categorical variables, if provided
  if (is.null(stratify_by)) {
    results_df <- collectResults_unstratified_calculate_univariate_correlations(results_df, cell_clock_df, all_clocks, cell_types, control_covariates)
    results_df <- results_df %>% select(-Categorical_Level, -Categorical_Variable)
  } else {
    results_df <- collectResults_stratified_calculate_univariate_correlations(results_df, cell_clock_df, all_clocks, cell_types, control_covariates, stratify_by)
  }
  
  # Save outputs
  save_outputs_calculate_univariate_correlations(study, output_dirs, results_df, save_results = save_results)
  
  return(results_df)
}