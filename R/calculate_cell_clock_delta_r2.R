#' Calculate Delta R-squared for Each Clock, With and Without Cells
#'
#' This function calculates the change in R-squared (delta R2) when adding
#' cell types to models predicting age-adjusted clock measures. It can stratify
#' results by categorical variables and control for specified covariates.
#'
#' @param data A data frame containing cell types, clock measures, age, and optional
#'   stratification variables and control covariates.
#' @param study A character string specifying the study name, used for naming output files.
#' @param cell_types A character vector of column names in `data` that contain cell type measurements.
#' @param all_clocks A character vector of column names in `data` that contain clock measurements.
#' @param age_column A character string specifying the name of the age column in `data`.
#' @param stratify_by A character vector of column names in `data` to use for stratification.
#'   Default is `NULL` (no stratification).
#' @param control_covariates A character vector of column names in `data` to use as control
#'   covariates in the regression models. Default is `NULL` (no control covariates).
#' @param output_dir A character string specifying the output directory. If `NULL`, the current
#'   working directory is used. Default is `NULL.`
#' @param save_results Logical; if `TRUE`, results are saved as a CSV file. Default is `FALSE.`
#'
#' @return A data frame with columns:
#'   \item{Clock_Variable}{Name of the clock measure}
#'   \item{Categorical_Variable}{Name of the stratification variable (if applicable)}
#'   \item{Categorical_Level}{Level of the stratification variable (if applicable)}
#'   \item{R2_no_cells}{R-squared value for the model without cell types}
#'   \item{R2_cells}{R-squared value for the model with cell types}
#'   \item{Delta_R2}{Difference between R2_cells and R2_no_cells}
#'   \item{n_obs}{Number of observations used in the model}
#'
#' @import dplyr tidyr
#' @importFrom broom glance
#' @importFrom rlang sym
#' @importFrom xfun dir_create
#'
#' @examples
#' \dontrun{
#' result <- calculate_cell_clock_delta_r2(
#'   data = my_data,
#'   study = "aging_study",
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
calculate_cell_clock_delta_r2<- function(data, study, cell_types, all_clocks, age_column, stratify_by = NULL, control_covariates = NULL, output_dir=NULL, save_results=FALSE) {
  
  # Error-Checking
  raise_errors_calculate_cell_clock_delta_r2(data, study, cell_types, all_clocks, age_column, stratify_by, control_covariates, output_dir, save_results)
  
  # Create output directory if it doesn't exist
  if (save_results) {
    output_dirs <- create_output_directories_calculate_cell_clock_delta_r2(output_dir=output_dir, save_results=save_results)
  } else {
    output_dirs <- NULL
  }
  
  # Combine All Columns
  all_columns <- c(all_clocks, cell_types)
  
  # Create standardized residuals for each column
  standardized_residuals_list <- lapply(all_columns, function(column) {
    regression_and_standardization_calculate_cell_clock_delta_r2(column, data, age_column, all_columns, control_covariates)
  })  
  
  # Name the columns of residualized outcomes
  names(standardized_residuals_list) <- paste0(all_columns, "_resids")
  
  # Add the standardized residuals to the data frame
  data <- cbind(data, as.data.frame(standardized_residuals_list))
  
  # Create an empty data frame to store the results
  results_df <- data.frame(
    Clock_Variable = character(),
    Categorical_Variable = character(),
    Categorical_Level = character(),
    R2_no_cells = numeric(),
    R2_cells = numeric(),
    Delta_R2 = numeric(),
    n_obs = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Append unstratified results
  results_df <- appendUnstratifiedResults_calculate_cell_clock_delta_r2(results_df, data, all_clocks, cell_types, 
                                          age_column, control_covariates=control_covariates, categorical_var=NULL, level=NULL)
  
  # Append stratified results, if categorical variables are provided
  if (!is.null(stratify_by)){
    results_df <- appendStratifiedResults_calculate_cell_clock_delta_r2(results_df, data, all_clocks, cell_types, 
                                          age_column, control_covariates=control_covariates, stratify_by)
  }
  
  # Save outputs
  save_outputs_calculate_cell_clock_delta_r2(study, output_dirs, results_df, save_results = save_results)
  
  
  return(results_df)
}