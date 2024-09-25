#' Calculate Associations Between Exposures and Age-Adjusted Cell Types and Clock Measures
#'
#' This function calculates associations between categorical exposures and
#' age-adjusted cell types or clock measures, optionally controlling for covariates.
#'
#' @param data A data frame containing cell types, clock measures, age, categorical
#'   variables (exposures), and optional control covariates.
#' @param study A character string specifying the study name, used for naming output files.
#' @param cell_types A character vector of column names in `data` that contain cell type measurements.
#' @param all_clocks A character vector of column names in `data` that contain clock measurements.
#' @param age_column A character string specifying the name of the age column in `data`.
#' @param categorical_variables A character vector of column names in `data` representing
#'   the categorical variables (exposures) to be analyzed.
#' @param control_covariates A character vector of column names in `data` to use as additional control
#'   covariates (Batch, Array, Plate, Row, Column etc.) in the regression models. Default is `NULL` (no additional control covariates).
#' @param output_dir A character string specifying the output directory. If `NULL`, the current
#'   working directory is used. Default is `NULL.`
#' @param save_results Logical; if `TRUE`, results are saved as a CSV file. Default is `FALSE.`
#'
#' @return A data frame with columns:
#'   \item{Cell_or_Clock}{Name of the age-adjusted cell type or clock measure}
#'   \item{Exposure}{Name of the exposure variable and its level}
#'   \item{Beta}{Regression coefficient}
#'   \item{Lower_CI}{Lower bound of the 95% confidence interval}
#'   \item{Upper_CI}{Upper bound of the 95% confidence interval}
#'
#' @import dplyr tidyr
#' @importFrom broom tidy
#' @importFrom rlang sym
#' @importFrom xfun dir_create
#'
#' @examples
#' \dontrun{
#' result <- exposure_cell_clock_associations(
#'   data = my_data,
#'   study = "aging_study",
#'   cell_types = c("Neutrophils", "Monocytes"),
#'   all_clocks = c("GrimAge", "PhenoAge"),
#'   age_column = "Age",
#'   categorical_variables = c("Smoking", "Alcohol"),
#'   control_covariates = c("Sex", "BMI"),
#'   save_results = TRUE
#' )
#' }
#'
#' @export
exposure_cell_clock_associations <- function(data, study, cell_types, all_clocks, age_column, categorical_variables, control_covariates = NULL, output_dir=NULL, save_results=FALSE) {
  
  # Raise Errors
  raise_errors_exposure_cell_clock_associations(data, study, cell_types, all_clocks, age_column, categorical_variables, control_covariates, output_dir, save_results)
  
  # Create output directory if it doesn't exist
  if (save_results) {
    output_dirs <- create_output_directories_exposure_cell_clock_associations(output_dir=output_dir, save_results=save_results)
  } else {
    output_dirs <- NULL
  }
  
  
  all_columns <- c(all_clocks, cell_types)
  
  # Create an empty data frame to store results
  results_df <- data.frame(Cell_or_Clock = character(),
                           Exposure = character(),
                           Beta = numeric(),
                           Lower_CI = numeric(),
                           Upper_CI = numeric(),
                           stringsAsFactors = FALSE)
  
  # Create standardized residuals for each column
  standardized_residuals_list <- lapply(all_columns, function(column) {
    regression_and_standardization_exposure_cell_clock_associations(column, data, age_column, all_columns, control_covariates)
  })  
  
  # Name the columns of residualized outcomes
  names(standardized_residuals_list) <- paste0(all_columns, "_resids")
  
  # Add the standardized residuals to the data frame
  data <- cbind(data, as.data.frame(standardized_residuals_list))
  
  # Iterate through categorical variables
  for (categorical_var in categorical_variables) {
    data_resids_nona <- data %>% dplyr::filter(!is.na(!!sym(categorical_var)))
    
    for (cell_clock_column in paste0(all_columns, "_resids")) {
      
      # Append results
      results_df <- append_results_exposure_cell_clock_associations(results_df, data, age_column, control_covariates, categorical_var, cell_clock_column)
      
    }
  }
  
  # Save outputs
  save_outputs_exposure_cell_clock_associations(study, output_dirs, results_df, save_results = save_results)
  
  
  return(results_df)
}
