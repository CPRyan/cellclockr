#' Calculate Delta R-squared for Continuous Variables
#'
#'
#' This function calculates R-squared values for the relationship between 
#' user-defined `control_covariates` and age-adjusted cell type and clock measures, 
#' both with and without the additional of user-defined `continuous_variables`.
#' It produces a table with the R2 values for a basic model (`control_covariates` only), a full model 
#' (`control_covariates` plus `continuous_variables`), and the difference between them (delta-R2).
#'
#' @param data A data frame containing cell types, clock measures, age, continuous variables,
#'   and control covariates.
#' @param study A character string specifying the study name, used for naming output files.
#' @param cell_types A character vector of column names in `data` that contain cell type measurements.
#' @param all_clocks A character vector of column names in `data` that contain clock measurements.
#' @param continuous_variables A character vector of column names in `data` representing
#'   the continuous variables to be analyzed.
#' @param control_covariates A character vector of column names in `data` to be used as control covariates.
#' @param age_column A character string specifying the name of the age column in `data`.
#' @param output_dir A character string specifying the output directory. If `NULL`, the current
#'   working directory is used. Default is `NULL.`
#' @param save_results Logical; if TRUE, results are saved as a CSV file. Default is FALSE.
#'
#' @return A data frame with columns:
#'   \item{Cell_or_Clock_Outcome}{Name of the age-adjusted cell type or clock measure}
#'   \item{Continuous_Variable}{Name of the continuous variable added to the model}
#'   \item{R2_no_continuous_variable}{R-squared value for the model without the continuous variable}
#'   \item{R2_w_continuous_variable}{R-squared value for the model with the continuous variable}
#'   \item{Delta_R2}{Difference between the two R-squared values}
#'   \item{n_obs}{Number of complete observations used in the model}
#'
#' @import dplyr tidyr
#' @importFrom broom glance
#' @importFrom xfun dir_create
#'
#' @examples
#' \dontrun{
#' result <- calculate_delta_r2_for_continuous(
#'   data = data,
#'   study = "my_study",
#'   cell_types = c("Neutrophils", "Monocytes"),
#'   all_clocks = c("PCGrimAge", "PCPhenoAge"),
#'   continuous_variables = c("BMI", "Cholesterol"),
#'   control_covariates = c("Sex", "Smoking"),
#'   age_column = "Age",
#'   save_results = TRUE
#' )
#' }
#'
#' @export

calculate_delta_r2_for_continuous <- function(data, study, cell_types, all_clocks, continuous_variables, control_covariates, age_column, output_dir=NULL, save_results=FALSE) {
  
  # Raise errors
  raise_errors_calculate_delta_r2_for_continuous(data, study, cell_types, all_clocks, continuous_variables, control_covariates, age_column, output_dir, save_results)
  
  # Set up the output directory
  if (save_results) {
    output_dirs <- create_output_directories_calculate_delta_r2_for_continuous(output_dir=output_dir, save_results=save_results)
  } else {
    output_dirs <- NULL
  }
  
  # Create all_columns vector
  all_columns <- c(all_clocks, cell_types)
  
  # Standardize residuals for each column
  standardized_residuals_list <- lapply(all_columns, function(column) {
    regression_and_standardization_calculate_delta_r2_for_continuous(column, data, age_column, all_columns, control_covariates)
  })
  
  # Name the columns of residualized outcomes
  names(standardized_residuals_list) <- paste0(all_columns, "_resids")
  
  # Add residuals to the original data frame
  data <- bind_cols(data, as_tibble(standardized_residuals_list))
  
  # Create an empty data frame to store the results
  results_df <- data.frame(Cell_or_Clock_Outcome = character(),
                           Continuous_Variable = character(),
                           R2_no_continuous_variable = numeric(),
                           R2_w_continuous_variable = numeric(),
                           Delta_R2 = numeric(),
                           n_obs = numeric(),
                           stringsAsFactors = FALSE)
  
  # For each of the cell and clock outcomes
  for (key_var in paste0(all_columns, "_resids")) {
    # Calculate and append results
    results_df <- append_results_calculate_delta_r2_for_continuous(results_df, data, key_var, control_covariates, continuous_variables)
  }
  
  # Save outputs
  save_outputs_calculate_delta_r2_for_continuous(study, output_dirs, results_df, save_results = save_results)
  
  return(results_df)
}