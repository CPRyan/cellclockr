#' Calculate Exposure-Clock Associations With and Without Cell Type Adjustment
#'
#' This function calculates associations between categorical exposures and
#' age-adjusted clock measures, both with and without adjusting for cell types.
#' It optionally controls for additional covariates.
#'
#' @param data A data frame containing clock measures, cell types, age, categorical
#'   variables (exposures), and optional control covariates.
#' @param study A character string specifying the study name, used for naming output files.
#' @param cell_types A character vector of column names in `data` that contain cell type measurements.
#' @param all_clocks A character vector of column names in `data` that contain clock measurements.
#' @param age_column A character string specifying the name of the age column in `data`.
#' @param categorical_variables A character vector of column names in `data` representing
#'   the categorical variables (exposures) to be analyzed.
#' @param control_covariates A character vector of column names in `data` to use as control
#'   covariates in the regression models. Default is `NULL` (no control covariates).
#' @param output_dir A character string specifying the output directory. If `NULL`, the current
#'   working directory is used. Default is `NULL`.
#' @param save_results Logical; if `TRUE`, results are saved as a CSV file. Default is `FALSE`.
#'
#' @return A data frame with columns:
#'   \item{Clock_Variable}{Name of the age-adjusted clock measure}
#'   \item{Exposure}{Name of the exposure variable and its level}
#'   \item{Beta_no_cells}{Regression coefficient without cell type adjustment}
#'   \item{Lower_CI_no_cells}{Lower bound of the 95% CI without cell type adjustment}
#'   \item{Upper_CI_no_cells}{Upper bound of the 95% CI without cell type adjustment}
#'   \item{Beta_cells}{Regression coefficient with cell type adjustment}
#'   \item{Lower_CI_cells}{Lower bound of the 95% CI with cell type adjustment}
#'   \item{Upper_CI_cells}{Upper bound of the 95% CI with cell type adjustment}
#'
#' @import dplyr tidyr
#' @importFrom broom tidy
#' @importFrom rlang sym
#' @importFrom xfun dir_create
#'
#' @examples
#' \dontrun{
#' result <- exposure_clock_association_with_without_cells(
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
exposure_clock_association_with_without_cells <- function(data, study, cell_types, all_clocks, age_column, categorical_variables, control_covariates = NULL, output_dir = NULL, save_results = FALSE) {
  
  # Raise Errors
  raise_errors_exposure_clock_association_with_without_cells(data, study, cell_types, all_clocks, age_column, categorical_variables, control_covariates, output_dir, save_results)
  
  # Create output directory if it doesn't exist
  if (save_results) {
    output_dirs <- create_output_directories_exposure_clock_association_with_without_cells(output_dir=output_dir, save_results=save_results)
  } else {
    output_dirs <- NULL
  }
  
  # Create list of standardized residuals
  all_columns <- c(all_clocks, cell_types)
  standardized_residuals_list <- lapply(all_columns, function(column) {
    regression_and_standardization_exposure_clock_association_with_without_cells(column, data, age_column, all_columns, control_covariates)
  })  
  
  # Name the columns of residualized outcomes and add to dataframe
  names(standardized_residuals_list) <- paste0(all_columns, "_resids")
  data <- cbind(data, as.data.frame(standardized_residuals_list))
  
  # Create an empty data frame to store the results
  results_df <- data.frame(Clock_Variable = character(),
                           Exposure = character(),
                           Beta_no_cells = numeric(),
                           Lower_CI_no_cells = numeric(),
                           Upper_CI_no_cells = numeric(),
                           Beta_cells = numeric(),
                           Lower_CI_cells = numeric(),
                           Upper_CI_cells = numeric(),
                           stringsAsFactors = FALSE)
  
  # Iterate through categorical variables
  for (categorical_var in categorical_variables) {
    # Filter out missing data for the categorical variable
    data_resids_nona <- data %>% dplyr::filter(!is.na(!!sym(categorical_var)))
    # Iterate through clock variables
    for (clock_var in paste0(all_clocks, "_resids")) {
      
      # Calculate and append Results
      results_df <- append_results_exposure_clock_association_with_without_cells(results_df, data_resids_nona, age_column, control_covariates, categorical_var, clock_var, cell_types)
      
    }
  }
  
  
  # Save outputs
  save_outputs_exposure_clock_association_with_without_cells(study, output_dirs, results_df, save_results = save_results)
  
  return(results_df)
}