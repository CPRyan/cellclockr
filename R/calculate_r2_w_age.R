#' Calculate R2 for Age, Accounting for Control Covariates
#'
#' This function calculates R-squared values for the relationship between age and 
#' various cell type and clock measures, both with and without control covariates. 
#' It produces a table with the R2 values for a basic model (age only), a full model 
#' (age plus control covariates), and the difference between them (delta-R2).
#'
#' @param data A data frame containing cell types, clock measures, age, and control covariates. 
#' @param study A string specifying the study name, used for naming output files.
#' @param cell_types A character vector of column names in `data` that contain cell type estimates. 
#' @param all_clocks A character vector of column names in `data` that contain clock estimates.
#' @param age_column A character string specifying the name of the age column in `data`. 
#' @param control_covariates A character vector of column names in `data` to use as 
#' additional control covariates (Batch, Array, Plate, Row, Column etc.) in the regression models.
#' @param output_dir A character string specifying the output directory. If `NULL`, the current 
#' working directory is used. Default is `NULL.` 
#' @param save_results Logical; if TRUE, results are saved as a CSV file. Default is FALSE.
#' 
#' @return A data frame with columns:
#'   \item{Cell_or_Clock_Outcome}{Name of the cell type or clock measure}
#'   \item{R2_with_age}{R-squared value for the model with age only}
#'   \item{R2_w_age_and_control_covariates}{R-squared value for the model with age and control covariates}
#'   \item{Delta_R2}{Difference between the two R-squared values}
#'   \item{n_obs}{Number of complete observations used in the model}
#' 
#' @import dplyr
#' @importFrom broom glance
#' @importFrom xfun dir_create
#'
#' @examples
#' \dontrun{
#' result <- calculate_r2_w_age(
#'   data = data,
#'   study = "my_study",
#'   cell_types = c("Neutrophils", "Monocytes"),
#'   all_clocks = c("PCGrimAge", "PCPhenoAge"),
#'   age_column = "Age",
#'   control_covariates = c("Sex", "BMI"),
#'   save_results = TRUE
#' )
#' }
#'
#' @export

calculate_r2_w_age <- function(data, study, cell_types, all_clocks, age_column, control_covariates, output_dir=NULL, save_results=FALSE) {
  
  raise_errors_calculate_r2_w_age(data, study, cell_types, all_clocks, age_column, control_covariates, output_dir, save_results)
  
  # Setting up the output directory
  if (save_results) {
    output_dirs <- create_output_directories_calculate_r2_w_age(output_dir=output_dir, save_results=save_results)
  } else {
    output_dirs <- NULL
  }
  
  # Create an empty data frame to store the results
  results_df <- data.frame(Cell_or_Clock_Outcome = character(),
                           R2_with_age = numeric(),
                           R2_w_age_and_control_covariates = numeric(),
                           Delta_R2 = numeric(),
                           n_obs = numeric(),
                           stringsAsFactors = FALSE)
  
  # Remove missing rows from key variables
  all_columns <- c(cell_types, all_clocks)
  data <- data %>% select(Age, all_of(all_columns),
                          any_of(control_covariates)) %>% na.omit()
  
  # For each of the cell and clock outcomes
  for (key_var in all_columns) {
    
    # Append results
    results_df <- append_results_calculate_r2_w_age(results_df, data, all_columns, key_var, age_column, control_covariates, all_clocks)
    
  }
  
  # Order the results
  results_df <- results_df %>%
    arrange(desc(Is_Clock), Cell_or_Clock_Outcome)%>%
    select(-Is_Clock)
  
  # Save outputs
  save_outputs_calculate_r2_w_age(study, output_dirs, results_df, save_results = save_results)
  
  return(results_df)
}