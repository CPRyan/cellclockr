#' Helper functions for exposure_clock_association_with_without_cells
#'
#' This script contains internal helper functions for the 
#' exposure_clock_association_with_without_cells main function. These functions
#' handle error checking, directory creation, data processing, and result compilation.
#'
#' @name exposure_clock_association_with_without_cells_helper
#' @keywords internal
NULL

# Error Checking
raise_errors_exposure_clock_association_with_without_cells <- function(data, study, cell_types, all_clocks, age_column, categorical_variables, control_covariates, output_dir, save_results) {
  
  tryCatch({
    stopifnot(
      "Data was not input" = !is.null(data),
      "'data' must be a dataframe" = is.data.frame(data),
      "Study name required" = !is.null(study),
      "'study' must be a single string" = is.character(study) && length(study) == 1,
      "Clock columns were not input" = !is.null(all_clocks),
      "'all_clocks' must be a non-empty character vector" = is.character(all_clocks) && length(all_clocks) > 0,
      "Cell columns were not input" = !is.null(cell_types),
      "'cell_types' must be a non-empty character vector" = is.character(cell_types) && length(cell_types) > 0,
      "'age_column' must be specified" = !is.null(age_column),
      "'age_column' must be a single string" = is.character(age_column) && length(age_column) == 1,
      "Categorical variables were not input" = !is.null(categorical_variables),
      "Categorical variables must be a non-empty character vector" = is.character(categorical_variables) && length(categorical_variables) > 0,
      "'control_covariates' must be a non-empty character vector or left as default (=NULL)" = 
        (is.null(control_covariates)) || 
        (is.character(control_covariates) && length(control_covariates) > 0),
      "'output_dir' must be a string or left as default (=NULL)" = is.null(output_dir) || (is.character(output_dir) && length(output_dir) == 1),
      "'save_results' must be a single logical value or left as default (=FALSE)" = is.logical(save_results) && length(save_results) == 1
    )
    
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) {
        stop(paste("The specified output directory does not exist:", output_dir))
      }
    }
    
  }, error = function(e) {
    stop(e$message, call. = FALSE)
  })
  
  
  if (is.character(control_covariates) && length(intersect(control_covariates, categorical_variables)) > 0) {
    overlapping_vars <- intersect(control_covariates, categorical_variables)
    stop(paste("Error: The following variables appear in both control_covariates and categorical_variables, which is not allowed:",
               paste(overlapping_vars, collapse = ", ")))
  }
  
  if (!(age_column %in% colnames(data))) {
    stop(paste("The provided age column is not found in the dataframe."))
  }
  
  if (age_column %in% control_covariates){
    stop(paste("'control_covariates' should not contain 'age_column'."))
  }
  
  if (age_column %in% categorical_variables){
    stop(paste("'categorical_variables' should not contain 'age_column'."))
  }
  
  if (!all(all_clocks %in% colnames(data))) {
    missing_cols <- all_clocks[!all_clocks %in% colnames(data)]
    stop(paste("The following clocks are not found in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!all(cell_types %in% colnames(data))) {
    missing_cols <- cell_types[!cell_types %in% colnames(data)]
    stop(paste("The following cell types are not found in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!all(categorical_variables %in% colnames(data))) {
    missing_cols <- categorical_variables[!categorical_variables %in% colnames(data)]
    stop(paste("The following categorical variables are not found in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!is.null(control_covariates) && !all(control_covariates %in% colnames(data))) {
    missing_cols <- control_covariates[!control_covariates %in% colnames(data)]
    stop(paste("The following control covariates are not found in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
  
  # Issue a warning if no control covariates are provided
  if (is.null(control_covariates)) {
    cat("\n\nProceeding to calculate exposure clock associations WITHOUT CONTROLLING for any technical covariates (e.g. Plate, Array, Batch)\n\n")
  } else {
    cat("\n\nProceeding to calculate exposure clock associations CONTROLLING for technical covariates (e.g. Plate, Array, Batch)\n\n")
  }
}


# Create Output Directories
create_output_directories_exposure_clock_association_with_without_cells <- function(output_dir = NULL, save_results = FALSE) {
  # If output_dir is NULL, use the current working directory
  base_dir <- if (is.null(output_dir)) getwd() else output_dir
  
  # Initialize empty list to store directory paths
  dirs <- list()
  
  if (save_results) {
    results_dir <- file.path(base_dir, "cellclockR_output", "Tables")
    xfun::dir_create(results_dir)
    dirs$results <- results_dir
  }
  
  # Return the list of created directories (will be empty if no directories were created)
  return(dirs)
  
}

# regression_and_standardization
regression_and_standardization_exposure_clock_association_with_without_cells <- function(column_name, data, age_column, all_columns, control_covariates) {
  data <- data %>% 
    select(all_of(age_column), all_of(all_columns), any_of(control_covariates)) %>% 
    na.omit()
  
  # Perform the linear regression
  regression_model <- lm(reformulate(age_column, response = column_name), data = data)
  
  # Get the residuals and standardize them
  residuals <- as.vector(scale(resid(regression_model)))
  
  # Return the standardized residuals as a named vector
  return(residuals)
}

# Calculate and append results to the dataframe
append_results_exposure_clock_association_with_without_cells <- function(results_df, data, age_column, control_covariates, categorical_var, clock_var, cell_columns){
  
  no_cell_model_formula <- reformulate(
    termlabels=c(age_column, control_covariates, categorical_var),
    response = c(clock_var)
  )
  
  no_cell_model <- lm(no_cell_model_formula, data=data)
  tidy_model_no_cells <- tidy(no_cell_model)
  exposure_coefficients_no_cells <- tidy_model_no_cells %>% dplyr::filter(grepl(categorical_var, term))
  se_no_cells <- exposure_coefficients_no_cells$std.error
  lower_ci_no_cells <- exposure_coefficients_no_cells$estimate - 1.96 * se_no_cells
  upper_ci_no_cells <- exposure_coefficients_no_cells$estimate + 1.96 * se_no_cells
  
  with_cell_model_formula <- reformulate(
    termlabels=c(age_column, control_covariates, categorical_var, paste0(cell_columns, "_resids")),
    response = c(clock_var)
  )
  
  model_with_cells <- lm(with_cell_model_formula, data=data)
  tidy_model_with_cells <- tidy(model_with_cells)
  exposure_coefficients_with_cells <- tidy_model_with_cells %>% dplyr::filter(grepl(categorical_var, term))
  se_with_cells <- exposure_coefficients_with_cells$std.error
  lower_ci_with_cells <- exposure_coefficients_with_cells$estimate - 1.96 * se_with_cells
  upper_ci_with_cells <- exposure_coefficients_with_cells$estimate + 1.96 * se_with_cells
  
  # Add the results to the results data frame
  results_df <- dplyr::bind_rows(results_df, data.frame(Clock_Variable = clock_var,
                                                        Exposure = exposure_coefficients_no_cells$term,
                                                        Beta_no_cells = exposure_coefficients_no_cells$estimate,
                                                        Lower_CI_no_cells = lower_ci_no_cells,
                                                        Upper_CI_no_cells = upper_ci_no_cells,
                                                        Beta_cells = exposure_coefficients_with_cells$estimate,
                                                        Lower_CI_cells = lower_ci_with_cells,
                                                        Upper_CI_cells = upper_ci_with_cells)) %>%
    mutate(across(where(is.numeric), round, 3))
  
  return(results_df)
  
}

# Save outputs
save_outputs_exposure_clock_association_with_without_cells <- function(study, output_dirs, result, save_results = FALSE) {
  if (is.null(output_dirs) || length(output_dirs) == 0) {
    return(invisible(NULL))
  }
  
  if (save_results) {
    if (!is.null(output_dirs$results)) {
      filename <- paste0(study, "_exposure_clock_association_with_without_cells.csv")
      full_path <- file.path(output_dirs$results, filename)
      write.csv(result, file = full_path, row.names = FALSE)
      cat(sprintf("Saved table: %s\n", full_path))
    } else {
      warning("Table directory not found in output_dirs. Table was not saved.")
    }
  }
}