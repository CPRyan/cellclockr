#' Helper functions for calculate_delta_r2_for_continuous
#'
#' This script contains internal helper functions for the 
#' calculate_delta_r2_for_continuous main function. These functions
#' handle error checking, directory creation, data processing, and result compilation.
#'
#' @name calculate_delta_r2_for_continuous_helper
#' @keywords internal
NULL

# The raise_errors() function takes all of the user input as input, and raises all common errors
raise_errors_calculate_delta_r2_for_continuous <- function(data, study, cell_types, all_clocks, continuous_variables, control_covariates, age_column, output_dir, save_results) {
  
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
      "Age column specification is required" = !is.null(age_column),
      "'age_column' must be a single string" = is.character(age_column) && length(age_column) == 1,
      "Continuous variables were not input" = !is.null(continuous_variables),
      "Continuous variables must be a non-empty character vector" = is.character(continuous_variables) && length(continuous_variables) > 0,
      "'control_covariates were not input" = !is.null(control_covariates),
      "'control_covariates' must be a non-empty character vector" = is.character(control_covariates) && length(control_covariates) > 0,
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
  
  if (is.character(control_covariates) && length(intersect(control_covariates, continuous_variables)) > 0) {
    overlapping_vars <- intersect(control_covariates, continuous_variables)
    stop(paste("Error: The following variables appear in both control_covariates and continuous_variables, which is not allowed:",
               paste(overlapping_vars, collapse = ", ")))
  }
  
  if (!(age_column %in% colnames(data))) {
    stop(paste("The provided age column is not found in the dataframe."))
  }
  
  if (age_column %in% continuous_variables){
    stop(paste("'continuous_variables' should not contain 'age_column'."))
  }
  
  if (age_column %in% control_covariates){
    stop(paste("'control_covariates' should not contain 'age_column'."))
  }
  
  if (!all(all_clocks %in% colnames(data))) {
    missing_cols <- all_clocks[!all_clocks %in% colnames(data)]
    stop(paste("The following clocks are not found in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!all(cell_types %in% colnames(data))) {
    missing_cols <- cell_types[!cell_types %in% colnames(data)]
    stop(paste("The following cell types are not found in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!all(continuous_variables %in% colnames(data))) {
    missing_cols <- continuous_variables[!continuous_variables %in% colnames(data)]
    stop(paste("The following continuous variables are not found in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!all(control_covariates %in% colnames(data))) {
    missing_cols <- control_covariates[!control_covariates %in% colnames(data)]
    stop(paste("The following control covariates are not found in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
}

## This creates the folder where the user specifies
## And allows the directory to be called later on
create_output_directories_calculate_delta_r2_for_continuous <- function(output_dir = NULL, save_results = FALSE) {
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


# Regression and Standardization Function, modified to be external
regression_and_standardization_calculate_delta_r2_for_continuous <- function(column_name, data, age_column, all_columns, control_covariates) {
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

# Calculate and append results
append_results_calculate_delta_r2_for_continuous <- function(results_df, data, key_var, control_covariates, continuous_variables){
  
  basic_model_formula <- reformulate(
    termlabels = c(control_covariates),
    response = c(key_var)
  )
  basic_model <- lm(basic_model_formula, data=data)
  basic_model_r2 <- glance(basic_model)$adj.r.squared
  
  # Then, fit a model with the continuous variable, looping over each one at a time
  for (continuous_var in continuous_variables) {
    
    full_model_formula <- reformulate(
      termlabels = c(continuous_var, control_covariates),
      response = c(key_var)
    )
    full_model <- lm(full_model_formula, data=data)
    full_model_r2 <- glance(full_model)$adj.r.squared
    nobs <- glance(full_model)$nobs
    
    results_df <- dplyr::bind_rows(results_df, data.frame(Cell_or_Clock_Outcome = key_var,
                                                          Continuous_Variable = continuous_var,
                                                          R2_no_continuous_variable = basic_model_r2,
                                                          R2_w_continuous_variable = full_model_r2,
                                                          Delta_R2 = full_model_r2 - basic_model_r2,
                                                          n_obs = nobs) %>%
                                     mutate(across(where(is.numeric), round, 3))
    )
  }
  return(results_df)
}

# Save Output
save_outputs_calculate_delta_r2_for_continuous <- function(study, output_dirs, result, save_results = FALSE) {
  if (is.null(output_dirs) || length(output_dirs) == 0) {
    return(invisible(NULL))
  }
  
  if (save_results) {
    if (!is.null(output_dirs$results)) {
      filename <- paste0(study, "_clock_cell_continuous_delta_R2.csv")
      full_path <- file.path(output_dirs$results, filename)
      write.csv(result, file = full_path, row.names = FALSE)
      cat(sprintf("Saved table: %s\n", full_path))
    } else {
      warning("Table directory not found in output_dirs. Table was not saved.")
    }
  }
}