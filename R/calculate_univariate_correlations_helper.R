#' Helper functions for calculate_univariate_correlations
#'
#' This script contains internal helper functions for the 
#' calculate_univariate_correlations main function. These functions
#' handle error checking, directory creation, data processing, and result compilation.
#'
#' @name calculate_univariate_correlations_helper
#' @keywords internal
NULL

# The raise_errors() function takes all of the user input as input, and raises all common errors
raise_errors_calculate_univariate_correlations <- function(data, study, cell_types, all_clocks, age_column, stratify_by, control_covariates, output_dir, save_results) {
  
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
      "Categorical variables ('stratify_by') must be a non-empty character vector or left as default (=NULL)" = is.null(stratify_by) || (is.character(stratify_by) && length(stratify_by) > 0),
      "'control_covariates' must be a non-empty character vector or left as default (=NULL)" = is.null(control_covariates) || (is.character(control_covariates) && length(control_covariates) > 0),
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
  
  
  if (is.character(control_covariates) && length(intersect(control_covariates, stratify_by)) > 0) {
    overlapping_vars <- intersect(control_covariates, stratify_by)
    stop(paste("Error: The following variables appear in both control_covariates and stratify_by, which is not allowed:",
               paste(overlapping_vars, collapse = ", ")))
  }
  
  if (!(age_column %in% colnames(data))) {
    stop(paste("The provided age column is not found in the dataframe."))
  }
  
  if (age_column %in% control_covariates){
    stop(paste("'control_covariates' should not contain 'age_column'."))
  }
  
  if (age_column %in% stratify_by){
    stop(paste("'stratify_by' should not contain 'age_column'."))
  }
  
  if (!all(all_clocks %in% colnames(data))) {
    missing_cols <- all_clocks[!all_clocks %in% colnames(data)]
    stop(paste("The following clocks are not found in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!all(cell_types %in% colnames(data))) {
    missing_cols <- cell_types[!cell_types %in% colnames(data)]
    stop(paste("The following cell types are not found in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!all(stratify_by %in% colnames(data))) {
    missing_cols <- stratify_by[!stratify_by %in% colnames(data)]
    stop(paste("The following categorical variables ('stratify_by') are not found in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!is.null(control_covariates) && !all(control_covariates %in% colnames(data))) {
    missing_cols <- control_covariates[!control_covariates %in% colnames(data)]
    stop(paste("The following control covariates are not found in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
  
  # Issue a warning if no control covariates are provided
  if (is.null(control_covariates)) {
    cat("\n\nProceeding to calculate stratified univariate correlations WITHOUT CONTROLLING for any technical covariates (e.g. Plate, Array, Batch)\n\n")
  } else {
    cat("\n\nProceeding to calculate stratified univariate correlations CONTROLLING for technical covariates (e.g. Plate, Array, Batch)\n\n")
  }
}



# Function to make output directories
create_output_directories_calculate_univariate_correlations <- function(output_dir = NULL, save_results = FALSE) {
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

# Regression and Standardization Function
regression_and_standardization_calculate_univariate_correlations <- function(column_name, data, age_column, all_columns, control_covariates) {
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

# Function to fit model and collect coefficients
# Used to collect both the stratified and unstratified results
fitModel_collectCoeffCIs_calculate_univariate_correlations <- function(data, cell_var, clock_var, control_covariates){
  # Fit a linear regression model
  model_formula <- reformulate(
    termlabels = c(cell_var, control_covariates),
    response = c(clock_var)
  )
  model <- lm(model_formula, data = data)
  
  # Extract coefficients and confidence intervals using broom
  tidy_model <- tidy(model)
  
  # Filter the coefficients for the cell variable
  cell_coefficients <- dplyr::filter(tidy_model, term == cell_var)
  
  # Calculate the standard error
  se <- cell_coefficients$std.error
  
  # Calculate the Wald confidence intervals
  cell_coefficients$lower_ci <- cell_coefficients$estimate - 1.96 * se
  cell_coefficients$upper_ci <- cell_coefficients$estimate + 1.96 * se
  
  return(cell_coefficients)
}


collectResults_stratified_calculate_univariate_correlations <- function(results_df, cell_clock_df, clock_columns, cell_columns, control_covariates, categorical_variables){
  
  # Iterate through categorical variables
  for (categorical_var in categorical_variables) {
    # Ensure that the categorical variable is treated as a factor
    cell_clock_df[[categorical_var]] <- as.factor(cell_clock_df[[categorical_var]])
    
    # Get unique non-NA levels of the categorical variable
    levels <- unique(na.omit(cell_clock_df[[categorical_var]]))
    
    # Skip categorical variables with fewer than two unique non-NA levels
    if (length(levels) < 2) {
      next
    }
    
    # Iterate through levels, clock variables, and cell variables
    for (level in levels) {
      for (clock_var in paste0(clock_columns, "_resids")) {
        for (cell_var in paste0(cell_columns, "_resids")) {
          
          # Subset the data for the current level of the categorical variable
          subset_data <- cell_clock_df %>% dplyr::filter(.data[[categorical_var]] == level)
          
          # Check if there are enough observations for the model
          if (nrow(subset_data) < 2) {
            next
          }
          
          # Fit a linear regression model for each clock variable, including all cell variables and covariates  
          cell_coefficients <- fitModel_collectCoeffCIs_calculate_univariate_correlations(subset_data, cell_var, clock_var, control_covariates)
          
          
          # Add the results to the results data frame
          results_df <- bind_rows(results_df, data.frame(Clock_Variable = clock_var,
                                                         Categorical_Variable = categorical_var,
                                                         Categorical_Level = level,
                                                         Cell_Variable = cell_var,
                                                         Beta = cell_coefficients$estimate,
                                                         Lower_CI = cell_coefficients$lower_ci,
                                                         Upper_CI = cell_coefficients$upper_ci)) %>%
            mutate(across(where(is.numeric), round, 3)) %>%
            arrange(desc(Clock_Variable), Cell_Variable, Categorical_Variable)
        }
      }
    }
  }
  
  return(results_df)
}

collectResults_unstratified_calculate_univariate_correlations <- function(results_df, cell_clock_df, clock_columns, cell_columns, control_covariates){
  
  # Iterate through clock_columns
  for (clock_var in paste0(clock_columns, "_resids")) {
    # Iterate through cell_columns
    for (cell_var in paste0(cell_columns, "_resids")) {
      
      cell_coefficients <- fitModel_collectCoeffCIs_calculate_univariate_correlations(cell_clock_df, cell_var, clock_var, control_covariates)
      
      # Add the results to the results data frame
      results_df <- bind_rows(results_df, data.frame(Clock_Variable = clock_var,
                                                     Cell_Variable = cell_var,
                                                     Beta = cell_coefficients$estimate,
                                                     Lower_CI = cell_coefficients$lower_ci,
                                                     Upper_CI = cell_coefficients$upper_ci)) %>%
        mutate(across(where(is.numeric), round, 3))
    }
  }
  return(results_df)
}

# Save outputs
save_outputs_calculate_univariate_correlations <- function(study, output_dirs, result, save_results = FALSE) {
  if (is.null(output_dirs) || length(output_dirs) == 0) {
    return(invisible(NULL))
  }
  
  if (save_results) {
    if (!is.null(output_dirs$results)) {
      filename <- paste0(study, "_clock_cell_univariate_correlations.csv")
      full_path <- file.path(output_dirs$results, filename)
      write.csv(result, file = full_path, row.names = FALSE)
      cat(sprintf("Saved table: %s\n", full_path))
    } else {
      warning("Table directory not found in output_dirs. Table was not saved.")
    }
  }
}
