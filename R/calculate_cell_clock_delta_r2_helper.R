#' Helper functions for calculate_cell_clock_delta_r2
#'
#' This script contains internal helper functions for the 
#' calculate_cell_clock_delta_r2 main function. These functions
#' handle error checking, directory creation, data processing, and result compilation.
#'
#' @name calculate_cell_clock_delta_r2_helper
#' @keywords internal
NULL

# Error-Checking
raise_errors_calculate_cell_clock_delta_r2 <- function(data, study, cell_types, all_clocks, age_column, stratify_by, control_covariates, output_dir, save_results) {
  
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
    cat("\n\nProceeding to calculate cell/clock R2 (for all and stratified) WITHOUT CONTROLLING for any technical covariates (e.g. Plate, Array, Batch)\n\n")
  } else {
    cat("\n\nProceeding to calculate cell/clock R2 (for all and stratified) CONTROLLING for technical covariates (e.g. Plate, Array, Batch)\n\n")
  }
}

# Create Output Directories
create_output_directories_calculate_cell_clock_delta_r2 <- function(output_dir = NULL, save_results = FALSE) {
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
regression_and_standardization_calculate_cell_clock_delta_r2 <- function(column_name, data, age_column, all_columns, control_covariates) {
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

# Fits model, appends unstratified results
appendUnstratifiedResults_calculate_cell_clock_delta_r2 <- function(results_df, data, clock_columns, cell_columns, 
                                      age_column, control_covariates=NULL, categorical_var=NULL, level=NULL){
  
  # Fit a linear regression model for each clock variable, including all cell variables and covariates
  for (clock_var in paste0(clock_columns, "_resids")) {
    
    # # Check if there are enough observations for the model
    # if (nrow(data) < 2) {
    #   next
    # }
    
    cell_columns <- as.character(cell_columns)
    if (is.null(control_covariates)) {
      termlabels_basic <- c(age_column)
      termlabels_full <- c(age_column, paste0(cell_columns, "_resids"))
    } else {
      termlabels_basic <- c(age_column, control_covariates)
      termlabels_full <- c(age_column, paste0(cell_columns, "_resids"), control_covariates)
    }
    
    # Formulate the models
    basic_formula <- reformulate(
      termlabels = termlabels_basic,
      response = clock_var
    )
    
    full_formula <- reformulate(
      termlabels = termlabels_full,
      response = clock_var
    )
    
    # Fit the model WITHOUT CELLS
    basic_model <- lm(basic_formula, data = data)
    
    # Extract R2 value
    basic_model_r2 <- glance(basic_model)$adj.r.squared
    
    # Fit the model WITH CELLS
    full_model <- lm(full_formula, data = data)
    
    # Extract R2 value
    full_model_r2 <- glance(full_model)$adj.r.squared
    nobs <- glance(full_model)$nobs
    
    if(is.null(categorical_var)) {
      categorical_var <- "Full Dataset"
      level <- "Full Dataset"
    }
    
    # Add the results to the results data frame
    results_df <- bind_rows(results_df, data.frame(
      Clock_Variable = clock_var,
      Categorical_Variable = categorical_var,
      Categorical_Level = level,
      R2_no_cells = basic_model_r2,
      R2_cells = full_model_r2,
      Delta_R2 = full_model_r2 - basic_model_r2,
      n_obs = nobs
    ) %>% mutate(across(where(is.numeric), round, 3)))
  }
  
  return(results_df)
  
}

# Fits model, appends stratified results
appendStratifiedResults_calculate_cell_clock_delta_r2 <- function(results_df, data, clock_columns, cell_columns, 
                                    age_column, control_covariates, categorical_variables){
  
  for (categorical_var in categorical_variables) {
    
    # Ensure that the categorical variable is treated as a factor
    data[[categorical_var]] <- as.factor(data[[categorical_var]])
    
    # Get unique non-NA levels of the categorical variable
    levels <- unique(na.omit(data[[categorical_var]]))
    
    # Skip categorical variables with fewer than two unique non-NA levels
    if (length(levels) < 2) {
      next
    }
    
    # Append results for each level of each categorical variable
    for (level in levels) {
      # Subset data
      subset_data <- data[data[[categorical_var]] == level, ]
      
      if (nrow(data) < 2) {
        next
      }

      # Append Results
      results_df <- appendUnstratifiedResults_calculate_cell_clock_delta_r2(results_df, subset_data, 
                                              clock_columns, cell_columns, age_column, 
                                              control_covariates, categorical_var, level)
    }
  }
  
  return(results_df)
  
}

# Save outputs
save_outputs_calculate_cell_clock_delta_r2 <- function(study, output_dirs, result, save_results = FALSE) {
  if (is.null(output_dirs) || length(output_dirs) == 0) {
    return(invisible(NULL))
  }
  
  if (save_results) {
    if (!is.null(output_dirs$results)) {
      filename <- paste0(study, "_cell_clock_delta_r2.csv")
      full_path <- file.path(output_dirs$results, filename)
      write.csv(result, file = full_path, row.names = FALSE)
      cat(sprintf("Saved table: %s\n", full_path))
    } else {
      warning("Table directory not found in output_dirs. Table was not saved.")
    }
  }
}
