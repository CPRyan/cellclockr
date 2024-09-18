#' Helper functions for calculate_r2_w_age
#'
#' This script contains internal helper functions for the 
#' calculate_r2_w_age main function. These functions
#' handle error checking, directory creation, data processing, and result compilation.
#'
#' @name calculate_r2_w_age_helper
#' @keywords internal
NULL

# The raise_errors() function takes all of the user input as input, and raises all common errors
raise_errors_calculate_r2_w_age <- function(data, study, cell_types, all_clocks, age_column, control_covariates, output_dir, save_results) {
  
  tryCatch({
    stopifnot(
      "Data was not input" = !is.null(data),
      "'data' must be a dataframe" = is.data.frame(data),
      "Study name required" = !is.null(study),
      "'study' must be a single string" = is.character(study) && length(study) == 1,
      "Age column specification is required" = !is.null(age_column),
      "'age_column' must be a single string" = is.character(age_column) && length(age_column) == 1,
      "Clock columns were not input" = !is.null(all_clocks),
      "'all_clocks' must be a non-empty character vector" = is.character(all_clocks) && length(all_clocks) > 0,
      "Cell columns were not input" = !is.null(cell_types),
      "'cell_types' must be a non-empty character vector" = is.character(cell_types) && length(cell_types) > 0,
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
  
  if (!(age_column %in% colnames(data))) {
    stop(paste("The provided age column is not found in the dataframe."))
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
  
  if (control_covariates != 1 && !all(control_covariates %in% colnames(data))) {
    missing_cols <- control_covariates[!control_covariates %in% colnames(data)]
    stop(paste("The following control covariates are not found in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
}

## This creates the folder where the user specifies
## The directories are stored to be called later on
create_output_directories_calculate_r2_w_age <- function(output_dir = NULL, save_results = FALSE) {
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

# Calculate and append Results
append_results_calculate_r2_w_age <- function(results_df, data, all_columns, key_var, age_column, control_covariates, all_clocks) {
    
  # Fit basic model, collect adj r2
  basic_age_model_formula <- reformulate(
    termlabels=c(age_column),
    response = c(key_var)
  )
  basic_age_model <- lm(basic_age_model_formula, data=data)
  basic_age_model_r2 <- glance(basic_age_model)$adj.r.squared
  
  # Fit full model, collect adj r2
  full_model_formula <- reformulate(
    termlabels=c(age_column, control_covariates),
    response = c(key_var)
  )
  full_model <- lm(full_model_formula, data=data)
  full_model_r2 <- glance(full_model)$adj.r.squared
  nobs <- glance(full_model)$nobs
  
  # Append Results
  results_df <- dplyr::bind_rows(results_df, data.frame(Cell_or_Clock_Outcome = key_var,
                                                        R2_with_age = basic_age_model_r2,
                                                        R2_w_age_and_control_covariates = full_model_r2,
                                                        Delta_R2 = full_model_r2 - basic_age_model_r2,
                                                        n_obs = nobs,
                                                        Is_Clock = key_var %in% all_clocks) %>%
                                   mutate(across(where(is.numeric), round, 3))
  )

  return(results_df)
}

# Saving Output
save_outputs_calculate_r2_w_age <- function(study, output_dirs, result, save_results = FALSE) {
  if (is.null(output_dirs) || length(output_dirs) == 0) {
    return(invisible(NULL))
  }
  
  if (save_results) {
    if (!is.null(output_dirs$results)) {
      filename <- paste0(study, "_age_correlations_w_wo_control_covar.csv")
      full_path <- file.path(output_dirs$results, filename)
      write.csv(result, file = full_path, row.names = FALSE)
      cat(sprintf("Saved table: %s\n", full_path))
    } else {
      warning("Table directory not found in output_dirs. Table was not saved.")
    }
  }
}
