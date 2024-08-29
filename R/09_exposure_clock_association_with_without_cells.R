#' Perform Stratified Regression with and without Cell Variables and Save Results
#'
#' This function performs regression analysis on clock variables with and without including cell variables,
#' stratified by categorical variables, and saves the results to a CSV file.
#'
#' @param data A data frame containing the data.
#' @param control_covariates A character vector of control covariates.
#' @param categorical_variables A character vector of categorical variables for stratification.
#' @param study A character string for naming the output file.
#'
#' @return None (saves results to a CSV file).
#' @export
exposure_clock_association_with_without_cells <- function(cell_clock_df, control_covariates = 1, categorical_variables, study) {
  # Define clock and cell columns
  cell_columns <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv", "Eos", "Mono",  "NK", "Treg", "Neu")
  clock_columns <- c("PCHorvath1", "PCPhenoAge", "PCGrimAge", "DunedinPACE")

  # Create necessary directories
  xfun::dir_create("cellclockR_output/Tables")

  #####################################################################################
  # Function to perform the regression and standardization
  #####################################################################################
  regression_and_standardization <- function(column_name, data) {
    data <-data %>% select(Age, all_of(all_columns),
                           any_of(control_covariates)) %>% na.omit()
    # Perform the linear regression
    regression_model <- lm(data[[column_name]] ~ Age, data = data)
    # Get the residuals and standardize them
    residuals <- as.vector(scale(resid(regression_model)))
    # Return the standardized residuals as a named vector
    return(residuals)
  }


  # Create standardized residuals for each column
  standardized_residuals_list <- lapply(all_columns, function(column) {
    regression_and_standardization(column, cell_clock_df)
  })

  # Name the columns of residualized outcomes
  names(standardized_residuals_list) <- paste0(all_columns, "_resids")

  # Add the standardized residuals to the data frame
  cell_clock_df <- cbind(cell_clock_df, as.data.frame(standardized_residuals_list))


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


  ############################
  # Run function
  ############################

  # Check if all specified columns are present in the data
  missing_cols <- setdiff(c(clock_columns, cell_columns, categorical_variables), names(cell_clock_df))
  if (length(missing_cols) > 0) {
    stop(paste("Error: The following columns are missing in the data:", paste(missing_cols, collapse = ", ")))
  }

  # Issue a warning if no control covariates are provided
  if (length(control_covariates) == 1 && control_covariates == 1) {
    cat("\n\nProceeding to calculate exposure clock associations WITHOUT CONTROLLING for any technical covariates (e.g. Plate, Array, Batch)\n\n")
  } else {
    cat("\n\nProceeding to calculate exposure clock associations CONTROLLING for technical covariates (e.g. Plate, Array, Batch)\n\n")
  }


  # Iterate through categorical variables to create subsets
  for (categorical_var in categorical_variables) {

    # Filter out missing data for the categorical variable
    data_resids_nona <- cell_clock_df %>% filter(!is.na(!!sym(categorical_var)))

    # Iterate through clock variables
    for (clock_var in paste0(clock_columns, "_resids")) {

      # Fit the model WITHOUT cells
      model_no_cells <- lm(paste0(clock_var, " ~ Age + ", paste(control_covariates, collapse = " + "), " + ", categorical_var),
                           data = data_resids_nona)

      # Extract coefficients and confidence intervals
      tidy_model_no_cells <- tidy(model_no_cells)
      exposure_coefficients_no_cells <- tidy_model_no_cells %>% filter(grepl(categorical_var, term))
      se_no_cells <- exposure_coefficients_no_cells$std.error
      lower_ci_no_cells <- exposure_coefficients_no_cells$estimate - 1.96 * se_no_cells
      upper_ci_no_cells <- exposure_coefficients_no_cells$estimate + 1.96 * se_no_cells

      # Fit the model WITH cells
      model_with_cells <- lm(paste0(clock_var, " ~ Age + ", paste(control_covariates, collapse = " + "), " + ",
                                    paste(paste0(cell_columns, "_resids"), collapse = " + "), " + ", categorical_var),
                             data = data_resids_nona)

      # Extract coefficients and confidence intervals
      tidy_model_with_cells <- tidy(model_with_cells)
      exposure_coefficients_with_cells <- tidy_model_with_cells %>% filter(grepl(categorical_var, term))
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
    }
  }

  # Generate the CSV filename dynamically
  exposure_clock_w_wo_cells <- paste0("cellclockR_output/Tables/", study, "_exposure_clock_association_with_without_cells.csv")

  # Save the summary data to CSV
  write_csv(x = results_df, file = exposure_clock_w_wo_cells)

  # Clean up
  rm(results_df)
}
