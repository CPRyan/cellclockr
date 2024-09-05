#' Perform Regression, Standardization, and Delta R2 Calculations
#'
#' This function performs the following steps:
#' 1. Creates standardized residuals for each specified column using an existing function.
#' 2. Runs linear regressions for each combination of clock and cell residuals, with specified covariates.
#' 3. Saves the regression results to a CSV file.
#'
#' @param data A data frame containing the data.
#' @param control_covariates A character vector of control covariates.
#' @param categorical_variables A character vector of categorical variables for regression.
#' @param study A character string for naming the output file.
#'
#' @return None (saves results to a CSV file).
#' @import tidyverse xfun patchwork broom tibble readr sjlabelled
#' @importFrom stats lm median na.omit quantile resid
#' @export
exposure_cell_clock_associations <- function(cell_clock_df, control_covariates = 1, categorical_variables, study) {

  # Create necessary directories
  xfun::dir_create("cellclockR_output/Tables")

  # Create a list of column names for cell and clock columns
  cell_columns <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv", "Eos", "Mono",  "NK", "Treg", "Neu")
  clock_columns <- c("PCHorvath1", "PCPhenoAge", "PCGrimAge", "DunedinPACE")
  all_columns <- c(clock_columns, cell_columns)


  # Create an empty data frame to store the results
  results_df <- data.frame(Cell_or_Clock = character(),
                           Exposure = character(),
                           Beta = numeric(),
                           Lower_CI = numeric(),
                           Upper_CI = numeric(),
                           stringsAsFactors = FALSE)

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


  # Check if all specified columns are present in the data
  missing_cols <- setdiff(c(clock_columns, cell_columns, categorical_variables), names(cell_clock_df))
  if (length(missing_cols) > 0) {
    stop(paste("Error: The following columns are missing in the data:", paste(missing_cols, collapse = ", ")))
  }

  # Issue a warning if no control covariates are provided
  if (length(control_covariates) == 1 && control_covariates == 1) {
    cat("\n\nProceeding to calculate stratified univariate correlations WITHOUT CONTROLLING for any technical covariates (e.g. Plate, Array, Batch)\n\n")
  } else {
    cat("\n\nProceeding to calculate stratified univariate correlations CONTROLLING for technical covariates (e.g. Plate, Array, Batch)\n\n")
  }


  # Iterate through categorical variables
  for (categorical_var in categorical_variables) {
    data_resids_nona <- cell_clock_df %>% filter(!is.na(!!sym(categorical_var)))

    # Iterate through the residualized columns
    for (cell_clock_column in paste0(all_columns, "_resids")) {

      # Fit a linear regression model
      model <- lm(paste0(cell_clock_column, " ~ ", "Age", " + ", paste(control_covariates, collapse = " + "), "+", categorical_var),
                  data = data_resids_nona)

      # Extract coefficients and confidence intervals using broom
      tidy_model <- tidy(model)

      # Filter the coefficients for the categorical variable
      exposure_coefficients <- tidy_model %>% filter(grepl(categorical_var, term))

      # Calculate the standard error
      se <- exposure_coefficients$std.error

      # Calculate the confidence intervals
      lower_ci <- exposure_coefficients$estimate - 1.96 * se
      upper_ci <- exposure_coefficients$estimate + 1.96 * se

      # Add the results to the results data frame
      results_df <- bind_rows(results_df, data.frame(Cell_or_Clock = cell_clock_column,
                                                     Exposure = exposure_coefficients$term,
                                                     Beta = exposure_coefficients$estimate,
                                                     Lower_CI = lower_ci,
                                                     Upper_CI = upper_ci)) %>%
        mutate(across(where(is.numeric), round, 3))
    }
  }

  # Generate the CSV filename dynamically
  exposure_cell_or_clock_assoc <- paste0("cellclockR_output/Tables/", study, "_exposure_cell_or_clock_associations.csv")

  # Save the summary data to CSV
  write_csv(x = results_df, file = exposure_cell_or_clock_assoc)

  # Clean up
  rm(results_df)
}

