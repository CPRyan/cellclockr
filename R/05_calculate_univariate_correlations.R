#' Run Regressions with Clocks and Cells and Save Results
#'
#' This function performs regressions with each clock variable against each cell variable,
#' including specified covariates, and calculates the confidence intervals for the cell estimates.
#'
#' @param cell_clock_df A data frame containing the data with cells, clocks, and other covariates.
#' @param control_covariates A character vector of covariates to include in the regression models.
#' @param study A character string specifying the study name, used for naming output files.
#'
#' @return A data frame with results including Beta estimates and confidence intervals.
#' @export



calculate_univariate_correlations <- function(cell_clock_df, control_covariates = 1, study) {

  # Ensure required libraries are available
  if (!requireNamespace("broom", quietly = TRUE)) {
    install.packages("broom")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  library(broom)
  library(dplyr)

  # Define cell and clock columns as in previous functions
  cell_columns <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv", "Eos", "Mono", "NK", "Treg", "Neu")
  clock_columns <- c("PCHorvath1", "PCPhenoAge", "PCGrimAge", "DunedinPACE")

  # Create necessary directories
  xfun::dir_create("cellclockR_output/Tables")

  # Define all columns for residuals
  all_columns <- c(clock_columns, cell_columns)


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

  ############################
  # Use lapply to run the regression and standardization for each column
  ############################
  standardized_residuals_list <- lapply(all_columns, function(column) {
    regression_and_standardization(column, cell_clock_df)
  })

  # Name the columns of residualized outcomes
  names(standardized_residuals_list) <- paste0(all_columns, "_resids")

  cell_clock_df <-bind_cols(cell_clock_df, as_tibble(standardized_residuals_list))

  ############################

  ############################
  # Create an empty data frame to store the results
  results_df <- data.frame(Clock_Variable = character(),
                           Cell_Variable = character(),
                           Beta = numeric(),
                           Lower_CI = numeric(),
                           Upper_CI = numeric(),
                           stringsAsFactors = FALSE)

  # Ensure all specified columns are present in the data
  missing_cols <- setdiff(c(clock_columns, cell_columns), names(cell_clock_df))
  if (length(missing_cols) > 0) {
    stop(paste("Error: The following columns are missing in the data:", paste(missing_cols, collapse = ", ")))
  }

  # Issue a warning if no control covariates are provided
  if (length(control_covariates) == 1 && control_covariates == 1) {
    cat("\n\nProceeding to calculate stratified univariate correlations WITHOUT CONTROLLING for any technical covariates (e.g. Plate, Array, Batch)\n\n")
  } else {
    cat("\n\nProceeding to calculate stratified univariate correlations CONTROLLING for technical covariates (e.g. Plate, Array, Batch)\n\n")
  }

  # Iterate through clock_columns
  for (clock_var in paste0(clock_columns, "_resids")) {
    # Iterate through cell_columns
    for (cell_var in paste0(cell_columns, "_resids")) {

      # Fit a linear regression model
      model <- lm(paste0(clock_var, " ~ ", cell_var, " + ", paste(control_covariates, collapse = " + ")), data = cell_clock_df)

      # Extract coefficients and confidence intervals using broom
      tidy_model <- tidy(model)

      # Filter the coefficients for the cell variable
      cell_coefficients <- filter(tidy_model, term == cell_var)

      # Calculate the standard error
      se <- cell_coefficients$std.error

      # Calculate the confidence intervals
      lower_ci <- cell_coefficients$estimate - 1.96 * se
      upper_ci <- cell_coefficients$estimate + 1.96 * se

      # Add the results to the results data frame
      results_df <- bind_rows(results_df, data.frame(Clock_Variable = clock_var,
                                                     Cell_Variable = cell_var,
                                                     Beta = cell_coefficients$estimate,
                                                     Lower_CI = lower_ci,
                                                     Upper_CI = upper_ci)) %>%
        mutate(across(where(is.numeric), round, 3))
    }
  }

  # Generate the CSV filename dynamically
  univariate_correlations_filename <- paste0("cellclockR_output/Tables/", study, "_clock_cell_univariate_correlations.csv")

  # Save the summary data to CSV
  write_csv(x = results_df, file = univariate_correlations_filename)

  return(results_df)
}
