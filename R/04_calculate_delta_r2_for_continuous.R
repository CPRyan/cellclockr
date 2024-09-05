#' Calculate Delta R2 for Continuous Variables
#'
#' This function fits a basic model that uses only user defined `control_covariates` (e.g. Array, Batch) to predict clock or cell outcomes and takes the R2. It then fits another model that adds one of the user defined `continuous_variables` to the basic model, and takes the R2 from this model. It produces a table with the basic model R2, the full model R2, and the difference between them (delta-R2) for each of the user defined `continuous_variables`.
#'
#'
#'
#' @param data A data frame containing the data with cells, clock, control covariates, and continuous variables of interest. Must include a variable `Age`, used for residualization of cells and clocks on age before model fitting.
#' @param continuous_variables A character vector of continuous variables to test for associations with cells and clocks. These variables should be present in `data`.
#' @param control_covariates A character vector of covariates to include in the regression models.
#' @param study A character string specifying the study name, used for naming output files. This argument is required.
#'
#' @return A table with the basic model R2, the full model R2, and the difference between them (delta-R2) for each of the user defined `continuous_variables`. Also includes number of observations.
#' @import tidyverse xfun patchwork broom tibble sjlabelled
#' @importFrom stats lm median na.omit quantile resid
#' @export
calculate_delta_r2_for_continuous <- function(data = data, study = your_study_name, continuous_variables = your_continuous_variables, control_covariates = 1) {

  # Check that "Age" is included in the data
  if (!"Age" %in% colnames(data)) {
    stop("Error: 'Age' must be present in the data. 'age' or 'subject_age' etc. will not work.")
  }

  # Check that all continuous variables are numeric
  if (!all(sapply(continuous_variables, function(var) is.numeric(data[[var]])))) {
    stop("Error: All continuous variables must be numeric.")
  }

  # Create necessary directories
  xfun::dir_create("cellclockR_output/Figures")
  xfun::dir_create("cellclockR_output/Tables")

  # Cell columns
  cell_columns <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv", "Eos", "Mono", "NK", "Treg", "Neu")
  # Clock columns
  clock_columns <- c("PCHorvath1", "PCPhenoAge", "PCGrimAge", "DunedinPACE")
  # All columns
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





  # Standardize residuals for each column
  standardized_residuals_list <- lapply(all_columns, function(column) {
    regression_and_standardization(column, data)
  })

  # Name the columns of residualized outcomes
  names(standardized_residuals_list) <- paste0(all_columns, "_resids")

  # Add residuals to the original data frame
  data <- bind_cols(data, as_tibble(standardized_residuals_list))

  # Create an empty data frame to store the results
  results_df <- data.frame(Cell_or_Clock_Outcome = character(),
                           Continuous_Variable = character(),
                           R2_no_continuous_variable = numeric(),
                           R2_w_continuous_variable = numeric(),
                           Delta_R2 = numeric(),
                           n_obs = numeric(),
                           stringsAsFactors = FALSE)

  # For each of the cell and clock outcomes
  for (key_var in paste0(all_columns, "_resids")) {
    # Fit the basic model without the continuous variable
    basic_model <- lm(paste0(key_var, " ~ ", paste(control_covariates, collapse = " + ")), data = data)
    basic_model_r2 <- glance(basic_model)$adj.r.squared

    # Then, fit a model with the continuous variable, looping over each one at a time
    for (continuous_var in continuous_variables) {
      # Fit the model with the continuous variable
      full_model <- lm(paste0(key_var, " ~ ", continuous_var, " + ", paste(control_covariates, collapse = " + ")), data = data)
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
  }

  # Generate the CSV filename dynamically
  continuous_correlations_filename <- paste0("cellclockR_output/Tables/", study, "_clock_cell_continuous_delta_R2.csv")

  # Save the summary data to CSV
  write_csv(x = results_df, file = continuous_correlations_filename)

  return(results_df)
}
