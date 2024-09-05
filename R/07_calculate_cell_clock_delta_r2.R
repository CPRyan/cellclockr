#' Calculate Delta R2 for Each Clock with and without Cells
#'
#' This function calculates the Delta R2 for each clock variable, comparing models with and without cell variables.
#' It also handles categorical variables by stratifying the analysis and calculating Delta R2 for each level of the categorical variable.
#'
#' @param data A data frame containing the data.
#' @param control_covariates A character vector of control covariates.
#' @param categorical_variables A character vector of categorical variables to stratify by.
#' @param study A character string specifying the study name, used for naming output files.
#'
#' @return A data frame with Delta R2 results and saves the results to a CSV file.
#' @import tidyverse xfun patchwork broom tibble readr sjlabelled
#' @importFrom stats lm median na.omit quantile resid
#' @export
calculate_cell_clock_delta_r2<- function(cell_clock_df, control_covariates = 1, categorical_variables, study) {

  # Create output directory if it does not exist
  output_dir <- "cellclockR_output/Tables/"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Define column names
  clock_columns <- c("PCHorvath1", "PCPhenoAge", "PCGrimAge", "DunedinPACE")
  cell_columns <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv", "Eos", "Mono", "NK", "Treg", "Neu")
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


  # Create standardized residuals for each column
  standardized_residuals_list <- lapply(all_columns, function(column) {
    regression_and_standardization(column, cell_clock_df)
  })

  # Name the columns of residualized outcomes
  names(standardized_residuals_list) <- paste0(all_columns, "_resids")

  # Add the standardized residuals to the data frame
  cell_clock_df <- cbind(cell_clock_df, as.data.frame(standardized_residuals_list))

  # Create an empty data frame to store the results
  results_df <- data.frame(
    Clock_Variable = character(),
    Categorical_Variable = character(),
    Categorical_Level = character(),
    R2_no_cells = numeric(),
    R2_cells = numeric(),
    Delta_R2 = numeric(),
    n_obs = numeric(),
    stringsAsFactors = FALSE
  )


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
    cat("\n\nProceeding to calculate cell/clock R2 (for all and stratified) WITHOUT CONTROLLING for any technical covariates (e.g. Plate, Array, Batch)\n\n")
  } else {
    cat("\n\nProceeding to calculate cell/clock R2 (for all and stratified) CONTROLLING for technical covariates (e.g. Plate, Array, Batch)\n\n")
  }


  # Basic model not stratified
  for (clock_var in paste0(clock_columns, "_resids")) {
    # Fit the model WITHOUT CELLS
    basic_model <- lm(paste0(clock_var, " ~ Age + ", paste(control_covariates, collapse = " + ")), data = cell_clock_df)

    # Extract R2 value
    basic_model_r2 <- glance(basic_model)$adj.r.squared

    # Fit the model WITH CELLS
    full_model <- lm(paste0(clock_var, " ~ Age + ", paste(paste0(cell_columns, "_resids"), collapse = " + "), " + ", paste(control_covariates, collapse = " + ")), data = cell_clock_df)

    # Extract R2 value
    full_model_r2 <- glance(full_model)$adj.r.squared
    nobs <- glance(full_model)$nobs

    # Add the results to the results data frame
    results_df <- bind_rows(results_df, data.frame(
      Clock_Variable = clock_var,
      Categorical_Variable = "Full Dataset",
      Categorical_Level = "Full Dataset",
      R2_no_cells = basic_model_r2,
      R2_cells = full_model_r2,
      Delta_R2 = full_model_r2 - basic_model_r2,
      n_obs = nobs
    ) %>% mutate(across(where(is.numeric), round, 3)))
  }

  # Stratified models
  for (categorical_var in categorical_variables) {
    # Ensure that the categorical variable is treated as a factor
    cell_clock_df[[categorical_var]] <- as.factor(cell_clock_df[[categorical_var]])

    # Get unique non-NA levels of the categorical variable
    levels <- unique(na.omit(cell_clock_df[[categorical_var]]))

    # Skip categorical variables with fewer than two unique non-NA levels
    if (length(levels) < 2) {
      next
    }

    # Iterate through levels
    for (level in levels) {
      # Fit a linear regression model for each clock variable, including all cell variables and covariates
      for (clock_var in paste0(clock_columns, "_resids")) {
        # Subset the data for the current level of the categorical variable
        subset_data <- cell_clock_df %>% filter(.data[[categorical_var]] == level)

        # Check if there are enough observations for the model
        if (nrow(subset_data) < 2) {
          next
        }

        # Fit the model WITHOUT CELLS
        basic_model <- lm(paste0(clock_var, " ~ Age + ", paste(control_covariates, collapse = " + ")), data = subset_data)

        # Extract R2 value
        basic_model_r2 <- glance(basic_model)$adj.r.squared

        # Fit the model WITH CELLS
        full_model <- lm(paste0(clock_var, " ~ Age + ", paste(paste0(cell_columns, "_resids"), collapse = " + "), " + ", paste(control_covariates, collapse = " + ")), data = subset_data)

        # Extract R2 value
        full_model_r2 <- glance(full_model)$adj.r.squared
        nobs <- glance(full_model)$nobs

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
    }
  }

  # Save results to CSV
  delta_R2_filename <- paste0(output_dir, study, "cell_clock_delta_r2.csv")
  write_csv(x = results_df, file = delta_R2_filename)

  # Clean up
  rm(results_df)
}
