#' Calculate R2 for Age, Accounting for Control Covariates
#'
#' This function fits a basic model correlates `Age` with each cell and clock measure. It then fits another model that adds user defined `control_covariates` to the basic model, and takes the R2 from this model. It produces a table with the basic `Age` model R2, the full model R2, and the difference between them (delta-R2).
#'
#'
#'
#' @param data A data frame containing the data with cells, clock, control covariates. Must include a variable `Age`.
#' #' @param study A character string specifying the study name, used for naming output files. This argument is required.
#' @param control_covariates A character vector of covariates to include in the regression models.
#' @return A table with the basic model R2 (including only Age), the full model R2 (that includes Age and user defined control covariates), and the difference between them (delta-R2). Also includes number of complete observations.
#' @export
calculate_r2_w_age <- function(data = data, study = your_study_name, control_covariates = 1) {
  # Ensure required libraries are available
  if (!requireNamespace("xfun", quietly = TRUE)) {
    install.packages("xfun")
  }
  if (!requireNamespace("broom", quietly = TRUE)) {
    install.packages("broom")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    install.packages("tibble")
  }
  library(xfun)
  library(broom)
  library(dplyr)
  library(tibble)

  # Check that "Age" is included in the data
  if (!"Age" %in% colnames(data)) {
    stop("Error: 'Age' must be present in the data. 'age' or 'subject_age' etc. will not work.")
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

  # Create an empty data frame to store the results
  results_df <- data.frame(Cell_or_Clock_Outcome = character(),
                           R2_with_age = numeric(),
                           R2_w_age_and_control_covariates = numeric(),
                           Delta_R2 = numeric(),
                           n_obs = numeric(),
                           stringsAsFactors = FALSE)

  # remove missing rows from key variables
  data <- data %>% select(Age, all_of(all_columns),
                    any_of(control_covariates)) %>% na.omit()

  # For each of the cell and clock outcomes
  for (key_var in all_columns) {
    # Fit the basic model without the continuous variable
    basic_age_model <- lm(paste0(key_var, " ~ ", "Age"), data = data)
    basic_age_model_r2 <- glance(basic_age_model)$adj.r.squared
    full_model <- lm(paste0(key_var, " ~ ", "Age", " + ", paste(control_covariates, collapse = " + ")), data = data)
    full_model_r2 <- glance(full_model)$adj.r.squared
    nobs <- glance(full_model)$nobs


    results_df <- dplyr::bind_rows(results_df, data.frame(Cell_or_Clock_Outcome = key_var,
                                                          R2_with_age = basic_age_model_r2,
                                                          R2_w_age_and_control_covariates = full_model_r2,
                                                          Delta_R2 = full_model_r2 - basic_age_model_r2,
                                                          n_obs = nobs) %>%
                                     mutate(across(where(is.numeric), round, 3))
      )
    }

  # Generate the CSV filename dynamically
  age_correlations_filename <- paste0("cellclockR_output/Tables/", study, "_age_correlations_w_wo_control_covar.csv")

  # Save the summary data to CSV
  write_csv(x = results_df, file =   age_correlations_filename)

  return(results_df)
}
