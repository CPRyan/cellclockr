# 1. stratified_boxplots_cells (01A)
# 2. stratified_boxplots_clocks (01A)
# 3. calculate_r2_w_age (new)
# 3. calculate_delta_r2_for_continuous (01B)
# 4. calculate_univariate_correlations (02A)
# 5. calculate_univariate_correlations_stratified (02B)
# 6. calculate_cell_clock_delta_r2 (03)
# 7. exposure_cell_clock_associations (04)
# exposure_cell_clock_associations_stratified (05)


pacman::p_load(tidyverse, readr)
clhns <-read_csv(here::here("../clhns_ics_complete_df.csv"))

# 1.
stratified_boxplots_cells(cell_clock_df = clhns, id = 'uncchdid', study = 'clhns', categorical_variables = c("Female", 'currently_pregnant'))

# 2.
stratified_boxplots_clocks(cell_clock_df = clhns, id = 'uncchdid', study = 'clhns', categorical_variables = c("Female", 'currently_pregnant'))

# 3.
calculate_r2_w_age(data = clhns, study = "clhns")

# 4.
calculate_delta_r2_for_continuous(cell_clock_df = clhns, continuous_variables = c("Age", "PCCystatinC", "PCDNAmTL"), study = "clhns", control_covariates = c("Array"))

# 5.
calculate_univariate_correlations(cell_clock_df = clhns, study = 'clhns', control_covariates = c('Sample_Plate', 'Array'))

# 6.
calculate_univariate_correlations_stratified(cell_clock_df = clhns, study = 'clhns', control_covariates = c('Sample_Plate', 'Array'), categorical_variables = c("Female", "currently_pregnant"))

# 7.
calculate_cell_clock_delta_r2(cell_clock_df = clhns, study = 'clhns', categorical_variables = c("Female","currently_pregnant"), control_covariates = c('Sample_Plate', 'Array'))

# 8.
exposure_cell_clock_associations(cell_clock_df = clhns, study = 'clhns', categorical_variables = c("Female","currently_pregnant"), control_covariates = c('Sample_Plate', 'Array'))

# 9.
exposure_clock_association_with_without_cells(cell_clock_df = clhns, study = 'clhns', categorical_variables = c("Female","currently_pregnant"))



