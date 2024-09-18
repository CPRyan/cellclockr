#' Generate Stratified Boxplots for Clocks
#'
#' This function creates boxplots for clocks and a summary table of boxplot
#' parameters for the entire dataset. If categorical variables are provided,
#' the output will be stratified by levels in each categorical variable.
#'
#' @param data A data frame containing clock data and optional categorical variables.
#' @param id A string specifying the column name for participant/subject IDs.
#' @param study A string specifying the study name, used for naming output files.
#' @param all_clocks A character vector of column names present in `data` that contain the clock data.
#' @param highlighted_clock_types A character vector of clocks to be highlighted
#' with an adjusted y-axis scale. Must be a subset of `all_clocks`. Default is NULL
#' where no clocks are highlighted.
#' @param categorical_variables A character vector of categorical variables for
#' stratification. These must be present in `data`. Default using `NULL` produces 
#' data without groupings.
#' @param colors A character vector of colors that can be used to replace 
#' the default color palette. If an insufficient number of colors are provided for the
#' levels of a categorical variable, default colors or a color ramp will be used.
#' @param output_dir A string specifying the output directory. If `NULL`, the current
#'   working directory is used.
#' @param save_plots Logical; if TRUE, plots are saved in 'cellclockR_output/Plots/'.
#' @param save_summaries Logical; if TRUE, summaries are saved in 'cellclockR_output/Summaries/'.
#'
#' @return A list containing the ggplot objects and summary data frames for each categorical variable.
#' 
#' @import dplyr tidyr ggplot2 patchwork
#' @importFrom xfun dir_create
#' 
#' @examples
#' \dontrun{
#' result <- stratified_boxplots_clocks(
#'   data = data,
#'   id = "participant_id",
#'   study = "my_study",
#'   all_clocks = c("DunedinPACE", "PCPhenoAge", "PCGrimAge"),
#'   highlighted_clocks = c("DunedinPACE")
#'   categorical_variables = c("Sex", "AgeGroup"),
#'   save_plots = TRUE,
#'   save_summaries = TRUE
#' )
#' }
#' 
#' @export

stratified_boxplots_clocks <- function(data, id, study, all_clocks, highlighted_clocks=NULL, categorical_variables = NULL, colors=NULL, output_dir=NULL, save_plots=FALSE, save_summaries=FALSE) {
  
  # Raise Errors
  raise_errors_stratified_boxplots_clocks(data, id, study, all_clocks, highlighted_clocks, categorical_variables, colors, output_dir, save_plots, save_summaries)
  
  # Create output directories, if required
  if (save_plots || save_summaries) {
    output_dirs <- create_output_directories_stratified_boxplots_clocks(output_dir, save_plots, save_summaries)
  } else {
    output_dirs <- NULL
  }
  
  # Assign color palette
  colors <- assign_color_palette_stratified_boxplots_clocks(categorical_variables, data, colors)
  
  # Transform data for plotting and summary generation
  long_data <- transform_and_subset_data_stratified_boxplots_clocks(data, categorical_variables, id, all_clocks, highlighted_clocks)
  highlighted_long <- long_data$highlighted_long
  non_highlighted_long <- long_data$non_highlighted_long
  all_long <- long_data$all_long # Used in the summaries
  
  # Create a list to store the plots and summaries
  plots <- list()
  summaries <- list()
  
  # Generate the overall plot and summary
  plots$overall <- generate_clock_boxplots_stratified_boxplots_clocks(highlighted_long, non_highlighted_long, colors, variable=NULL) +
    plot_annotation(title = "Overall Clock Estimates")
  summaries$overall <- generate_summary_stratified_boxplots_clocks(all_long)
  
  # Generate boxplots and summaries for each categorical variable
  if (!is.null(categorical_variables) && length(categorical_variables) > 0) {
    for (variable in categorical_variables) {
      plots[[variable]] <- generate_clock_boxplots_stratified_boxplots_clocks(highlighted_long, non_highlighted_long, colors, variable)
      summaries[[variable]] <- generate_summary_stratified_boxplots_clocks(all_long, variable)
    }
  }
  
  # Display plots
  print(plots)
  
  # Save outputs
  save_outputs_stratified_boxplots_clocks(study, output_dirs, plots, summaries, save_plots = save_plots, save_summaries = save_summaries)
  
}