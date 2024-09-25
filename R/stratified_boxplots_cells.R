#' Generate Stratified Boxplots for Cell Types
#'
#' This function creates boxplots for cell types and a summary table of boxplot
#' parameters for the entire dataset. If categorical variables are provided,
#' the output will be stratified by levels in each categorical variable.
#'
#' @param data A data frame containing cell type data and optional categorical variables.
#' @param id A string specifying the column name for participant/subject IDs.
#' @param study A string specifying the study name, used for naming output files.
#' @param cell_types A character vector of column names present in `data` that contain the cell type data.
#' @param highlighted_cell_types A character vector of cell types to be highlighted
#' with an adjusted y-axis scale. Must be a subset of `cell_types`. Default is NULL
#' where no cell types are highlighted. Neutrophil cell counts are typically highlighted 
#' as they greatly exceed those of other cell types.
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
#' result <- stratified_boxplots_cells(
#'   data = data,
#'   id = "participant_id",
#'   study = "my_study",
#'   cell_types = c("Neutrophils", "Basophils", "Monocytes"),
#'   highlighted_cell_types = c("Neutrophils")
#'   categorical_variables = c("Sex", "AgeGroup"),
#'   save_plots = TRUE,
#'   save_summaries = TRUE
#' )
#' }
#' 
#' @export
stratified_boxplots_cells <- function(data, id, study, cell_types, highlighted_cell_types=NULL, categorical_variables = NULL, colors=NULL, output_dir=NULL, save_plots=FALSE, save_summaries=FALSE) {
  
  # Raise errors
  raise_errors_stratified_boxplots_cells(data, id, study, cell_types, highlighted_cell_types, categorical_variables, colors, output_dir, save_plots, save_summaries)
  
  # Set output directory as specified
  if (save_plots || save_summaries) {
    output_dirs <- create_output_directories_stratified_boxplots_cells(output_dir, save_plots, save_summaries)
  } else {
    output_dirs <- NULL
  }
  
  # Set colors as specified, ensure the color palette has sufficient coverage
  colors <- assign_color_palette_stratified_boxplots_cells(categorical_variables, data, colors)
  
  # Transform and subset data for the highlighted and non-highlighted cells
  long_data <- transform_and_subset_data_stratified_boxplots_cells(data, categorical_variables, id, cell_types, highlighted_cell_types)
  highlighted_long <- long_data$highlighted_long
  non_highlighted_long <- long_data$non_highlighted_long
  all_long <- long_data$all_long # Used in the summaries
  
  plots <- list()
  summaries <- list()
  
  # Generate the overall plot and summary
  plots$overall <- generate_cell_boxplots_stratified_boxplots_cells(highlighted_long, non_highlighted_long, colors, variable=NULL) +
    plot_annotation(title = "Overall Cell Type Distribution")
  summaries$overall <- generate_summary_stratified_boxplots_cells(all_long)
  
  # Generate boxplots and summaries for each categorical variable
  if (!is.null(categorical_variables) && length(categorical_variables) > 0) {
    for (variable in categorical_variables) {
      plots[[variable]] <- generate_cell_boxplots_stratified_boxplots_cells(highlighted_long, non_highlighted_long, colors, variable)
      summaries[[variable]] <- generate_summary_stratified_boxplots_cells(all_long, variable)
    }
  }
  
  # Display plots
  print(plots)
  
  # Save outputs, dependent on user choice
  save_outputs_stratified_boxplots_cells(study, output_dirs, plots, summaries, save_plots, save_summaries)
  
  return(list(plots = plots, summaries = summaries))
}