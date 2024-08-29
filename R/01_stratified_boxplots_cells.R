#' Stratified Boxplots for Cells
#'
#' This function generates boxplots for cells and a summary table of boxplot parameters for the whole dataset. If categorical variables/factors are provided under `categorical_variables`, output will also be provided stratified for levels in each categorical variable.
#'
#' @param data A data frame containing the data with cells and any categorical variables that will be used to create groups for the boxplots.
#' @param categorical_variables A character vector of categorical variables on which to stratify the data for boxplots. These variables must be present in `data`. Default using `NULL` produces plots data without groupings.
#' @param id A character string specifying the column name in `data` that contains participant or subject IDs (e.g. "id", "participant", etc.). This argument is required.
#' @param study A character string specifying the study name, used for naming output files (e.g. "framingham", "whi", etc.). This argument is required.
#'
#' @return A list of ggplot objects containing the generated boxplots for each categorical variable. The function also saves summary data frames as CSV files in the "Cells_Clocks_Output/Tables/" directory.
#' @export
stratified_boxplots_cells <- function(data = your_data, id = your_id, study = your_study_name, categorical_variables = NULL) {

  # Check if ID is provided
  if (is.null(id)) {
    stop("Error: Participant/Subject ID required.")
  }

  # Check if Study name is provided
  if (is.null(study)) {
    stop("Error: Study name required.")
  }

  # Load necessary libraries
  if (!requireNamespace("xfun", quietly = TRUE)) {
    install.packages("xfun")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  library(xfun)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)

  # Create necessary directories
  xfun::dir_create("cellclockR_output/Figures")
  xfun::dir_create("cellclockR_output/Tables")

  # Set default color palette
  cols <- c("#875692FF", "#F38400FF", "#A1CAF1FF", "#BE0032FF", "#C2B280FF",  "#008856FF",
            "#E68FACFF", "#0067A5FF", "#F99379FF", "#604E97FF", "#F6A600FF", "#B3446CFF",
            "#DCD300FF", "#882D17FF", "#8DB600FF", "#654522FF", "#E25822FF", "#2B3D26FF")

  # Add "Null" to the data to get a Null group
  categorical_variables_w_null <- c(categorical_variables, "Null")

  # Organize the Data
  all_long <- data %>%
    mutate(Null = "Null") %>%
    mutate(across(all_of(categorical_variables_w_null), as.factor)) %>%
    select(all_of(id), all_of(categorical_variables_w_null),
           Bas, Bmem, Bnv, CD4mem, CD4nv, CD8mem, CD8nv, Eos, Mono, NK, Treg, Neu) %>%
    pivot_longer(cols = c(Bas, Bmem, Bnv, CD4mem, CD4nv, CD8mem, CD8nv, Eos, Mono, NK, Treg, Neu),
                 names_to = "Cell Type", values_to = "Cell Count")

  # Create a list to store the plots
  plots <- list()

  # Generate and save boxplots for each categorical variable
  for (variable in categorical_variables_w_null) {
    # Split data into non-neutrophils and neutrophils
    non_neu_long <- all_long %>%
      filter(`Cell Type` != "Neu", !is.na(!!sym(variable)))

    neu_long <- all_long %>%
      filter(`Cell Type` == "Neu", !is.na(!!sym(variable)))

    # Extract levels for the variable
    non_neu_levels <- levels(factor(non_neu_long[[variable]]))
    neu_levels <- levels(factor(neu_long[[variable]]))

    # Plot non-neutrophil data
    non_neu_plot <- non_neu_long %>%
      ggplot() +
      stat_boxplot(aes(x = `Cell Type`, y = `Cell Count`,
                       color = !!sym(variable)),
                   geom = "errorbar", lwd = 1) +
      geom_boxplot(aes(x = `Cell Type`, y = `Cell Count`,
                       color = !!sym(variable)),
                   outlier.shape = NA, lwd = 1) +
      theme_bw() +
      labs(x = "Cell Types", y = "Cell Counts") +
      theme(legend.position = "top") +
      scale_x_discrete(labels = levels(as_factor(non_neu_long$`Cell Type`))) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.6, hjust = 0.5)) +
      theme(axis.title.x = element_blank()) +
      scale_color_manual(values = cols)

    # Plot neutrophil data
    neu_plot <- neu_long %>%
      ggplot() +
      stat_boxplot(aes(x = `Cell Type`, y = `Cell Count`,
                       color = !!sym(variable)),
                   geom = "errorbar", lwd = 1) +
      geom_boxplot(aes(x = `Cell Type`, y = `Cell Count`,
                       color = !!sym(variable)),
                   outlier.shape = NA, lwd = 1) +
      theme_bw() +
      labs(x = "Cell Types", y = "Cell Counts") +
      theme(legend.position = "top") +
      scale_x_discrete(labels = "Neutrophils") +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.6, hjust = 0.5)) +
      theme(axis.title.x = element_blank()) +
      scale_color_manual(values = cols)

    # Combine the plots
    full_plot <- non_neu_plot + neu_plot + theme(legend.position = "none") + plot_layout(ncol = 2, widths = c(9.3, 1))

    # Save the plot using the variable name as part of the filename
    filename <- paste0(study, "_cells_", variable, "_boxplots.png")
    ggsave(filename = paste0("cellclockR_output/Figures/", filename), plot = full_plot)

    # Extract and save summary data
    non_neu_summary <- non_neu_long %>%
      group_by(`Cell Type`, !!sym(variable)) %>%
      summarise(ymin = min(`Cell Count`),
                lower = quantile(`Cell Count`, 0.25),
                middle = median(`Cell Count`),
                upper = quantile(`Cell Count`, 0.75),
                ymax = max(`Cell Count`), .groups = 'drop')

    neu_summary <- neu_long %>%
      group_by(`Cell Type`, !!sym(variable)) %>%
      summarise(ymin = min(`Cell Count`),
                lower = quantile(`Cell Count`, 0.25),
                middle = median(`Cell Count`),
                upper = quantile(`Cell Count`, 0.75),
                ymax = max(`Cell Count`), .groups = 'drop') %>%
      mutate(`Cell Type` = "Neutrophils")

    # Combine summary data
    combined_summary <- bind_rows(non_neu_summary, neu_summary)

    # Generate the CSV filename dynamically
    summary_filename <- paste0("cellclockR_output/Tables/", study, "_cells_", variable, "_summary.csv")

    # Save the summary data to CSV
    write_csv(x = combined_summary %>% mutate(across(where(is.numeric), round, 3)), file = summary_filename)

    # Store the plot in the list
    plots[[variable]] <- full_plot
  }

  return(plots)
}
