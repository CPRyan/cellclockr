#' Stratified Boxplots for Clocks
#'
#' This function generates boxplots for clocks and a summary table of boxplot parameters for the whole dataset. If categorical variables/factors are provided under `categorical_variables`, output will also be provided stratified for levels in each categorical variable.
#'
#' @param data A data frame containing the data with clocks and any categorical variables that will be used to create groups for the boxplots.
#' @param categorical_variables A character vector of categorical variables on which to stratify the data for boxplots. These variables must be present in `data`. Default using `NULL` produces plots data without groupings.
#' @param id A character string specifying the column name in `data` that contains participant or subject IDs (e.g. "id", "participant", etc.). This argument is required.
#' @param study A character string specifying the study name, used for naming output files (e.g. "framingham", "whi", etc.). This argument is required.
#'
#' @return A ggplot object containing the generated boxplots. The function also saves a summary data frame as a CSV file in the "Cells_Clocks_Output/Tables/" directory.
#' @export
stratified_boxplots_clocks <- function(data = your_data, id = your_id, study = your_study_name, categorical_variables = NULL) {

  # Check for required parameters
  if (missing(id) || is.null(id)) {
    stop("Error: Participant/Subject ID required.")
  }
  if (missing(study) || is.null(study)) {
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
           PCHorvath1, PCPhenoAge, PCGrimAge, DunedinPACE) %>%
    pivot_longer(cols = c(PCHorvath1, PCPhenoAge, PCGrimAge, DunedinPACE),
                 names_to = "Clock Type", values_to = "Clock Estimate")

  # Create a list to store the plots
  plots <- list()

  # Generate and save boxplots for each categorical variable
  for (variable in categorical_variables_w_null) {
    non_pace_long <- all_long %>%
      filter(`Clock Type` != "DunedinPACE", !is.na(!!sym(variable)))

    pace_long <- all_long %>%
      filter(`Clock Type` == "DunedinPACE", !is.na(!!sym(variable)))

    # Plot non-DunedinPACE data
    non_pace_plot <- non_pace_long %>%
      ggplot() +
      stat_boxplot(aes(x = `Clock Type`, y = `Clock Estimate`,
                       color = !!sym(variable)),
                   geom = "errorbar", lwd = 1) +
      geom_boxplot(aes(x = `Clock Type`, y = `Clock Estimate`,
                       color = !!sym(variable)),
                   outlier.shape = NA, lwd = 1) +
      theme_bw() +
      labs(x = "Clock", y = "Clock Estimate") +
      theme(legend.position = "top") +
      scale_x_discrete(labels = levels(as_factor(non_pace_long$`Clock Type`))) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.6, hjust = 0.5)) +
      theme(axis.title.x = element_blank()) +
      scale_color_manual(values = cols)

    # Plot DunedinPACE data
    pace_plot <- pace_long %>%
      ggplot() +
      stat_boxplot(aes(x = `Clock Type`, y = `Clock Estimate`,
                       color = !!sym(variable)),
                   geom = "errorbar", lwd = 1) +
      geom_boxplot(aes(x = `Clock Type`, y = `Clock Estimate`,
                       color = !!sym(variable)),
                   outlier.shape = NA, lwd = 1) +
      theme_bw() +
      labs(x = "Clock", y = "Clock Estimate") +
      theme(legend.position = "top") +
      scale_x_discrete(labels = "DunedinPACE") +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.6, hjust = 0.5)) +
      theme(axis.title.x = element_blank()) +
      scale_color_manual(values = cols)

    # Combine the plots
    full_plot <- non_pace_plot + pace_plot + theme(legend.position = "none") + plot_layout(ncol = 2, widths = c(2.9, 1))

    # Save the plot using the variable name as part of the filename
    filename <- paste0(study, "_clocks_", variable, "_boxplots.png")
    ggsave(filename = paste0("cellclockR_output/Figures/", filename), plot = full_plot)

    # Extract and save summary data
    non_pace_summary <- non_pace_long %>%
      group_by(`Clock Type`, !!sym(variable)) %>%
      summarise(ymin = min(`Clock Estimate`),
                lower = quantile(`Clock Estimate`, 0.25),
                middle = median(`Clock Estimate`),
                upper = quantile(`Clock Estimate`, 0.75),
                ymax = max(`Clock Estimate`), .groups = 'drop')

    pace_summary <- pace_long %>%
      group_by(`Clock Type`, !!sym(variable)) %>%
      summarise(ymin = min(`Clock Estimate`),
                lower = quantile(`Clock Estimate`, 0.25),
                middle = median(`Clock Estimate`),
                upper = quantile(`Clock Estimate`, 0.75),
                ymax = max(`Clock Estimate`), .groups = 'drop')
    # Combine summary data
    combined_summary <- bind_rows(non_pace_summary, pace_summary)


    # Generate the CSV filename dynamically
    summary_filename <- paste0("cellclockR_output/Tables/", study, "_clocks_", variable, "_summary.csv")

    # Save the summary data to CSV
    write_csv(x = combined_summary %>% mutate(across(where(is.numeric), round, 3)), file = summary_filename)

    # Store the plot in the list
    plots[[variable]] <- full_plot
  }

  return(plots)
}
