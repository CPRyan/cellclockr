#' Helper functions for stratified_boxplots_cells
#'
#' This script contains internal helper functions for the 
#' stratified_boxplots_cells main function. These functions
#' handle error checking, directory creation, data processing, and result compilation.
#'
#' @name stratified_boxplots_cells_helper
#' @keywords internal
NULL

# The raise_errors() function takes all of the user input as input, and raises all common errors
raise_errors_stratified_boxplots_cells <- function(data, id, study, cell_types, highlighted_cell_types, categorical_variables, colors, output_dir, save_plots, save_summaries) {
  
  tryCatch({
    stopifnot(
      "Data was not input" = !is.null(data),
      "'data' must be a dataframe" = is.data.frame(data),
      "Participant/Subject ID required" = !is.null(id),
      "ID column not found in the dataframe" = id %in% colnames(data),
      "Study name required" = !is.null(study),
      "'study' must be a single string" = is.character(study) && length(study) == 1,
      "Cell type columns were not input" = !is.null(cell_types),
      "'cell_types' must be a non-empty character vector" = is.character(cell_types) && length(cell_types) > 0,
      "highlighted_cell_types must be NULL or a non-empty character vector" = 
        is.null(highlighted_cell_types) || (is.character(highlighted_cell_types) && length(highlighted_cell_types) > 0),
      "'categorical_variables' must be NULL or a non-empty character vector" = 
        is.null(categorical_variables) || (is.character(categorical_variables) && length(categorical_variables) > 0 ),
      "'colors' must be a non-empty character vector or left as default (=NULL)" = is.null(colors) || (is.character(colors) && length(colors) > 0),
      "'output_dir' must be a string or left as default (=NULL)" = is.null(output_dir) || (is.character(output_dir) && length(output_dir) == 1),
      "'save_plots' must be a single logical value or left as default (=FALSE)" = is.logical(save_plots) && length(save_plots) == 1,
      "'save_summaries' must be a single logical value or left as default (=FALSE)" = is.logical(save_summaries) && length(save_summaries) == 1
    )
    
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) {
        stop(paste("The specified output directory does not exist:", output_dir))
      }
    }
    
  }, error = function(e) {
    stop(e$message, call. = FALSE)
  })
  
  if (!is.null(highlighted_cell_types) && identical(sort(unique(cell_types)), sort(unique(highlighted_cell_types)))) {
    stop("'highlighted_cell_types' must not be identical to 'cell_types'. Highlighted cell types should be a subset of cell types.")
  }
  
  if (!all(cell_types %in% colnames(data))) {
    missing_cols <- cell_types[!cell_types %in% colnames(data)]
    stop(paste("The following cell types are not found in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!is.null(highlighted_cell_types) && !all(highlighted_cell_types %in% cell_types)) {
    invalid_cells <- highlighted_cell_types[!highlighted_cell_types %in% cell_types]
    stop(paste("The following highlighted cell types are not in 'cell_types':", paste(invalid_cells, collapse = ", ")))
  }
  
  if (!is.null(categorical_variables) && !all(categorical_variables %in% colnames(data))) {
    missing_cols <- categorical_variables[!categorical_variables %in% colnames(data)]
    stop(paste("The following categorical variables are not found in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
}

## This function creates the folder where the user specifies
## It also returns the new directories so that they can be called later
create_output_directories_stratified_boxplots_cells <- function(output_dir = NULL, save_plots = FALSE, save_summaries = FALSE) {
  # If output_dir is NULL, use the current working directory
  base_dir <- if (is.null(output_dir)) getwd() else output_dir
  
  # Initialize empty list to store directory paths
  dirs <- list()
  
  # Create the full paths and directories only if the corresponding save option is TRUE
  if (save_plots) {
    figures_dir <- file.path(base_dir, "cellclockR_output", "Plots")
    xfun::dir_create(figures_dir)
    dirs$figures <- figures_dir
  }
  
  if (save_summaries) {
    tables_dir <- file.path(base_dir, "cellclockR_output", "Summaries")
    xfun::dir_create(tables_dir)
    dirs$tables <- tables_dir
  }
  
  # Return the list of created directories (will be empty if no directories were created)
  return(dirs)
}

## This function loads the default colors, and changes them to provided colors if they're specified
## It checks all of the categorical variables to see how many levels they have.
## If there are fewer colors than the maximum number of levels in a categorical variable
## It returns a warning and generates a color palette of that exact size
assign_color_palette_stratified_boxplots_cells <- function(categorical_variables, data, colors = NULL) {
  default_colors <- c("#875692FF", "#F38400FF", "#A1CAF1FF", "#BE0032FF", "#C2B280FF", "#008856FF",
                      "#E68FACFF", "#0067A5FF", "#F99379FF", "#604E97FF", "#F6A600FF", "#B3446CFF",
                      "#DCD300FF", "#882D17FF", "#8DB600FF", "#654522FF", "#E25822FF", "#2B3D26FF")
  
  # If no categorical variables, return the first default color or user-provided color
  if (is.null(categorical_variables)) {
    return(if(is.null(colors)) default_colors[1] else colors[1])
  }
  
  # Determine the number of colors needed
  max_levels <- max(sapply(categorical_variables, function(x) length(unique(data[[x]]))))
  
  if (!is.null(colors)) {
    if (length(colors) < max_levels) {
      warning("Not enough user-provided colors for all levels. Switching to default colors.")
      colors <- default_colors
    }
  }
  
  if (is.null(colors) || length(colors) < max_levels) {
    colors <- default_colors
    if (length(colors) < max_levels) {
      warning("Not enough default colors for all levels. Generating a larger color palette.")
      colors <- colorRampPalette(colors)(max_levels)
    }
  }
  
  return(colors)
}

# Prepares the data for plot and summary generation
transform_and_subset_data_stratified_boxplots_cells <- function(data, categorical_variables, id, cell_types, highlighted_cell_types) {
  
  # Transform and subset data
  all_long <- data %>%
    # Convert categorical variables to factors if they exist
    {if (!is.null(categorical_variables)) 
      mutate(., across(all_of(categorical_variables), as.factor)) 
      else .} %>%
    # Select relevant columns
    select(all_of(id), if(!is.null(categorical_variables)) all_of(categorical_variables), all_of(cell_types)) %>%
    # Pivot to long format
    pivot_longer(
      cols = all_of(cell_types),
      names_to = "Cell Type",
      values_to = "Cell Count"
    )
  
  # Handle case when highlighted_cell_types is NULL
  if (is.null(highlighted_cell_types) || length(highlighted_cell_types) == 0) {
    non_highlighted_long <- all_long
    highlighted_long <- all_long[0, ]  # Empty dataframe with same structure as all_long
  } else {
    non_highlighted_cell_types <- setdiff(cell_types, highlighted_cell_types)
    
    non_highlighted_long <- all_long %>%
      dplyr::filter(`Cell Type` %in% non_highlighted_cell_types)
    
    highlighted_long <- all_long %>%
      dplyr::filter(`Cell Type` %in% highlighted_cell_types)
  }
  
  return(list(all_long = all_long, 
              highlighted_long = highlighted_long, 
              non_highlighted_long = non_highlighted_long))
}

# Generates the stratified boxplots
generate_cell_boxplots_stratified_boxplots_cells <- function(highlighted_long, non_highlighted_long, colors, variable=NULL){
  
  
  # Define the correct order of cell types
  all_cell_types <- unique(c(highlighted_long$`Cell Type`, non_highlighted_long$`Cell Type`))
  
  # Separate cell types for highlighted and non-highlighted data
  highlighted_cell_types <- unique(highlighted_long$`Cell Type`)
  non_highlighted_cell_types <- setdiff(all_cell_types, highlighted_cell_types)
  
  generate_boxplot <- function(long_data, variable = NULL, colors, is_highlighted = FALSE) {
    
    # Ensure Cell Type is treated as a factor with consistent levels
    cell_types <- if(is_highlighted) highlighted_cell_types else non_highlighted_cell_types
    
    # Ensure Cell Type is treated as a factor with appropriate levels
    long_data$`Cell Type` <- factor(long_data$`Cell Type`, levels = cell_types)
    
    if (!is.null(variable)) {
      plot <- ggplot(long_data, aes(x = `Cell Type`, y = `Cell Count`, color = !!sym(variable)))
    } else {
      plot <- ggplot(long_data, aes(x = `Cell Type`, y = `Cell Count`))
    }
    
    plot <- plot +
      theme_bw() +
      labs(x = if(is_highlighted) "" else "Cell Types",
           y = if(is_highlighted) "" else "Cell Counts") +
      theme(legend.position = if(!is_highlighted || is.null(variable)) "top" else "none",
            axis.text.x = element_text(angle = 45, vjust = 0.6, hjust = 0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = if(is_highlighted) 0 else 12),
            axis.text.y = element_text(size = 8)) +
      scale_x_discrete(limits = cell_types)
    
    if (!is.null(variable)) {
      plot <- plot +
        stat_boxplot(geom = "errorbar", lwd = 1) +
        geom_boxplot(outlier.shape = NA, lwd = 1) +
        scale_color_manual(values = colors)
    } else {
      plot <- plot +
        stat_boxplot(geom = "errorbar", lwd = 1, color = colors[1]) +
        geom_boxplot(outlier.shape = NA, lwd = 1, color = colors[1])
    }
    
    return(plot)
  }
  
  # Set factor levels appropriately for each dataset
  highlighted_long$`Cell Type` <- factor(highlighted_long$`Cell Type`, levels = highlighted_cell_types)
  non_highlighted_long$`Cell Type` <- factor(non_highlighted_long$`Cell Type`, levels = non_highlighted_cell_types)
  
  non_highlighted_plot <- generate_boxplot(non_highlighted_long, variable, colors, is_highlighted = FALSE)
  
  # Check if highlighted_long is empty
  if (nrow(highlighted_long) > 0) {
    
    highlighted_plot <- generate_boxplot(highlighted_long, variable, colors, is_highlighted = TRUE)
    
    # Calculate proportions based on the number of unique cell types
    total_cell_types <- length(unique(c(highlighted_long$`Cell Type`, non_highlighted_long$`Cell Type`)))
    highlighted_proportion <- length(unique(highlighted_long$`Cell Type`)) / total_cell_types
    non_highlighted_proportion <- length(unique(non_highlighted_long$`Cell Type`)) / total_cell_types
    
    # Ensure there's always some space for non-highlighted plot
    overall_plot <- non_highlighted_plot + highlighted_plot +
      plot_layout(ncol = 2, widths = c(non_highlighted_proportion, highlighted_proportion))
    
  } else {
    
    overall_plot <- non_highlighted_plot
  }
  
  return(overall_plot)
}


# Function to generate summaries
# Generates differently based on whether or not a variable is provided
generate_summary_stratified_boxplots_cells <- function(data, variable = NULL) {
  if (is.null(variable)) {
    summary <- data %>%
      group_by(`Cell Type`) %>%
      summarise(ymin = min(`Cell Count`),
                lower = quantile(`Cell Count`, 0.25),
                middle = median(`Cell Count`),
                upper = quantile(`Cell Count`, 0.75),
                ymax = max(`Cell Count`), 
                .groups = 'drop')
  } else {
    summary <- data %>%
      group_by(`Cell Type`, !!sym(variable)) %>%
      summarise(ymin = min(`Cell Count`),
                lower = quantile(`Cell Count`, 0.25),
                middle = median(`Cell Count`),
                upper = quantile(`Cell Count`, 0.75),
                ymax = max(`Cell Count`), 
                .groups = 'drop')
  }
  
  return(summary)
}

# Save outputs
save_outputs_stratified_boxplots_cells <- function(study, output_dirs, plots, summaries, save_plots = FALSE, save_summaries = FALSE) {
  if (is.null(output_dirs) || length(output_dirs) == 0) {
    return(invisible(NULL))
  }
  
  if (save_plots) {
    if (!is.null(output_dirs$figures)) {
      for (plot_name in names(plots)) {
        filename <- paste0(study, "_cells_", plot_name, "_boxplots.png")
        full_path <- file.path(output_dirs$figures, filename)
        ggsave(filename = full_path, plot = plots[[plot_name]], width = 12, height = 8, dpi = 300)
        cat(sprintf("Saved plot: %s\n", full_path))
      }
    } else {
      warning("Figures directory not found in output_dirs. Plots were not saved.")
    }
  }
  
  if (save_summaries) {
    if (!is.null(output_dirs$tables)) {
      for (summary_name in names(summaries)) {
        filename <- paste0(study, "_cells_", summary_name, "_summary.csv")
        full_path <- file.path(output_dirs$tables, filename)
        write.csv(summaries[[summary_name]], file = full_path, row.names = FALSE)
        cat(sprintf("Saved summary: %s\n", full_path))
      }
    } else {
      warning("Tables directory not found in output_dirs. Summaries were not saved.")
    }
  }
}