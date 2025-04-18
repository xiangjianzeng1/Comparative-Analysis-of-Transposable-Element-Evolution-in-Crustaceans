##################################################
library(tidyverse)
library(cowplot)

# Read data
all_data1 <- list.files(pattern = "*.merged_100windows.txt") %>% 
  map_df(~read_tsv(.x) %>% 
           group_by(species) %>% 
           mutate(window_id = row_number()))

# Create plots for each species
plot_list <- map(unique(all_data1$species), ~{
  current_species <- .x
  sp_data <- filter(all_data1, species == current_species)
  
  # Calculate scaling factors
  max_gene <- max(sp_data$gene_count)
  max_te <- max(sp_data$te_count)
  ratio <- max_gene / max_te
  
  ggplot(sp_data, aes(x = window_id)) +
    # Trend lines
    geom_line(aes(y = gene_count, color = "Gene Count"), size = 0.6) +
    geom_line(aes(y = te_count * ratio, color = "TE Count"), size = 0.6) +
    
    # Axis settings
    scale_x_continuous(
      name = "Window Sequence",
      breaks = seq(0, max(sp_data$window_id), by = 20),  # Show tick every 20 windows
      expand = c(0.02, 0.02)  # Reduce margins
    ) +
    scale_y_continuous(
      name = "Gene Count",
      sec.axis = sec_axis(~./ratio, name = "TE Count"),
      breaks = scales::pretty_breaks(n = 4),  # Auto generate 4 major ticks
      expand = c(0.05, 0.05)
    ) +
    
    # Color settings
    scale_color_manual(
      values = c("Gene Count" = "#1b9e77", "TE Count" = "#d95f02")
    ) +
    
    # Theme settings
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      panel.border = element_rect(color = "grey70", fill = NA, size = 0.3),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(margin = margin(t = 5)),  # Top margin for X-axis
      axis.title.y.left = element_text(color = "#1b9e77", margin = margin(r = 6)),
      axis.title.y.right = element_text(color = "#d95f02", margin = margin(l = 6)),
      axis.ticks = element_line(color = "grey60", size = 0.3),
      axis.text = element_text(color = "grey40", size = 7),
      legend.position = "none",
      plot.margin = unit(c(5,5,5,5), "mm")
    ) +
    labs(title = current_species)
})

# Combine plots
final_plot <- plot_grid(
  plotlist = plot_list,
  ncol = 2,
  align = "hv",  # Ensure axis alignment
  axis = "tblr",
  labels = "AUTO",
  label_size = 10,
  label_x = 0.1,
  label_y = 0.95
) 
final_plot
# Save high-quality image
ggsave("gene_te_final_plot.pdf", final_plot, width = 28, height = 24, units = "cm")

# --------------------------
# Data Reading and Preprocessing Module
# --------------------------

# Verify file existence
data_files <- list.files(pattern = "*merged_100windows.txt$")
if(length(data_files) == 0) stop("No input files found. Please confirm:
1. Correct working directory? (Current: ", getwd(), ")
2. Files contain 'merged_100windows.txt' suffix")

# Safely read data
all_data <- map_df(data_files, ~{
  tryCatch({
    df <- read_tsv(.x, show_col_types = FALSE) %>%
      group_by(species) %>%
      mutate(window_id = row_number())  # Add sequential window IDs
    
    # Data validation
    required_cols <- c("chr", "start", "end", "gene_count", "te_count", "species")
    if(!all(required_cols %in% names(df))) {
      stop("File ", .x, " missing required columns: ", 
           paste(setdiff(required_cols, names(df)), collapse = ", "))
    }
    df
  }, error = function(e) {
    message("Error processing file ", .x, ": ", e$message)
    NULL
  })
})

# Validate data content
if(nrow(all_data) == 0) stop("All files failed to load. Please check file formats")
print(str(all_data))  # Print data structure

# --------------------------
# Plot Generation Module
# --------------------------

# Create output directory
output_dir <- "Species_Plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define standardized plotting function
create_species_plot <- function(sp) {
  sp_data <- all_data %>%
    filter(species == sp) %>%
    arrange(window_id)  # Ensure window order
  
  # Calculate scaling factors
  gene_max <- max(sp_data$gene_count, na.rm = TRUE)
  te_max <- max(sp_data$te_count, na.rm = TRUE)
  ratio <- ifelse(te_max == 0, 1, gene_max / te_max)  # Prevent division by zero
  
  # Core plotting logic
  ggplot(sp_data, aes(x = window_id)) +
    geom_line(aes(y = gene_count, color = "Gene"), linewidth = 0.8) +
    geom_line(aes(y = te_count * ratio, color = "TE"), linewidth = 0.8) +
    scale_color_manual(
      name = "",
      values = c(Gene = "#d95f02", TE = "#1b9e77"),
      labels = c(
        Gene = paste("Gene Count (max", gene_max, ")"), 
        TE = paste("TE Count (max", te_max, ")")
      )
    ) +
    scale_x_continuous(
      name = "Window Sequence",
      breaks = scales::breaks_pretty(6),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      name = "Gene Count",
      sec.axis = sec_axis(~./ratio, name = "TE Count"),
      breaks = scales::breaks_pretty(5)
    ) +
    labs(
      title = paste("Species:", sp),
      subtitle = paste("Total Windows:", nrow(sp_data)),
      caption = paste("Data source:", paste(data_files, collapse = "; "))
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      axis.title.y.left = element_text(color = "#d95f02", margin = margin(r = 10)),
      axis.title.y.right = element_text(color = "#1b9e77", margin = margin(l = 10)),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.border = element_rect(fill = NA, color = "grey70"),
      legend.position = "bottom",
      plot.margin = margin(15, 15, 15, 15)
    )
}

# --------------------------
# Batch Output Module
# --------------------------

# Get valid species list
valid_species <- unique(all_data$species)
if(length(valid_species) == 0) stop("No valid species found in the data")

# Generate and save plots
walk(valid_species, ~{
  sp <- .x
  message("Processing species: ", sp)
  
  # Generate plot
  plot <- tryCatch({
    create_species_plot(sp)
  }, error = function(e) {
    message("Plot generation failed for ", sp, ": ", e$message)
    NULL
  })
  
  # Safely save output
  if(!is.null(plot)) {
    filename <- file.path(output_dir, paste0(gsub("[^[:alnum:]]", "_", sp), ".pdf"))
    
    ggsave(
      filename = filename,
      plot = plot,
      device = cairo_pdf,
      width = 20, 
      height = 15,
      units = "cm",
      dpi = 300
    )
    message("Successfully saved: ", filename)
  }
})

message("Processing complete! Output directory: ", normalizePath(output_dir))
