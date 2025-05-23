################################################################################
# 05. Plot read depth for each genome assembly and save all plots in a grid arrange.
################################################################################

# Clear global environment
rm(list = ls())

# Load required libraries
library(tidyr)
library(dplyr)
library(Biostrings)
library(ggplot2)
library(gridExtra)

# Set working directory
setwd("xxxxxxxx")

# List all FASTA files starting with "scaffolds_"
fasta_files <- list.files(pattern = "^scaffolds_.*\\.fasta$")

# Function to process each FASTA file and generate plots
process_genome <- function(fasta_file) {
  # Read the FASTA file
  sequences <- readDNAStringSet(fasta_file)
  
  # Extract sequence names (headers)
  sequence_names <- names(sequences)
  
  # Split the sequence names by "_"
  split_names <- strsplit(sequence_names, "_")
  
  # Convert the list to a dataframe
  df <- do.call(rbind, lapply(split_names, function(x) data.frame(t(x), stringsAsFactors = FALSE)))
  
  # Extract desired columns
  new_df <- df %>% select(X2, X4, X6)
  
  # Rename columns
  colnames(new_df) <- c("scaffold", "length", "coverage")
  new_df$scaffold <- as.numeric(new_df$scaffold)
  new_df$coverage <- as.numeric(new_df$coverage)
  new_df$length <- as.numeric(new_df$length)
  
  # Filter data to exclude outliers and very small scaffolds
  filtered_df <- new_df[new_df$coverage >= 8 & new_df$coverage <= 200, ]
  filtered_df <- filtered_df[filtered_df$length >= 200, ]
  
  # Calculate the peak bin for coverage
  coverage_hist_data <- ggplot_build(ggplot(filtered_df, aes(x = coverage)) +
                                       geom_histogram(binwidth = 1))$data[[1]]
  peak_bin <- coverage_hist_data$x[which.max(coverage_hist_data$count)]
  
  # Generate the genome coverage distribution plot
  cvg_plot <- ggplot(filtered_df, aes(x = coverage)) +
    geom_histogram(binwidth = 1, fill = "black", color = "black") +
    labs(title = paste("Genome Coverage Distribution:", fasta_file), 
         x = "Coverage", y = "Frequency") +
    theme_minimal() +
    geom_vline(aes(xintercept = peak_bin), color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text", x = peak_bin, y = Inf, label = paste("Peak Bin:", round(peak_bin, 2)), vjust = 1, color = "black") +
    theme(plot.title = element_text(hjust = 0.5, size = 11))
  
  return(cvg_plot)
}

# Initialize a list to store all plots
plot_list <- list()

# Loop through all FASTA files, generate plots, and save them individually as TIFF files
for (fasta_file in fasta_files) {
  plot <- process_genome(fasta_file)
  plot_list[[fasta_file]] <- plot
  
  # Save each individual plot as a TIFF file
  tiff_filename <- paste0("Coverage_Distribution_", gsub(".fasta", "", fasta_file), ".tiff")
  tiff(tiff_filename, width = 800, height = 600, res = 120)
  print(plot)
  dev.off()
}

# Combine plots into grids of three columns and save as TIFF files
plots_per_page <- 9 # Adjust to fit the number of plots per grid
num_pages <- ceiling(length(plot_list) / plots_per_page)

for (page in seq_len(num_pages)) {
  start_idx <- (page - 1) * plots_per_page + 1
  end_idx <- min(page * plots_per_page, length(plot_list))
  grid_plots <- plot_list[start_idx:end_idx]
  
  # Save the grid as a TIFF file
  tiff_filename <- paste0("Coverage_Distribution_Grid_Page_", page, ".tiff")
  tiff(tiff_filename, width = 2400, height = 1800, res = 120) # Larger size for grid
  do.call(grid.arrange, c(grid_plots, ncol = 3))
  dev.off()
  
  # Save the grid as a PNG file
  png_filename <- paste0("Coverage_Distribution_Grid_Page_", page, ".png")
  png(png_filename, width = 2400, height = 1800, res = 120) # Larger size for grid
  do.call(grid.arrange, c(grid_plots, ncol = 3))
  dev.off()
  
}

# Message to indicate the process is complete
cat("All individual and grid plots saved as TIFF and PNG files.")
