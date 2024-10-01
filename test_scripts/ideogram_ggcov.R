# Load required libraries
library(tidyverse)
library(ggcoverage)
library(ggplot2)
library(cowplot)

# Tidy BED data
tidy_bed_data <- function(bed_file) {
  # Load and tidy the BED file
  bed <- read.delim(bed_file, header = FALSE) %>%
    rename(chr = V1, start = V2, end = V3, gene = V4, coverage = V5) %>%
    arrange(chr, start, end)

  # Convert chromosome column to factor with numeric ordering
  bed$chr <- factor(bed$chr, levels = str_sort(unique(bed$chr), numeric = TRUE))

  # Reformat the data to match the ggcoverage structure
  ggbed <- bed %>%
    rename(seqnames = chr, score = coverage, label = gene) %>%
    mutate(width = end - start, Type = seqnames, Group = seqnames)

  # Split the data by chromosome
  ggbed_list <- split(ggbed, ggbed$seqnames)

  return(ggbed_list)
}

# Generate a coverage plot with the ideogram:
ggcov <- function(df) {
  ggcoverage(
    data = df,
    color = "red",
    mark.region = select(df, start, end, label),
    mark.color = "lavenderblush3",
    mark.label.size = 3,
    range.position = "out"
  ) +
    geom_ideogram(
      genome = "hg38",
      plot.space = 0,
      highlight.centromere = TRUE
    )
}

# Generate all plots on a single page and save as PDF
generate_coverage_plots <- function(ggbed_list, output_file) {
  # Generate a list of all coverage plots
  plot_list <- lapply(ggbed_list, ggcov)

  # Create a grid layout for all the plots in one page
  grid <- plot_grid(plotlist = plot_list, ncol = 2)  # Adjust ncol for your layout

  # Save the plot grid to a PDF
  pdf(output_file, width = 25, height = 18)
  print(grid)
  dev.off()
}

# Usage example:
bed_file <- ""     # Path to bed file
output_pdf <- "ggcov_ideograms.pdf"

# Tidy the BED data
ggbed_list <- tidy_bed_data(bed_file)

# Generate and save coverage plots
generate_coverage_plots(ggbed_list, output_pdf)



