library(ggplot2)
library(tidyverse)

# Reorder the data in function of chromosome and genomic position:
tidy_bed_data <- function(input_bed) {
  bed <- read.delim(input_bed, header = FALSE) %>%
    rename(chr = V1, start = V2, end = V3, gene = V4, coverage = V5) %>%
    arrange(chr, start, end)

  bed$chr <- factor(bed$chr, levels = str_sort(unique(bed$chr), numeric = TRUE))
  bed$gene <- factor(bed$gene, levels = bed$gene)

  return(bed)
}

# Function to generate the plot
generate_plot <- function(bed, output_pdf) {
  # Calculate maximum coverage and ticks for the y-axis
  max_coverage <- ceiling(max(bed$coverage) / 10) * 10
  axis_ticks <- seq(0, max_coverage, length.out = 5)

  # Generate the plot and save as PDF
  pdf(output_pdf, width = 20, height = 10)
  print( ggplot(bed, aes(x = gene, y = coverage)) +
    geom_bar(stat = "identity", fill = "turquoise4") +
    facet_wrap(~ chr, nrow = 2, scales = "free_x") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8.5),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 18),
      strip.text.x = element_text(size = 18),
      strip.background = element_rect(fill = NA),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
      panel.background = element_blank(),
      panel.grid.major = element_line(colour = "grey")
    ) +
    scale_y_continuous(limits = c(0, max(axis_ticks)), breaks = axis_ticks) +
    coord_cartesian(expand = FALSE) +
    labs(y = "Mean coverage") )
  dev.off()
}

# Usage
input_bed <- ""     #Change to bed file
output_pdf <- "coverage_adaptive_sampling.pdf"

# Call functions
bed_data <- tidy_bed_data(input_bed)
generate_plot(bed_data, output_pdf)
