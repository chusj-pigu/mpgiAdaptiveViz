# Load required libraries
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggrepel)

# Function to tidy cytoband and bed data
tidy_data <- function(cyto_file, bed_file) {
  # Load hg38 cytobands file
  cyto1 <- read.table(file = cyto_file, header = FALSE, sep = "\t")
  colnames(cyto1) <- c("Chrom", "Start", "End", "CytobandName", "Stain")

  # Color mapping based on cyto Stain
  colors <- data.frame(
    Stain = unique(cyto1$Stain),
    colors = c("white", "grey87", "grey75", "grey50", "black", "red", "grey93", "cornflowerblue")
  )

  # Tidy cytoband data
  cyto1 <- cyto1 %>%
    select(-Start) %>%
    filter(nchar(Chrom) < 6) %>%
    arrange(Chrom, End) %>%
    mutate(pos = End) %>%
    filter(Chrom != "chrM") %>%
    left_join(colors)

  cyto_list <- split(cyto1, cyto1$Chrom)

  # Adjust for spacing in cytoband positions
  for (i in seq_along(cyto_list)) {
    for (j in 2:nrow(cyto_list[[i]])) {
      cyto_list[[i]][j, 5] <- cyto_list[[i]][j, 2] - cyto_list[[i]][j - 1, 2]
    }
    cyto_list[[i]]$CytobandName <- factor(cyto_list[[i]]$CytobandName, levels = cyto_list[[i]]$CytobandName)
  }

  # Load bed file and reorder
  bed <- read.delim(bed_file, header = FALSE) %>%
    rename(chr = V1, start = V2, end = V3, gene = V4, coverage = V5) %>%
    arrange(chr, start, end)

  bed$chr <- factor(bed$chr, levels = str_sort(unique(bed$chr), numeric = TRUE))
  bed$gene <- factor(bed$gene, levels = bed$gene)

  bed_list <- split(bed, bed$chr)

  # Re-order cytoband list according to chromosome
  cyto_list <- cyto_list[names(bed_list)]

  return(list(cyto_list = cyto_list, bed_list = bed_list, max_coverage = ceiling(max(bed$coverage) / 10) * 10))
}

# Function to plot coverage and ideogram for a given chromosome
plotCovIdeo <- function(cyto, cov, max_coverage, axis_ticks) {
  # Ideogram plot
  ideogram <- ggplot(cyto, aes(fill = CytobandName, y = pos, x = Chrom)) +
    geom_bar(position = "stack", stat = "identity", color = "black") +
    theme_classic() +
    theme(
      legend.position = 'none',
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_text(color = "white")
    ) +
    coord_flip(ylim = c(0, max(cyto$End))) +
    scale_fill_manual(values = cyto$colors, breaks = rev(levels(cyto$CytobandName))) +
    scale_y_continuous(limits = c(0, max(cyto$End)))

  # Coverage plot
  cov_plot <- ggplot(cov, aes(xmin = start, xmax = end, ymin = 0, ymax = coverage)) +
    geom_rect(color = "turquoise4", fill = "turquoise4") +
    theme_classic() +
    ggtitle(unique(cov$chr)) +
    theme(
      axis.title = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
    ) +
    geom_text_repel(aes(x = start, y = coverage, label = gene)) +
    scale_x_continuous(limits = c(0, max(cyto$End))) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_coverage), breaks = axis_ticks)

  # Combine plots
  grid <- plot_grid(cov_plot, ideogram, ncol = 1, align = "h", rel_heights = c(4, 1))
  return(grid)
}

# Function to generate all plots and save them as PDF
generate_plots <- function(cyto_list, bed_list, max_coverage, output_file) {
  axis_ticks <- seq(0, max_coverage, length.out = 5)
  plots <- list()

  for (i in seq_along(bed_list)) {
    plots[[i]] <- plotCovIdeo(cyto_list[[i]], bed_list[[i]], max_coverage, axis_ticks)
  }

  # Save plots in a PDF file
  pdf(output_file, width = 20, height = 10)
  num_plots <- length(plots)
  for (i in seq(1, num_plots, by = 6)) {
    print(do.call(plot_grid, plots[i:min(i + 5, num_plots)]))
  }
  dev.off()
}

# Usage
cyto_file <- "cyto.txt"
bed_file <- ""    # Path to bed file
output_pdf <- "ideogram_coverage.pdf"

# Tidy data
tidied_data <- tidy_data(cyto_file, bed_file)

# Generate plots
generate_plots(tidied_data$cyto_list, tidied_data$bed_list, tidied_data$max_coverage, output_pdf)
