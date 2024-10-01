#' CoveragePlot R6 Class
#'
#' @importFrom R6 R6Class
#' @importFrom dplyr rename arrange
#' @importFrom stringr str_sort
#' @importFrom ggplot2 ggplot aes geom_bar facet_wrap theme element_text element_blank element_rect scale_y_continuous coord_cartesian labs
#'
#' @return
#' @export
#'
#' @examples
#' @name CoveragePlot
CoveragePlot <- R6Class("CoveragePlot",
                        public = list(
                          #' @field bedFile The path to the bed file
                          bedFile = NULL,
                          #' Initialize the CoveragePlot object
                          #' @description
                          #' Initialize the CoveragePlot object
                          #' with the specified bedfile
                          initialize = function(bedFile) {
                            stopifnot(ncol(bedFile) < 5)
                            self@bedFile <- bedFile
                          }

                          #' Tidy the bed file as a data frame
                          #' @description
                          #' This method renames the columns and reorder the rows
                          #' according to chromosome and genomic position in preparation to plot
                          #' @return NULL
                          tidy_bed = function(bedFile) {
                            self@bedDf <- read.delim(self@bedFile, header = FALSE) %>%
                              rename(chr = V1, start = V2, end = V3, gene = V4, coverage = V5) %>%
                              arrange(chr, start, end)
                            self@bedDf$chr <- factor(self@bedDf$chr, levels = str_sort(unique(self@bedDf$chr), numeric = TRUE))
                            self@bedDf$gene <- factor(self@bedDf$gene, levels = self@bedDf$gene)

                            self@max_coverage <- ceiling(max(self@bedDf$coverage) / 10) * 10
                            self@axis_ticks <- seq(0, self@max_coverage, length.out = 5)
                          }
                          #' Make a grid barplot for each chromosome
                          #' @description
                          #' This method creates a bar plot representing
                          #' mean gene coverage in each chromosome
                          #' @param output_pdf Path to the output PDF file, default is "~/coverage_adaptive_sampling.pdf"
                          plot = function(output_pdf = "~/coverage_adaptive_sampling.pdf") {
                            pdf(output_pdf, width = 20, height = 10)
                            print( ggplot(self@bedDf, aes(x = gene, y = coverage)) +
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
                                     scale_y_continuous(limits = c(0, max(self@axis_ticks)), breaks = self@axis_ticks) +
                                     coord_cartesian(expand = FALSE) +
                                     labs(y = "Mean coverage") )
                            dev.off()
                          }
                        ))
