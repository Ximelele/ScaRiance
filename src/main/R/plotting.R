library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

#' Create a haplotype‐BAF scatter plot
#'
#' @param df A data.frame with columns: chrom.position, baf (0–1), and optionally group (“blue” or “red”).
#' @param chrom Character, chromosome label for the title.
#' @param xlim Numeric(2), min and max position (bp).
#' @param ylim Numeric(2), y‐axis limits (default c(-0.1,1.1)).
#' @param colors Named vector of length 2, e.g. c(blue="#1874CD", red="#CD5555").
#' @param point_size Numeric, point size.
#' @param cytoband_data Optional data.frame with columns: chrom, start, end, gieStain.
#'        Will be shaded in the background.
#' @param show_x_axis Logical, whether to show x‐axis text & ticks.
#' @return A ggplot2 object.
#' @export
create_haplotype_plot <- function(
  df,
  chrom,
  xlim,
  ylim = c(-0.1, 1.1),
  colors = c(blue = "deepskyblue4", red = "brown2"),
  point_size = 1.5,
  cytoband_data = NULL,
  show_x_axis = FALSE
) {
  # reshape to long form
  df_long <- df %>%
    select(chrom.position, blue = points.blue, red = points.red) %>%
    pivot_longer(-chrom.position, names_to = "group", values_to = "baf")

  p <- ggplot()

  # optional cytoband shading
  if (!is.null(cytoband_data)) {
    cyto_chr <- cytoband_data %>% filter(chrom == paste0("chr", chrom))
    p <- p +
      geom_rect(
        data = cyto_chr,
        aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = gieStain),
        inherit.aes = FALSE, alpha = 0.15
      ) +
      scale_fill_manual(
        values = c(gpos = "grey80", gneg = "white", stalk = "grey90"),
        guide = "none"
      )
  }

  p +
    geom_point(
      data = df_long,
      aes(x = chrom.position, y = baf, color = group),
      size = point_size
    ) +
    scale_color_manual(
      values = colors,
      labels = c(blue = "BAF", red = "1 – BAF"),
      name = NULL
    ) +
    scale_x_continuous(
      limits = xlim,
      expand = expansion(0.02),
      labels = label_number(scale = 1e-6, suffix = " Mb")
    ) +
    scale_y_continuous(limits = ylim, breaks = seq(0, 1, by = 0.2)) +
    labs(
      title = paste("Chromosome", chrom),
      x = "Position",
      y = "BAF"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 10) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      axis.title    = element_text(face = "bold", size = 8),
      axis.text     = element_text(size = 6, face = "bold"),
      axis.text.x   = if (show_x_axis) element_text() else element_blank(),
      axis.ticks.x  = if (show_x_axis) element_line() else element_blank(),
      plot.margin   = margin(t = 20, r = 20, b = 40, l = 20)
    )
}

#' Read haplotype BAF, generate and save plot
#'
#' @param baf_file Path to tab‐delimited BAF file with header, must contain Position and BAF.
#' @param out_file Path to save the PNG.
#' @param chrom Chromosome label (numeric or character).
#' @param cytoband_data Optional cytoband data.frame.
#' @export
plot_haplotype_data <- function(
  baf_file,
  out_file,
  chrom,
  cytoband_data = NULL
) {
  df <- read.table(baf_file, header = TRUE, sep = "\t") %>%
    rename(points.blue = 3) %>%           # assume col3 is BAF
    mutate(points.red = 1 - points.blue)

  # dynamic limits & point size
  if (nrow(df) > 0) {
    xlim <- range(df$Position, na.rm = TRUE)
    psize <- max(0.5, min(3, 1000 / nrow(df)))
  } else {
    xlim <- c(1, 2)
    psize <- 1.5
  }

  p <- create_haplotype_plot(
    df = df %>% rename(chrom.position = Position),
    chrom = chrom,
    xlim = xlim,
    colors = c(blue = "#1874CD", red = "#CD5555"),
    point_size = psize,
    cytoband_data = cytoband_data,
    show_x_axis = TRUE
  )

  # save with ragg for crispness (optional)
  ggsave(
    filename = out_file,
    plot     = p,
    width    = 20,
    height   = 5,
    dpi      = 500,
    device   = ragg::agg_png
  )
}
x