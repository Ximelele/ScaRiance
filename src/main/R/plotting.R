library(ggplot2)
library(dplyr)
library(readr)
library(viridis)
library(gtools)

create_segmented_plot <- function(BAFoutputchr, bkps_chrom = NULL, samplename, chr, output_png) {
  # Prepare data
  BAFoutputchr$PositionMb <- BAFoutputchr$Position / 1e6
  if (!is.null(bkps_chrom)) {
    bkps_chrom$PositionMb <- bkps_chrom$position / 1e6
  }

  # Create the plot
  p <- ggplot(BAFoutputchr, aes(x = PositionMb)) +
    geom_point(aes(y = BAF), color = "brown3", size = 0.5, alpha = 0.5) +
    geom_point(aes(y = tempBAFsegm), color = "forestgreen", size = 0.3, alpha = 0.7) +
    scale_x_continuous(limits = c(min(BAFoutputchr$PositionMb), max(BAFoutputchr$PositionMb)), expand = c(0.02, 0.02)) +
    scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, by = 0.2)) +
    labs(
      title = paste0(samplename, ", chromosome ", chr),
      x = "Position (Mb)",
      y = "BAF (phased)"
    ) +
    coord_cartesian(clip = "off") +
    # Add fixed aspect ratio to maintain proportions
    theme_classic() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, color = "black", face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10, face = "bold"),
      plot.margin = margin(t = 20, r = 20, b = 40, l = 20),
      # Keep x-axis text visible
      axis.text.x = element_text(size = 10, face = "bold")
    )
    # Set a fixed aspect ratio to prevent squashing


  # Add vertical lines for breakpoints if provided
  if (!is.null(bkps_chrom)) {
    p <- p + geom_vline(data = bkps_chrom, aes(xintercept = PositionMb),
                        color = "black", linetype = "dashed", alpha = 0.6)
  }

  # Save the plot with improved height - KEY CHANGES HERE
  output_path <- paste0(output_png, "RAFseg_chr", chr, ".png")
  # Increase height to at least 8 inches, and maintain wider width for chromosome
  ggsave(output_path, plot = p, width = 20, height = 5, dpi = 500)
}

create.haplotype.plot <- function(chrom.position, points.blue, points.red, x.min, x.max, title, xlab, ylab, point.size = 1, cytoband_data, alpha = 0.7) {  # Added alpha parameter with default value
  data <- data.frame(
    chrom.position = chrom.position,
    points.light = points.blue,
    points.dark = points.red
  )

  plot <- ggplot(data) +
    geom_point(aes(x = chrom.position, y = points.blue), color = "deepskyblue4", size = point.size, alpha = alpha) +  # Added alpha here
    geom_point(aes(x = chrom.position, y = points.red), color = "brown2", size = point.size, alpha = alpha) +  # Added alpha here
    scale_x_continuous(limits = c(x.min, x.max), expand = c(0.02, 0.02)) +
    scale_y_continuous(limits = c(-0.1, 1.1), breaks = seq(0, 1, by = 0.2)) +
    labs(title = title, x = xlab, y = ylab) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5, color = "black", face = "bold"),
      axis.title.x = element_text(size = 8, face = "bold"),
      axis.title.y = element_text(size = 8, face = "bold"),
      axis.text = element_text(size = 6, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 20, r = 20, b = 40, l = 20)
    )

  return(plot)
}

plot.haplotype.data <- function(haplotyped.baf.file, imageFileName, chrom ) {
  mut_data <- read.table(haplotyped.baf.file, sep = "\t", header = TRUE)

  # Read and filter cytoband data for the relevant chromosome
  # cyto_data <- read_cytoband_data(cytoband_file, paste0("chr", chrom))

  # Explicitly set x_min to 0 and use the maximum end of cytoband data for x_max
  if (nrow(mut_data) > 0) {
    x_min = min(mut_data$Position,na.rm=T)
    x_max = max(mut_data$Position,na.rm=T)
  } else {
    x_min = 1
    x_max = 2
  }

  # Dynamically calculate point size
  point.size <- if (nrow(mut_data) > 0) max(0.5, min(3, 1000 / nrow(mut_data))) else 1

  haplotype_plot <- create.haplotype.plot(
    chrom.position = mut_data$Position,
    points.blue = mut_data[, 3],
    points.red = 1 - mut_data[, 3],
    x.min = x_min,
    x.max = x_max,
    title = paste("Chromosome", chrom, sep = " "),
    xlab = "pos",
    ylab = "BAF",
    point.size = point.size,
    alpha = 0.5  # Set transparency to 0.5 - adjust as needed
  )

  # Save the plot as a PNG file
  ggsave(filename = imageFileName, plot = haplotype_plot, width = 20, height = 5, dpi = 500)
}



library(ggplot2)

create_baf_plot <- function(BAFoutputchr, samplename, chr, output_png, bkps_chrom = NULL) {
  # Prepare data
  BAFoutputchr$PositionMb <- BAFoutputchr$Position / 1e6
  BAFoutputchr$plotColor <- ifelse(BAFoutputchr$tempBAFsegm > 0.5, "firebrick2", "dodgerblue3")
  BAFoutputchr$BAFseg_flipped <- 1 - BAFoutputchr$BAFseg

  if (!is.null(bkps_chrom)) {
    bkps_chrom$PositionMb <- bkps_chrom$position / 1e6
  }

  # Create plot
  p <- ggplot(BAFoutputchr, aes(x = PositionMb)) +
    geom_point(aes(y = BAF, color = plotColor), size = 0.5, alpha = 0.7) +
    geom_point(aes(y = BAFseg), color = "orange2", size = 0.3, alpha = 0.7) +
    geom_point(aes(y = BAFseg_flipped), color = "purple3", size = 0.3, alpha = 0.7) +
    scale_color_manual(values = c("red" = "brown3", "blue" = "deepskyblue4")) +
    scale_x_continuous(limits = c(min(BAFoutputchr$PositionMb), max(BAFoutputchr$PositionMb)), expand = c(0.02, 0.02)) +
    scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, by = 0.2)) +
    labs(
      title = paste0(samplename, ", chromosome ", chr),
      x = "Position (Mb)",
      y = "BAF (phased)"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5, color = "black", face = "bold"),
      axis.title.x = element_text(size = 8, face = "bold"),
      axis.title.y = element_text(size = 8, face = "bold"),
      axis.text = element_text(size = 6, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 20, b = 40, l = 20)
    )

  # Add vertical lines for breakpoints if provided
  if (!is.null(bkps_chrom)) {
    p <- p + geom_vline(data = bkps_chrom, aes(xintercept = PositionMb),
                        color = "black", linetype = "dashed", alpha = 0.6)
  }

  # Save the plot
  output_path <- paste0(output_png, "segment_chr", chr, ".png")
  ggsave(output_path, plot = p, width = 20, height = 5, dpi = 500)
}






create_cnv_plot <- function(cnv_file, baf_file, output_png) {
  # 1. Load CNV
  cnv <- read_tsv(cnv_file,
                  col_types = cols(
                    Chromosome = col_character(),
                    Start.Pos = col_double(),
                    End.Pos = col_double(),
                    BAF = col_double(),
                    Total.Copy.Number = col_double()
                  )) %>%
    mutate(
      Chr = Chromosome,
      Start = Start.Pos,
      End = End.Pos,
      StartMb = Start.Pos / 1e6,
      EndMb = End.Pos / 1e6
    )

  # 2. Load BAF
  baf_raw <- read_tsv(baf_file)
  baf_col <- setdiff(names(baf_raw), c("Chromosome", "Position"))[1]
  baf <- baf_raw %>%
    rename(Chr = Chromosome, BAF = !!baf_col) %>%
    mutate(
      Start = Position,
      End = Position,
      PositionMb = Position / 1e6
    )

  # 3. Chromosome ordering
  all_chr_levels <- mixedsort(unique(c(cnv$Chr, baf$Chr)))
  cnv$Chr <- factor(cnv$Chr, levels = all_chr_levels)
  baf$Chr <- factor(baf$Chr, levels = all_chr_levels)

  # 4. Plot
  p <- ggplot() +
    # Raw BAF dots
    geom_point(data = baf, aes(x = PositionMb, y = BAF),
               color = "lavender", size = 0.3, alpha = 0.3) +
    # Colored BAF segment line (same y = BAF)
    geom_segment(data = cnv,
                 aes(x = StartMb, xend = EndMb,
                     y = BAF, yend = BAF,
                     color = Total.Copy.Number),
                 size = 1.2, alpha = 0.95) +
    facet_wrap(~ Chr, scales = "free_x", ncol = 4) +
    scale_color_viridis(
      option = "plasma", direction = -1,
      name = "Total CN",
      limits = c(floor(min(cnv$Total.Copy.Number)), ceiling(max(cnv$Total.Copy.Number))),
      guide = guide_colorbar(title.position = "bottom", title.hjust = 0.5,
                             barwidth = 15, barheight = 0.5)
    ) +
    labs(
      title = "CNV per chromsome"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(t = 10, b = 10)),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 9),
      strip.text = element_text(size = 10, face = "bold"),
      panel.spacing = unit(2, "lines"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.position = "bottom",
      legend.box = "horizontal"
    )

  # 5. Save
  ggsave(
    filename = paste0(output_png, "cnv_per_chrom.png"),
    plot = p,
    width = 14,
    height = 16,
    dpi = 300
  )
}






