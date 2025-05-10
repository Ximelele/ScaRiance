library(ggplot2)  # Note: the correct package name is ggplot2, not ggplot

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



create.baf.plot = function(chrom.position, points.red.blue, plot.red, points.darkred, points.darkblue, x.min, x.max, title, xlab, ylab, prior_bkps_pos=NULL) {
  par(mar = c(5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2, cex.lab = 2)
  plot(c(x.min,x.max), c(0,1), pch=".", type = "n", main=title, xlab=xlab, ylab=ylab)
  points(chrom.position, points.red.blue, pch=".", col=ifelse(plot.red, "red", "blue"), cex=2)
  points(chrom.position, points.darkred, pch=19, cex=0.5, col="darkred")
  points(chrom.position, points.darkblue, pch=19, cex=0.5, col="darkblue")
  if (!is.null(prior_bkps_pos)) {
    for (i in 1:length(prior_bkps_pos)) {
      abline(v=prior_bkps_pos[i])
    }
  }
}



