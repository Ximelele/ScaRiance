source("/app/.ScaRiance/src/main/R/fastPCF.R")
source("/app/.ScaRiance/src/main/R/plotting.R")
library("readr")


adjustSegmValues <- function(baf_chrom) {
  segments <- rle(baf_chrom$BAFseg)
  for (i in 1:length(segments$lengths)) {
    end <- cumsum(segments$lengths[1:i])
    end <- end[length(end)]
    start <- (end - segments$lengths[i]) + 1

    baf_chrom$BAFseg[start:end] <- median(baf_chrom$BAFphased[start:end])
  }
  return(baf_chrom)
}



segment.baf.phased = function(samplename, inputfile, outputfile, prior_breakpoints_file = NULL, gamma = 10, phasegamma = 3, kmin = 3, phasekmin = 3, no_segmentation = F, calc_seg_baf_option = 3,output_png) {
  # Function that takes SNPs that belong to a single segment and looks for big holes between
  # each pair of SNPs. If there is a big hole it will add another breakpoint to the breakpoints data.frame
  addin_bigholes = function(breakpoints, positions, chrom, startpos, maxsnpdist) {
    # If there is a big hole (i.e. centromere), add it in as a separate set of breakpoints

    # Get the chromosome coordinate right before a big hole
    bigholes = which(diff(positions) >= maxsnpdist)
    if (length(bigholes) > 0) {
      for (endindex in bigholes) {
        breakpoints = rbind(breakpoints,
                            data.frame(chrom = chrom, start = startpos, end = positions[endindex]))
        startpos = positions[endindex + 1]
      }
    }
    return(list(breakpoints = breakpoints, startpos = startpos))
  }

  # Helper function that creates segment breakpoints from SV calls
  # @param bkps_chrom Breakpoints for a single chromosome
  # @param BAFrawchr Raw BAF values of germline heterozygous SNPs on a single chromosome
  # @param addin_bigholes Flag whether bog holes in data are to be added as breakpoints
  # @return A data.frame with chrom, start and end columns
  # @author sd11
  bkps_to_presegment_breakpoints = function(chrom, bkps_chrom, BAFrawchr, addin_bigholes) {
    maxsnpdist = 3000000

    bkps_breakpoints = bkps_chrom$position

    # If there are no prior breakpoints, we cannot insert any
    if (length(bkps_breakpoints) > 0) {
      breakpoints = data.frame()

      # check which comes first, the breakpoint or the first SNP
      if (BAFrawchr$Position[1] < bkps_breakpoints[1]) {
        startpos = BAFrawchr$Position[1]
        startfromsv = 1 # We're starting from SNP data, so the first SV should be added first
      } else {
        startpos = bkps_breakpoints[1]
        startfromsv = 2 # We've just added the first SV, don't use it again
      }

      for (svposition in bkps_breakpoints[startfromsv:length(bkps_breakpoints)]) {
        selectedsnps = BAFrawchr$Position >= startpos & BAFrawchr$Position <= svposition
        if (sum(selectedsnps, na.rm = T) > 0) {

          if (addin_bigholes) {
            # If there is a big hole (i.e. centromere), add it in as a separate set of breakpoints
            res = addin_bigholes(breakpoints, BAFrawchr$Position[selectedsnps], chrom, startpos, maxsnpdist)
            breakpoints = res$breakpoints
            startpos = res$startpos
          }

          endindex = max(which(selectedsnps))
          breakpoints = rbind(breakpoints, data.frame(chrom = chrom, start = startpos, end = BAFrawchr$Position[endindex]))
          # Previous SV is the new starting point for the next segment
          startpos = BAFrawchr$Position[endindex + 1]
        }
      }

      # Add the remainder of the chromosome, if available
      if (BAFrawchr$Position[nrow(BAFrawchr)] > bkps_breakpoints[length(bkps_breakpoints)]) {
        endindex = nrow(BAFrawchr)
        breakpoints = rbind(breakpoints, data.frame(chrom = chrom, start = startpos, end = BAFrawchr$Position[endindex]))
      }
    } else {
      # There are no SVs, so create one big segment
      print("No prior breakpoints found")
      startpos = BAFrawchr$Position[1]
      breakpoints = data.frame()

      if (addin_bigholes) {
        # If there is a big hole (i.e. centromere), add it in as a separate set of breakpoints
        res = addin_bigholes(breakpoints, BAFrawchr$Position, chrom, startpos, maxsnpdist = maxsnpdist)
        breakpoints = res$breakpoints
        startpos = res$startpos
      }

      breakpoints = rbind(breakpoints, data.frame(chrom = chrom, start = startpos, end = BAFrawchr$Position[nrow(BAFrawchr)]))
    }
    return(breakpoints)
  }

  # Run PCF on presegmented data
  # @param BAFrawchr Raw BAF for this chromosome
  # @param presegment_chrom_start
  # @param presegment_chrom_end
  # @param phasekmin
  # @param phasegamma
  # @param kmin
  # @param gamma
  # @param no_segmentation Do not perform segmentation. This step will switch the haplotype blocks, but then just takes the mean BAFphased as BAFsegm
  # @return A data.frame with columns Chromosome,Position,BAF,BAFphased,BAFseg
  run_pcf = function(BAFrawchr, presegment_chrom_start, presegment_chrom_end, phasekmin, phasegamma, kmin, gamma, no_segmentation = F) {
    row.indices = which(BAFrawchr$Position >= presegment_chrom_start &
                          BAFrawchr$Position <= presegment_chrom_end)

    BAF = BAFrawchr[row.indices, 2]
    pos = BAFrawchr[row.indices, 1]


    sdev <- getMad(ifelse(BAF < 0.5, BAF, 1 - BAF), k = 25)
    # Standard deviation is not defined for a single value
    if (is.na(sdev)) {
      sdev = 0
    }
    #DCW 250314
    #for cell lines, sdev goes to zero in regions of LOH, which causes problems.
    #0.09 is around the value expected for a binomial distribution around 0.5 with depth 30
    if (sdev < 0.09) {
      sdev = 0.09
    }

    print(paste("BAFlen=", length(BAF), sep = ""))
    if (length(BAF) < 50) {
      BAFsegm = rep(mean(BAF), length(BAF))
    }else {
      res = selectFastPcf(BAF, phasekmin, phasegamma * sdev, T)
      BAFsegm = res$yhat
    }

    BAFphased = ifelse(BAFsegm > 0.5, BAF, 1 - BAF)

    if (length(BAFphased) < 50 | no_segmentation) {
      BAFphseg = rep(mean(BAFphased), length(BAFphased))
    }else {
      res = selectFastPcf(BAFphased, kmin, gamma * sdev, T)
      BAFphseg = res$yhat
    }

    if (length(BAF) > 0) {

      #
      # Note: When adding options, also add to merge_segments
      #

      # Recalculate the BAF of each segment, if required
      if (calc_seg_baf_option == 1) {
        # Adjust the segment BAF to not take the mean as that is sensitive to improperly phased segments
        BAFphseg = adjustSegmValues(data.frame(BAFphased = BAFphased, BAFseg = BAFphseg))$BAFseg
      } else if (calc_seg_baf_option == 2) {
        # Don't do anything, the BAF is already the mean
      } else if (calc_seg_baf_option == 3) {
        # Take the median, unless the median is exactly 0 or 1. At the extreme
        # there is no difference between lets say 40 and 41 copies and BB cannot
        # fit a copy number state. The mean is less prone to become exactly 0 or 1
        # but the median is generally a better estimate that is less sensitive to
        # how well the haplotypes have been reconstructed
        BAFphseg_median = adjustSegmValues(data.frame(BAFphased = BAFphased, BAFseg = BAFphseg))$BAFseg
        BAFphseg = ifelse(BAFphseg_median %in% c(0, 1), BAFphseg, BAFphseg_median)
        # if (BAFphseg_median!=0 & BAFphseg_median!=1) {
        #   BAFphseg = BAFphseg_median
        # }
      } else {
        warning("Supplied calc_seg_baf_option to segment.baf.phased not valid, using mean BAF by default")
      }
    }

    return(data.frame(Chromosome = rep(chr, length(row.indices)),
                      Position = BAFrawchr[row.indices, 1],
                      BAF = BAF,
                      BAFphased = BAFphased,
                      BAFseg = BAFphseg,
                      tempBAFsegm = BAFsegm)) # Keep track of BAFsegm for the plot below
  }

  BAFraw = as.data.frame(readr::read_tsv(file = inputfile, col_names = TRUE, col_types = "cin"))
  if (!is.null(prior_breakpoints_file)) { bkps = read.table(prior_breakpoints_file, header = T, stringsAsFactors = F) } else { bkps = NULL }

  BAFoutput = NULL
  for (chr in unique(BAFraw[, 1])) {
    print(paste0("Segmenting ", chr))
    BAFrawchr = BAFraw[BAFraw[, 1] == chr, c(2, 3)]

    BAFrawchr = BAFrawchr[!is.na(BAFrawchr[, 2]),]
    if (!is.null(bkps)) {
      bkps_chrom = bkps[bkps$chromosome == chr,]
    } else {
      bkps_chrom = data.frame(chromosome = character(), position = numeric())
    }

    breakpoints_chrom = bkps_to_presegment_breakpoints(chr, bkps_chrom, BAFrawchr, addin_bigholes = T)
    BAFoutputchr = NULL

    for (r in 1:nrow(breakpoints_chrom)) {
      BAFoutput_preseg = run_pcf(BAFrawchr, breakpoints_chrom$start[r], breakpoints_chrom$end[r], phasekmin, phasegamma, kmin, gamma, no_segmentation)
      BAFoutputchr = rbind(BAFoutputchr, BAFoutput_preseg)
    }
    # output_png
    # png(filename = paste("/app/ScalaBattenberg/Lynch.1121.03.N.bam/Plots/",samplename,"_RAFseg_chr",chr,".png",sep=""), width = 2000, height = 1000, res = 200)
    create_segmented_plot(BAFoutputchr=BAFoutputchr,bkps_chrom = bkps_chrom,chr = chr,samplename=samplename,output_png=output_png)
    # dev.off()

    # png(filename = paste("/app/ScalaBattenberg/Lynch.1121.03.N.bam/Plots/",samplename,"_segment_chr",chr,".png",sep=""), width = 2000, height = 1000, res = 200)
    create_baf_plot(
      BAFoutputchr = BAFoutputchr,
      samplename = samplename,
      chr = chr,
      output_png = output_png,
      bkps_chrom = bkps_chrom
    )




    BAFoutputchr$BAFphased = ifelse(BAFoutputchr$tempBAFsegm > 0.5, BAFoutputchr$BAF, 1 - BAFoutputchr$BAF)
    # Remove the temp BAFsegm values as they are only needed for plotting
    BAFoutput = rbind(BAFoutput, BAFoutputchr[, c(1:5)])
  }
  colnames(BAFoutput) = c("Chromosome", "Position", "BAF", "BAFphased", "BAFseg")
  print("Printing segmentation")
  print(BAFoutput)
  write.table(BAFoutput, outputfile, sep = "\t", row.names = F, col.names = T, quote = F)
}

