source("/app/ScalaBattenberg/src/main/R/clonal_ascat.R")
source("/app/ScalaBattenberg/src/main/R/orderEdges.R")
library("ASCAT")
library("assertthat")

#' Fit copy number
#'
#' Function that will fit a clonal copy number profile to segmented data. It first
#' matches the raw LogR with the segmented BAF to create segmented LogR. Then ASCAT
#' is run to obtain a clonal copy number profile. Beyond logRsegmented it produces
#' the rho_and_psi file and the cellularity_ploidy file.
#' @param samplename Samplename used to name the segmented logr output file
#' @param outputfile.prefix Prefix used for all output file names, except logRsegmented
#' @param inputfile.baf.segmented Filename that points to the BAF segmented data
#' @param inputfile.baf Filename that points to the raw BAF data
#' @param inputfile.logr Filename that points to the raw LogR data
#' @param dist_choice The distance metric that is used internally to rank clonal copy number solutions
#' @param ascat_dist_choice The distance metric used to obtain an initial cellularity and ploidy estimate
#' @param min.ploidy The minimum ploidy to consider (Default 1.6)
#' @param max.ploidy The maximum ploidy to consider (Default 4.8)
#' @param min.rho The minimum cellularity to consider (Default 0.1)
#' @param max.rho The maximum cellularity to consider (Default 1.0)
#' @param min.goodness The minimum goodness of fit for a solution to have to be considered (Default 63)
#' @param uninformative_BAF_threshold The threshold beyond which BAF becomes uninformative (Default 0.51)
#' @param gamma_param Technology parameter, compaction of Log R profiles. Expected decrease in case of deletion in diploid sample, 100 "\%" aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays (Default 1)
#' @param use_preset_rho_psi Boolean whether to use user specified rho and psi values (Default F)
#' @param preset_rho A user specified rho to fit a copy number profile to (Default NA)
#' @param preset_psi A user specified psi to fit a copy number profile to (Default NA)
#' @param read_depth Legacy parameter that is no longer used (Default 30)
#' @param analysis A String representing the type of analysis to be run, this determines whether the distance figure is produced (Default paired)
#' @author dw9, sd11
#' @export
fit.copy.number = function(samplename, outputfile.prefix, inputfile.baf.segmented, inputfile.baf, inputfile.logr, dist_choice = 0, ascat_dist_choice = 1, min.ploidy = 1.6, max.ploidy = 8.1, min.rho = 0.1, max.rho = 1.0, min.goodness = 30, uninformative_BAF_threshold = 0.51, gamma_param = 1, use_preset_rho_psi = F, preset_rho = NA, preset_psi = NA, read_depth = 30, log_segment_file, ploting_prefix) {


  assert_that(file.exists(inputfile.baf.segmented), msg = "Baf segment file not found")
  assert_that(file.exists(inputfile.baf), msg = "Baf file not found")
  assert_that(file.exists(inputfile.logr), msg = "LogR file not found")

  # Check for enough options supplied for rho and psi
  if ((max.ploidy - min.ploidy) < 0.05) {
    stop(paste("Supplied ploidy range must be larger than 0.05: ", min.ploidy, "-", max.ploidy, sep = ""))
  }
  if ((max.rho - min.rho) < 0.01) {
    stop(paste("Supplied rho range must be larger than 0.01: ", min.rho, "-", max.rho, sep = ""))
  }

  # Read in the required data

  segmented.BAF.data = as.data.frame(readr::read_tsv(file = inputfile.baf.segmented, col_names = T, col_types = "cinnn"))

  raw.BAF.data = as.data.frame(readr::read_tsv(file = inputfile.baf, col_names = T, col_types = "cin"))

  raw.logR.data = as.data.frame(readr::read_tsv(file = inputfile.logr, col_names = T, col_types = "cin"))


  identifiers = paste(segmented.BAF.data[, 1], segmented.BAF.data[, 2], sep = "_")
  dups = which(duplicated(identifiers))
  if (length(dups) > 0) {
    segmented.BAF.data = segmented.BAF.data[-dups,]
    identifiers = identifiers[-dups]
  }
  rownames(segmented.BAF.data) = identifiers

  # Drop NAs
  raw.BAF.data = raw.BAF.data[!is.na(raw.BAF.data[, 3]),]
  raw.logR.data = raw.logR.data[!is.na(raw.logR.data[, 3]),]


  BAF.data = list()
  logR.data = list()
  segmented.logR.data = list()
  matched.segmented.BAF.data = list()
  gsubchr = function(chr) gsub("chr", "", as.character(chr))

  chr.names = gsubchr(unique(segmented.BAF.data[, 1]))

  segmented.BAF.data$Chromosome = gsubchr(segmented.BAF.data$Chromosome)
  raw.BAF.data$Chromosome = gsubchr(raw.BAF.data$Chromosome)
  raw.logR.data$Chromosome = gsubchr(raw.logR.data$Chromosome)

  baf_segmented_split = split(segmented.BAF.data, f = segmented.BAF.data$Chromosome)
  baf_split = split(raw.BAF.data, f = raw.BAF.data$Chromosome)
  logr_split = split(raw.logR.data, f = raw.logR.data$Chromosome)

  # For each chromosome
  for (chr in chr.names) {
    chr.BAF.data = baf_split[[chr]]

    # Skip the rest if there is no data for this chromosome
    if (nrow(chr.BAF.data) == 0) { next }
    # Match segments with chromosome position
    chr.segmented.BAF.data = baf_segmented_split[[chr]]
    indices = match(chr.segmented.BAF.data[, 2], chr.BAF.data$Position)

    if (sum(is.na(indices)) == length(indices) | length(indices) == 0) {
      next
    }

    # Drop NAs here too
    chr.segmented.BAF.data = chr.segmented.BAF.data[!is.na(indices),]

    # Append the segmented data
    matched.segmented.BAF.data[[chr]] = chr.segmented.BAF.data
    BAF.data[[chr]] = chr.BAF.data[indices[!is.na(indices)],]

    # Append raw LogR
    chr.logR.data = logr_split[[chr]]
    indices = match(chr.segmented.BAF.data[, 2], chr.logR.data$Position)
    logR.data[[chr]] = chr.logR.data[indices[!is.na(indices)],]
    chr.segmented.logR.data = chr.logR.data[indices[!is.na(indices)],]

    # Append segmented LogR
    segs = rle(chr.segmented.BAF.data[, 5])$lengths
    cum.segs = c(0, cumsum(segs))
    for (s in 1:length(segs)) {
      chr.segmented.logR.data[(cum.segs[s] + 1):cum.segs[s + 1], 3] = mean(chr.segmented.logR.data[(cum.segs[s] + 1):cum.segs[s + 1], 3], na.rm = T)
    }
    segmented.logR.data[[chr]] = chr.segmented.logR.data
  }

  # Sync the dataframes
  selection = c()
  for (chrom in chr.names) {
    matched.segmented.BAF.data.chr = matched.segmented.BAF.data[[chrom]] #matched.segmented.BAF.data[matched.segmented.BAF.data[,1]==chrom,]
    logR.data.chr = logR.data[[chrom]] #logR.data[logR.data[,1]==chrom,]

    selection = matched.segmented.BAF.data.chr[, 2] %in% logR.data.chr[, 2]
    matched.segmented.BAF.data[[chrom]] = matched.segmented.BAF.data.chr[selection,]
    segmented.logR.data[[chrom]] = segmented.logR.data[[chrom]][selection,]
  }

  # Combine the split data frames into a single for the subsequent steps
  matched.segmented.BAF.data = do.call(rbind, matched.segmented.BAF.data)
  segmented.logR.data = do.call(rbind, segmented.logR.data)
  BAF.data = do.call(rbind, BAF.data)
  logR.data = do.call(rbind, logR.data)
  names(matched.segmented.BAF.data)[5] = samplename

  # write out the segmented logR data
  row.names(segmented.logR.data) = row.names(matched.segmented.BAF.data)
  row.names(logR.data) = row.names(matched.segmented.BAF.data)

  # log segment file from scala
  # log_segment_file
  write.table(segmented.logR.data, log_segment_file, sep = "\t", quote = F, col.names = F, row.names = F)

  # Prepare the data for going into the runASCAT functions
  segBAF = 1 - matched.segmented.BAF.data[, 5]
  segLogR = segmented.logR.data[, 3]
  logR = logR.data[, 3]
  names(segBAF) = rownames(matched.segmented.BAF.data)
  names(segLogR) = rownames(matched.segmented.BAF.data)
  names(logR) = rownames(matched.segmented.BAF.data)

  chr.segs = NULL
  for (ch in 1:length(chr.names)) {
    chr.segs[[ch]] = which(logR.data[, 1] == chr.names[ch])
  }


  distance.outfile = paste(ploting_prefix, "distance.png", sep = "", collapse = "") # kjd 20-2-2014
  copynumberprofile.outfile = paste(ploting_prefix, "copynumberprofile.png", sep = "", collapse = "") # kjd 20-2-2014
  nonroundedprofile.outfile = paste(ploting_prefix, "nonroundedprofile.png", sep = "", collapse = "") # kjd 20-2-2014
  cnaStatusFile = paste(outputfile.prefix, "copynumber_solution_status.txt", sep = "", collapse = "")

  ascat_optimum_pair = runASCAT(logR, 1 - BAF.data[, 3], segLogR, segBAF, chr.segs, ascat_dist_choice, distance.outfile, copynumberprofile.outfile, nonroundedprofile.outfile, cnaStatusFile = cnaStatusFile, gamma = gamma_param, allow100percent = T, reliabilityFile = NA, min.ploidy = min.ploidy, max.ploidy = max.ploidy, min.rho = min.rho, max.rho = max.rho, min.goodness, chr.names = chr.names, analysis = "paired") # kjd 4-2-2014


  distance.outfile = paste(ploting_prefix, "second_distance.png", sep = "", collapse = "") # kjd 20-2-2014
  copynumberprofile.outfile = paste(ploting_prefix, "second_copynumberprofile.png", sep = "", collapse = "") # kjd 20-2-2014
  nonroundedprofile.outfile = paste(ploting_prefix, "second_nonroundedprofile.png", sep = "", collapse = "") # kjd 20-2-2014

  # All is set up, now run ASCAT to obtain a clonal copynumber profile
  out = run_clonal_ASCAT(logR, 1 - BAF.data[, 3], segLogR, segBAF, chr.segs, matched.segmented.BAF.data, ascat_optimum_pair, dist_choice, distance.outfile, copynumberprofile.outfile, nonroundedprofile.outfile, gamma_param = gamma_param, read_depth, uninformative_BAF_threshold, allow100percent = T, reliabilityFile = NA, psi_min_initial = min.ploidy, psi_max_initial = max.ploidy, rho_min_initial = min.rho, rho_max_initial = max.rho, chr.names = chr.names) # kjd 21-2-2014

  ascat_optimum_pair_fraction_of_genome = out$output_optimum_pair_without_ref
  ascat_optimum_pair_ref_seg = out$output_optimum_pair
  is.ref.better = out$is.ref.better

  # Save rho, psi and ploidy for future reference
  rho_psi_output = data.frame(rho = c(ascat_optimum_pair$rho, ascat_optimum_pair_fraction_of_genome$rho, ascat_optimum_pair_ref_seg$rho), psi = c(ascat_optimum_pair$psi, ascat_optimum_pair_fraction_of_genome$psi, ascat_optimum_pair_ref_seg$psi), ploidy = c(ascat_optimum_pair$ploidy, ascat_optimum_pair_fraction_of_genome$ploidy, ascat_optimum_pair_ref_seg$ploidy), distance = c(NA, out$distance_without_ref, out$distance), is.best = c(NA, !is.ref.better, is.ref.better), row.names = c("ASCAT", "FRAC_GENOME", "REF_SEG"))
  write.table(rho_psi_output, paste(outputfile.prefix, "rho_and_psi.txt", sep = ""), quote = F, sep = "\t")
}

#' Fit subclonal copy number
#'
#' This function fits a subclonal copy number profile where a clonal profile is unlikely.
#' It goes over each segment of a clonal copy number profile and does a simple t-test. If the
#' test is significant it is unlikely that the data can be explained by a single copy number
#' state. We therefore fit a second state, i.e. there are two cellular populations with each
#' a different state: Subclonal copy number.
#' @param sample.name Name of the sample, used in figures
#' @param baf.segmented.file String that points to a file with segmented BAF output
#' @param logr.file String that points to the raw LogR file to be used in the subclonal copy number figures
#' @param rho.psi.file String pointing to the rho_and_psi file generated by \code{fit.copy.number}
#' @param output.file Filename of the file where the final copy number fit will be written to
#' @param output.figures.prefix Prefix of the filenames for the chromosome specific copy number figures
#' @param output.gw.figures.prefix Prefix of the filenames for the genome wide copy number figures
#' @param chr_names Vector of allowed chromosome names
#' @param masking_output_file Filename of where the masking details need to be written. Masking is performed to remove very high copy number state segments
#' @param max_allowed_state The maximum CN state allowed (Default 100)
#' @param prior_breakpoints_file A two column file with prior breakpoints, possibly from structural variants. This file must contain two columns: chromosome and position. These are used when making the figures
#' @param gamma Technology specific scaling parameter for LogR (Default 1)
#' @param segmentation.gamma Legacy parameter that is no longer used (Default NA)
#' @param siglevel Threshold under which a p-value becomes significant. When it is significant a second copy number state will be fitted (Default 0.05)
#' @param maxdist Slack in BAF space to allow a segment to be off it's optimum before becoming significant. A segment becomes significant very quickly when a breakpoint is missed, this parameter alleviates the effect (Default 0.01)
#' @param noperms The number of permutations to be run when bootstrapping the confidence intervals on the copy number state of each segment (Default 1000)
#' @param seed Seed to set when performing bootstrapping (Default: Current time)
#' @param calc_seg_baf_option Various options to recalculate the BAF of a segment. Options are: 1 - median, 2 - mean, 3 - ifelse median==0|1, mean, median. (Default: 3)
#' @author dw9, sd11
#' @export

callSubclones = function(sample.name, baf.segmented.file, logr.file, rho.psi.file, output.file, output.figures.prefix, output.gw.figures.prefix, chr_names, masking_output_file, max_allowed_state = 250, prior_breakpoints_file = NULL, gamma = 1, segmentation.gamma = NA, siglevel = 0.05, maxdist = 0.01, noperms = 1000, seed = as.integer(Sys.time()), calc_seg_baf_option = 3) {


  set.seed(seed)

  # Load rho/psi/goodness of fit
  res = load.rho.psi.file(rho.psi.file)
  rho = res$rho
  psit = res$psit
  psi = rho * psit + 2 * (1 - rho) # psi of all cells


  BAFvals = as.data.frame(readr::read_tsv(file = baf.segmented.file, col_names = T, col_types = "cinnn"))
  if (colnames(BAFvals)[1] == "X" || colnames(BAFvals)[1] == "chrX") {
    BAFvals = BAFvals[, -1]
  }


  # Load the raw LogR data

  LogRvals = as.data.frame(readr::read_tsv(file = logr.file, col_names = T, col_types = "cin"))
  if (colnames(LogRvals)[1] == "X" || colnames(LogRvals)[1] == "chrX") {
    # If there were rownames, then delete this column. Should not be an issue with new BB runs
    LogRvals = LogRvals[, -1]
  }


  ctrans = c(1:length(chr_names))
  names(ctrans) = chr_names
  ctrans.logR = c(1:length(chr_names))
  names(ctrans.logR) = chr_names


  ################################################################################################
  # Determine copy number for each segment
  ################################################################################################
  res = determine_copynumber(BAFvals, LogRvals, rho, psi, gamma, ctrans, ctrans.logR, maxdist, siglevel, noperms)
  subcloneres = res$subcloneres
  write.table(subcloneres, gsub(".txt", "_1.txt", output.file), quote = F, col.names = T, row.names = F, sep = "\t")

  # Scan the segments for cases that should be merged
  res = merge_segments(subcloneres, BAFvals, LogRvals, rho, psi, gamma, calc_seg_baf_option)
  BAFvals = res$bafsegmented

  res = determine_copynumber(BAFvals, LogRvals, rho, psi, gamma, ctrans, ctrans.logR, maxdist, siglevel, noperms)
  subcloneres = res$subcloneres

  # Scan for very high copy number segments and set those to NA - This is in part an artifact of small segments
  res = mask_high_cn_segments(subcloneres, BAFvals, max_allowed_state)
  subcloneres = res$subclones

  # Write the masking details to file
  masking_details = data.frame(samplename = sample.name, masked_count = res$masked_count, masked_size = res$masked_size, max_allowed_state = max_allowed_state)
  write.table(masking_details, file = masking_output_file, quote = F, col.names = T, row.names = F, sep = "\t")

  write.table(subcloneres[, c(1:3, 8:13)], output.file, quote = F, col.names = T, row.names = F, sep = "\t")
  write.table(subcloneres, gsub(".txt", "_extended.txt", output.file), quote = F, col.names = T, row.names = F, sep = "\t")


  subclones = as.data.frame(subcloneres)
  subclones[, 2:ncol(subclones)] = sapply(2:ncol(subclones), function(x) { as.numeric(as.character(subclones[, x])) })

  # Recalculate the ploidy based on the actual fit
  seg_length = floor((subclones$endpos - subclones$startpos) / 1000)
  is_subclonal_maj = abs(subclones$nMaj1_A - subclones$nMaj2_A) > 0
  is_subclonal_min = abs(subclones$nMin1_A - subclones$nMin2_A) > 0
  is_subclonal_maj[is.na(is_subclonal_maj)] = F
  is_subclonal_min[is.na(is_subclonal_min)] = F
  segment_states_min = subclones$nMin1_A * ifelse(is_subclonal_min, subclones$frac1_A, 1) + ifelse(is_subclonal_min, subclones$nMin2_A, 0) * ifelse(is_subclonal_min, subclones$frac2_A, 0)
  segment_states_maj = subclones$nMaj1_A * ifelse(is_subclonal_maj, subclones$frac1_A, 1) + ifelse(is_subclonal_maj, subclones$nMaj2_A, 0) * ifelse(is_subclonal_maj, subclones$frac2_A, 0)
  ploidy = sum((segment_states_min + segment_states_maj) * seg_length, na.rm = T) / sum(seg_length, na.rm = T)

  # Create user friendly cellularity and ploidy output file
  cellularity_ploidy_output = data.frame(purity = c(rho), ploidy = c(ploidy), psi = c(psit))
  cellularity_file = gsub("_copynumber.txt", "_purity_ploidy.txt", output.file) # NAP: updated the name of the output file, consistent with new title
  write.table(cellularity_ploidy_output, cellularity_file, quote = F, sep = "\t", row.names = F)
}

#' Given all the determined values make a copy number call for each segment
#'
#' @param BAFvals BAFsegmented data.frame with 5 columns
#' @param LogRvals Raw logR values in data.frame with 3 columns
#' @param rho Optimal rho value, the choosen cellularity
#' @param psi Optimal psi value, the choosen ploidy
#' @param gamma Platform gamma parameter
#' @param ctrans Named vector of chromosome names
#' @param ctrans.logR Named vector of chromosome names
#' @param maxdist Max distance a segment is tolerated to be not considered for subclonal copy number
#' @param siglevel Level at which a segment can become significantly different from the nearest clonal state
#' @param noperms Number of bootstrap permutations
#' @return A data.frame with copy number determined for each segment
#' @author dw9
#' @noRd
determine_copynumber = function(BAFvals, LogRvals, rho, psi, gamma, ctrans, ctrans.logR, maxdist, siglevel) {
  BAFphased = BAFvals[, 4]
  BAFseg = BAFvals[, 5]
  BAFpos = as.vector(ctrans[as.vector(BAFvals[, 1])] * 1000000000 + BAFvals[, 2])
  LogRpos = as.vector(ctrans.logR[as.vector(LogRvals[, 1])] * 1000000000 + LogRvals[, 2])

  switchpoints = c(0, which(BAFseg[-1] != BAFseg[-(length(BAFseg))] | BAFvals[-1, 1] != BAFvals[-nrow(BAFvals), 1]), length(BAFseg))
  BAFlevels = BAFseg[switchpoints[-1]]

  pval = NULL
  BAFpvals = vector(length = length(BAFseg))
  subcloneres = NULL

  for (i in 1:length(BAFlevels)) {

    l = BAFlevels[i]

    # Make sure that BAF>=0.5, otherwise nMajor and nMinor may be the wrong way around
    l = max(l, 1 - l)

    BAFke = BAFphased[(switchpoints[i] + 1):switchpoints[i + 1]]

    startpos = min(BAFpos[(switchpoints[i] + 1):switchpoints[i + 1]])
    endpos = max(BAFpos[(switchpoints[i] + 1):switchpoints[i + 1]])

    # Assuming all SNPs in this segment are on the same chromosome
    chrom = BAFvals[(switchpoints[i] + 1):switchpoints[i + 1],]$Chromosome[1]
    LogR = mean(LogRvals[LogRpos >= startpos &
                           LogRpos <= endpos &
                           !is.infinite(LogRvals[, 3]), 3], na.rm = T)

    # if we don't have a value for LogR, fill in 0
    if (is.na(LogR)) {
      LogR = 0
    }
    nMajor = (rho - 1 + l * psi * 2^(LogR / gamma)) / rho
    nMinor = (rho - 1 + (1 - l) * psi * 2^(LogR / gamma)) / rho

    # Increase nMajor and nMinor together, to avoid impossible combinations (with negative subclonal fractions)
    if (nMinor < 0) {
      if (l == 1) {
        # Avoid calling infinite copy number
        nMajor = 1000
      } else {
        nMajor = nMajor + l * (0.01 - nMinor) / (1 - l)
      }
      nMinor = 0.01
    }

    # Note that these are sorted in the order of ascending BAF:
    nMaj = c(floor(nMajor), ceiling(nMajor), floor(nMajor), ceiling(nMajor))
    nMin = c(ceiling(nMinor), ceiling(nMinor), floor(nMinor), floor(nMinor))
    x = floor(nMinor)
    y = floor(nMajor)

    # Total copy number, to determine priority options
    ntot = nMajor + nMinor

    levels = (1 - rho + rho * nMaj) / (2 - 2 * rho + rho * (nMaj + nMin))
    # Problem if rho=1 and nMaj=0 and nMin=0
    levels[nMaj == 0 & nMin == 0] = 0.5

    #DCW - just test corners on the nearest edge to determine clonality
    # If the segment is called as subclonal, this is the edge that will be used to determine the subclonal proportions that are reported first
    all.edges = orderEdges(levels, l, ntot, x, y)
    nMaj.test = all.edges[1, c(1, 3)]
    nMin.test = all.edges[1, c(2, 4)]
    test.levels = (1 - rho + rho * nMaj.test) / (2 - 2 * rho + rho * (nMaj.test + nMin.test))
    whichclosestlevel.test = which.min(abs(test.levels - l))

    # Test whether a segment should be subclonal
    if (is.na(sd(BAFke)) || sd(BAFke) == 0) {
      pval[i] = 0 # problem caused by segments with constant BAF (usually 1 or 2)
    } else {
      pval[i] = t.test(BAFke, alternative = "two.sided", mu = test.levels[whichclosestlevel.test])$p.value
    }
    if (abs(l - test.levels[whichclosestlevel.test]) < maxdist) {
      pval[i] = 1
    }

    #DCW 240314
    BAFpvals[(switchpoints[i] + 1):switchpoints[i + 1]] = pval[i]

    # If the difference is significant, call subclonal level
    if (pval[i] <= siglevel) {

      all.edges = orderEdges(levels, l, ntot, x, y)
      # Switch order, so that negative copy numbers are at the end
      na.indices = which(is.na(rowSums(all.edges)))
      if (length(na.indices) > 0) {
        all.edges = rbind(all.edges[-na.indices,], all.edges[na.indices,])
      }
      nMaj1 = all.edges[, 1]
      nMin1 = all.edges[, 2]


      subcloneres = rbind(subcloneres, c(chrom, startpos - floor(startpos / 1000000000) * 1000000000,
                                         endpos - floor(endpos / 1000000000) * 1000000000, l, pval[i], LogR, ntot,
                                         nMaj1[1], nMin1[1]))
    }else {
      #if called as clonal, use the best corner from the nearest edge
      subcloneres = rbind(subcloneres, c(chrom, startpos - floor(startpos / 1000000000) * 1000000000,
                                         endpos - floor(endpos / 1000000000) * 1000000000, l, pval[i], LogR, ntot,
                                         nMaj.test[whichclosestlevel.test], nMin.test[whichclosestlevel.test], 1, rep(NA, 57)))

    }
  }
  colnames(subcloneres) = c("Chromosome", "Start.Pos", "End.Pos", "BAF", "pval", "LogR", "Total.Copy.Number",
                            "Major.Copy.Number", "Minor.Copy.Number")
  subcloneres = as.data.frame(subcloneres)
  for (i in 2:ncol(subcloneres)) {
    subcloneres[, i] = as.numeric(as.character(subcloneres[, i]))
  }
  return(list(subcloneres = subcloneres, BAFpvals = BAFpvals))
}


#' Merge copy number segments
#'
#' Merges segments if there is not enough evidence for them to be separate. Two adjacent segments are merged
#' when they are either fit with the same clonal copy number state or when their BAF is not significantly different
#' and their logR puts them in the same square.
#' @param subclones A completely fit copy number profile in Battenberg output format
#' @param bafsegmented A BAFsegmented data.frame with the 5 columns that corresponds to the subclones file
#' @param logR The raw logR data
#' @param rho The rho estimate that the profile was fit with
#' @param psi the psi estimate that the profile was fit with
#' @param platform_gamma The gamma parameter for this platform
#' @param calc_seg_baf_option Various options to recalculate the BAF of a segment. Options are: 1 - median, 2 - mean, 3 - ifelse median== 0|1, mean, median. (Default: 3)
#' @param verbose A boolean to show merging operations (Default: FALSE)
#' @return A list with two fields: bafsegmented and subclones. The subclones field contains a data.frame in
#' Battenberg output format with the merged segments. The bafsegmented field contains the BAFsegmented data
#' corresponding to the provided subclones data.frame.
#' @author sd11, tl
#' @noRd
merge_segments = function(subclones, bafsegmented, logR, rho, psi, platform_gamma, calc_seg_baf_option = 3, verbose = F) {

  calc_nmin = function(rho, psi, baf, logr, platform_gamma) {
    return((rho -
      1 -
      (baf - 1) *
        2^(logr / platform_gamma) *
        ((1 - rho) * 2 + rho * psi)) / rho)
  }

  calc_nmaj = function(rho, psi, baf, logr, platform_gamma) {
    return((rho - 1 + baf *
      2^(logr / platform_gamma) *
      ((1 - rho) * 2 + rho * psi)) / rho)
  }

  # Convert DF into GRanges objects
  df2gr = function(DF, chr, pos1, pos2) {
    return(GenomicRanges::makeGRangesFromDataFrame(df = DF,
                                                   keep.extra.columns = T,
                                                   ignore.strand = T,
                                                   seqinfo = NULL,
                                                   seqnames.field = chr,
                                                   start.field = pos1,
                                                   end.field = pos2,
                                                   starts.in.df.are.0based = F))
  }

  # Function called when two segments have not been merged so there is no need to recheck those again
  updateNeighbour = function(subclones, INDEX, INDEX_N) {
    if (INDEX_N > INDEX) {
      subclones$Next_checked[INDEX] = T
      subclones$Prev_checked[INDEX_N] = T
    } else {
      subclones$Prev_checked[INDEX] = T
      subclones$Next_checked[INDEX_N] = T
    }
    return(subclones)
  }

  # Function called when two segments have been merged so we need to recheck its two neighbours
  updateAround = function(subclones, INDEX) {
    if (INDEX > 1) {
      subclones$Prev_checked[INDEX] = F
      subclones$Next_checked[INDEX - 1] = F
    } else {
      subclones$Prev_checked[INDEX] = T
    }
    if (INDEX < length(subclones)) {
      subclones$Next_checked[INDEX] = F
      subclones$Prev_checked[INDEX + 1] = F
    } else {
      subclones$Next_checked[INDEX] = T
    }
    return(subclones)
  }

  # Function called to test whether two segments must be checked
  checkStatus = function(subclones, INDEX, INDEX_N) {
    if (INDEX_N > INDEX) {
      # Largest segment (INDEX_N) is after smallest one (INDEX)
      stopifnot(subclones$Next_checked[INDEX] == subclones$Prev_checked[INDEX_N])
      if (subclones$Next_checked[INDEX] && subclones$Prev_checked[INDEX_N]) {
        return(T)
      } else {
        return(F)
      }
    } else {
      # Largest segment (INDEX_N) is before smallest one (INDEX)
      stopifnot(subclones$Prev_checked[INDEX] == subclones$Next_checked[INDEX_N])
      if (subclones$Prev_checked[INDEX] && subclones$Next_checked[INDEX_N]) {
        return(T)
      } else {
        return(F)
      }
    }
  }

  # Function to merge two segments
  merge_seg = function(subclones, bafsegmented, logR, INDEX, INDEX_N, calc_seg_baf_option) {
    # Update start/end information
    if (INDEX_N < INDEX) {
      GenomicRanges::end(subclones[INDEX_N]) = GenomicRanges::end(subclones[INDEX])
    } else {
      GenomicRanges::start(subclones[INDEX_N]) = GenomicRanges::start(subclones[INDEX])
    }
    # Remove segment
    subclones = subclones[-INDEX]
    if (INDEX_N < INDEX) INDEX = INDEX - 1
    # Reset neighbour checking
    subclones = updateAround(subclones, INDEX)
    if (calc_seg_baf_option == 1) {
      # This uses median
      NEW_BAF = median(bafsegmented$BAFphased[GenomicRanges::findOverlaps(subclones[INDEX], bafsegmented)@to], na.rm = T)
    } else if (calc_seg_baf_option == 2) {
      # This uses mean
      NEW_BAF = mean(bafsegmented$BAFphased[GenomicRanges::findOverlaps(subclones[INDEX], bafsegmented)@to], na.rm = T)
    } else if (calc_seg_baf_option == 3) {
      # We'll prefer the median BAF as a segment summary
      # but change to the mean when the median is extreme
      # as at 0 or 1 the BAF is uninformative for the fitting
      median_BAF = median(bafsegmented$BAFphased[GenomicRanges::findOverlaps(subclones[INDEX], bafsegmented)@to], na.rm = T)
      mean_BAF = mean(bafsegmented$BAFphased[GenomicRanges::findOverlaps(subclones[INDEX], bafsegmented)@to], na.rm = T)
      if (median_BAF != 0 && median_BAF != 1) {
        NEW_BAF = median_BAF
      } else {
        NEW_BAF = mean_BAF
      }
      rm(median_BAF, mean_BAF)
    }
    # Update both BAF and logR information
    subclones[INDEX]$BAF = NEW_BAF
    INDEX_logR = GenomicRanges::findOverlaps(subclones[INDEX], logR)@to
    if (length(INDEX_logR) == 0) {
      subclones[INDEX]$LogR = 0
    } else {
      subclones[INDEX]$LogR = mean(logR$logR[INDEX_logR], na.rm = T)
    }
    rm(INDEX_logR)
    # Update segmented baf
    bafsegmented$BAFseg[GenomicRanges::findOverlaps(subclones[INDEX], bafsegmented)@to] = NEW_BAF
    # TODO Update the logRseg as well
    # Reset IDs
    subclones$ID = 1:length(subclones)
    return(list(subclones = subclones, bafsegmented = bafsegmented))
  }

  requireNamespace("GenomicRanges")
  if (!(calc_seg_baf_option %in% 1:3)) calc_seg_baf_option = 3
  # Convert DFs into GRanges objects
  if (verbose) print('Convert DFs into GRanges objects')
  subclones = df2gr(subclones, 'chr', 'startpos', 'endpos')
  bafsegmented = df2gr(bafsegmented, 'Chromosome', 'Position', 'Position')
  logR = df2gr(logR, 'Chromosome', 'Position', 'Position')
  names(GenomicRanges::mcols(logR)) = 'logR'
  # Split GRanges objects by chromosomes
  chr_names = GenomicRanges::seqnames(GenomicRanges::seqinfo(bafsegmented))
  subclones = lapply(chr_names, function(x) subclones[GenomicRanges::seqnames(subclones) == x])
  bafsegmented = lapply(chr_names, function(x) bafsegmented[GenomicRanges::seqnames(bafsegmented) == x])
  logR = lapply(chr_names, function(x) logR[GenomicRanges::seqnames(logR) == x])
  stopifnot(all(sapply(subclones, length) > 0) &&
              all(sapply(bafsegmented, length) > 0) &&
              all(sapply(logR, length) > 0))
  names(subclones) = chr_names
  names(bafsegmented) = chr_names
  names(logR) = chr_names
  # For each chromosome
  for (CHR in chr_names) {
    if (verbose) print(paste0('Merging segments within: ', CHR))
    # Define ID, Prev_checked and Next_checked to help processing data
    subclones[[CHR]]$ID = 1:length(subclones[[CHR]])
    subclones[[CHR]]$Prev_checked = F
    subclones[[CHR]]$Next_checked = F
    subclones[[CHR]]$Prev_checked[1] = T
    subclones[[CHR]]$Next_checked[length(subclones[[CHR]])] = T
    # Pick all possible IDs
    IDs = subclones[[CHR]]$ID
    while (length(IDs) != 0) {
      # Amongst all IDs, select the ones that must be checked
      IDs = subclones[[CHR]]$ID[which(!subclones[[CHR]]$Prev_checked | !subclones[[CHR]]$Next_checked)]
      if (length(IDs) == 0) break
      # Amongst all of those, select the smallest one
      INDEX = IDs[which.min(GenomicRanges::width(subclones[[CHR]][which(subclones[[CHR]]$ID %in% IDs)]))]
      # Select neighbours (two or one if segments is first or last)
      if (INDEX == 1) {
        Neighbours = order(GenomicRanges::distance(subclones[[CHR]][INDEX], subclones[[CHR]][INDEX + 1]))
        names(Neighbours) = INDEX + 1
      } else if (INDEX == length(subclones[[CHR]])) {
        Neighbours = order(GenomicRanges::distance(subclones[[CHR]][INDEX], subclones[[CHR]][INDEX - 1]))
        names(Neighbours) = INDEX - 1
      } else {
        Neighbours = order(GenomicRanges::distance(subclones[[CHR]][INDEX], subclones[[CHR]][INDEX + c(-1, 1)]))
        names(Neighbours) = INDEX + c(-1, 1)
      }
      if (verbose) print(paste0('Working on segment: ', INDEX, ' (', subclones[[CHR]][INDEX], ')'))
      # For each neighbour
      for (i in Neighbours) {
        INDEX_N = as.numeric(names(Neighbours[i]))
        if (verbose) print(paste0('Checking neighbour: ', INDEX_N, ' (', subclones[[CHR]][INDEX_N], '; distance=', GenomicRanges::distance(subclones[[CHR]][INDEX], subclones[[CHR]][INDEX_N]), ')'))
        # Test whether seg and neighbour (INDEX and INDEX_N) have already been checked
        if (checkStatus(subclones[[CHR]], INDEX, INDEX_N)) { if (verbose) { print('Already checked') }; next }
        # Test whether seg and neighbour are far away from each other
        if (GenomicRanges::distance(subclones[[CHR]][INDEX], subclones[[CHR]][INDEX_N]) > 3e6) {
          if (verbose) print('Distance > 3Mb - do not merge')
          subclones[[CHR]] = updateNeighbour(subclones[[CHR]], INDEX, INDEX_N)
        } else {
          # Test whether seg and neighbour have the same clonal CN solution
          if (subclones[[CHR]]$nMaj1_A[INDEX] == subclones[[CHR]]$nMaj1_A[INDEX_N] &&
            subclones[[CHR]]$nMin1_A[INDEX] == subclones[[CHR]]$nMin1_A[INDEX_N] &&
            subclones[[CHR]]$frac1_A[INDEX] == 1 &&
            subclones[[CHR]]$frac1_A[INDEX_N] == 1) {
            if (verbose) print('Same clonal CN solution - merge')
            res = merge_seg(subclones[[CHR]], bafsegmented[[CHR]], logR[[CHR]], INDEX, INDEX_N, calc_seg_baf_option)
            subclones[[CHR]] = res$subclones
            bafsegmented[[CHR]] = res$bafsegmented
            rm(res)
            break
          } else {
            # Test whether seg and neighbour have different BAF/logR distributions
            if (verbose) print('Different CN solutions: check BAF and logR')
            nmin_curr = round(calc_nmin(rho, psi, subclones[[CHR]]$BAF[INDEX], subclones[[CHR]]$LogR[INDEX], platform_gamma))
            nmaj_curr = round(calc_nmaj(rho, psi, subclones[[CHR]]$BAF[INDEX], subclones[[CHR]]$LogR[INDEX], platform_gamma))
            nmin_other = round(calc_nmin(rho, psi, subclones[[CHR]]$BAF[INDEX_N], subclones[[CHR]]$LogR[INDEX_N], platform_gamma))
            nmaj_other = round(calc_nmaj(rho, psi, subclones[[CHR]]$BAF[INDEX_N], subclones[[CHR]]$LogR[INDEX_N], platform_gamma))
            if (nmin_curr == nmin_other || nmaj_curr == nmaj_other) {
              # Test whether there are more than 10 values to check significance
              if (sum(!is.na(logR[[CHR]]$logR[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX], logR[[CHR]])@to])) > 10 &&
                sum(!is.na(logR[[CHR]]$logR[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX_N], logR[[CHR]])@to])) > 10 &&
                sum(!is.na(bafsegmented[[CHR]]$BAFphased[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX], bafsegmented[[CHR]])@to])) > 10 &&
                sum(!is.na(bafsegmented[[CHR]]$BAFphased[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX_N], bafsegmented[[CHR]])@to])) > 10) {
                logr_significant = t.test(logR[[CHR]]$logR[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX], logR[[CHR]])@to],
                                          logR[[CHR]]$logR[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX_N], logR[[CHR]])@to])$p.value < 0.05
                baf_significant = t.test(bafsegmented[[CHR]]$BAFphased[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX], bafsegmented[[CHR]])@to],
                                         bafsegmented[[CHR]]$BAFphased[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX_N], bafsegmented[[CHR]])@to])$p.value < 0.05
                if ((!logr_significant) && (!baf_significant)) {
                  if (verbose) print('No significant difference - merge')
                  res = merge_seg(subclones[[CHR]], bafsegmented[[CHR]], logR[[CHR]], INDEX, INDEX_N, calc_seg_baf_option)
                  subclones[[CHR]] = res$subclones
                  bafsegmented[[CHR]] = res$bafsegmented
                  rm(res)
                  break
                } else {
                  if (verbose) print('Significant difference - do not merge')
                  subclones[[CHR]] = updateNeighbour(subclones[[CHR]], INDEX, INDEX_N)
                }
              } else {
                if (verbose) print('Too few values - do not merge')
                subclones[[CHR]] = updateNeighbour(subclones[[CHR]], INDEX, INDEX_N)
              }
            } else {
              if (verbose) print('Different squares - do not merge')
              subclones[[CHR]] = updateNeighbour(subclones[[CHR]], INDEX, INDEX_N)
            }
          }
        }
      }; rm(i)
    }
  }; rm(CHR)
  if (verbose) print('Convert GRanges objects into DFs')
  bafsegmented = data.frame(Reduce(c, bafsegmented), stringsAsFactors = F)[, -c(3:5)]
  bafsegmented$seqnames = as.character(bafsegmented$seqnames)
  colnames(bafsegmented)[1:2] = c('Chromosome', 'Position')
  subclones = data.frame(Reduce(c, subclones), stringsAsFactors = F)[, -c(4:5)]
  subclones$seqnames = as.character(subclones$seqnames)
  colnames(subclones)[1:3] = c('Chromosome', 'Start.Pos', 'End.Pos')
  subclones$ID = NULL
  subclones$Prev_checked = NULL
  subclones$Next_checked = NULL
  return(list(bafsegmented = bafsegmented, subclones = subclones))
}

#' Mask segments that have a too high CN state
#' @param subclones Subclones output data
#' @param bafsegmented BAFsegmented data
#' @param max_allowed_state The maximum state allowed before overruling takes place
#' @return A list with the masked subclones, bafsegmented and the number of segments masked and their total genome size
#' @author sd11
mask_high_cn_segments = function(subclones, bafsegmented, max_allowed_state) {
  count = 0
  masked_size = 0
  for (i in 1:nrow(subclones)) {
    if (subclones$nMaj1_A[i] > max_allowed_state | subclones$nMin1_A[i] > max_allowed_state) {
      # Mask this segment

      subclones[i, "Major.Copy.Number"] = NA
      subclones[i, "Minor.Copy.Number"] = NA
      # Mask the BAFsegmented
      bafsegmented[subclones$chr[i] == bafsegmented$Chromosome &
                     subclones$startpos[i] < bafsegmented$Position &
                     subclones$endpos[i] >= bafsegmented$Position, c("BAFseg")] = NA
      count = count + 1
      masked_size = masked_size + (subclones$endpos[i] - subclones$startpos[i])
    }
  }
  return(list(subclones = subclones, bafsegmented = bafsegmented, masked_count = count, masked_size = masked_size))
}


#' Plot the copy number genome wide in two different ways. This creates the Battenberg average
#' profile where subclonal copy number is represented as a mixture of two different states and
#' the Battenberg subclones profile where subclonal copy number is plotted as two different
#' separate states. The thickness of the line represents the fraction of tumour cells carying
#' the particular state.
#' @noRd
plot.gw.subclonal.cn = function(subclones, BAFvals, rho, ploidy, goodness, output.gw.figures.prefix, chr.names, tumourname) {
  # Map start and end of each segment into the BAF values. The plot uses the index of this BAF table as x-axis
  pos_min = array(NA, nrow(subclones))
  pos_max = array(NA, nrow(subclones))
  for (i in 1:nrow(subclones)) {
    segm_chr = subclones$chr[i] == BAFvals$Chromosome &
      subclones$startpos[i] < BAFvals$Position &
      subclones$endpos[i] >= BAFvals$Position
    pos_min[i] = min(which(segm_chr))
    pos_max[i] = max(which(segm_chr))
  }

  # For those segments that are subclonal, Obtain the second state.
  is_subclonal = which(subclones$frac1_A < 1)
  subcl_min = array(NA, length(is_subclonal))
  subcl_max = array(NA, length(is_subclonal))
  for (i in 1:length(is_subclonal)) {
    segment_index = is_subclonal[i]
    segm_chr = subclones$chr[segment_index] == BAFvals$Chromosome &
      subclones$startpos[segment_index] < BAFvals$Position &
      subclones$endpos[segment_index] >= BAFvals$Position
    subcl_min[i] = min(which(segm_chr))
    subcl_max[i] = max(which(segm_chr))
  }

  # Determine whether it's the major or the minor allele that is represented by two states
  is_subclonal_maj = abs(subclones$nMaj1_A - subclones$nMaj2_A) > 0
  is_subclonal_min = abs(subclones$nMin1_A - subclones$nMin2_A) > 0
  is_subclonal_maj[is.na(is_subclonal_maj)] = F
  is_subclonal_min[is.na(is_subclonal_min)] = F

  # BB represents subclonal CN as a mixture of two CN states. Calculate this mixture for both minor allele and total CN.
  #segment_states_min = subclones$nMin1_A * ifelse(is_subclonal_min, subclones$frac1_A, 1)  + ifelse(is_subclonal_min, subclones$nMin2_A, 0) * ifelse(is_subclonal_min, subclones$frac2_A, 0)
  #segment_states_tot = (subclones$nMaj1_A+subclones$nMin1_A) * ifelse(is_subclonal_maj, subclones$frac1_A, 1) + ifelse(is_subclonal_maj, subclones$nMaj2_A+subclones$nMin2_A, 0) * ifelse(is_subclonal_maj, subclones$frac2_A, 0)

  segment_states_min = subclones$nMin1_A * ifelse(is_subclonal_min, subclones$frac1_A, 1) + ifelse(is_subclonal_min, subclones$nMin2_A, 0) * ifelse(is_subclonal_min, subclones$frac2_A, 0)
  segment_states_maj = subclones$nMaj1_A * ifelse(is_subclonal_maj, subclones$frac1_A, 1) + ifelse(is_subclonal_maj, subclones$nMaj2_A, 0) * ifelse(is_subclonal_maj, subclones$frac2_A, 0)
  segment_states_tot = segment_states_maj + segment_states_min

  # Determine which SNPs are on which chromosome, to be used as a proxy for chromosome size in the plots
  chr.segs = lapply(1:length(chr.names), function(ch) { which(BAFvals$Chromosome == chr.names[ch]) })

  # Plot subclonal copy number as mixtures of two states
  png(filename = paste(output.gw.figures.prefix, "_average.png", sep = ""), width = 2000, height = 500, res = 200)
  create.bb.plot.average(bafsegmented = BAFvals,
                         ploidy = ploidy,
                         rho = rho,
                         goodnessOfFit = goodness,
                         pos_min = pos_min,
                         pos_max = pos_max,
                         segment_states_min = segment_states_min,
                         segment_states_tot = segment_states_tot,
                         chr.segs = chr.segs,
                         chr.names = chr.names,
                         tumourname = tumourname)
  dev.off()

  # Plot subclonal copy number as two separate states
  png(filename = paste(output.gw.figures.prefix, "_subclones.png", sep = ""), width = 2000, height = 500, res = 200)
  create.bb.plot.subclones(bafsegmented = BAFvals,
                           subclones = subclones,
                           ploidy = ploidy,
                           rho = rho,
                           goodnessOfFit = goodness,
                           pos_min = pos_min,
                           pos_max = pos_max,
                           subcl_min = subcl_min,
                           subcl_max = subcl_max,
                           is_subclonal = is_subclonal,
                           is_subclonal_maj = is_subclonal_maj,
                           is_subclonal_min = is_subclonal_min,
                           chr.segs = chr.segs,
                           chr.names = chr.names,
                           tumourname = tumourname)
  dev.off()
}

#' Load the rho and psi estimates from a file.
#' @noRd
load.rho.psi.file = function(rho.psi.file) {
  rho_psi_info = read.table(rho.psi.file, header = T, sep = "\t", stringsAsFactors = F)
  # Always use best solution from grid search - reference segment sometimes gives strange results
  rho = rho_psi_info$rho[rownames(rho_psi_info) == "FRAC_GENOME"] # rho = tumour percentage (called tp in previous versions)
  psit = rho_psi_info$psi[rownames(rho_psi_info) == "FRAC_GENOME"] # psi of tumour cells
  goodness = rho_psi_info$distance[rownames(rho_psi_info) == "FRAC_GENOME"] # goodness of fit
  return(list(rho = rho, psit = psit, goodness = goodness))
}

#' Collapse a BAFsegmented file into segment start and end points
#'
#' This function looks through the BAFsegmented for stretches of equal
#' BAFseg and records the start and end coordinates in a data.frame
#' @param bafsegmented The BAFsegmented output from segmentation
#' @return A data.frame with columns chromosome, start and end
#' @author sd11
#' @noRd
collapse_bafsegmented_to_segments = function(bafsegmented) {
  segments_collapsed = data.frame()
  for (chrom in unique(bafsegmented$Chromosome)) {
    bafsegmented_chrom = bafsegmented[bafsegmented$Chromosome == chrom,]
    segments = rle(bafsegmented_chrom$BAFseg)
    startpoint = 1
    for (i in 1:length(segments$lengths)) {
      endpoint = startpoint + segments$lengths[i] - 1
      segments_collapsed = rbind(segments_collapsed,
                                 data.frame(chromosome = chrom, start = bafsegmented_chrom$Position[startpoint], end = bafsegmented_chrom$Position[endpoint]))
      startpoint = endpoint + 1
    }
  }
  return(segments_collapsed)
}

#' Function to make additional figures
#'
#' @param samplename Name of the sample for the plot title
#' @param logr_file File containing all logR data
#' @param bafsegmented_file File containing the BAFsegmented data
#' @param logrsegmented_file File with the logRsegmented data
#' @param allelecounts_file Optional file with raw allele counts (Default: NULL)
#' @author sd11
#' @export
make_posthoc_plots = function(samplename, logr_file, bafsegmented_file, logrsegmented_file, allelecounts_file = NULL) {
  # Make some post-hoc plots
  logr = Battenberg::read_table_generic(logr_file)
  bafsegmented = as.data.frame(Battenberg::read_table_generic(bafsegmented_file))
  logrsegmented = as.data.frame(Battenberg::read_table_generic(logrsegmented_file, header = F))
  colnames(logrsegmented) = c("Chromosome", "Position", "logRseg")
  outputfile = paste0(samplename, "_alleleratio.png")
  allele_ratio_plot(samplename = samplename, logr = logr, bafsegmented = bafsegmented, logrsegmented = logrsegmented, outputfile = outputfile, max.plot.cn = 8)

  if (!is.null(allelecounts_file)) {
    allelecounts = as.data.frame(Battenberg::read_table_generic(allelecounts_file))
    outputfile = paste0(samplename, "_coverage.png")
    coverage_plot(samplename, allelecounts, outputfile)
  }
}

