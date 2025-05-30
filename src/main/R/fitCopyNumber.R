source("/app/.ScaRiance/src/main/R/clonal_ascat.R")
source("/app/.ScaRiance/src/main/R/orderEdges.R")
library("ASCAT")
library("assertthat")


fit.copy.number = function(samplename, outputfile.prefix, inputfile.baf.segmented, inputfile.baf, inputfile.logr, dist_choice = 0, ascat_dist_choice = 1, min.ploidy = 1.6, max.ploidy = 4.8, min.rho = 0.1, max.rho = 1.0, min.goodness = 63, uninformative_BAF_threshold = 0.51, gamma_param = 1, read_depth = 30, log_segment_file, ploting_prefix) {


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


  distance.outfile = paste(ploting_prefix, "distance.png", sep = "", collapse = "")
  copynumberprofile.outfile = paste(ploting_prefix, "copynumberprofile.png", sep = "", collapse = "")
  nonroundedprofile.outfile = paste(ploting_prefix, "nonroundedprofile.png", sep = "", collapse = "")
  cnaStatusFile = paste(outputfile.prefix, "copynumber_solution_status.txt", sep = "", collapse = "")

  ascat_optimum_pair = runASCAT(logR, 1 - BAF.data[, 3], segLogR, segBAF, chr.segs, ascat_dist_choice, distance.outfile, copynumberprofile.outfile, nonroundedprofile.outfile, cnaStatusFile = cnaStatusFile, gamma = gamma_param, allow100percent = T, reliabilityFile = NA, min.ploidy = min.ploidy, max.ploidy = max.ploidy, min.rho = min.rho, max.rho = max.rho, min.goodness, chr.names = chr.names, analysis = "paired")

  # All is set up, now run ASCAT to obtain a clonal copynumber profile
  out = run_clonal_ASCAT(logR, 1 - BAF.data[, 3], segLogR, segBAF, chr.segs, matched.segmented.BAF.data, ascat_optimum_pair, dist_choice, NA, NA, NA, gamma_param = gamma_param, read_depth, uninformative_BAF_threshold, allow100percent = T, reliabilityFile = NA, psi_min_initial = min.ploidy, psi_max_initial = max.ploidy, rho_min_initial = min.rho, rho_max_initial = max.rho, chr.names = chr.names)

  ascat_optimum_pair_fraction_of_genome = out$output_optimum_pair_without_ref
  ascat_optimum_pair_ref_seg = out$output_optimum_pair
  is.ref.better = out$is.ref.better

  # Save rho, psi and ploidy for future reference
  rho_psi_output = data.frame(rho = c(ascat_optimum_pair$rho, ascat_optimum_pair_fraction_of_genome$rho, ascat_optimum_pair_ref_seg$rho), psi = c(ascat_optimum_pair$psi, ascat_optimum_pair_fraction_of_genome$psi, ascat_optimum_pair_ref_seg$psi), ploidy = c(ascat_optimum_pair$ploidy, ascat_optimum_pair_fraction_of_genome$ploidy, ascat_optimum_pair_ref_seg$ploidy), distance = c(NA, out$distance_without_ref, out$distance), is.best = c(NA, !is.ref.better, is.ref.better), row.names = c("ASCAT", "FRAC_GENOME", "REF_SEG"))
  print("fit.copy.number")
  print(rho_psi_output)
  write.table(rho_psi_output, paste(outputfile.prefix, "rho_and_psi.txt", sep = ""), quote = F, sep = "\t")
}

callSubclones = function(baf.segmented.file, logr.file, rho.psi.file, output.file, chr_names, gamma = 1, siglevel = 0.05, maxdist = 0.01, seed = as.integer(Sys.time())) {


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
    LogRvals = LogRvals[, -1]
  }


  ctrans = c(1:length(chr_names))
  names(ctrans) = chr_names
  ctrans.logR = c(1:length(chr_names))
  names(ctrans.logR) = chr_names


  ################################################################################################
  # Determine copy number for each segment
  ################################################################################################
  res = determine_copynumber(BAFvals, LogRvals, rho, psi, gamma, ctrans, ctrans.logR, maxdist, siglevel)
  subcloneres = res$subcloneres

  write.table(subcloneres, output.file, quote = F, col.names = T, row.names = F, sep = "\t")

}

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


    # If the segment is called as subclonal, this is the edge that will be used to determine the subclonal proportions that are reported first
    all.edges = orderEdges(levels, l, ntot, x, y)
    nMaj.test = all.edges[1, 1]
    nMin.test = all.edges[1, 2]
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


#' Load the rho and psi estimates from a file.
#' @noRd
load.rho.psi.file = function(rho.psi.file) {
  rho_psi_info = read.table(rho.psi.file, header = T, sep = "\t", stringsAsFactors = F)
  rho = rho_psi_info$rho[rownames(rho_psi_info) == "FRAC_GENOME"] # rho = tumour percentage (called tp in previous versions)
  psit = rho_psi_info$psi[rownames(rho_psi_info) == "FRAC_GENOME"] # psi of tumour cells
  goodness = rho_psi_info$distance[rownames(rho_psi_info) == "FRAC_GENOME"] # goodness of fit
  return(list(rho = rho, psit = psit, goodness = goodness))
}




