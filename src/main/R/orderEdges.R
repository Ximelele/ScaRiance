#' Convenience function that orders edges or squares
#' @author dw9, kd7
#' @noRd
orderEdges <- function(levels, observedBAF, totalCopyNumber, minorCNV, majorCNV) {
  result <- matrix(nrow = 6, ncol = 2)
  colnames(result) <- c("Major.Copy.Number", "Minor.Copy.Number")
  isLogRCase <- (totalCopyNumber < minorCNV + majorCNV + 1)

  if (observedBAF > levels[3]) {
    if (isLogRCase) {
      result[, "Major.Copy.Number"] <- c(majorCNV, majorCNV - 1, majorCNV, majorCNV + 1, majorCNV + 1, majorCNV + 1)
      result[, "Minor.Copy.Number"] <- c(minorCNV, minorCNV, minorCNV, minorCNV, minorCNV - 1, minorCNV)
    }else {
      result[, "Major.Copy.Number"] <- c(majorCNV + 1, majorCNV + 1, majorCNV + 1, majorCNV, majorCNV - 1, majorCNV)
      result[, "Minor.Copy.Number"] <- c(minorCNV, minorCNV - 1, minorCNV, minorCNV, minorCNV, minorCNV)
    }
  }
  else if (observedBAF > levels[2]) {
    if (isLogRCase) {
      result[, "Major.Copy.Number"] <- c(majorCNV, majorCNV, majorCNV, majorCNV + 1, majorCNV + 1, majorCNV + 1)
      result[, "Minor.Copy.Number"] <- c(minorCNV, minorCNV - 1, minorCNV, minorCNV, minorCNV - 1, minorCNV)
    } else {
      result[, "Major.Copy.Number"] <- c(majorCNV + 1, majorCNV + 1, majorCNV + 1, majorCNV, majorCNV, majorCNV)
      result[, "Minor.Copy.Number"] <- c(minorCNV, minorCNV - 1, minorCNV, minorCNV, minorCNV - 1, minorCNV)
    }
  }
  else {
    # Case 2b
    if (isLogRCase) {
      result[, "Major.Copy.Number"] <- c(majorCNV, majorCNV, majorCNV, majorCNV, majorCNV - 1, majorCNV)
      result[, "Minor.Copy.Number"] <- c(minorCNV, minorCNV - 1, minorCNV, minorCNV + 1, minorCNV + 1, minorCNV + 1)
    } else {
      result[, "Major.Copy.Number"] <- c(majorCNV, majorCNV - 1, majorCNV, majorCNV, majorCNV, majorCNV)
      result[, "Minor.Copy.Number"] <- c(minorCNV + 1, minorCNV + 1, minorCNV + 1, minorCNV, minorCNV - 1, minorCNV)
    }
  }

  has_negative <- apply(result, 1, function(row) any(row < 0))
  if (any(has_negative)) {
    result[has_negative,] <- NA
  }

  return(result)

}


          #' Function that fetches the nearest edge for a given a rho, psi, BAF and major and minor allele
          #' that corresponds to a certain mix of two copy number states. It first identifies the nearest edge
          #' and then just compares the vertices at the end of this edge to find the best corner.
          #' @author dw9, kd7
          #' @noRd
GetNearestCorners_bestOption <- function(rho, BAFreq, nMajor, nMinor) {
  nMaj = c(floor(nMajor), ceiling(nMajor), floor(nMajor), ceiling(nMajor))
  nMin = c(ceiling(nMinor), ceiling(nMinor), floor(nMinor), floor(nMinor))
  x = floor(nMinor)
  y = floor(nMajor)

  # total copy number, to determine priority options
  ntot = nMajor + nMinor

  BAF_levels = (1 - rho + rho * nMaj) / (2 - 2 * rho + rho * (nMaj + nMin))
  #problem if rho=1 and nMaj=0 and nMin=0
  BAF_levels[nMaj == 0 & nMin == 0] = 0.5

  nMaj1 = NULL
  nMin1 = NULL


  # case 1 or 2a:
  if (BAFreq > BAF_levels[3]) { #DCW
    #LogR criterion: ntot < x+y+1
    if (ntot < x + y + 1) {
      # take the six options, sorted according to LogR priority (3+3) + simplicity (1+2+1+2)
      nMaj1 = y
      nMin1 = x
    }
    else {
      nMaj1 = y + 1
      nMin1 = x
    }
  }
    # case 2c:
    #else if( is.finite(BAF_levels[2]) &&  (BAFreq>BAF_levels[2]) ) { # kjd 14-2-2014
  else if (BAFreq > BAF_levels[2]) { #DCW
    if (ntot < x + y + 1) {
      nMaj1 = y
      nMin1 = x
    }
    else {
      nMaj1 = y + 1
      nMin1 = x

    }
  }
    # case 2b:
  else {
    if (ntot < x + y + 1) {
      nMaj1 = y
      nMin1 = x

    }
    else {
      nMaj1 = y
      nMin1 = x + 1

    }
  }

  nMaj_vect = c(nMaj1)
  nMin_vect = c(nMin1)

  nearest_segment = list(nMaj = nMaj_vect, nMin = nMin_vect)

  return(nearest_segment)
}
