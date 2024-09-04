
#' Title
#'
#' @param reads
#' @param peaks
#' @param singleEnd
#'
#' @return
#' @export
#'
#' @examples
frip <- function(reads, peaks, singleEnd=T) {
  reads <- BamFile(reads)
  # find reads in peaks
  overlaps <- summarizeOverlaps(
    peaks,
    reads,
    mode="IntersectionNotEmpty",
    ignore.strand=T,
    singleEnd=singleEnd,
    count.mapped.reads=T
  )

  # sum up all reads in peaks and divide by all mapped reads to get FRiP
  readsInPeaks <- sum(assay(overlaps))
  if (class(reads) %in% c("BamViews", "BamFile")) {
    allReads <- colData(overlaps)$mapped
  } else {
    allReads <- colData(overlaps)$records
  }
  result <- readsInPeaks/allReads

  return(result)
}
