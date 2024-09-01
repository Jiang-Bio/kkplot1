#' Title remove_gr_grangse
#'
#' @param b Grange object
#' @param minoverlap Overlap at least a few bases
#'
#' @return Granges
#' @export
#'
#' @examples remove_gr_grangse(GRANGE,minoverlap=4)
remove_gr_grangse <- function(b,minoverlap=4){
  overlaps <- findOverlaps(b,minoverlap = minoverlap, drop.self = TRUE, drop.redundant = TRUE)
  keep <- setdiff(seq_along(b), unique(subjectHits(overlaps)))
  gr_no_overlaps <- b[keep]
  gr_no_overlaps
}
