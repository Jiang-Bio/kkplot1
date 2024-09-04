
#' Title
#'
#' @param PWMatrix
#' @param DNAStringSet
#' @param shuffle_num
#' @param pseudo
#' @param p.cutoff
#'
#' @return
#' @export
#'
#' @examples
enrich_in_DNAstringSet <- function(PWMatrix,DNAStringSet,shuffle_num=10,
                                   pseudo=20,p.cutoff=1e-05){
  match <- motifmatchr::matchMotifs(pwms = PWMatrix, DNAStringSet,
                                    out = "matches",bg ="subject",p.cutoff = p.cutoff)
  match <- match@assays@data@listData$motifMatches %>% colSums()
  match_s <- lapply(1:shuffle_num, function(x){
    df <- motifmatchr::matchMotifs(pwms =PWMatrix,subject = DNAStringSet %>% universalmotif::shuffle_sequences(),
                                   out ="matches",p.cutoff=p.cutoff,bg="subject")
    match_s <- df@assays@data@listData$motifMatches %>% colSums()
    match_s})%>% largeListToDf() %>% rowMeans()
  result <- (match+pseudo )/(match_s+pseudo)
  result
}
