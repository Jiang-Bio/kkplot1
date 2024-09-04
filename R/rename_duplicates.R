#' Title
#'
#' @param strings
#'
#' @return
#' @export
#'
#' @examples
rename_duplicates <- function(strings) {
  counts <- table(strings)
  seen <- setNames(rep(1, length(counts)), names(counts))
  new_strings <- character(length(strings))

  for (i in seq_along(strings)) {
    current_string <- strings[i]
    if (counts[current_string] > 1) {
      seen_count <- seen[current_string]
      new_strings[i] <- paste(current_string, seen_count, sep = "_")
      seen[current_string] <- seen_count + 1
    } else {
      new_strings[i] <- paste(current_string, 1, sep = "_")
    }
  }

  return(new_strings)
}
