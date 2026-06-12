#' @importFrom dplyr .data
NULL

#' @noRd
roundAny <- function(x, accuracy){
  round(x / accuracy) * accuracy
}
