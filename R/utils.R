#' @noRd
roundAny <- function(x, accuracy){
  round(x / accuracy) * accuracy
}
