# script to make colour palette
#' @title make colour palette
#' @description assign palette to use palette_discrete <- palette()
#'
palette <- function() {
  palette_discrete <- wesanderson::wes_palette(names(wes_palettes)[8], type = "discrete")
  return(palette_discrete)
}
