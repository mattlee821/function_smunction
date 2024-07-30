# script to make colour palette
#' @title make colour palettes
#' @description assign palette to use palettes <- palette()
#' @returns returns a list of three colour palettes: 2 discrete and 1 continous
#' @export
#' @importFrom wesanderson wes_palette
#' @importFrom grDevices colorRampPalette
palette <- function() {
  palette_discrete1 <- wesanderson::wes_palette(name = "Royal1")
  palette_discrete2 <- wesanderson::wes_palette(name = "Darjeeling1")
  palette_continuous <- structure(.Data = grDevices::colorRampPalette(c(wesanderson::wes_palette(name = "Royal1")[1],
                                                                        wesanderson::wes_palette(name = "Royal1")[2]))(100),
                                  class = "palette", name = "Royal1 continuous")
  # create a list to store all palettes
  palettes <- list(
    discrete1 = palette_discrete1,
    discrete2 = palette_discrete2,
    continuous = palette_continuous
  )

  # return list
  return(palettes)

}
