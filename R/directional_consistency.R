# function to make a data_frame of directional consistency results

#' @title Workout consistency of effect estimates
#' @description Before using, your data needs to be in wide format, e.g., rows are exposures and columns are the betas for those IDs you want to check consistency for
#' Convert to wide: `data <- as.data.frame(pivot_wider(data, names_from = group, values_from = beta))`
#' @param data your data frame
#'

directional_consistency <- function(data) {

  # make col1 rownames
  rownames(data) <- data[,1]
  data <- data[,c(2:ncol(data))]

  # assign directions
  data$direction <- sapply(1:nrow(data), function(x) ifelse(all(sign(data[x,]) == 1), 1, ifelse(all(sign(data[x,]) == -1), 2, 0)))

  # convert directions to labels
  data <- data %>%
  plyr::mutate(direction_group = dplyr::case_when(direction == 1 ~ "positive",
                                                  direction == 2 ~ "negative",
                                                  direction == 0 ~ "inconsistent"))

  # make factor and order
  data$direction_group <- factor(data$direction_group,
                                levels = c("positive", "negative", "inconsistent"),ordered = TRUE)

  # extract cols of interest
  data <- data[,c("direction", "direction_group")]
  return(data)

}
