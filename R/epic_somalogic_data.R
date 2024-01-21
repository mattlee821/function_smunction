#'  Create a list of dataframes with sex-specific and subtype data for Epic SomaLogic data.
#' @description
#' This function is a master function for 'epic_somalogic_create_sex_specific_dataframes()'
#' and 'epic_somalogic_create_subtypes()'. It takes a dataframe, splits it by
#' 'ID_col', and creates a list of sex-specific and subtype dataframes
#' for each 'ID_col' group. 'ID_col' is the filename (e.g., cancer-transformed.txt)
#'
#' @param data A dataframe containing the Epic SomaLogic data.
#' @return A list of dataframes with sex-specific and subtype data.
#'
#' @export
epic_somalogic_create_subtypes <- function(data) {
  # Function to create subtypes within a dataframe
  # Subset for cases and non-cases
  data_cancer <- subset(data, cncr_mal_clrt == 1)
  data_cancer <- droplevels(data_cancer)
  data_control <- subset(data, cncr_mal_clrt == 0)
  data_control <- droplevels(data_control)
  data <- bind_rows(data_cancer, data_control)
  data <- droplevels(data)

  # Subset for subtypes - Colon
  data_cancer_colon <- subset(data, cncr_mal_clrt_colon == 1)
  data_colon <- bind_rows(data_cancer_colon, data_control)
  data_colon <- droplevels(data_colon)

  # Subset for subtypes - Rectum
  data_cancer_rectum <- subset(data, cncr_mal_clrt_rectum == 1)
  data_rectum <- bind_rows(data_cancer_rectum, data_control)
  data_rectum <- droplevels(data_rectum)

  # Subset for subtypes - Early Onset
  data_cancer_eo <- subset(data, ageevent < 55)
  data_control_eo <- subset(data_control, age < 55)
  data_control_eo <- data_control_eo %>%
    mutate(age_exit_cancer_1st = ifelse(age_exit_cancer_1st > 55, 55, age_exit_cancer_1st))
  data_eo <- bind_rows(data_cancer_eo, data_control_eo)
  data_eo <- droplevels(data_eo)

  # Update the 'ID' column
  data$ID <- paste0(data$ID, ";cancer")
  data_colon$ID <- paste0(data_colon$ID, ";colon")
  data_rectum$ID <- paste0(data_rectum$ID, ";rectum")
  data_eo$ID <- paste0(data_eo$ID, ";earlyonset")

  # Return the list of subtypes
  return(list(
    cancer = data,
    colon = data_colon,
    rectum = data_rectum,
    earlyonset = data_eo
  ))
}

#' Create sex-specific dataframes for a given 'ID_col' group for Epic SomaLogic data
#' @description
#' This function takes a subset of data for a specific 'ID_col' group and creates
#' sex-specific dataframes - sex-combined, male, and female.
#'
#' @param data_group A subset of the original dataframe for a specific 'ID_col' group.
#' @param ID_col The unique identifier for the 'ID_col' group.
#' @return A list containing sex-specific dataframes.
#'
#' @export
epic_somalogic_create_sex_specific_dataframes <- function(data_group, ID_col) {
  # Create sex-combined and sex-specific dataframes for a given 'ID_col' group
  sex_combined <- data_group
  sex_combined$ID <- paste0(data_group[[ID_col]][1], ";combined")

  male <- data_group[data_group$sex == 2, ]
  male$ID <- paste0(data_group[[ID_col]][1], ";male")

  female <- data_group[data_group$sex == 1, ]
  female$ID <- paste0(data_group[[ID_col]][1], ";female")

  return(list(
    sex_combined = sex_combined,
    male = male,
    female = female
  ))
}

#'  Create a list of dataframes with sex-specific and subtype data for Epic SomaLogic data.
#' @description
#' This function is a master function for 'epic_somalogic_create_sex_specific_dataframes()'
#' and 'epic_somalogic_create_subtypes()'. It takes a dataframe, splits it by
#' 'ID_col', and creates a list of sex-specific and subtype dataframes
#' for each 'ID_col' group. 'ID_col' is the filename (e.g., cancer-transformed.txt)
#'
#' @param data A dataframe containing the Epic SomaLogic data.
#' @param ID_col The column containing the unique identifier for the 'data_id' group.
#' @return A list of dataframes with sex-specific and subtype data.
#'
#' @export
epic_somalogic_create_subtype_dataframes <- function(data, ID_col) {
  # Create a list to store the dataframes
  list_data <- list()

  # Split the data by 'ID_col'
  data_split <- split(data, data[[ID_col]])

  # Iterate over each 'ID_col' group
  for (ID_val in names(data_split)) {
    # Access the group dataframe
    data_group <- data_split[[ID_val]]

    # Create sex-combined and sex-specific dataframes
    sex_dataframes <- epic_somalogic_create_sex_specific_dataframes(data_group, ID_col)
    list_data[[ID_val]] <- sex_dataframes

    # Access each sex-specific dataframe
    for (data_type in names(sex_dataframes)) {
      # Access the dataframe for the current type
      current_data <- sex_dataframes[[data_type]]

      # Create subtypes within the dataframe
      subtypes <- epic_somalogic_create_subtypes(current_data)

      # Store the updated dataframes in the list
      list_data[[ID_val]][[paste0(data_type, "_cancer")]] <- subtypes$cancer
      list_data[[ID_val]][[paste0(data_type, "_colon")]] <- subtypes$colon
      list_data[[ID_val]][[paste0(data_type, "_rectum")]] <- subtypes$rectum
      list_data[[ID_val]][[paste0(data_type, "_eo")]] <- subtypes$earlyonset
    }
  }

  return(list_data)
}
