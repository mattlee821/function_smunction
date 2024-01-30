#'  `epic_somalogic_create_subtypes()`: Create a list of dataframes with sex-specific and subtype data for Epic SomaLogic data.
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
  data_cancer_eo <- subset(data_cancer, ageevent < 55)
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

#' `epic_somalogic_create_sex_specific_dataframes()`: create sex-specific dataframes for a given 'ID_col' group for Epic SomaLogic data
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
  combined <- data_group
  combined$ID <- paste0(data_group[[ID_col]][1], ";combined")

  male <- data_group[data_group$sex == 2, ]
  male$ID <- paste0(data_group[[ID_col]][1], ";male")

  female <- data_group[data_group$sex == 1, ]
  female$ID <- paste0(data_group[[ID_col]][1], ";female")

  return(list(
    combined = combined,
    male = male,
    female = female
  ))
}

#'  `epic_somalogic_create_subtype_dataframes()`: create sex-specific and subtype data for Epic SomaLogic data.
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

#' `epic_somalogic_exclusions()`: process exclusions from `epic_somalogic_create_subtype_dataframes()`
#' @description
#' This function takes a nested list of data frames and a nested list of exclusions, and applies
#' the exclusion logic to remove specified rows (samples) and columns (features) from the data frames.
#'
#' @param list_data A list of data frames.
#' @param exclusions A list of exclusions for each data frame.
#' @return The modified list_data after applying exclusions.
#'
#' @export
epic_somalogic_exclusions <- function(list_data, exclusions) {
# Iterate over the main list
for (i in seq_along(list_data)) {
  # Get the name of the current list
  current_list_name <- names(list_data)[i]

  # Check if the current list is in the exclusions list
  if (current_list_name %in% names(exclusions)) {
    # Iterate over the sublists in exclusions
    for (j in seq_along(exclusions[[current_list_name]])) {
      # Get the name of the current sublist
      current_sublist_name <- names(exclusions[[current_list_name]])[j]

      # Iterate over the data frames in list_data
      for (k in seq_along(list_data[[i]])) {
        # Get the name of the current data frame
        current_df_name <- names(list_data[[i]])[k]

        # Check if the current data frame is in the exclusion sublist
        if (current_df_name %in% current_sublist_name) {
          # Convert to data.table
          dt <- as.data.table(list_data[[i]][[k]])

          # 1. Print the name of the dataframe
          cat("# Processing:", current_list_name, "/", current_sublist_name, "/", "\n")

          # 2. Print the number of rows and columns containing "seq"
          num_rows_before <- nrow(dt)
          num_cols_seq_before <- sum(grepl("^seq", colnames(dt)))
          cat("## Number of samples before:", num_rows_before, "\n")
          cat("## Number of features before:", num_cols_seq_before, "\n")

          # Get the exclusion vectors
          id_samples_exclude <- exclusions[[current_list_name]][[current_sublist_name]]$id_samples_exclude
          id_features_exclude <- exclusions[[current_list_name]][[current_sublist_name]]$id_features_exclude

          # 3. Print the number of rows to exclude and the values of "id_samples_exclude"
          num_rows_to_exclude <- sum(dt$idepic %in% id_samples_exclude)
          cat("## Number of samples to exclude:", num_rows_to_exclude, "\n")
          cat("## samples to exclude:", paste(id_samples_exclude, collapse = ", "), "\n")

          # 4. Print the number of columns to exclude and the values of "id_features_exclude"
          num_cols_to_exclude <- sum(colnames(dt) %in% id_features_exclude)
          cat("## Number of features to exclude:", num_cols_to_exclude, "\n")
          cat("## features to exclude:", paste(id_features_exclude, collapse = ", "), "\n")

          # Identify rows to keep
          rows_to_keep <- !(dt$idepic %in% id_samples_exclude)

          # Identify columns to keep
          cols_to_keep <- !(colnames(dt) %in% id_features_exclude)

          # Remove specified rows and columns from the data frame
          list_data[[i]][[k]] <- dt[rows_to_keep, cols_to_keep, with = FALSE]

          # 5. Print the number of rows and the number of columns containing "seq" after making the exclusions
          num_rows_after <- nrow(list_data[[i]][[k]])
          num_cols_seq_after <- sum(grepl("^seq", colnames(list_data[[i]][[k]])))
          cat("## Number of samples after:", num_rows_after, "\n")
          cat("## Number of features after:", num_cols_seq_after, "\n")
        }
      }
    }
  }
}

return(list_data)
}


