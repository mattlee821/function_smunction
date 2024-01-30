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

#' `epic_somalogic_models()`: a list of coxph models
#' @description
#' This function defines three Cox proportional hazards models including
#'
#' @return A list containing three Cox proportional hazards models
#'
epic_somalogic_models <- function() {
  # Define the list of models
  models <- list(
    model_1 = function(df, exposure_var) {
      coxph(Surv(age, ageevent, cncr_mal_clrt_colon) ~ get(exposure_var) +
              cluster(idepic),
            data = df)
    },
    model_2 = function(df, exposure_var) {
      coxph(Surv(age, ageevent, cncr_mal_clrt_colon) ~ get(exposure_var) +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },
    model_3 = function(df, exposure_var) {
      coxph(Surv(age, ageevent, cncr_mal_clrt_colon) ~ get(exposure_var) +
              bmi_adj + alc_re + smoke_stat +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    }
  )

  return(models)
}

#' `epic_somalogic_analysis()`: Perform analysis on EPIC somalogic data
#' @description
#' This function takes a processed list from `epic_somalogic_process_data()` and
#' performs Cox proportional hazards analysis using predefined models from
#' `epic_somalogic_models()`. It iterates over each combination of data frames
#' and exposures, fits the models, and stores the results in a data frame.
#' It prints the number of unique exposures and data frame combinations for sanity
#'
#' @param list_processed A list from `epic_somalogic_process_data()`
#' @return A data frame containing the results of Cox proportional hazards analysis
#'
epic_somalogic_analysis <- function(list_processed) {
  # Initialize an empty data frame to store results
  data_result <- data.frame(
    exposure = character(),
    data1 = character(),
    data2 = character(),
    data3 = character(),
    data4 = character(),
    model_type = character(),
    n = numeric(),
    nevent = numeric(),
    coef = numeric(),
    coef_exp = numeric(),
    se = numeric(),
    ci_lower_exp = numeric(),
    ci_upper_exp = numeric(),
    se_robust = numeric(),
    pval = numeric(),
    degrees_freedom = numeric(),
    concordance = numeric(),
    concordance_se = numeric(),
    likelihood_ratio_test = numeric(),
    likelihood_ratio_test_pval = numeric(),
    wald_test = numeric(),
    wald_test_pval = numeric(),
    score_test = numeric(),
    score_test_pval = numeric(),
    robust_score_test = numeric(),
    robust_score_test_pval = numeric(),
    exclusion_missing = numeric()
  )

  # Define a list of models
  models <- epic_somalogic_models()

  # Counters for unique values
  unique_exposures <- vector()
  unique_data1 <- vector()
  unique_data2 <- vector()
  unique_data3 <- vector()
  unique_data4 <- vector()
  unique_data5 <- vector()

  # Iterate over each list within list_processed
  for (i in seq_along(list_processed)) {
    list1_name <- names(list_processed)[i]
    cat("Processing list", list1_name, "\n")
    list_processed_i <- list_processed[[i]]

    for (j in seq_along(list_processed_i)) {
      list2_name <- names(list_processed_i)[j]
      cat("  -", list2_name, "\n")
      list_processed_j <- list_processed_i[[j]]

      for (k in seq_along(list_processed_j)) {
        list3_name <- names(list_processed_j)[k]
        cat("  -", list3_name, "\n")
        list_processed_k <- list_processed_j[[k]]

        for (l in seq_along(list_processed_k)) {
          data_name <- names(list_processed_k)[l]
          cat("  -", data_name, "\n")
          data <- list_processed_k[[l]]

          LABEL <- paste0(list1_name, "_", list2_name, "_", list3_name, "_", data_name)

          for (exposure_var in names(data)[grepl("seq.10000", names(data))]) {
            for (model_name in names(models)) {
              cat(paste0("      + running:", model_name, "; exposure:", exposure_var, "; data:", LABEL, "\n"))
              model_result <- models[[model_name]](data, exposure_var)

              data_result <- rbind(data_result, data.frame(
                exposure = exposure_var,
                data1 = list1_name,
                data2 = list2_name,
                data3 = list3_name,
                data4 = data_name,
                model_type = as.numeric(gsub("model_", "", model_name)),
                n = model_result$n,
                nevent = model_result$nevent,
                coef = summary(model_result)$coef[1, 1],
                coef_exp = summary(model_result)$conf.int[1, 1],
                se = summary(model_result)$coef[1, 3],
                ci_lower_exp = summary(model_result)$conf.int[1, 3],
                ci_upper_exp = summary(model_result)$conf.int[1, 4],
                se_robust = summary(model_result)$coef[1, 4],
                pval = summary(model_result)$coef[1, 6],
                degrees_freedom = summary(model_result)$logtest[["df"]],
                concordance = summary(model_result)$concordance[["C"]],
                concordance_se = summary(model_result)$concordance[["se(C)"]],
                likelihood_ratio_test = summary(model_result)$logtest[["test"]],
                likelihood_ratio_test_pval = summary(model_result)$logtest[["pvalue"]],
                wald_test = summary(model_result)$waldtest[["test"]],
                wald_test_pval = summary(model_result)$waldtest[["pvalue"]],
                score_test = summary(model_result)$sctest[["test"]],
                score_test_pval = summary(model_result)$sctest[["pvalue"]],
                robust_score_test = summary(model_result)$robscore[["test"]],
                robust_score_test_pval = summary(model_result)$robscore[["pvalue"]],
                exclusion_missing = length(model_result$na.action)
              ))

              # Track unique values
              unique_exposures <- unique(c(unique_exposures, exposure_var))
              unique_data1 <- unique(c(unique_data1, list1_name))
              unique_data2 <- unique(c(unique_data2, list2_name))
              unique_data3 <- unique(c(unique_data3, list3_name))
              unique_data4 <- unique(c(unique_data4, data_name))
            }
          }
        }
      }
    }
  }
  # END: print info
  cat("# END")
  cat("## N exposures:", length(unique_exposures), "\n")
  cat("## N data1:", length(unique_data1), "\n")
  cat("## N data2:", length(unique_data2), "\n")
  cat("## N data3:", length(unique_data3), "\n")
  cat("## N data4:", length(unique_data4), "\n")
  cat("## N analyses for 1 model:", length(unique_exposures) * length(unique_data1) * length(unique_data2) * length(unique_data3) * length(unique_data4), "\n")
  # Now data_result contains the desired output
  return(data_result)
}

