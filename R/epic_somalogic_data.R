#' `epic_somalogic_format_data()`: process standard variables
#' @description
#' This function processes EPIC-SOMALOGIC data by calculating follow-up time,
#' converting certain categorical variables to numeric, and reordering factor levels.
#'
#' @param df A dataframe containing EPIC-SOMALOGIC data.
#'
#' @return A processed dataframe with calculated follow-up time, converted categorical variables, and reordered factor levels.
#'
#' @details
#' The function calculates the follow-up time by subtracting the age at cancer onset
#' from the current age. It converts the "cvd_t2d_coh" variable to numeric, where "No"
#' is converted to 0 and "Yes" to 1. It creates age groups with a 5-year interval
#' and reorders factor levels for variables "smoke_stat", "pa_index", and "l_school".
#'
#'
#' @export
epic_somalogic_format_data <- function(df) {
  df$followup <- df$ageevent - df$age
  df <- df %>%
    mutate(cvd_t2d_coh = ifelse(cvd_t2d_coh == "No", 0, 1)) %>%
    mutate(age_group = as.integer(factor(cut(age,
                                             breaks = seq(0, max(age) + 5, by = 5),
                                             right = FALSE,
                                             labels = seq(0, max(age), by = 5)))))  # Update labels
  df$smoke_stat <- factor(df$smoke_stat, levels = c("Never", "Former", "Smoker"))
  df$pa_index <- factor(df$pa_index, levels = c("Inactive", "Moderately inactive", "Moderately active", "Active", "Missing"))
  df$l_school <- factor(df$l_school, levels = c("None", "Primary school completed", "Secondary school", "Technical/professional school", "Longer education (incl. University deg.)", "Not specified"))
  return(df)
}

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
epic_somalogic_create_subtypes <- function(df) {
  # Function to create subtypes within a dataframe
  # Subset for cases and non-cases
  data_cancer <- subset(df, cncr_mal_clrt == 1)
  data_cancer$indevent <- 1
  data_control <- subset(df, cncr_mal_clrt == 0)
  data_control$indevent <- 0
  data_overall <- bind_rows(data_cancer, data_control)
  data_overall <- droplevels(data_overall)

  # Subset for subtypes - Colon
  data_cancer_colon <- subset(df, cncr_mal_clrt_colon == 1)
  data_cancer_colon$indevent <- 1
  data_control_colon <- subset(df, cncr_mal_clrt_colon == 0)
  data_control_colon$indevent <- 0
  data_colon <- bind_rows(data_cancer_colon, data_control_colon)
  data_colon <- droplevels(data_colon)

  # Subset for subtypes - Rectum
  data_cancer_rectum <- subset(df, cncr_mal_clrt_rectum == 1)
  data_cancer_rectum$indevent <- 1
  data_control_rectum <- subset(df, cncr_mal_clrt_rectum == 0)
  data_control_rectum$indevent <- 0
  data_rectum <- bind_rows(data_cancer_rectum, data_control_rectum)
  data_rectum <- droplevels(data_rectum)

  # Subset for subtypes - Early Onset
  data_cancer_eo <- subset(data_cancer, ageevent < 55)
  data_control_eo <- subset(data_control, age < 55)
  data_control_eo <- data_control_eo %>%
    mutate(ageevent = ifelse(ageevent > 55, 55, ageevent))
  data_eo <- bind_rows(data_cancer_eo, data_control_eo)
  data_eo <- droplevels(data_eo)

  # Update the 'ID' column
  data_overall$ID <- paste0(data_overall$ID, ";overall")
  data_colon$ID <- paste0(data_colon$ID, ";colon")
  data_rectum$ID <- paste0(data_rectum$ID, ";rectum")
  data_eo$ID <- paste0(data_eo$ID, ";earlyonset")

  # Return the list of subtypes
  return(list(
    overall = data_overall,
    colon = data_colon,
    rectum = data_rectum,
    earlyonset = data_eo
  ))
}

#' `epic_somalogic_create_sex_specific_dataframes()`: create sex-specific dataframes for a given 'ID_col' group for Epic SomaLogic data
#' @description
#' This function takes a subset of data for a specific 'ID_col' group and creates
#' sex-specific dataframes - sex-combined, male, and female.  Must use `lapply(list_data, epic_somalogic_create_sex)`
#' if you want to run this over a list of dataframes.
#'
#' @param data_group A subset of the original dataframe for a specific 'ID_col' group.
#' @return A list containing sex-specific dataframes.
#'
#' @export
epic_somalogic_create_sex <- function(df) {
  combined <- df
  female <- split(df, df$sex)[[1]]
  male <- split(df, df$sex)[[2]]
  return(list(combined = combined, female = female, male = male))
}

#' `epic_somalogic_prentice_weighting()`: Prentice weighting function for EPIC somalogic analysis
#' @description
#' This function applies Prentice weighting to the provided data based on follow-up years.
#' It filters out cases where the age at event is less than the age plus follow-up years,
#' and then creates different subsets of the data based on disease status and cohort membership.
#'
#' @param data A data frame containing the EPIC somalogic data
#' @param followup_years The number of follow-up years used for weighting
#' @return A processed data frame with Prentice weighting applied
epic_somalogic_prentice_weighting <- function(df, followup_years) {
  df <- df %>% filter(ageevent > age + followup_years) %>% mutate(age = age + followup_years)
  cat("* N samples before Prentice weighting (follow-up:", followup_years, "years):", nrow(df), "\n")
  CasesOutsideSubCohort <- df %>% filter(cvd_t2d_coh == 0 & indevent == 1) %>% mutate(age = ageevent - 1e-4)
  CasesInSubcohort <- df %>% filter(cvd_t2d_coh == 1 & indevent == 1)
  ControlsSubcohort <- df %>% filter(cvd_t2d_coh == 1 & indevent == 0)
  cat("* N samples after Prentice weighting (follow-up:", followup_years, "years):", nrow(bind_rows(CasesOutsideSubCohort, ControlsSubcohort, ControlsSubcohort)), "\n")
  df <- bind_rows(CasesOutsideSubCohort, CasesInSubcohort, ControlsSubcohort)
  return(df)
}

#' `epic_somalogic_create_analysis()`: process EPIC somalogic data for analysis
#' @description
#' This function preprocesses EPIC somalogic data for analysis by applying
#' Prentice weighting using `epic_somalogic_process_data()`.
#'
#' @param list_df A nested list containing EPIC somalogic data frames
#' @return A processed nested list of EPIC somalogic data frames for analysis
epic_somalogic_create_analysis <- function(list_df) {
  list_processed <- list()  # Initialize the list to store processed data

  # Loop through each nested list in the main list
  for (i in seq_along(list_df)) {
    df_name <- names(list_df)[i]  # Get the name of the nested list

    cat("# Processing nested list:", df_name, "\n")

    nested_list <- list_df[[i]]
    processed_nested_list <- list()  # Initialize a nested list to store processed dataframes

    # Loop through each dataframe in the nested list
    for (j in seq_along(nested_list)) {
      df_name_nested <- names(nested_list)[j]  # Get the name of the dataframe
      df <- nested_list[[j]]  # Get the dataframe

      cat("## Processing dataframe:", df_name_nested, "\n")

      # Prentice weighting
      processed_df <- list(
        analysis = epic_somalogic_prentice_weighting(df, 0),
        followup_excl2year = epic_somalogic_prentice_weighting(df, 2),
        followup_excl5year = epic_somalogic_prentice_weighting(df, 5)
      )

      processed_nested_list[[df_name_nested]] <- processed_df
    }

    # Store the processed nested list in the main processed list with its name
    list_processed[[df_name]] <- processed_nested_list
  }

  return(list_processed)
}

#' `epic_somalogic_make_datacomplete()`: process missing values in EPIC somalogic data
#' @description
#' This function converts missing values in EPIC somalogic data to -9 and
#' calculates the percentage of missing values for each variable.
#'
#' @param df_complete A complete data frame containing EPIC somalogic data
#' @param followup_years The number of follow-up years for the analysis
#' @return A processed data frame with missing values replaced by -9
epic_somalogic_make_datacomplete <- function(df_complete, followup_years) {
  cat("## Converting missing values for data_analysis to -9 (follow-up:", followup_years, "years)", "\n")
  df_complete$bmi_missing <- ifelse(is.na(df_complete$bmi_c), 1, 0)
  df_complete$bmi_c[is.na(df_complete$bmi_c)] <- -9
  cat("* bmi_c:", round(table(df_complete$bmi_missing)[2] / nrow(df_complete) * 100, 2), "% missing", "\n")
  df_complete$alc_re_missing <- ifelse(is.na(df_complete$alc_re), 1, 0)
  df_complete$alc_re[is.na(df_complete$alc_re)] <- -9
  cat("* alc_re:", round(table(df_complete$alc_re_missing)[2] / nrow(df_complete) * 100, 2), "% missing", "\n")
  df_complete$pa_mets_missing <- ifelse(is.na(df_complete$pa_mets), 1, 0)
  df_complete$pa_mets[is.na(df_complete$pa_mets)] <- -9
  cat("* pa_mets:", round(table(df_complete$pa_mets_missing)[2] / nrow(df_complete) * 100, 2), "% missing", "\n")
  df_complete$QE_ENERGY_missing <- ifelse(is.na(df_complete$QE_ENERGY), 1, 0)
  df_complete$QE_ENERGY[is.na(df_complete$QE_ENERGY)] <- -9
  cat("* QE_ENERGY:", round(table(df_complete$QE_ENERGY_missing)[2] / nrow(df_complete) * 100, 2), "% missing", "\n")
  df_complete$QGE0701_missing <- ifelse(is.na(df_complete$QGE0701), 1, 0)
  df_complete$QGE0701[is.na(df_complete$QGE0701)] <- -9
  cat("* QGE0701:", round(table(df_complete$QGE0701_missing)[2] / nrow(df_complete) * 100, 2), "% missing", "\n")
  df_complete$QGE0704_missing <- ifelse(is.na(df_complete$QGE0704), 1, 0)
  df_complete$QGE0704[is.na(df_complete$QGE0704)] <- -9
  cat("* QGE0704:", round(table(df_complete$QGE0704_missing)[2] / nrow(df_complete) * 100, 2), "% missing", "\n")
  df_complete$QE_FIBT_missing <- ifelse(is.na(df_complete$QE_FIBT), 1, 0)
  df_complete$QE_FIBT[is.na(df_complete$QE_FIBT)] <- -9
  cat("* QE_FIBT:", round(table(df_complete$QE_FIBT_missing)[2] / nrow(df_complete) * 100, 2), "% missing", "\n")
  df_complete$QE_CA_missing <- ifelse(is.na(df_complete$QE_CA), 1, 0)
  df_complete$QE_CA[is.na(df_complete$QE_CA)] <- -9
  cat("* QE_CA:", round(table(df_complete$QE_CA_missing)[2] / nrow(df_complete) * 100, 2), "% missing", "\n")

  return(df_complete)
}

#' `epic_somalogic_process_data()`: creates nested raw and complete dataframes
#' @description
#' the function applies `epic_somalogic_create_datacomplete()` to a nested list and creates
#' a nested list of raw and complete data frames
#'
#' @param list_df A nested list containing EPIC-SOMALOGIC dataframes.
#'
#' @return A nested list of processed EPIC-SOMALOGIC dataframes with calculated follow-up time,
#' converted categorical variables, and reordered factor levels.
#'
#' @export
epic_somalogic_process_data <- function(list_df) {
  list_processed <- list()  # Initialize the list to store processed data

  # Loop through each nested list in the main list
  for (i in seq_along(list_df)) {
    name_nested_list <- names(list_df)[i]  # Get the name of the nested list
    cat("# Processing list:", name_nested_list, "\n")
    nested_list <- list_df[[i]]
    processed_nested_list <- list()  # Initialize a nested list to store processed dataframes

    # Loop through each dataframe in the nested list
    for (j in seq_along(nested_list)) {
      name_nested_list1 <- names(nested_list)[j]  # Get the name of the nested list
      cat("## Processing dataframe:", name_nested_list1, "\n")
      nested_list1 <- nested_list[[j]]

      # Make list: processed data
      processed_df <- list(
        data_analysis = list(
          data_raw = nested_list1$analysis,
          data_complete = epic_somalogic_make_datacomplete(nested_list1$analysis, 0)
        ),
        followup_excl2year = list(
          data_raw = nested_list1$followup_excl2year,
          data_complete = epic_somalogic_make_datacomplete(nested_list1$followup_excl2year, 2)
        ),
        followup_excl5year = list(
          data_raw = nested_list1$followup_excl5year,
          data_complete = epic_somalogic_make_datacomplete(nested_list1$followup_excl5year, 5)
        )
      )

      processed_nested_list[[name_nested_list1]] <- processed_df
    }

    # Store the processed nested list in the main processed list with its name
    list_processed[[name_nested_list]] <- processed_nested_list
  }

  # Return the processed data list
  return(list_processed)
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
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },
    model_1_extra = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) +
              fasting_c_misscat + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_cos_4_pi + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_sin_4_pi +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },

    model_2 = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) +
              bmi_c + alc_re + smoke_stat + pa_index + l_school +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },
    model_2_extra = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) +
              bmi_c + alc_re + smoke_stat + pa_index + l_school +
              fasting_c_misscat + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_cos_4_pi + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_sin_4_pi +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },

    model_3 = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) +
              bmi_c + alc_re + smoke_stat + pa_index + l_school +
              QE_ENERGY + QGE0701 + QGE0704 + QE_FIBT + QE_CA +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },
    model_3_extra = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) +
              bmi_c + alc_re + smoke_stat + pa_index + l_school +
              QE_ENERGY + QGE0701 + QGE0704 + QE_FIBT + QE_CA +
              fasting_c_misscat + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_cos_4_pi + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_sin_4_pi +
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

          for (exposure_var in names(data)[grepl("seq", names(data))]) {
            if (any(!is.na(data[[exposure_var]]))) {
              for (model_name in names(models)) {
                cat(paste0("      + running: ", model_name, "; exposure:", exposure_var, "; data:", LABEL, "\n"))
                model_result <- models[[model_name]](data, exposure_var)

                tryCatch({
                  coef_summary <- summary(model_result)$coef
                  pval <- if (nrow(coef_summary) >= 1 && ncol(coef_summary) >= 6) {
                    coef_summary[1, 6]
                  } else {
                    NA
                  }
                }, error = function(e) {
                  print(paste("Error occurred:", e))
                  pval <- NA
                })

                # Construct data frame
                data_result <- rbind(data_result, data.frame(
                  exposure = exposure_var,
                  data1 = list1_name,
                  data2 = list2_name,
                  data3 = list3_name,
                  data4 = data_name,
                  model_type = model_name,
                  n = model_result$n,
                  nevent = model_result$nevent,
                  coef = ifelse(is.numeric(summary(model_result)$coef[1,1]), summary(model_result)$coef[1,1], NA),
                  coef_exp = ifelse(is.numeric(summary(model_result)$conf.int[1,1]), summary(model_result)$conf.int[1,1], NA),
                  se = ifelse(is.numeric(summary(model_result)$coef[1,3]), summary(model_result)$coef[1,3], NA),
                  ci_lower_exp = ifelse(is.numeric(summary(model_result)$conf.int[1,3]), summary(model_result)$conf.int[1,3], NA),
                  ci_upper_exp = ifelse(is.numeric(summary(model_result)$conf.int[1,4]), summary(model_result)$conf.int[1,4], NA),
                  se_robust = ifelse(is.numeric(summary(model_result)$coef[1,4]), summary(model_result)$coef[1,4], NA),
                  pval = pval,
                  degrees_freedom = ifelse(is.numeric(summary(model_result)$logtest[["df"]]), summary(model_result)$logtest[["df"]], NA),
                  concordance = ifelse(is.numeric(summary(model_result)$concordance[["C"]]), summary(model_result)$concordance[["C"]], NA),
                  concordance_se = ifelse(is.numeric(summary(model_result)$concordance[["se(C)"]]), summary(model_result)$concordance[["se(C)"]], NA),
                  likelihood_ratio_test = ifelse(is.numeric(summary(model_result)$logtest[["test"]]), summary(model_result)$logtest[["test"]], NA),
                  likelihood_ratio_test_pval = ifelse(is.numeric(summary(model_result)$logtest[["pvalue"]]), summary(model_result)$logtest[["pvalue"]], NA),
                  wald_test = ifelse(is.numeric(summary(model_result)$waldtest[["test"]]), summary(model_result)$waldtest[["test"]], NA),
                  wald_test_pval = ifelse(is.numeric(summary(model_result)$waldtest[["pvalue"]]), summary(model_result)$waldtest[["pvalue"]], NA),
                  score_test = ifelse(is.numeric(summary(model_result)$sctest[["test"]]), summary(model_result)$sctest[["test"]], NA),
                  score_test_pval = ifelse(is.numeric(summary(model_result)$sctest[["pvalue"]]), summary(model_result)$sctest[["pvalue"]], NA),
                  robust_score_test = ifelse(is.numeric(summary(model_result)$robscore[["test"]]), summary(model_result)$robscore[["test"]], NA),
                  robust_score_test_pval = ifelse(is.numeric(summary(model_result)$robscore[["pvalue"]]), summary(model_result)$robscore[["pvalue"]], NA),
                  exclusion_missing = ifelse(is.numeric(length(model_result$na.action)), length(model_result$na.action), NA)
                ))
              }
            } else {
              # Add a row of NA values if exposure is completely missing
              data_result <- rbind(data_result, data.frame(
                exposure = exposure_var,
                data1 = list1_name,
                data2 = list2_name,
                data3 = list3_name,
                data4 = data_name,
                model_type = NA,
                n = NA,
                nevent = NA,
                coef = NA,
                coef_exp = NA,
                se = NA,
                ci_lower_exp = NA,
                ci_upper_exp = NA,
                se_robust = NA,
                pval = NA,
                degrees_freedom = NA,
                concordance = NA,
                concordance_se = NA,
                likelihood_ratio_test = NA,
                likelihood_ratio_test_pval = NA,
                wald_test = NA,
                wald_test_pval = NA,
                score_test = NA,
                score_test_pval = NA,
                robust_score_test = NA,
                robust_score_test_pval = NA,
                exclusion_missing = NA
              ))
            }
          }
        }
      }
    }
  }
  # END: print info
  cat("# END\n")
  # Now data_result contains the desired output
  return(data_result)
}


#' `epic_somalogic_models_heterogeneity_sex()`: a list of coxph models
#' @description
#' This function defines three Cox proportional hazards models including
#'
#' @return A list containing three Cox proportional hazards models
#'
epic_somalogic_models_heterogeneity_sex <- function() {
  # Define the list of models
  models <- list(
    model_1 = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },
    model_1_heterogeneity = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) + get(exposure_var):sex +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },
    model_1_extra = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) +
              fasting_c_misscat + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_cos_4_pi + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_sin_4_pi +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },
    model_1_heterogeneity_extra = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) + get(exposure_var):sex +
              fasting_c_misscat + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_cos_4_pi + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_sin_4_pi +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },

    model_2 = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) +
              bmi_c + alc_re + smoke_stat + pa_index + l_school +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },
    model_2_heterogeneity = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) + get(exposure_var):sex +
              bmi_c + alc_re + smoke_stat + pa_index + l_school +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },
    model_2_extra = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) +
              bmi_c + alc_re + smoke_stat + pa_index + l_school +
              fasting_c_misscat + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_cos_4_pi + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_sin_4_pi +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },
    model_2_heterogeneity_extra = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) + get(exposure_var):sex +
              bmi_c + alc_re + smoke_stat + pa_index + l_school +
              fasting_c_misscat + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_cos_4_pi + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_sin_4_pi +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },

    model_3 = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) +
              bmi_c + alc_re + smoke_stat + pa_index + l_school +
              QE_ENERGY + QGE0701 + QGE0704 + QE_FIBT + QE_CA +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },
    model_3_heterogeneity = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) + get(exposure_var):sex +
              bmi_c + alc_re + smoke_stat + pa_index + l_school +
              QE_ENERGY + QGE0701 + QGE0704 + QE_FIBT + QE_CA +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },
    model_3_extra = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) +
              bmi_c + alc_re + smoke_stat + pa_index + l_school +
              QE_ENERGY + QGE0701 + QGE0704 + QE_FIBT + QE_CA +
              fasting_c_misscat + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_cos_4_pi + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_sin_4_pi +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    },
    model_3_heterogeneity_extra = function(df, exposure_var) {
      coxph(Surv(age, ageevent, indevent) ~ get(exposure_var) + get(exposure_var):sex +
              bmi_c + alc_re + smoke_stat + pa_index + l_school +
              QE_ENERGY + QGE0701 + QGE0704 + QE_FIBT + QE_CA +
              fasting_c_misscat + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_cos_4_pi + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_sin_4_pi +
              strata(sex, center, age_group) +
              cluster(idepic),
            data = df)
    }

  )

  return(models)
}


#' `epic_somalogic_analysis_heterogeneity()`: Perform analysis on EPIC somalogic data
#' @description
#' This function takes a processed list from `epic_somalogic_process_data()` and
#' performs Cox proportional hazards analysis using predefined models from
#' `epic_somalogic_models_heterogeneity_sex()`. It iterates over each combination of data frames
#' and exposures, fits the models, and stores the results in a data frame.
#' It prints the number of unique exposures and data frame combinations for sanity
#'
#' @param list_processed A list from `epic_somalogic_process_data()`
#' @return A data frame containing the results of Cox proportional hazards analysis
#'
epic_somalogic_analysis_heterogeneity <- function(list_processed) {
  # Initialize an empty data frame to store results
  data_result <- data.frame(
    exposure = character(),
    data1 = character(),
    data2 = character(),
    data3 = character(),
    data4 = character(),
    heterogeneity_pval_model1 = numeric(),
    heterogeneity_pval_model2 = numeric(),
    heterogeneity_pval_model3 = numeric(),
    heterogeneity_pval_model1_extra = numeric(),
    heterogeneity_pval_model2_extra = numeric(),
    heterogeneity_pval_model3_extra = numeric()
  )

  # Define a list of models
  models <- epic_somalogic_models_heterogeneity_sex()

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

          for (exposure_var in names(data)[grepl("seq", names(data))]) {
            if (any(!is.na(data[[exposure_var]]))) {
              cat(paste0("      + running: heterogeneity; exposure: ", exposure_var, "; data: ", LABEL, "\n"))

              # run normal and interaction models within likelihood ratio test and extract p-value for heterogeneity
              lrtest_m1 <- lmtest::lrtest(models[[1]](data, exposure_var),
                                          models[[2]](data, exposure_var))$Pr[2]
              lrtest_m2 <- lmtest::lrtest(models[[5]](data, exposure_var),
                                          models[[6]](data, exposure_var))$Pr[2]
              lrtest_m3 <- lmtest::lrtest(models[[9]](data, exposure_var),
                                          models[[10]](data, exposure_var))$Pr[2]

              lrtest_m1_extra <- lmtest::lrtest(models[[3]](data, exposure_var),
                                                models[[4]](data, exposure_var))$Pr[2]
              lrtest_m2_extra <- lmtest::lrtest(models[[7]](data, exposure_var),
                                                models[[8]](data, exposure_var))$Pr[2]
              lrtest_m3_extra <- lmtest::lrtest(models[[11]](data, exposure_var),
                                                models[[12]](data, exposure_var))$Pr[2]

              # Construct data frame
              data_result <- rbind(data_result, data.frame(
                exposure = exposure_var,
                data1 = list1_name,
                data2 = list2_name,
                data3 = list3_name,
                data4 = data_name,
                heterogeneity_pval_model1 = lrtest_m1,
                heterogeneity_pval_model2 = lrtest_m2,
                heterogeneity_pval_model3 = lrtest_m3,
                heterogeneity_pval_model1_extra = lrtest_m1_extra,
                heterogeneity_pval_model2_extra = lrtest_m2_extra,
                heterogeneity_pval_model3_extra = lrtest_m3_extra
              ))
            } else {
              # Add a row of NA values if exposure is completely missing
              data_result <- rbind(data_result, data.frame(
                exposure = exposure_var,
                data1 = list1_name,
                data2 = list2_name,
                data3 = list3_name,
                data4 = data_name,
                heterogeneity_pval_model1 = NA,
                heterogeneity_pval_model2 = NA,
                heterogeneity_pval_model3 = NA,
                heterogeneity_pval_model1_extra = NA,
                heterogeneity_pval_model2_extra = NA,
                heterogeneity_pval_model3_extra = NA
              ))
            }
          }
        }
      }
    }
  }

  # END: print info
  cat("# END\n")

  # Now data_result contains the desired output
  return(data_result)
}


#' Perform enrichment analysis
#'
#' This function conducts enrichment analysis on the provided data using the enricher function.
#'
#' @param data A data frame containing the data for analysis.
#' @param sig_direction The direction of significance (e.g., "up" or "down").
#' @param group_across The column name representing the grouping variable across which enrichment analysis is performed.
#' @param group_within The column name representing the grouping variable within which enrichment analysis is performed.
#' @param universe_geneinfo A vector containing gene information for the entire universe.
#' @param database A data frame containing the database for enrichment analysis.
#' @param minGSSize The minimum gene set size.
#' @param maxGSSize The maximum gene set size.
#' @param pvalueCutoff The p-value cutoff for significance.
#' @param qvalueCutoff The q-value cutoff for significance.
#' @param file_path The file path where results and plots will be saved.
#' @return NULL
#' @export
#'
enrichment_analysis <- function(data, sig_direction, group_across, group_within,
                                universe_geneinfo, database,
                                minGSSize, maxGSSize, pvalueCutoff, qvalueCutoff, file_path) {
  for (i in unique(data[[group_across]])) {
    # Run enrichment analysis using the enricher function
    results_enricher <- data %>%
      filter(.data[[group_across]] == i) %>%
      group_by(.data[[group_within]]) %>%
      do({
        dat <- .
        pull(dat, gene) %>%
          enricher(universe = universe_geneinfo,
                   TERM2GENE = database,
                   minGSSize = minGSSize,
                   maxGSSize = maxGSSize,
                   pvalueCutoff = pvalueCutoff,
                   qvalueCutoff = qvalueCutoff) %>%
          as_tibble()
      })

    # Write results to file
    write.table(results_enricher, file = paste0(file_path, "results_", sig_direction, "_across-", i, "_within-", group_within, ".txt"),
                row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

    # Process and store the results
    results_enricher_wide_sig <- results_enricher %>%
      filter(pvalue < 0.05) %>%
      select(all_of(group_within), ID, pvalue) %>%
      mutate(`log10(p-adjust)` = -log10(pvalue)) %>%
      select(-pvalue) %>%
      spread(all_of(group_within), `log10(p-adjust)`, fill = 0) %>%
      column_to_rownames("ID")

    order <- results_enricher_wide_sig %>%
      dist(method = "euclidean") %>%
      hclust(method = "ward.D2") %>%
      with(labels[order])

    # Generate bubble heatmap and store it in the results list
    plot_i <- results_enricher %>%
      filter(pvalue < 0.05) %>%
      mutate(Term = factor(ID, levels = order),
             Cancer = factor(all_of(group_within))) %>%
      ggplot(aes(Term, Cancer, size = Count, color = -log10(pvalue))) +
      coord_flip() +
      geom_point() +
      scale_color_gradient(low = "grey", high = "red") +
      theme_bw() +
      theme(axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8,
                                       angle = 90,
                                       vjust = 0.5,
                                       hjust = 1))

    tiff(paste0(file_path, "plot-enrichment_", sig_direction, "_across-", i, "_within-", group_within, ".tiff"),
         height = 50, width = 20, units = "cm", res = 300, compression = "lzw")
    print(plot_i)
    dev.off()
  }
}
