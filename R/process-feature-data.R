#' Process feature data
#'
#' This function processes feature data given specified metadata.
#' It supports exclusion of features with extreme missingness,
#' various imputation methods, transformation, plate correction, centering,
#' and case-control data handling.
#'
#' @param data A data frame with the first column as sample IDs and remaining columns containing feature values, where feature IDs are the column names.
#' @param data_meta_features A data frame containing metadata for the features, including information such as limit of detection and missingness percentage.
#' @param data_meta_samples A data frame containing metadata for the samples, including information necessary for plate correction and case-control analysis.
#' @param col_samples A string specifying the column name in \code{data} that contains sample IDs (e.g., \code{"Idepic_Bio"}).
#' @param col_features A string specifying the column name in \code{data} that contains feature IDs, which should match the column names of \code{data} (e.g., \code{"UNIPROT"}).
#' @param save A logical for whether you want to save the feature data, plots, and exclusion info. Default is \code{FALSE}/.
#' @param path_out A string specifying the output directory where the processed data will be saved.
#' @param path_outliers A string specifying the output directory where outlier information will be saved.
#' @param exclusion_extreme_feature A logical flag indicating whether to exclude features with extreme missingness. Default is \code{FALSE}.
#' @param col_missing_feature A string specifying the column name in \code{data_meta_features} that contains the percentage of missing values for each feature (e.g., \code{"missing_pct"}).
#' @param missing_pct_feature A numeric value specifying the threshold percentage for missingness above which features will be excluded (e.g., \code{0.9}).
#' @param exclusion_extreme_sample A logical flag indicating whether to exclude samples with extreme missingness. Default is \code{FALSE}.
#' @param col_missing_sample A string specifying the column name in \code{data_meta_samples} that contains the percentage of missing values for each sample (e.g., \code{"missing_pct_samples"}).
#' @param missing_pct_sample A numeric value specifying the threshold percentage for missingness above which samples will be excluded (e.g., \code{0.9}).
#' @param imputation A logical flag indicating whether imputation should be performed. Default is \code{FALSE}.
#' @param imputation_method A string specifying the method to use for imputation. Options include \code{"LOD"}, \code{"1/5th"}, \code{"KNN"}, \code{"ppca"}, \code{"median"}, \code{"mean"}, \code{"rf"}, and \code{"left-censored"}.
#' @param col_LOD A string specifying the column name in \code{data_meta_features} that contains the limit of detection (LOD) values, required if \code{imputation_method} is \code{"LOD"}.
#' @param transformation A logical flag indicating whether transformation should be performed. Default is \code{FALSE}.
#' @param transformation_method A string specifying the method to use for transformation. Options include \code{"InvRank"}, \code{"Log10"}, \code{"Log10Capped"}, and \code{"Log10ExclExtremes"}.
#' @param outlier A logical flag indicating whether outlier exclusion should be performed across features and samples. Default is \code{FALSE}.
#' @param plate_correction A logical flag indicating whether plate correction should be performed. Default is \code{FALSE}.
#' @param cols_listRandom A string or vector specifying columns in \code{data_meta_samples} to be treated as random effects in plate correction (e.g., \code{"batch_plate"}).
#' @param cols_listFixedToKeep A vector specifying columns in \code{data_meta_samples} to be treated as fixed effects in plate correction and retained in the model (e.g., \code{c("Center", "Country")}).
#' @param cols_listFixedToRemove A vector specifying columns in \code{data_meta_samples} to be treated as fixed effects in plate correction and removed from the model. Default is \code{NULL}.
#' @param col_HeteroSked A string specifying the column in \code{data_meta_samples} to be used for heteroskedasticity correction.
#' @param case_control A logical flag indicating whether the data are case-control and if matched samples should be handled accordingly. Default is \code{FALSE}.
#' @param col_case_control A string specifying the column name in \code{data_meta_samples} that contains case-control matching information (e.g., \code{"Match_Caseset"}).
#' @param centre_scale A logical flag indicating whether to center and scale the data. Default is \code{FALSE}.
#'
#' @return The processed data is returned and saved as a `.rds` file to the specified output directory.
#'
#' @details
#' This function performs several data processing steps (in order):
#' \itemize{
#'   \item Excludes features with extreme missingness based on a specified threshold.
#'   \item Excludes samples with extreme missingness based on a specified threshold.
#'   \item Imputes missing values using various methods.
#'   \item Transforms the data using specified methods.
#'   \item Excludes outlying samples and features using PCA and LOF.
#'   \item Handles case-control data to ensure matched samples are treated appropriately.
#'   \item Corrects for plate effects using specified random and fixed effects.
#'   \item Centers and scales the data if \code{centre_scale} is \code{TRUE}.
#' }
#' @export

process_data <- function(
    data,
    data_meta_features,
    data_meta_samples,
    col_samples,
    col_features,
    save = FALSE,
    path_out,
    path_outliers,
    exclusion_extreme_feature = FALSE, col_missing_feature = NULL, missing_pct_feature = NULL,
    exclusion_extreme_sample = FALSE, col_missing_sample = NULL, missing_pct_sample = NULL,
    imputation = FALSE, imputation_method = NULL, col_LOD = NULL,
    transformation = FALSE, transformation_method = NULL,
    outlier = FALSE,
    plate_correction = FALSE, cols_listRandom = NULL, cols_listFixedToKeep = NULL, cols_listFixedToRemove = NULL, col_HeteroSked = NULL,
    centre_scale = FALSE,
    case_control = FALSE, col_case_control = NULL
) {

  # df and VAR checks ====
  ## Check if specified columns exist in the data frames
  process_data_data_check(data = df,
                          data_meta_features = df_features,
                          data_meta_samples = df_samples,
                          col_samples = col_samples,
                          col_features = col_features,
                          exclusion_extreme_feature = exclusion_extreme_feature,
                          exclusion_extreme_sample = exclusion_extreme_sample,
                          col_missing_feature = col_missing_feature,
                          missing_pct_feature = missing_pct_feature,
                          col_missing_sample = col_missing_sample,
                          missing_pct_sample = missing_pct_sample)

  # make df ====
  df <- df %>%
    tibble::column_to_rownames(col_samples)
  df_samples <- tibble::as_tibble(df_samples)
  df_features <- tibble::as_tibble(df_features)

  # make labels ====
  labels <- process_data_generate_labels(exclusion_extreme_feature = exclusion_extreme_feature,
                                         missing_pct_feature = missing_pct_feature,
                                         exclusion_extreme_sample = exclusion_extreme_sample,
                                         missing_pct_sample = missing_pct_sample,
                                         imputation = imputation,
                                         imputation_method = imputation_method,
                                         transformation = transformation,
                                         transformation_method = transformation_method,
                                         outlier = outlier,
                                         plate_correction = plate_correction,
                                         centre_scale = centre_scale)

  LABEL_exclusion_extreme_feature = labels$LABEL_exclusion_extreme_feature
  LABEL_exclusion_extreme_sample = labels$LABEL_exclusion_extreme_sample
  LABEL_imputation = labels$LABEL_imputation
  LABEL_transformation = labels$LABEL_transformation
  LABEL_outlier = labels$LABEL_outlier
  LABEL_plate_correction = labels$LABEL_plate_correction
  LABEL_centre_scale = labels$LABEL_centre_scale

  # Extreme exclusion features ====
  if (exclusion_extreme_feature) {
    list_exclude_features <- process_data_exclude_features(df = df,
                                                           df_features = df_features,
                                                           missing_pct_feature = missing_pct_feature,
                                                           col_missing_feature = col_missing_feature,
                                                           col_features = col_features)
    df <- list_exclude_features$df
    df_features <- list_exclude_features$df_features
    id_features_exclude <- list_exclude_features$id_features_exclude
  }else{
    id_features_exclude <- NULL
  }

  # Extreme exclusion samples ====
  if (exclusion_extreme_sample) {
    list_exclude_samples <- process_data_exclude_samples(df = df,
                                                         df_samples = df_samples,
                                                         missing_pct_sample = missing_pct_sample,
                                                         col_missing_sample = col_missing_sample,
                                                         col_samples = col_samples)
    df <- list_exclude_samples$df
    df_samples <- list_exclude_samples$df_samples
    id_samples_exclude <- list_exclude_samples$id_samples_exclude
  }else{
    id_samples_exclude <- NULL
  }

  # imputation ====
  if (imputation) {
    df <- impute_data(df = df,
                              df_features = df_features,
                              imputation_method = imputation_method,
                              col_LOD = col_LOD,
                              col_features = col_features)
  }

  # transformation ====
  if (transformation) {
    list_transformation <- transform_data(df = df,
                                          transformation_method = transformation_method)
    df <- list_transformation$df
    LABEL_transformation <- list_transformation$transformation_method
    LABEL_centre_scale <- list_transformation$centre_scale
  }

  # outlier exclusion ====
  if(outlier){
    list_outliers <- outlier_pca_lof(df = df,
                                     col_samples= col_samples,
                                     df_samples = df_samples)
    df <- list_outliers$df
    df_samples <- list_outliers$df_samples
    plot_samples_outlier <- list_outliers$plot_samples_outlier
    id_samples_outlier <- list_outliers$id_samples_outlier
  }else{
    plot_samples_outlier <- NULL
    id_samples_outlier <- NULL
  }


  # case-control data ====
  if (case_control) {
    list_casecontrol <- filter_case_control(df = df,
                                            df_samples = df_samples,
                                            col_case_control = col_case_control,
                                            col_samples = col_samples)
    df <- list_casecontrol$df
    df_samples <- list_casecontrol$df_samples
    id_samples_casecontrol <- list_casecontrol$id_samples_casecontrol
  }else{
    id_samples_casecontrol <- NULL
  }

  # plate correction ====
  if (plate_correction) {
    cat("# Plate correction \n")
    # Check required columns for plate correction
    if (is.null(cols_listRandom) || is.null(cols_listFixedToKeep)) {
      stop("* Both 'cols_listRandom' and 'cols_listFixedToKeep' must be provided for plate correction.")
    }

    # make list for mixed model function
    list <- df %>%
      tibble::rownames_to_column(var = col_samples) %>%
      {
        list(data_features = .,
             data_samples = df_samples,
             data_meta_features = df_features %>%
               dplyr::rename(Name = tidyselect::all_of(col_features)))
      }

    ## transformation
    list <- functions::normalization_residualMixedModels(list = list,
                                                         forIdentifier = col_samples, # sample ID
                                                         listRandom = cols_listRandom, # variables to model as random effects; effects will be removed
                                                         listFixedToKeep = cols_listFixedToKeep, # variables to model as fixed effects; effects will be kept
                                                         listFixedToRemove = cols_listFixedToRemove, # variables to model as fixed effects; effects will be removed
                                                         HeteroSked = col_HeteroSked # variable for which heteroscedasticity will be accounted for
    )
    df <- list$data %>%
      dplyr::select(-tidyselect::all_of(c(cols_listRandom, cols_listFixedToKeep))) %>%
      tibble::column_to_rownames(col_samples)
  }

  # centre and scale ====
  if (centre_scale) {
    cat("# Centre and scale \n")
      df <- df %>%
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ scale(.)))
      LABEL_centre_scale <- TRUE
      cat("## data centred and scaled \n")
  }

  # LABEL ====
  LABEL <- paste0("exclusion-feature-", LABEL_exclusion_extreme_feature,
                  "_exclusion-sample-", LABEL_exclusion_extreme_sample,
                  "_imputation-", LABEL_imputation,
                  "_transformation-", LABEL_transformation,
                  "_outlier-", LABEL_outlier,
                  "_platecorrection-", LABEL_plate_correction,
                  "_centre-scale-", LABEL_centre_scale)

  # create id exclusions ====
  id_exclusions <- list(
    if (exists("id_features_exclude")) id_features_exclude else NULL,
    if (exists("id_samples_exclude")) id_samples_exclude else NULL,
    if (exists("id_samples_outlier")) id_samples_outlier else NULL,
    if (exists("id_samples_casecontrol")) id_samples_casecontrol else NULL
  )
  names(id_exclusions) <- c("id_features_exclude", "id_samples_exclude", "id_samples_outlier", "id_samples_casecontrol")

  cat("### total feature(s) excluded =", length(id_features_exclude), "\n")

  total_excluded_samples <- 0
  if (exists("id_samples_exclude")) {
    total_excluded_samples <- total_excluded_samples + length(id_samples_exclude)
  }
  if (exists("id_samples_outlier")) {
    total_excluded_samples <- total_excluded_samples + length(id_samples_outlier)
  }
  if (exists("id_samples_casecontrol")) {
    total_excluded_samples <- total_excluded_samples + length(id_samples_casecontrol)
  }
  cat("### total sample(s) excluded =", total_excluded_samples, "\n")

  # save ====
  if(save){
    cat("# Saving to: ", path_out, "\n")
    ## save outlier figure ====
    if(LABEL_outlier){
      cat("## saving outlier plots \n")
      tiff(paste0(path_outliers, "plot-outliers_", LABEL, ".tiff"), width = 1500, height = 1500, units = "px")
      print(cowplot::plot_grid(
        GGally::ggmatrix_gtable(plot_samples),
        GGally::ggmatrix_gtable(plot_features),
        nrow = 2))
      dev.off()
    }
    ## save excluded IDs ====
    cat("## saving exclusion info \n")
    saveRDS(object = id_exclusions, file = paste0(path_outliers, "id-exclusions_", LABEL, ".rds"))

    ## save feature data ====
    cat("## saving feature data \n")
    df_obj <- df %>%
      tibble::rownames_to_column(col_samples)
    saveRDS(object = df_obj, file = paste0(path_out, "data-features_", LABEL, ".rds"))
  }
  # return ====
  if(LABEL_outlier) {
    cat("# returning a list of data, outlier plot, and exclusion IDs")
    df <- df %>%
      tibble::rownames_to_column(col_samples)
    return(list(df = df, df_samples = df_samples, df_features = df_features, plot_samples_outlier = plot_samples_outlier, id_exclusions = id_exclusions))
  } else {
    cat("# returning a list of data and exclusion IDs")
    df <- df %>%
      tibble::rownames_to_column(col_samples)
    return(list(df = df, df_samples = df_samples, df_features = df_features, id_exclusions = id_exclusions))
  }
}

#' Normalize feature Data Using Residual Mixed Models
#'
#' This function normalizes feature data by removing unwanted effects using mixed models.
#' It accounts for random effects and fixed effects specified by the user, and optionally
#' corrects for heteroskedasticity. Residuals are calculated with `FUNnormalization_residualMixedModels()`.
#' The method is details in [https://www.mdpi.com/2218-1989/11/9/631](https://www.mdpi.com/2218-1989/11/9/631)
#'
#' @param list A list containing the following elements:
#' \describe{
#'   \item{data_features}{A data frame or tibble of dimensions n x (K+p) with:
#'     \describe{
#'       \item{n}{Number of observations.}
#'       \item{p}{Number of features.}
#'       \item{K}{Number of variables used for unique identification of individuals.}
#'     }}
#'   \item{data_samples}{A data frame or tibble of dimensions n x (K+d) with:
#'     \describe{
#'       \item{n}{Number of observations.}
#'       \item{K}{Number of unique identifiers.}
#'       \item{d}{Additional variables useful in final analysis (e.g., country, age, BMI).}
#'     }}
#'   \item{data_meta_features}{A p x 3 matrix indicating each feature's Name, Class, and Type.}
#' }
#' @param forIdentifier A character vector of strings indicating the names of variables used for unique identification of individuals.
#' @param listRandom A character vector of strings containing variable names modeled as random effects to be removed. If not NULL, should be either of length 1 or contain nested variables.
#' @param listFixedToKeep A character vector of strings containing variable names modeled as fixed effects to be kept.
#' @param listFixedToRemove A character vector of strings containing variable names modeled as fixed effects to be removed.
#' @param HeteroSked A string or NULL. If not NULL, the name of the variable for which heteroskedasticity will be accounted for. Must be included in `listRandom`.
#'
#' @return A list with:
#'   \item{data}{A tibble with unwanted variation removed.}
#'   \item{data_samples}{The input 'data_samples' data frame, ordered by IdentifierPipeline.}
#'   \item{data_meta_features}{The input 'data_meta_features' matrix.}
#'
#' @export
normalization_residualMixedModels <- function(list,
                                              forIdentifier = c("Study", "Batch", "CaseSet", "Idepic"),
                                              listRandom = NULL,
                                              listFixedToKeep = NULL,
                                              listFixedToRemove = NULL,
                                              HeteroSked = NULL) {

  # Convert list of data frames to tibble and unite identifiers ====
  data_features <- tibble::as_tibble(list$data_features) %>%
    tidyr::unite(IdentifierPipeline, tidyselect::all_of(forIdentifier), sep = "", remove = FALSE) %>%
    dplyr::arrange(IdentifierPipeline)
  data_samples <- tibble::as_tibble(list$data_samples) %>%
    tidyr::unite(IdentifierPipeline, tidyselect::all_of(forIdentifier), sep = "", remove = FALSE) %>%
    dplyr::arrange(IdentifierPipeline)
  data_meta_features <- tibble::as_tibble(list$data_meta_features)

  # Check for consistency between data_features and data_samples ====
  if (sum(data_features$IdentifierPipeline == data_samples$IdentifierPipeline) < nrow(data_features)) {
    stop("col_samples should be the same in data_features and data_samples")
  }

  # Merge data_features with data_samples and select relevant columns ====
  var.context <- unique(c("IdentifierPipeline", forIdentifier, listRandom, listFixedToKeep, listFixedToRemove))
  data <- dplyr::left_join(data_features,
                           data_samples %>%
                             dplyr::select(tidyselect::all_of(var.context)) %>%
                             dplyr::select(-tidyselect::any_of(forIdentifier)),
                           by = "IdentifierPipeline") %>%
    dplyr::select(tidyselect::all_of(var.context), tidyselect::all_of(colnames(list$data_features)[which(colnames(list$data_features) %in% data_meta_features$Name)]))

  # Prepare data for modeling
  data <- reshape2::melt(data,
                         id.vars = which(colnames(data) %in% var.context),
                         value.name = "Y",
                         variable.name = "id_feature")
  data <- data %>%
    dplyr::arrange(dplyr::across(c("id_feature", tidyselect::all_of(var.context))))
  data.context <- data[1:table(data$id_feature)[[1]], 1:length(var.context)]

  # Remove unwanted effects using mixed models
  id_feature <- unique(data$id_feature)
  for (i in id_feature) {
    cat(i, "\n")
    df <- data %>%
      dplyr::filter(id_feature == i)
    res <- FUNnormalization_residualMixedModels(df = df,
                                                listRandom = listRandom,
                                                listFixedToKeep = listFixedToKeep,
                                                listFixedToRemove = listFixedToRemove,
                                                HeteroSked = HeteroSked,
                                                i = i)
    data.context <- cbind(data.context, res)
  }

  # Finalize output
  colnames(data.context) <- c(var.context, as.character(id_feature))
  data.context <- tibble::as_tibble(data.context) %>%
    dplyr::arrange(IdentifierPipeline) %>%
    dplyr::select(-IdentifierPipeline)
  data_samples <- data_samples %>%
    dplyr::arrange(IdentifierPipeline) %>%
    dplyr::select(-IdentifierPipeline)

  return(list(data = data.context, data_samples = data_samples, data_meta_features = data_meta_features))
}

#' Compute Residuals Using Mixed Models
#'
#' This function calculates residuals from a mixed model to remove unwanted effects.
#'
#' @param df A data frame with feature measurements and relevant variables.
#' @param listRandom A character vector of strings specifying random effect variables.
#' @param listFixedToKeep A character vector of strings specifying fixed effect variables to keep.
#' @param listFixedToRemove A character vector of strings specifying fixed effect variables to remove.
#' @param HeteroSked A string or NULL. If not NULL, specifies the variable for heteroskedasticity correction.
#' @param i A string specifying the current feature being processed.
#'
#' @return A numeric vector of residuals for the given feature.
#'
#' @export
FUNnormalization_residualMixedModels <- function(df,
                                                 listRandom,
                                                 listFixedToKeep,
                                                 listFixedToRemove,
                                                 HeteroSked,
                                                 i) {

  listFixedInit <- c(listFixedToKeep, listFixedToRemove)
  listFixed <- if (is.null(listFixedInit)) "1" else listFixedInit
  y <- df$Y

  # Prepare data for model fitting
  df_lmer0 <- cbind.data.frame(y = y, df[, c(listRandom, listFixedInit)])
  colnames(df_lmer0)[-1] <- c(listRandom, listFixedInit)
  indNoNAs <- which(rowSums(is.na(df_lmer0)) == 0)
  df_lmer <- df_lmer0[indNoNAs, ]
  Y1 <- rep(NA, length(y))

  if (!is.null(listFixedToKeep)) {
    df_temp <- data.frame(df_lmer[, listFixedToKeep])
    colnames(df_temp) <- listFixedToKeep
    XFixedToKeep <- stats::model.matrix(~ ., data = df_temp)
  }

  if (is.null(HeteroSked)) {
    # Fit mixed model without heteroskedasticity correction
    textformulaMixed <- paste("y ~", paste0(listFixed, collapse = "+"),
                              paste0(" + (1|", paste0(listRandom, collapse = ") + (1|"), ")"))
    tochecklmer <- tryCatch(lme4::lmer(stats::as.formula(textformulaMixed), data = df_lmer,
                                       control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
                            error = function(e) NULL)

    if (!is.null(tochecklmer)) {
      if (!is.null(listFixedToKeep)) {
        Y1[indNoNAs] <- XFixedToKeep %*% matrix(lme4::fixef(tochecklmer)[names(lme4::fixef(tochecklmer)) %in% colnames(XFixedToKeep)], ncol = 1) + stats::resid(tochecklmer)
      } else {
        Y1[indNoNAs] <- mean(y) + stats::resid(tochecklmer)
      }

      # Check for collinearity (VIF)
      cat("## Collinearity (VIF) for:", i, "\n")
      vif_values <- car::vif(tochecklmer)
      print(vif_values)

      # Check the random effects variance
      cat("## Random effects variance for:", i, "\n")
      rand_var <- lme4::VarCorr(tochecklmer)
      print(rand_var)

      # Check for singularity
      cat("## Model for", i, "is singular:", lme4::isSingular(tochecklmer), "\n")

    } else {
      cat("## Model did not converge for:", i, "\n")
    }
  } else {
    # Fit mixed model with heteroskedasticity correction
    listrandomeffect <- lapply(listRandom, function(x) stats::as.formula(paste0("~1|", x)))
    forweights <- stats::as.formula(paste0("~1|", HeteroSked))
    textformulaFixedHeteroSked <- paste("y ~", paste0(listFixed, collapse = "+"))

    tochecklme <- tryCatch(nlme::lme(stats::as.formula(textformulaFixedHeteroSked), random = listrandomeffect, weights = nlme::varIdent(form = forweights),
                                     data = df_lmer,
                                     control = list(maxIter = 900, msMaxIter = 900, niterEM = 900, msMaxEval = 900, opt = "optim")),
                           error = function(e) NULL)

    if (!is.null(tochecklme)) {
      forrescale <- sd(stats::residuals(tochecklme, type = "pearson")) / sd(stats::residuals(tochecklme))
      if (!is.null(listFixedToKeep)) {
        Y1[indNoNAs] <- XFixedToKeep %*% matrix(nlme::fixed.effects(tochecklme)[names(nlme::fixed.effects(tochecklme)) %in% colnames(XFixedToKeep)], ncol = 1) +
          stats::residuals(tochecklme, type = "pearson") / forrescale
      } else {
        Y1[indNoNAs] <- mean(y) + stats::residuals(tochecklme, type = "pearson") / forrescale
      }

      # Check for collinearity (VIF)
      cat("## Collinearity (VIF) for:", i, "\n")
      vif_values <- car::vif(tochecklme)
      print(vif_values)

      # Check the random effects variance
      cat("## Random effects variance for:", i, "\n")
      rand_var <- nlme::VarCorr(tochecklme)
      print(rand_var)

      # Check for singularity (only applies to lmer models, so not applicable here)
      cat("## Singular check not applicable to nlme models for:", i, "\n")
    } else {
      cat("## Model did not converge for:", i, "\n")
    }
  }

  return(Y1)
}

#' `process_data()` function: check required columns
#'
#' This function checks for the presence of required columns
#' It verifies that necessary columns are present based on specified exclusion flags.
#'
#' @param data Data frame containing sample data, which must include `col_samples`.
#' @param data_meta_features Data frame containing feature metadata, which may need to include
#'        `col_features`, `col_missing_feature`, and `missing_pct_feature` depending on parameters.
#' @param data_meta_samples Data frame containing sample metadata, which may need to include
#'        `col_missing_sample` and `missing_pct_sample` depending on parameters.
#' @param col_samples Character vector specifying the column name(s) in `data` representing sample IDs.
#' @param col_features Character string specifying the column name in `data_meta_features` representing
#'        feature IDs (optional, only checked if provided).
#' @param exclusion_extreme_feature Logical; if `TRUE`, both `col_missing_feature` and `missing_pct_feature`
#'        must be provided and exist in `data_meta_features`.
#' @param exclusion_extreme_sample Logical; if `TRUE`, both `col_missing_sample` and `missing_pct_sample`
#'        must be provided and exist in `data_meta_samples`.
#' @param col_missing_feature Character string specifying the column name in `data_meta_features` for
#'        missing data in features (required if `exclusion_extreme_feature` is `TRUE`).
#' @param missing_pct_feature Numeric value indicating the acceptable percentage of missing data for features
#'        (required if `exclusion_extreme_feature` is `TRUE`).
#' @param col_missing_sample Character string specifying the column name in `data_meta_samples` for
#'        missing data in samples (required if `exclusion_extreme_sample` is `TRUE`).
#' @param missing_pct_sample Numeric value indicating the acceptable percentage of missing data for samples
#'        (required if `exclusion_extreme_sample` is `TRUE`).
#'
#' @return Logical; `TRUE` if all required columns are present, otherwise an error is raised.
#'
process_data_data_check <- function(data, data_meta_features, data_meta_samples,
                                    col_samples, col_features = NULL,
                                    exclusion_extreme_feature = FALSE,
                                    exclusion_extreme_sample = FALSE,
                                    col_missing_feature = NULL, missing_pct_feature = NULL,
                                    col_missing_sample = NULL, missing_pct_sample = NULL) {

  # Check if required sample columns exist in data
  required_columns_data <- c(col_samples)
  if (!all(required_columns_data %in% colnames(data))) {
    stop("Sample ID column is missing.")
  }

  # Check if the col_features column exists in data_meta_features
  if (!is.null(col_features) && !col_features %in% colnames(data_meta_features)) {
    stop("The specified 'col_features' does not exist in 'data_meta_features'.")
  }

  # Check if metadata columns exist for feature exclusions
  if (exclusion_extreme_feature) {
    if (is.null(col_missing_feature) || is.null(missing_pct_feature)) {
      stop("Both 'col_missing_feature' and 'missing_pct_feature' must be provided when 'exclusion_extreme_feature' is TRUE.")
    }
    if (!col_missing_feature %in% colnames(data_meta_features)) {
      stop("The specified 'col_missing_feature' does not exist in 'data_meta_features'.")
    }
  }

  # Check if metadata columns exist for sample exclusions
  if (exclusion_extreme_sample) {
    if (is.null(col_missing_sample) || is.null(missing_pct_sample)) {
      stop("Both 'col_missing_sample' and 'missing_pct_sample' must be provided when 'exclusion_extreme_sample' is TRUE.")
    }
    if (!col_missing_sample %in% colnames(data_meta_samples)) {
      stop("The specified 'col_missing_sample' does not exist in 'data_meta_samples'.")
    }
  }

  # If all checks pass, return TRUE to indicate success
  return(TRUE)
}

#' `process_data()` function: generate labels
#'
#' This function creates labels for saving based on several data processing options
#' (e.g., exclusion, imputation, transformation, outlier handling, etc.).
#' It returns a list of character labels that indicate the specified settings,
#' which can be used for file naming or tracking parameter choices.
#'
#' @param exclusion_extreme_feature Logical, whether extreme feature exclusion is applied.
#' @param missing_pct_feature Numeric, the percentage threshold for excluding features based on missing data (only used if `exclusion_extreme_feature` is TRUE).
#' @param exclusion_extreme_sample Logical, whether extreme sample exclusion is applied.
#' @param missing_pct_sample Numeric, the percentage threshold for excluding samples based on missing data (only used if `exclusion_extreme_sample` is TRUE).
#' @param imputation Logical, whether imputation is applied.
#' @param imputation_method Character, the method used for imputation (only used if `imputation` is TRUE).
#' @param transformation Logical, whether data transformation is applied.
#' @param transformation_method Character, the method used for transformation (only used if `transformation` is TRUE).
#' @param outlier Logical, whether outlier handling is applied.
#' @param plate_correction Logical, whether plate correction is applied.
#' @param centre_scale Logical, whether centering and scaling are applied.
#'
#' @return A list of character labels indicating the settings for each parameter.
#'
#' @export
process_data_generate_labels <- function(exclusion_extreme_feature,
                                         missing_pct_feature = NULL,
                                         exclusion_extreme_sample,
                                         missing_pct_sample = NULL,
                                         imputation,
                                         imputation_method = NULL,
                                         transformation,
                                         transformation_method = NULL,
                                         outlier,
                                         plate_correction,
                                         centre_scale) {

  # Create label for extreme feature exclusion
  LABEL_exclusion_extreme_feature <- if (exclusion_extreme_feature) {
    missing_pct_feature <- as.numeric(missing_pct_feature)
    as.character(missing_pct_feature)
  } else {
    "FALSE"
  }

  # Create label for extreme sample exclusion
  LABEL_exclusion_extreme_sample <- if (exclusion_extreme_sample) {
    missing_pct_sample <- as.numeric(missing_pct_sample)
    as.character(missing_pct_sample)
  } else {
    "FALSE"
  }

  # Create label for imputation method
  LABEL_imputation <- if (imputation) {
    as.character(imputation_method)
  } else {
    "FALSE"
  }

  # Create label for transformation method
  LABEL_transformation <- if (transformation) {
    as.character(transformation_method)
  } else {
    "FALSE"
  }

  # Create label for outlier handling
  LABEL_outlier <- if (outlier) "TRUE" else "FALSE"

  # Create label for plate correction
  LABEL_plate_correction <- if (plate_correction) "TRUE" else "FALSE"

  # Create label for centering and scaling
  LABEL_centre_scale <- if (centre_scale) "TRUE" else "FALSE"

  # Return labels as a list
  return(list(
    LABEL_exclusion_extreme_feature = LABEL_exclusion_extreme_feature,
    LABEL_exclusion_extreme_sample = LABEL_exclusion_extreme_sample,
    LABEL_imputation = LABEL_imputation,
    LABEL_transformation = LABEL_transformation,
    LABEL_outlier = LABEL_outlier,
    LABEL_plate_correction = LABEL_plate_correction,
    LABEL_centre_scale = LABEL_centre_scale
  ))
}

#' Exclude features with X misingness
#'
#' This function excludes features from a dataset (`df`) based on a specified
#' threshold of missing data percentage in a meta data file.
#'
#' @param df A data frame containing the main data with features to be filtered.
#' @param df_features A data frame containing metadata for features, including missingness percentages.
#' @param missing_pct_feature Numeric, the threshold (between 0 and 1) for feature missing data exclusion (e.g., 0.2 for 20%).
#' @param col_missing_feature Character, the column name in `df_features` indicating the percentage of missing data per feature.
#' @param col_features Character, the column name in `df_features` with feature IDs corresponding to columns in `df`.
#'
#' @return A list containing:
#' \item{df}{The filtered data frame with excluded features removed.}
#' \item{df_features}{The filtered feature metadata with excluded features removed.}
#'
#' @export
process_data_exclude_features <- function(df, df_features,
                                          missing_pct_feature = NULL,
                                          col_missing_feature,
                                          col_features) {

  # Initialize the excluded features vector
  excluded_features <- character(0)

  # Check if extreme exclusion of features is enabled
    VAR_missing_pct <- missing_pct_feature
    cat(paste0("# Exclusion features: excluding features with more than ", VAR_missing_pct * 100, "% missingness \n"))

    # Identify features to exclude based on missing data threshold
    excluded_features <- df_features %>%
      dplyr::filter(.data[[col_missing_feature]] > VAR_missing_pct) %>%
      dplyr::pull(col_features)

    # Remove identified features from main data
    df <- df %>%
      dplyr::select(-tidyselect::all_of(excluded_features))

    # Filter feature metadata to match remaining features in df
    df_features <- df_features[df_features[[col_features]] %in% colnames(df), ]

    # Report excluded features
    cat(paste0("## Exclusion features: excluded ", length(excluded_features), " feature(s) \n"))

  # Return the filtered data, feature metadata, and excluded features
  return(list(df = df, df_features = df_features, id_features_exclude = excluded_features))
}

#' Exclude samples with X missingness
#'
#' This function excludes samples from a dataset (`df`) based on a specified
#' threshold of missing data percentage in a meta data file.
#'
#' @param df A data frame containing the main data with samples to be filtered.
#' @param df_samples A data frame containing metadata for samples, including missingness percentages.
#' @param missing_pct_sample Numeric, the threshold (between 0 and 1) for sample missing data exclusion (e.g., 0.2 for 20%).
#' @param col_missing_sample Character, the column name in `df_samples` indicating the percentage of missing data per sample.
#' @param col_samples Character, the column name in `df_samples` with sample IDs corresponding to row names in `df`.
#'
#' @return A list containing:
#' \item{df}{The filtered data frame with excluded samples removed.}
#' \item{df_samples}{The filtered sample metadata with excluded samples removed.}
#' \item{excluded_samples}{A vector of the names of samples that were excluded.}
#'
#' @export
process_data_exclude_samples <- function(df, df_samples,
                                         missing_pct_sample = NULL,
                                         col_missing_sample,
                                         col_samples) {

  # Initialize the excluded samples vector
  excluded_samples <- character(0)

  # Check if extreme exclusion of samples is enabled
  VAR_missing_pct <- missing_pct_sample
  cat(paste0("# Exclusion samples: excluding samples with more than ", VAR_missing_pct * 100, "% missingness \n"))

  # Identify samples to exclude based on missing data threshold
  excluded_samples <- df_samples %>%
    dplyr::filter(.data[[col_missing_sample]] > VAR_missing_pct) %>%
    dplyr::pull(col_samples)

  # Remove identified samples from main data
  df <- df[!rownames(df) %in% excluded_samples, ]

  # Filter sample metadata to match remaining samples in df
  df_samples <- df_samples %>%
    dplyr::filter(!(!!rlang::sym(col_samples) %in% excluded_samples))

  # Report excluded samples
  cat(paste0("## Exclusion samples: excluded ", length(excluded_samples), " sample(s) \n"))

  # Return the filtered data, sample metadata, and excluded samples
  return(list(df = df, df_samples = df_samples, id_samples_exclude = excluded_samples))
}

#' Impute missing data in a data frame
#'
#' This function imputes missing data in a data frame (`df`) using the specified
#' `imputation_method`. Supported methods include "LOD", "1/5th", "KNN", "ppca", "median", "mean", "rf", and "left-censored".
#'
#' @param df A data frame with missing values to be imputed.
#' @param df_features A data frame containing feature metadata, required for "LOD" imputation.
#' @param imputation_method A character string specifying the imputation method. Valid options are:
#'   "LOD", "1/5th", "KNN", "ppca", "median", "mean", "rf", "left-censored".
#' @param col_LOD (Optional) Character name of the column in `df_features` containing the limit of detection (LOD) values, required for "LOD" imputation.
#' @param col_features (Optional) Character name of the column in `df_features` containing feature names, required for "LOD" imputation.
#'
#' @return A data frame (`df`) with imputed values.
#'
#' @export
impute_data <- function(df, df_features = NULL,
                                imputation_method = "mean",
                                col_LOD = NULL,
                                col_features = NULL) {

  # Define valid imputation methods
  valid_imputation_methods <- c("LOD", "1/5th", "KNN", "ppca", "median", "mean", "rf", "left-censored")

  # Validate imputation method
  if (is.null(imputation_method) || !(imputation_method %in% valid_imputation_methods)) {
    stop(paste("* Invalid imputation method. Choose from:", paste(valid_imputation_methods, collapse = ", ")))
  }

  # Validate LOD-specific parameters
  if (imputation_method == "LOD" && is.null(col_LOD)) {
    stop("* The 'col_LOD' parameter must be provided when using 'LOD' imputation method.")
  }

  # Perform imputation based on selected method
  switch(imputation_method,
         "LOD" = {
           cat("## Imputation using LOD \n")
           data_LOD <- stats::setNames(df_features[[col_LOD]], df_features[[col_features]])
           replace_with_lod <- function(df, lod_vector) {
             for (col in names(lod_vector)) {
               if (col %in% names(df)) {
                 df[[col]] <- ifelse(is.na(df[[col]]) | is.infinite(df[[col]]), lod_vector[[col]], df[[col]])
               }
             }
             return(df)
           }
           df <- df %>% replace_with_lod(data_LOD) %>%
             tibble::as_tibble() %>%
             tibble::add_column(id = rownames(df), .before = 1) %>%
             tibble::column_to_rownames(var = "id")
         },
         "1/5th" = {
           cat("## Imputation using 1/5th lowest detected value \n")
           df <- df %>%
             dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ ifelse(is.na(.), 1/5 * min(., na.rm = TRUE), .))) %>%
             tibble::as_tibble() %>%
             tibble::add_column(id = rownames(df), .before = 1) %>%
             tibble::column_to_rownames(var = "id")
         },
         "KNN" = {
           cat("## Imputation using KNN \n")
           df <- df %>%
             as.matrix() %>% # Convert to matrix
             t() %>% # Transpose the matrix
             impute::impute.knn(colmax = 1) %>% # Perform KNN imputation
             .$data %>% # Extract the imputed data
             t() %>%
             tibble::as_tibble() %>%
             tibble::add_column(id = rownames(df), .before = 1) %>%
             tibble::column_to_rownames(var = "id")
         },
         "ppca" = {
           cat("## Imputation using probabilistic PCA \n")
           df <- df %>%
             as.matrix() %>%
             pcaMethods::pca(nPcs = 3, method = "ppca") %>%
             pcaMethods::completeObs() %>%
             tibble::as_tibble() %>%
             tibble::add_column(id = rownames(df), .before = 1) %>%
             tibble::column_to_rownames(var = "id")
         },
         "median" = {
           cat("## Imputation using median \n")
           df <- df %>%
             as.matrix() %>%
             missMethods::impute_median(type = "columnwise") %>%
             tibble::as_tibble() %>%
             tibble::add_column(id = rownames(df), .before = 1) %>%
             tibble::column_to_rownames(var = "id")
         },
         "mean" = {
           cat("## Imputation using mean \n")
           df <- df %>%
             as.matrix() %>%
             missMethods::impute_mean(type = "columnwise") %>%
             tibble::as_tibble() %>%
             tibble::add_column(id = rownames(df), .before = 1) %>%
             tibble::column_to_rownames(var = "id")
         },
         "rf" = {
           cat("## Imputation using random forest \n")
           cl <- parallel::makeCluster(7)
           doParallel::registerDoParallel(cl)
           df <- df %>%
             as.matrix() %>%
             missForest::missForest(parallelize = 'variables', verbose = TRUE) %>%
             .$ximp %>% # Extract the imputed data
             tibble::as_tibble() %>%
             tibble::add_column(id = rownames(df), .before = 1) %>%
             tibble::column_to_rownames(var = "id")
           parallel::stopCluster(cl)
         },
         "left-censored" = {
           cat("## Imputation using left-censored \n")
           df <- df %>%
             as.matrix() %>%
             imputeLCMD::impute.MAR.MNAR(model.selector = imputeLCMD::model.Selector(.), method.MNAR = 'QRILC') %>%
             tibble::as_tibble() %>%
             tibble::add_column(id = rownames(df), .before = 1) %>%
             tibble::column_to_rownames(var = "id")
         }
  )

  return(df)
}

#' Apply Data Transformation to Data Frame
#'
#' @description This function applies various data transformations to a numeric data frame based on a specified transformation method.
#' Supported methods include "InvRank", "Log10", "Log10Capped", and "Log10ExclExtremes".
#'
#' @param df A data frame containing numeric data for transformation.
#' @param transformation_method A string specifying the transformation method to apply.
#' Supported values are "InvRank", "Log10", "Log10Capped", and "Log10ExclExtremes".
#'
#' @return A list containing:
#' \item{df}{The transformed data frame with the applied method.}
#' \item{transformation_method}{The label of the applied transformation method.}
#' \item{centre_scale}{A label indicating whether the data was centered and scaled (TRUE or NA).}
#' If conditions are not met for the specified transformation method, the original data frame will be returned in the list.
transform_data <- function(df,
                           transformation_method) {
  cat("# Transformation \n")
  valid_transformation_methods <- c("InvRank", "Log10", "Log10Capped", "Log10ExclExtremes")

  # Validate transformation method
  if (is.null(transformation_method) || !(transformation_method %in% valid_transformation_methods)) {
    stop(paste("* Invalid transformation method. Choose from:", paste(valid_transformation_methods, collapse = ", ")))
  }

  # Define transformation functions based on method
  ## InvRank
  if (transformation_method == "InvRank") {
    cat("## transformation using InvRank \n")
    fun_scale <- function(.x) {
      stats::qnorm((rank(.x, na.last = "keep") - 0.5) / sum(!is.na(.x))) # rank-inverse normalization
    }
    LABEL_transformation <- "InvRank"
    LABEL_centre_scale <- "FALSE"

    ## Log10
  } else if (transformation_method == "Log10") {
    if (any(df < 0, na.rm = TRUE)) {
      cat("* log10 not possible as negative values present. Transformation will be skipped.\n")
      return(list(df = df, transformation_method = FALSE, centre_scale = FALSE))
    } else {
      cat("## transformation using Log10 \n")
      fun_scale <- function(.x) scale(log10(.x)) # log-transform and scale
      LABEL_transformation <- "Log10"
      LABEL_centre_scale <- "TRUE"
    }

    ## Log10Capped
  } else if (transformation_method == "Log10Capped") {
    if (any(df < 0, na.rm = TRUE)) {
      cat("* Log10Capped not possible as negative values present. Transformation will be skipped.\n")
      return(list(df = df, transformation_method = FALSE, centre_scale = FALSE))
    } else {
      cat("## transformation using Log10Capped \n")
      fun_scale <- function(.x) {
        out <- scale(log10(.x))
        pmin(pmax(out, -5), 5) # cap values to [-5, 5]
      }
      LABEL_transformation <- "Log10Capped"
      LABEL_centre_scale <- "TRUE"
    }

    ## Log10ExclExtremes
  } else if (transformation_method == "Log10ExclExtremes") {
    if (any(df < 0, na.rm = TRUE)) {
      cat("* Log10ExclExtremes not possible as negative values present. Transformation will be skipped.\n")
      return(list(df = df, transformation_method = FALSE, centre_scale = FALSE))
    } else {
      cat("## transformation using Log10ExclExtremes \n")
      fun_scale <- function(.x) {
        out <- scale(log10(.x))
        ifelse(abs(out) <= 5, out, NA) # exclude extremes outside [-5, 5]
      }
      LABEL_transformation <- "Log10ExclExtremes"
      LABEL_centre_scale <- "TRUE"
    }
  }

  # Apply transformation if function defined
  if (exists("fun_scale")) {
    df <- df %>%
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), fun_scale))
  }

  return(list(df = df, transformation_method = LABEL_transformation, centre_scale = LABEL_centre_scale))
}

#' Perform PCA and Identify Outliers
#'
#' This function performs Principal Component Analysis (PCA) on the provided
#' data frame and identifies outliers based on the Local Outlier Factor (LOF).
#' It generates plots for sample outliers and excludes identified outliers from
#' the original data frame.
#'
#' @param df A data frame containing the features to analyze. Rows represent samples
#'           and columns represent features.
#' @param col_samples A string specifying the column name in `df_samples` to use for sample filtering.
#' @param col_features A string specifying the column name in `df_features` to use for feature filtering.
#' @param df_samples A data frame containing sample data for filtering.
#' @param df_features A data frame containing feature data for filtering.
#'
#' @return A list containing the filtered data frame, the sample outlier plot,
#'         and a vector of IDs of excluded samples.
#' @export
#'
#' @examples
#' \dontrun{
#'   result <- outlier_pca_lof(df, col_samples = "sample_col",
#'                               col_features = "feature_col",
#'                               df_samples = df_samples,
#'                               df_features = df_features)
#'   print(result$filtered_df)
#'   print(result$plot_samples)
#'   print(result$excluded_samples)
#' }
outlier_pca_lof <- function(df,
                            col_samples,
                            col_features,
                            df_samples,
                            df_features) {
  cat("# Outlier exclusion using PCA and LOF \n")

  # Check for missing values in the dataframe
  if (any(is.na(df))) {
    cat("* There are missing values in the data; PCA is not possible with missing data\n")
    cat("* Outlier exclusions will be skipped and the label will be changed to reflect this\n")
    return(list(df = df, plot_samples_outlier = NULL, id_samples_outlier = NULL))
  } else {
    # Step 1: PCA sample-wise
    pca_samples <- prcomp(x = df, rank. = 10)
    data_samples_analysis <- tibble::as_tibble(pca_samples$x) %>%
      tibble::add_column(id = rownames(df), .before = 1) %>%
      tibble::column_to_rownames(var = "id")

    # Step 2: Sample Outlier Detection
    llof_samples <- bigutilsr::LOF(data_samples_analysis) # Compute distances using local outlier factor
    outliers_samples <- which(llof_samples > bigutilsr::tukey_mc_up(llof_samples)) # Identify outlier threshold
    id_samples_outlier <- rownames(data_samples_analysis[outliers_samples, ])
    cat(paste0("## Excluded ", length(id_samples_outlier), " sample(s)\n"))

    # Step 3: Plot for sample outliers
    cat("## Creating sample outlier plot\n")
    data_samples_analysis$llof <- llof_samples
    plot_samples_outlier <- GGally::ggpairs(data_samples_analysis,
                                            upper = list(continuous = "blank"),
                                            diag = list(continuous = "blankDiag"),
                                            lower = list(continuous = "points"),
                                            columns = 1:10,
                                            mapping = ggplot2::aes(color = llof),
                                            legend = 11) +
      ggplot2::scale_color_viridis_c(direction = -1, option = "F") +
      ggplot2::labs(color = paste("Local outlier factor:", round(bigutilsr::tukey_mc_up(llof_samples), 2)),
                    title = paste0("Samples to exclude: ", length(outliers_samples))) +
      ggplot2::theme_minimal() +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey"),
                     legend.position = "bottom",
                     panel.spacing = grid::unit(1, "lines"))

    # Step 4: Exclude Outliers
    if (length(outliers_samples) > 0) {
      df <- df[!rownames(df) %in% id_samples_outlier, ] # Filter data
      df_samples <- df_samples %>%
        dplyr::filter(!(!!rlang::sym(col_samples) %in% id_samples_outlier)) # Filter sample data
    }

    return(list(df = df, df_samples = df_samples, plot_samples_outlier = plot_samples_outlier, id_samples_outlier = id_samples_outlier))
  }
}

#' Filter Samples Based on Case Control Criteria
#'
#' This function filters the samples in the provided data frame based on case control criteria.
#' It excludes individuals without a matched caseset and those matched casesets with
#' more than two individuals. It also provides information on excluded samples.
#'
#' @param df A data frame containing feature data where rows represent samples and columns represent features.
#' @param df_samples A data frame containing sample metadata with at least the column specified by `col_case_control`.
#' @param col_case_control A string specifying the column name in `df_samples` that indicates case control groups.
#' @param col_samples A string specifying the column name in `df_samples` to use for sample filtering (e.g., ID column).
#'
#' @return A list containing:
#'   - `filtered_df`: The filtered feature data frame.
#'   - `filtered_samples`: The filtered sample data frame.
#'   - `excluded_ids`: A vector of IDs of excluded samples.
#' @export
filter_case_control <- function(df,
                                df_samples,
                                col_case_control,
                                col_samples) {
  cat("# Filter for case control status \n")

  # Check if col_case_control is provided and exists in df_samples
  if (is.null(col_case_control)) {
    stop("* The 'col_case_control' parameter must be provided when 'case_control' is TRUE.")
  }
  if (!col_case_control %in% colnames(df_samples)) {
    stop("* The specified 'col_case_control' does not exist in 'df_samples'.")
  }


  # identify IDs with a single match
  id_exclusion_matchcaseset <- df_samples %>%
    dplyr::group_by(!!rlang::sym(col_case_control)) %>%
    dplyr::count() %>%
    dplyr::filter(n != 2) %>%
    dplyr::pull(!!rlang::sym(col_case_control))

  id_exclusion_matchcaseset <- df_samples %>%
    dplyr::filter(!!rlang::sym(col_case_control) %in% id_exclusion_matchcaseset) %>%
    dplyr::pull(col_samples)

  cat("## Excluded", length(id_exclusion_matchcaseset), "samples due to already excluded matched caseset \n")

  # filter df
  df <- df %>%
    filter(!rownames(df) %in% id_exclusion_matchcaseset)

  # filter df_samples
  df_samples <- df_samples %>%
    dplyr::filter(!(!!rlang::sym(col_samples) %in% id_exclusion_matchcaseset))

  ## Filter feature data
  df <- df[rownames(df) %in% df_samples[[col_samples]], ]

  return(list(df = df, df_samples = df_samples, id_samples_casecontrol = id_exclusion_matchcaseset))
}

#' Calculate Intraclass Correlation Coefficient (ICC)
#'
#' This function calculates the Intraclass Correlation Coefficient (ICC) for a set of samples and quality controls
#' using a linear mixed-effects model.
#' From: https://pmc.ncbi.nlm.nih.gov/articles/PMC6570933/
#' Based on: https://github.com/courtneyschiffman/Metabolomics-Filtering/blob/master/ICC.R
#'
#' @param df A data frame where rows = features and columns = samples
#' column names should be sample IDs; row names should be feature IDs
#' @param id_samples A vector of column names or indices representing the sample columns in `df`.
#' @param id_qc A vector of column names or indices representing the quality control columns in `df`.
#'
#' @return A named numeric vector of ICC values for each row in `df`.
#' @export
calculate_ICC <- function(df, id_samples, id_qc) {
  vars1 <- foreach::foreach(i = 1:nrow(df),
                            .packages = 'nlme',
                            .combine = 'rbind') %dopar% {
    reps <- factor(c(1:length(id_samples),
                     rep(length(id_samples) + 1,
                         length(id_qc))))

    # Create a data frame
    data <- data.frame(y = c(as.numeric(df[i, id_samples]),
                             as.numeric(df[i, id_qc])),
                       reps = reps)

    # Fit the linear mixed-effects model
    mm <- nlme::lme(y ~ 1, random = ~1 | reps,
                    data = data,
                    na.action = stats::na.omit)

    # Return variance components
    return(as.numeric(nlme::VarCorr(mm)[1:2]))
  }

  # Calculate ICC
  ICC <- apply(vars1, 1, function(x) x[1] / sum(x))
  names(ICC) <- rownames(df)

  return(ICC)
}
