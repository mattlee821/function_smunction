#' Process and Normalize feature data
#'
#' This function processes a given dataset with specified metadata features and samples. It supports exclusion of features with extreme missingness, various imputation methods, transformation, plate correction, and case-control data handling.
#'
#' @param data A data frame with col1 sample ID and remaining columns of features with feature IDs as colnames
#' @param data_meta_features A data frame containing metadata for the features
#' @param data_meta_samples A data frame containing metadata for the samples
#' @param col_samples A string specifying the column name in \code{data} that contains sample IDs. For example, \code{"Idepic_Bio"}.
#' @param col_features A string specifying the column name in \code{data} that contains feature IDs (should be the same values as the colnames of \code{data}. For example, \code{"UNIPROT"}.
#' @param path_out A string specifying the output directory where the processed data will be saved
#' @param path_outliers A string specifying the output directory where outlier information will be saved
#' @param exclusion_extreme A logical flag to indicate if features with extreme missingness should be excluded. Default is \code{FALSE}.
#' @param col_missing A string specifying the column name in \code{data_meta_features} that contains the percentage of missing values for each feature. For example, \code{"missing_pct"}.
#' @param missing_pct A numeric value specifying the threshold percentage for missingness above which features will be excluded. For example, \code{0.9}.
#' @param imputation A logical flag to indicate if imputation should be performed. Default is \code{FALSE}.
#' @param imputation_method A string specifying the method to use for imputation. Options are \code{"LOD"}, \code{"1/5th"}, \code{"KNN"}, \code{"ppca"}, \code{"median"}, \code{"mean"}, \code{"rf"}, and \code{"left-censored"}.
#' @param col_LOD A string specifying the column name in \code{data_meta_features} that contains the limit of detection (LOD) values, required if \code{imputation_method} is \code{"LOD"}.
#' @param transformation A logical flag to indicate if transformation should be performed. Default is \code{FALSE}.
#' @param transformation_method A string specifying the method to use for transformation. Options are \code{"InvRank"}, \code{"Log10"}, \code{"Log10Capped"}, and \code{"Log10ExclExtremes"}.
#' @param outlier A logical flag to indicate if outlier exclusion should be performed across features and samples. Default is \code{FALSE}.
#' @param plate_correction A logical flag to indicate if plate correction should be performed. Default is \code{FALSE}.
#' @param cols_listRandom A string or vector specifying columns in \code{data_meta_samples} to be treated as random effects in plate correction. For example, \code{"batch_plate"}.
#' @param cols_listFixedToKeep A vector specifying columns in \code{data_meta_samples} to be treated as fixed effects in plate correction and kept in the model. For example, \code{c("Center", "Country")}.
#' @param cols_listFixedToRemove NOT IMPLEMENTED YET: A vector specifying columns in \code{data_meta_samples} to be treated as fixed effects in plate correction and removed from the model.
#' @param col_HeteroSked A string specifying the column in \code{data_meta_samples} to be used for heteroskedasticity correction.
#' @param case_control A logical flag to indicate if the data are case-control and matched samples should be handled accordingly. Default is \code{FALSE}.
#' @param col_case_control A string specifying the column name in \code{data_meta_samples} that contains case-control matching information. For example, \code{"Match_Caseset"}.
#'
#' @return The processed and normalized data is saved to the specified output directory.
#'
#' @details
#' This function performs several data processing steps:
#' \itemize{
#'   \item Excludes features with extreme missingness based on a specified threshold.
#'   \item Imputes missing values using various methods.
#'   \item Transforms the data using specified methods.
#'   \item Corrects for plate effects using specified random and fixed effects.
#'   \item Handles case-control data to ensure matched samples are treated appropriately.
#' }
#'
#' @import dplyr
#' @import tibble
#' @import impute
#' @import pcaMethods
#' @import missMethods
#' @import missForest
#' @import imputeLCMD
#' @import bigutilsr
#' @import GGally
#' @import viridis
#' @import cowplot
#' @import parallel
#' @export
process_data <- function(
    data,
    data_meta_features,
    data_meta_samples,
    col_samples,
    col_features,
    path_out,
    path_outliers,
    exclusion_extreme = FALSE, col_missing = NULL, missing_pct = NULL,
    imputation = FALSE, imputation_method = NULL, col_LOD = NULL,
    transformation = FALSE, transformation_method = NULL,
    outlier = FALSE,
    plate_correction = FALSE, cols_listRandom = NULL, cols_listFixedToKeep = NULL, cols_listFixedToRemove = NULL, col_HeteroSked = NULL,
    case_control = FALSE, col_case_control = NULL
) {

  # Create label for saving ====
  if(exclusion_extreme == FALSE){
    LABEL_exclusion_extreme <- as.character("FALSE")
  } else{
    LABEL_exclusion_extreme <- as.character(missing_pct)
  }
  if(imputation == FALSE){
    LABEL_imputation <- as.character("FALSE")
  } else{
    LABEL_imputation <- as.character(imputation_method)
  }
  if(transformation == FALSE){
    LABEL_transformation <- as.character("FALSE")
  } else{
    LABEL_transformation <- as.character(transformation_method)
  }
  if(plate_correction == FALSE){
    LABEL_plate_correction <- as.character("FALSE")
  } else{
    LABEL_plate_correction <- as.character("TRUE")
  }
  LABEL <- paste0("exclusion-", LABEL_exclusion_extreme,
                  "_imputation-", LABEL_imputation,
                  "_transformation-", LABEL_transformation,
                  "_platecorrection-", LABEL_plate_correction)

  # Prepare data ====
  df <- data %>%
    column_to_rownames(col_samples)
  df_samples <- as_tibble(data_meta_samples)
  df_features <- as_tibble(data_meta_features)

  # Extreme exclusion ====
  if (exclusion_extreme) {
    VAR_missing_pct <- missing_pct
    cat(paste0("excluding features with more than ", VAR_missing_pct*100, "% missingness \n"))
    ## remove features with >X% missing
    id_features <- df_features %>%
      filter(.data[[col_missing]] > VAR_missing_pct) %>%
      pull(col_features)

    df <- df %>%
      select(-all_of(id_features))
    ## filter feature data
    df_features <- df_features[df_features[[col_features]] %in% colnames(df), ]
  }

  # imputation: replace missing values with X ====
  if (imputation) {
    ## LOD ====
    if (imputation_method == "LOD") {
      cat("imputing using LOD \n")
      if(imputation_method == "LOD"){
        data_LOD <- setNames(df_features[[col_LOD]], df_features[[col_features]])
      }
      replace_with_lod <- function(df, lod_vector) {
        for (col in names(lod_vector)) {
          if (col %in% names(df)) {
            df[[col]] <- ifelse(is.na(df[[col]]) | is.infinite(df[[col]]), lod_vector[[col]], df[[col]])
          }
        }
        return(df)
      }
      df <- df %>%
        replace_with_lod(data_LOD) %>%
        as_tibble() %>%
        add_column(id = rownames(df), .before = 1) %>%
        column_to_rownames(var = "id")
    }

    ## 1/5th lowest detected value for each feature ====
    if (imputation_method == "1/5th") {
      cat("imputing using 1/5th lowest detected value \n")
      df <- df %>%
        mutate(across(where(is.numeric), ~ ifelse(is.na(.), 1/5 * min(., na.rm = TRUE), .))) %>%
        as_tibble() %>%
        add_column(id = rownames(df), .before = 1) %>%
        column_to_rownames(var = "id")
    }

    ## KNN ====
    if (imputation_method == "KNN") {
      cat("imputing using KNN \n")
      df <- df %>%
        as.matrix() %>% # Convert to matrix
        t() %>% # Transpose the matrix
        impute::impute.knn(colmax = 1) %>% # Perform KNN imputation
        .$data %>% # Extract the imputed data
        t() %>%
        as_tibble() %>%
        add_column(id = rownames(df), .before = 1) %>%
        column_to_rownames(var = "id")
    }

    ## probabilistic pca ====
    if (imputation_method == "ppca") {
      cat("imputing using probabilistic PCA \n")
      df <- df %>%
        as.matrix() %>%
        pcaMethods::pca(nPcs = 3, method = "ppca") %>%
        pcaMethods::completeObs() %>%
        as_tibble() %>%
        add_column(id = rownames(df), .before = 1) %>%
        column_to_rownames(var = "id")
    }

    ## median ====
    if (imputation_method == "median") {
      cat("imputing using median \n")
      df <- df %>%
        as.matrix() %>%
        missMethods::impute_median(type = "columnwise") %>%
        as_tibble() %>%
        add_column(id = rownames(df), .before = 1) %>%
        column_to_rownames(var = "id")
    }

    ## mean ====
    if (imputation_method == "mean") {
      cat("imputing using mean \n")
      df <- df %>%
        as.matrix() %>%
        missMethods::impute_mean(type = "columnwise") %>%
        as_tibble() %>%
        add_column(id = rownames(df), .before = 1) %>%
        column_to_rownames(var = "id")
    }

    ## random forest ====
    if (imputation_method == "rf") {
      cat("imputing using random forest \n")
      cl <- parallel::makeCluster(7)
      doParallel::registerDoParallel(cl)
      df <- df %>%
        as.matrix() %>%
        missForest::missForest(parallelize = 'variables', verbose = TRUE) %>%
        .$ximp %>% # Extract the imputed data
        as_tibble() %>%
        add_column(id = rownames(df), .before = 1) %>%
        column_to_rownames(var = "id")
      parallel::stopCluster(cl)
    }

    ## left-censored ====
    if (imputation_method == "left-censored") {
      cat("imputing using left-censored \n")
      df <- df %>%
        as.matrix() %>%
        imputeLCMD::impute.MAR.MNAR(model.selector = imputeLCMD::model.Selector(.), method.MNAR = 'QRILC') %>%
        as_tibble() %>%
        add_column(id = rownames(df), .before = 1) %>%
        column_to_rownames(var = "id")
    }
  }

  # transformation ====
  if (transformation) {
    # transformation method ====
    if (transformation_method == "InvRank") {
      fun_scale <- function(.x) {
        out <- qnorm((rank(.x, na.last = "keep") - 0.5) / sum(!is.na(.x))) # rank-inverse normalization
      }
    }
    if (transformation_method == "Log10") {
      fun_scale <- function(.x) { # log-transform + center/scale vector
        out <- log10(.x) # log-transform
        out <- out - mean(out) # center
        out <- out / sd(out) # scale to unit-variance
      }
    }
    if (transformation_method == "Log10Capped") {
      fun_scale <- function(.x) { # log-transform + center/scale vector
        out <- log10(.x) # log-transform
        out <- out - mean(out) # center
        out <- out / sd(out) # scale to unit-variance
        out <- pmin(pmax(out, -5), 5) # now we replace any values out of the [-5, 5] range by +/- 5.
      }
    }
    if (transformation_method == "Log10ExclExtremes") {
      fun_scale <- function(.x) { # log-transform + center/scale vector
        out <- log10(.x) # log-transform
        out <- out - mean(out) # center
        out <- out / sd(out) # scale to unit-variance
        out <- ifelse(abs(out) <= 5, out, NA) # now we replace any values out of the [-5, 5] range by +/- 5.
      }
    }

    # transform ====
    df <- df %>%
      mutate(across(everything(), fun_scale))
  }
  # outlier exclusion ====
  if(outlier){
    # Outlier identification
    ## based on: https://privefl.github.io/blog/detecting-outlier-samples-in-pca/
    # step 1: PCA sample wise and feature wise ====
    pca_samples <- prcomp(x = df, rank. = 10)
    data_samples_analysis <- pca_samples$x %>%
      as_tibble() %>%
      add_column(id = rownames(df), .before = 1) %>%
      column_to_rownames(var = "id")

    data_features <- t(df)
    pca_features <- prcomp(x = data_features, rank. = 10)
    data_features_analysis <- as.data.frame(pca_features$x)

    # step 2: sample outliers ====
    llof <- LOF(data_samples_analysis) # compute distances using local outlier factor
    outliers_samples <- which(llof > tukey_mc_up(llof)) # identify outlier threshold using tukeys rule
    id_samples_exclude <- rownames(data_samples_analysis[outliers_samples, ])
    # plot PCs with lof highlighting
    cat("making sample outlier plot \n")
    data_samples_analysis$llof <- llof
    plot_samples <- GGally::ggpairs(data_samples_analysis,
                                    upper = list(continuous = "blank"),
                                    diag = list(continuous = "blankDiag"),
                                    lower = list(continuous = "points"),
                                    columns = 1:10,
                                    mapping = aes(color = llof),
                                    legend = 11) +
      scale_color_viridis_c(direction = -1, option = "F") +
      labs(color = paste("Local outlier factor:", round(tukey_mc_up(llof), 2)), title = paste0("Samples to exclude: ",  length(outliers_samples))) +
      theme_minimal() +
      theme(strip.background = element_rect(fill = "grey")) +
      theme(legend.position = "bottom") +
      theme(panel.spacing = unit(1, "lines"))


    # step 3: feature outliers ====
    llof <- LOF(data_features_analysis) # compute distances using local outlier factor
    outliers_features <- which(llof > tukey_mc_up(llof)) # identify outlier threshold using tukeys rule
    id_features_exclude <- rownames(data_features_analysis[outliers_features, ])
    # plot PCs with lof highlighting
    cat("making feature outlier plot \n")
    data_features_analysis$llof <- llof
    plot_features <- GGally::ggpairs(data_features_analysis,
                                     upper = list(continuous = "blank"),
                                     diag = list(continuous = "blankDiag"),
                                     lower = list(continuous = "points"),
                                     columns = 1:10,
                                     mapping = aes(color = llof),
                                     legend = 11) +
      scale_color_viridis_c(direction = -1, option = "F") +
      labs(color = paste("Local outlier factor:", round(tukey_mc_up(llof), 2)), title = paste0("Features to exclude: ",  length(outliers_features))) +
      theme_minimal() +
      theme(strip.background = element_rect(fill = "grey")) +
      theme(legend.position = "bottom") +
      theme(panel.spacing = unit(1, "lines"))

    # step 4: exclusion ====
    ## save figure
    cat("saving outlier plots \n")
    tiff(paste0(path_outliers, "plot-outliers_", LABEL, ".tiff"), width = 1000, height = 1000, units = "px")
    print(cowplot::plot_grid(
      GGally::ggmatrix_gtable(plot_samples),
      GGally::ggmatrix_gtable(plot_features),
      nrow = 2))
    dev.off()
    ## make exclusions
    if(length(outliers_samples) > 0) {df <- df[-outliers_samples, ]}
    if(length(outliers_features) > 0) {df <- df[, -outliers_features]}
    ## save excluded IDs
    cat("saving outlier info \n")
    id_outliers <- list(features_extreme = id_features, features_lof = id_features_exclude, samples_lof = id_samples_exclude)
    saveRDS(object = id_outliers, file = paste0(path_outliers, "id-outliers_", LABEL, ".rds"))
  }

  # case-control data ====
  if(case_control){
    ## re-filter feature data based on matched cases
    ## filter meta data for included samples
    col_case_control <- sym(col_case_control)
    df_samples <- df_samples %>%
      filter(.data[[col_samples]] %in% rownames(df))

    ## exclude individuals without a matched caseset
    df_samples <- df_samples %>%
      group_by(!!col_case_control) %>%
      filter(n() >= 2) %>%
      ungroup()

    ## exclude matched casesets with more than two people
    df_samples <- df_samples %>%
      group_by(!!col_case_control) %>%
      filter(n() <= 2) %>%
      ungroup()

    ## filter feature data
    df <- df[rownames(df) %in% df_samples[[col_samples]], ]
  }

  # normalisation and plate correction ====
  if(plate_correction){
    ## convert rownames to column1
    df <- df %>%
      rownames_to_column(var = col_samples)
    ## make df
    df <- list(data_features = df,
               data_samples = df_samples,
               data_meta_features = df_features %>%
                 rename(Name = all_of(col_features)))
    ## transformation
    df <- normalization_residualMixedModels(df,
                                            forIdentifier = col_samples, # sample ID
                                            listRandom = cols_listRandom, # variables to model as random effects; effects will be removed
                                            listFixedToKeep = cols_listFixedToKeep, # variables to model as fixed effects; effects will be kept
                                            listFixedToRemove = cols_listFixedToRemove, # variables to model as fixed effects; effects will be removed
                                            HeteroSked = col_HeteroSked # variable for which heteroscedasticity will be accounted for
    )
    df <- df$data %>%
      select(-all_of(c(cols_listRandom, cols_listFixedToKeep)))
    saveRDS(object = df, file = paste0(path_out, "data-features_", LABEL, ".rds"))
  } else {
    saveRDS(object = df, file = paste0(path_out, "data-features_", LABEL, ".rds"))
  }
  return(df)
}

#' Normalize feature Data Using Residual Mixed Models
#'
#' This function normalizes feature data by removing unwanted effects using mixed models.
#' It accounts for random effects and fixed effects specified by the user, and optionally
#' corrects for heteroskedasticity.
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
#' @import tidyverse reshape2 hglm nlme lme4
#' @export
normalization_residualMixedModels <- function(list,
                                              forIdentifier = c("Study", "Batch", "CaseSet", "Idepic"),
                                              listRandom = NULL,
                                              listFixedToKeep = NULL,
                                              listFixedToRemove = NULL,
                                              HeteroSked = NULL) {

  # Convert list of data frames to tibble and unite identifiers ====
  data_features <- as_tibble(list$data_features) %>%
    tidyr::unite(IdentifierPipeline, all_of(forIdentifier), sep = "", remove = FALSE) %>%
    arrange(IdentifierPipeline)
  data_samples <- as_tibble(list$data_samples) %>%
    tidyr::unite(IdentifierPipeline, all_of(forIdentifier), sep = "", remove = FALSE) %>%
    arrange(IdentifierPipeline)
  data_meta_features <- as_tibble(list$data_meta_features)

  # Check for consistency between data_features and data_samples ====
  if (sum(data_features$IdentifierPipeline == data_samples$IdentifierPipeline) < nrow(data_features)) {
    stop("Idepic should be the same in data_features and data_samples")
  }

  # Merge data_features with data_samples and select relevant columns ====
  var.context <- unique(c("IdentifierPipeline", forIdentifier, listRandom, listFixedToKeep, listFixedToRemove))
  data <- left_join(data_features,
                    data_samples %>%
                      dplyr::select(all_of(var.context)) %>%
                      dplyr::select(-any_of(forIdentifier)),
                    by = "IdentifierPipeline") %>%
    dplyr::select(all_of(var.context), all_of(colnames(list$data_features)[which(colnames(list$data_features) %in% data_meta_features$Name)]))

  # Prepare data for modeling
  data <- reshape2::melt(data,
                         id.vars = which(colnames(data) %in% var.context),
                         value.name = "Y",
                         variable.name = "id_feature")
  data <- data %>%
    arrange(across(c("id_feature", all_of(var.context))))
  data.context <- data[1:table(data$id_feature)[[1]], 1:length(var.context)]

  # Remove unwanted effects using mixed models
  id_feature <- unique(data$id_feature)
  for (i in id_feature) {
    cat(i, "\n")
    df <- data %>%
      filter(id_feature == i)
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
  data.context <- as_tibble(data.context) %>%
    arrange(IdentifierPipeline) %>%
    dplyr::select(-IdentifierPipeline)
  data_samples <- data_samples %>%
    arrange(IdentifierPipeline) %>%
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
#' @import lme4 nlme
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
    tochecklmer <- tryCatch(lmer(as.formula(textformulaMixed), data = df_lmer,
                                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
                            error = function(e) NULL)
    if (!is.null(tochecklmer)) {
      if (!is.null(listFixedToKeep)) {
        Y1[indNoNAs] <- XFixedToKeep %*% matrix(fixef(tochecklmer)[names(fixef(tochecklmer)) %in% colnames(XFixedToKeep)], ncol = 1) + stats::resid(tochecklmer)
      } else {
        Y1[indNoNAs] <- mean(y) + stats::resid(tochecklmer)
      }
    }
  } else {
    # Fit mixed model with heteroskedasticity correction
    listrandomeffect <- lapply(listRandom, function(x) as.formula(paste0("~1|", x)))
    forweights <- as.formula(paste0("~1|", HeteroSked))
    textformulaFixedHeteroSked <- paste("y ~", paste0(listFixed, collapse = "+"))

    tochecklme <- tryCatch(lme(as.formula(textformulaFixedHeteroSked), random = listrandomeffect, weights = varIdent(form = forweights),
                               data = df_lmer,
                               control = list(maxIter = 900, msMaxIter = 900, niterEM = 900, msMaxEval = 900, opt = "optim")),
                           error = function(e) NULL)

    if (!is.null(tochecklme)) {
      forrescale <- sd(stats::residuals(tochecklme, type = "pearson")) / sd(stats::residuals(tochecklme))
      if (!is.null(listFixedToKeep)) {
        Y1[indNoNAs] <- XFixedToKeep %*% matrix(fixed.effects(tochecklme)[names(fixed.effects(tochecklme)) %in% colnames(XFixedToKeep)], ncol = 1) +
          stats::residuals(tochecklme, type = "pearson") / forrescale
      } else {
        Y1[indNoNAs] <- mean(y) + stats::residuals(tochecklme, type = "pearson") / forrescale
      }
    }
  }

  return(Y1)
}
