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
  required_columns_data <- c(col_samples)
  if (!all(required_columns_data %in% colnames(data))) {
    stop("sample ID col is missing.")
  }

  ## Check if metadata columns exist
  if (exclusion_extreme_feature) {
    if (is.null(col_missing_feature) || is.null(missing_pct_feature)) {
      stop("Both 'col_missing_feature' and 'missing_pct_feature' must be provided when 'exclusion_extreme_feature' is TRUE.")
    }
    if (!col_missing_feature %in% colnames(data_meta_features)) {
      stop("The specified 'col_missing_feature' does not exist in 'data_meta_features'.")
    }
  }

  if (exclusion_extreme_sample) {
    if (is.null(col_missing_sample) || is.null(missing_pct_sample)) {
      stop("Both 'col_missing_sample' and 'missing_pct_sample' must be provided when 'exclusion_extreme_sample' is TRUE.")
    }
    if (!col_missing_sample %in% colnames(data_meta_samples)) {
      stop("The specified 'col_missing_sample' does not exist in 'data_meta_samples'.")
    }
  }

  # Create label for saving ====
  if(exclusion_extreme_feature == FALSE){
    LABEL_exclusion_extreme_feature <- as.character("FALSE")
  } else{
    LABEL_exclusion_extreme_feature <- as.character(missing_pct_feature)
    missing_pct_feature <- as.numeric(missing_pct_feature)
  }
  if(exclusion_extreme_sample == FALSE){
    LABEL_exclusion_extreme_sample <- as.character("FALSE")
  } else{
    LABEL_exclusion_extreme_sample <- as.character(missing_pct_sample)
    missing_pct_sample <- as.numeric(missing_pct_sample)
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
  if(outlier == FALSE){
    LABEL_outlier <- as.character("FALSE")
  } else{
    LABEL_outlier <- as.character("TRUE")
  }
  if(plate_correction == FALSE){
    LABEL_plate_correction <- as.character("FALSE")
  } else{
    LABEL_plate_correction <- as.character("TRUE")
  }
  if(centre_scale == FALSE){
    LABEL_centre_scale <- as.character("FALSE")
  } else{
    LABEL_centre_scale <- as.character("TRUE")
  }

  # Prepare data ====
  df <- data %>%
    tibble::column_to_rownames(col_samples)
  df_samples <- tibble::as_tibble(data_meta_samples)
  df_features <- tibble::as_tibble(data_meta_features)

  # Extreme exclusion features ====
  if (exclusion_extreme_feature) {
    VAR_missing_pct <- missing_pct_feature
    cat(paste0("# Exclusion features: excluding features with more than ", VAR_missing_pct*100, "% missingness \n"))
    ## remove features with >X% missing
    id_exclusion_extreme_feature <- df_features %>%
      dplyr::filter(.data[[col_missing_feature]] > VAR_missing_pct) %>%
      dplyr::pull(col_features)
    # filter data
    df <- df %>%
      dplyr::select(-tidyselect::all_of(id_exclusion_extreme_feature))
    ## filter feature data
    df_features <- df_features[df_features[[col_features]] %in% colnames(df), ]
    cat(paste0("## Exclusion features: excluded ", length(id_exclusion_extreme_feature), " feature(s) \n"))
  }

  # Extreme exclusion samples ====
  if (exclusion_extreme_sample) {
    VAR_missing_pct <- missing_pct_sample
    cat(paste0("# Exclusion samples: excluding samples with more than ", VAR_missing_pct*100, "% missingness \n"))
    ## remove samples with >X% missing
    id_exclusion_extreme_sample <- df_samples %>%
      dplyr::filter(.data[[col_missing_sample]] > VAR_missing_pct) %>%
      dplyr::pull(col_samples)
    # filter data
    df <- df[!rownames(df) %in% id_exclusion_extreme_sample, ]
    ## filter sample data
    df_samples <- df_samples %>%
      filter(!(!!sym(col_samples) %in% id_exclusion_extreme_sample))
    cat(paste0("## Exclusion samples: excluded ", length(id_exclusion_extreme_sample), " sample(s) \n"))
  }

  # imputation ====
  if (imputation) {
    cat("# Imputation \n ")
    valid_imputation_methods <- c("LOD", "1/5th", "KNN", "ppca", "median", "mean", "rf", "left-censored")
    if (is.null(imputation_method) || !(imputation_method %in% valid_imputation_methods)) {
      stop(paste("* Invalid imputation method. Choose from:", paste(valid_imputation_methods, collapse = ", ")))
    }
    if (imputation_method == "LOD" && is.null(col_LOD)) {
      stop("* The 'col_LOD' parameter must be provided when using 'LOD' imputation method.")
    }

    ## LOD ====
    if (imputation_method == "LOD") {
      cat("## imputation using LOD \n")
      if(imputation_method == "LOD"){
        data_LOD <- stats::setNames(df_features[[col_LOD]], df_features[[col_features]])
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
        tibble::as_tibble() %>%
        tibble::add_column(id = rownames(df), .before = 1) %>%
        tibble::column_to_rownames(var = "id")
    }

    ## 1/5th lowest detected value for each feature ====
    if (imputation_method == "1/5th") {
      cat("## imputation using 1/5th lowest detected value \n")
      df <- df %>%
        dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ ifelse(is.na(.), 1/5 * min(., na.rm = TRUE), .))) %>%
        tibble::as_tibble() %>%
        tibble::add_column(id = rownames(df), .before = 1) %>%
        tibble::column_to_rownames(var = "id")
    }

    ## KNN ====
    if (imputation_method == "KNN") {
      cat("## imputation using KNN \n")
      df <- df %>%
        as.matrix() %>% # Convert to matrix
        t() %>% # Transpose the matrix
        impute::impute.knn(colmax = 1) %>% # Perform KNN imputation
        .$data %>% # Extract the imputed data
        t() %>%
        tibble::as_tibble() %>%
        tibble::add_column(id = rownames(df), .before = 1) %>%
        tibble::column_to_rownames(var = "id")
    }

    ## probabilistic pca ====
    if (imputation_method == "ppca") {
      cat("## imputation using probabilistic PCA \n")
      df <- df %>%
        as.matrix() %>%
        pcaMethods::pca(nPcs = 3, method = "ppca") %>%
        pcaMethods::completeObs() %>%
        tibble::as_tibble() %>%
        tibble::add_column(id = rownames(df), .before = 1) %>%
        tibble::column_to_rownames(var = "id")
    }

    ## median ====
    if (imputation_method == "median") {
      cat("## imputation using median \n")
      df <- df %>%
        as.matrix() %>%
        missMethods::impute_median(type = "columnwise") %>%
        tibble::as_tibble() %>%
        tibble::add_column(id = rownames(df), .before = 1) %>%
        tibble::column_to_rownames(var = "id")
    }

    ## mean ====
    if (imputation_method == "mean") {
      cat("## imputation using mean \n")
      df <- df %>%
        as.matrix() %>%
        missMethods::impute_mean(type = "columnwise") %>%
        tibble::as_tibble() %>%
        tibble::add_column(id = rownames(df), .before = 1) %>%
        tibble::column_to_rownames(var = "id")
    }

    ## random forest ====
    if (imputation_method == "rf") {
      cat("## imputation using random forest \n")
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
    }

    ## left-censored ====
    if (imputation_method == "left-censored") {
      cat("## imputation using left-censored \n")
      df <- df %>%
        as.matrix() %>%
        imputeLCMD::impute.MAR.MNAR(model.selector = imputeLCMD::model.Selector(.), method.MNAR = 'QRILC') %>%
        tibble::as_tibble() %>%
        tibble::add_column(id = rownames(df), .before = 1) %>%
        tibble::column_to_rownames(var = "id")
    }
  }

  # transformation ====
  if (transformation) {
    cat("# Transformation \n")
    valid_transformation_methods <- c("InvRank", "Log10", "Log10Capped", "Log10ExclExtremes")

    if (is.null(transformation_method) || !(transformation_method %in% valid_transformation_methods)) {
      stop(paste("* Invalid transformation method. Choose from:", paste(valid_transformation_methods, collapse = ", ")))
    }

    ## InvRank ====
    if (transformation_method == "InvRank") {
      cat("## transformation using InvRank \n")
      fun_scale <- function(.x) {
        out <- stats::qnorm((rank(.x, na.last = "keep") - 0.5) / sum(!is.na(.x))) # rank-inverse normalization
      }
    }

    ## Log10 ====
    if (transformation_method == "Log10") {
      if (any(df < 0, na.rm = TRUE)) {
        cat("* log10 not possible as negative values present. Transformation will be skipped and label set to FALSE\n")
        LABEL_transformation <- "FALSE"
      } else {
        cat("## transformation using Log10 \n")
        fun_scale <- function(.x) {
          out <- log10(.x) # log-transform
          out <- scale(out) # center and scale to unit-variance
        }
        LABEL_transformation <- "Log10"
        LABEL_centre_scale <- "TRUE"
      }
    }

    ## Log10Capped ====
    if (transformation_method == "Log10Capped") {
      if (any(df < 0, na.rm = TRUE)) {
        cat("* Log10Capped not possible as negative values present. Transformation will be skipped and label set to FALSE\n")
        LABEL_transformation <- "FALSE"
      } else {
        cat("## transformation using Log10Capped \n")
        fun_scale <- function(.x) {
          out <- log10(.x) # log-transform
          out <- scale(out) # center and scale to unit-variance
          out <- pmin(pmax(out, -5), 5) # cap values to [-5, 5]
        }
        LABEL_transformation <- "Log10Capped"
        LABEL_centre_scale <- "TRUE"
      }
    }

    ## Log10ExclExtremes ====
    if (transformation_method == "Log10ExclExtremes") {
      if (any(df < 0, na.rm = TRUE)) {
        cat("* Log10ExclExtremes not possible as negative values present. Transformation will be skipped and label set to FALSE\n")
        LABEL_transformation <- "FALSE"
      } else {
        cat("## transformation using Log10ExclExtremes \n")
        fun_scale <- function(.x) {
          out <- log10(.x) # log-transform
          out <- scale(out) # center and scale to unit-variance
          out <- ifelse(abs(out) <= 5, out, NA) # exclude extremes outside [-5, 5]
        }
        LABEL_transformation <- "Log10ExclExtremes"
        LABEL_centre_scale <- "TRUE"

      }
    }

    ## transform ====
    if (exists("fun_scale")) {
      df <- df %>%
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric), fun_scale))
    }
  }

  # outlier exclusion ====
  if(outlier){
    cat(paste0("# Outlier exclusion \n"))

    if (any(is.na(df))) {
      cat("* There are missing values in the data; PCA is not possible with missing data\n")
      cat("* Outlier exclusions will be skipped and the label will be changed to reflect this")
      LABEL_outlier <- FALSE

    } else {

      # Outlier identification
      ## based on: https://privefl.github.io/blog/detecting-outlier-samples-in-pca/
      # step 1: PCA sample wise and feature wise ====
      pca_samples <- prcomp(x = df, rank. = 10)
      data_samples_analysis <- pca_samples$x %>%
        tibble::as_tibble() %>%
        tibble::add_column(id = rownames(df), .before = 1) %>%
        tibble::column_to_rownames(var = "id")

      data_features <- t(df)
      pca_features <- stats::prcomp(x = data_features, rank. = 10)
      data_features_analysis <- as.data.frame(pca_features$x)

      # step 2: sample outliers ====
      llof <- bigutilsr::LOF(data_samples_analysis) # compute distances using local outlier factor
      outliers_samples <- which(llof > bigutilsr::tukey_mc_up(llof)) # identify outlier threshold using tukeys rule
      id_exclusion_outlier_sample <- rownames(data_samples_analysis[outliers_samples, ])
      cat(paste0("## excluded ", length(id_exclusion_outlier_sample), " sample(s) \n"))
      # plot PCs with lof highlighting
      cat("## making sample outlier plot \n")
      data_samples_analysis$llof <- llof
      plot_samples <- GGally::ggpairs(data_samples_analysis,
                                      upper = list(continuous = "blank"),
                                      diag = list(continuous = "blankDiag"),
                                      lower = list(continuous = "points"),
                                      columns = 1:10,
                                      mapping = ggplot2::aes(color = llof),
                                      legend = 11) +
        ggplot2::scale_color_viridis_c(direction = -1, option = "F") +
        ggplot2::labs(color = paste("Local outlier factor:", round(bigutilsr::tukey_mc_up(llof), 2)), title = paste0("Samples to exclude: ",  length(outliers_samples))) +
        ggplot2::theme_minimal() +
        ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey")) +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::theme(panel.spacing = grid::unit(1, "lines"))


      # step 3: feature outliers ====
      llof <- bigutilsr::LOF(data_features_analysis) # compute distances using local outlier factor
      outliers_features <- which(llof > bigutilsr::tukey_mc_up(llof)) # identify outlier threshold using tukeys rule
      id_exclusion_outlier_feature <- rownames(data_features_analysis[outliers_features, ])
      cat(paste0("## excluded ", length(id_exclusion_outlier_feature), " feature(s) \n"))
      # plot PCs with lof highlighting
      cat("## making feature outlier plot \n")
      data_features_analysis$llof <- llof
      plot_features <- GGally::ggpairs(data_features_analysis,
                                       upper = list(continuous = "blank"),
                                       diag = list(continuous = "blankDiag"),
                                       lower = list(continuous = "points"),
                                       columns = 1:10,
                                       mapping = ggplot2::aes(color = llof),
                                       legend = 11) +
        ggplot2::scale_color_viridis_c(direction = -1, option = "F") +
        ggplot2::labs(color = paste("Local outlier factor:", round(bigutilsr::tukey_mc_up(llof), 2)), title = paste0("Features to exclude: ",  length(outliers_features))) +
        ggplot2::theme_minimal() +
        ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey")) +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::theme(panel.spacing = grid::unit(1, "lines"))

      # step 4: exclusion ====
      if(length(outliers_samples) > 0) {
        # filter data
        df <- df[!rownames(df) %in% id_exclusion_outlier_sample, ]
        ## filter sample data
        df_samples <- df_samples %>%
          filter(!(!!sym(col_samples) %in% id_exclusion_outlier_sample))
      }

      if(length(outliers_features) > 0) {
        # filter data
        df <- df %>%
          dplyr::select(-tidyselect::all_of(id_exclusion_outlier_feature))
        ## filter feature data
        df_features <- df_features[df_features[[col_features]] %in% colnames(df), ]
      }

    }
  }

  # case-control data ====
  if (case_control) {
    cat("# Case control \n")
    if (is.null(col_case_control)) {
      stop("* The 'col_case_control' parameter must be provided when 'case_control' is TRUE.")
    }
    if (!col_case_control %in% colnames(data_meta_samples)) {
      stop("* The specified 'col_case_control' does not exist in 'data_meta_samples'.")
    }    ## re-filter feature data based on matched cases
    id_exclusion_matchcaseset <- df_samples %>%
      dplyr::group_by(!!sym(col_case_control)) %>%  # Group by the Match_Caseset
      dplyr::summarize(count = n(), .groups = 'drop') %>%
      dplyr::filter(count < 2 | count > 2) %>%  # Keep those with less than or more than 2
      dplyr::inner_join(df_samples, by = col_case_control) %>%  # Join back to df_samples
      dplyr::pull(!!sym(col_samples))  # Pull the Idepic_Bio values to exclude

    cat("## excluded", length(id_exclusion_matchcaseset), "samples due to already excluded matched caseset \n")
    # col_case_control <- rlang::sym(col_case_control)

    ## exclude individuals without a matched caseset
    df_samples <- df_samples %>%
      dplyr::group_by(!!sym(col_case_control)) %>%
      dplyr::filter(n() >= 2) %>%
      dplyr::ungroup()

    ## exclude matched casesets with more than two people
    df_samples <- df_samples %>%
      dplyr::group_by(!!sym(col_case_control)) %>%
      dplyr::filter(n() <= 2) %>%
      dplyr::ungroup()

    ## filter feature data
    df <- df[rownames(df) %in% df_samples[[col_samples]], ]
  }

  # plate correction ====
  if (plate_correction) {
    cat("# Plate correction \n")
    # Check required columns for plate correction
    if (is.null(cols_listRandom) || is.null(cols_listFixedToKeep)) {
      stop("* Both 'cols_listRandom' and 'cols_listFixedToKeep' must be provided for plate correction.")
    }
    ## convert rownames to column1
    df <- df %>%
      tibble::rownames_to_column(var = col_samples)
    ## make df
    df <- list(data_features = df,
               data_samples = df_samples,
               data_meta_features = df_features %>%
                 dplyr::rename(Name = tidyselect::all_of(col_features)))
    ## transformation
    df <- normalization_residualMixedModels(df,
                                            forIdentifier = col_samples, # sample ID
                                            listRandom = cols_listRandom, # variables to model as random effects; effects will be removed
                                            listFixedToKeep = cols_listFixedToKeep, # variables to model as fixed effects; effects will be kept
                                            listFixedToRemove = cols_listFixedToRemove, # variables to model as fixed effects; effects will be removed
                                            HeteroSked = col_HeteroSked # variable for which heteroscedasticity will be accounted for
    )
    df <- df$data %>%
      dplyr::select(-tidyselect::all_of(c(cols_listRandom, cols_listFixedToKeep)))
  }

  # centre and scale ====
  if (centre_scale) {
    cat("# Centre and scale \n")

    # Check if LABEL_transformation indicates transformation was applied
    if (LABEL_transformation == "TRUE" && grepl("log10", transformation_method, ignore.case = TRUE)) {
      cat("* Log10 transformation was applied and included centering and scaling. Skipping this centering and scaling step.\n")
    } else {
      df <- df %>%
        mutate(across(where(is.numeric), ~ scale(.)))
      cat("## data centred and scaled \n")

    }
  }

  # LABEL ====
  LABEL <- paste0("exclusion-feature-", LABEL_exclusion_extreme_feature,
                  "_exclusion-sample-", LABEL_exclusion_extreme_sample,
                  "_imputation-", LABEL_imputation,
                  "_transformation-", LABEL_transformation,
                  "_outlier-", LABEL_outlier,
                  "_platecorrection-", LABEL_plate_correction,
                  "_centre-scale-", LABEL_centre_scale)

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
    id_outliers <- list(
      if (exists("id_exclusion_extreme_feature")) id_exclusion_extreme_feature else NULL,
      if (exists("id_exclusion_extreme_sample")) id_exclusion_extreme_sample else NULL,
      if (exists("id_exclusion_outlier_feature")) id_exclusion_outlier_feature else NULL,
      if (exists("id_exclusion_outlier_sample")) id_exclusion_outlier_sample else NULL,
      if (exists("id_exclusion_matchcaseset")) id_exclusion_matchcaseset else NULL
    )
    names(id_outliers) <- c("id_exclusion_extreme_feature", "id_exclusion_extreme_sample", "id_exclusion_outlier_feature", "id_exclusion_outlier_sample", "id_exclusion_matchcaseset")

    ### check and print exclusions
    total_excluded_features <- 0
    if (exists("id_exclusion_extreme_feature")) {
      total_excluded_features <- total_excluded_features + length(id_exclusion_extreme_feature)
    }
    if (exists("id_exclusion_outlier_feature")) {
      total_excluded_features <- total_excluded_features + length(id_exclusion_outlier_feature)
    }
    cat("### total feature(s) excluded =", total_excluded_features, "\n")


    total_excluded_samples <- 0
    if (exists("id_exclusion_extreme_sample")) {
      total_excluded_samples <- total_excluded_samples + length(id_exclusion_extreme_sample)
    }
    if (exists("id_exclusion_outlier_sample")) {
      total_excluded_samples <- total_excluded_samples + length(id_exclusion_outlier_sample)
    }
    if (exists("id_exclusion_matchcaseset")) {
      total_excluded_samples <- total_excluded_samples + length(id_exclusion_matchcaseset)
    }
    cat("### total sample(s) excluded =", total_excluded_samples, "\n")

    ### save excluded info
    saveRDS(object = id_outliers, file = paste0(path_outliers, "id-exclusions-outliers_", LABEL, ".rds"))

    ## save feature data ====
    cat("## saving feature data \n")
    saveRDS(object = df, file = paste0(path_out, "data-features_", LABEL, ".rds"))
  }
  # return ====
  if(outlier) {
    cat("# returning a list of feature data and outlier plots")
    return(list(df = df, plot_samples = plot_samples, plot_features = plot_features))
  } else {
    cat("# returning feature data")
    return(df)
  }
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
    stop("Idepic should be the same in data_features and data_samples")
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
    }
  }

  return(Y1)
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
