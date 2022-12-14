# function to run all possible models for a given set of covariates in a data frame and identfiy the 'best' model using akaike information criterion

#' Provide akaike information criterion for all possible covariate combinations
#'
#' @param data your data frame
#' @param exposure name of the exposure/predictor/independent variable column
#' @param outcome name of the outcome/dependent variable column
#' @param covariates name(s) of the covariate column(s)
#' @param match_ID name of the column with your case/control match ID
#' @param model name of the model you want to run: "lm" = stats::lm(); "clogit" = survival::clogit()
#' @returns an object from AICcmodavg::aictab()
#' @examples
#' model_test(data = my_data,
#' exposure = "BMI",
#' outcome = "breast_cancer",
#' covariates = c("age", "sex"),
#' match_ID = "match_caseset",
#' model = "clogit")
#'

model_test <- function(data, exposure, outcome, covariates,  match_ID = NA, model) {
  # lm() linear model ====
  if (model == "lm") {
    formulas <- unlist(lapply(1:1, function(n) {
      utils::combn(covariates, n,
        FUN = function(row) {
          paste0(
            outcome, " ~ ", exposure,
            " + ", paste0(row, collapse = " + ")
          )
        }
      )
    }))
    models <- lapply(formulas, function(frml) stats::lm(frml, data = data))
    names(models) <- formulas # rename models using formula

    ## model selection paramaters
    model_selection <- AICcmodavg::aictab(cand.set = models)
  }

  # clogit() logistic model ====
  else if (model == "clogit") {
    formulas <- unlist(lapply(1:1, function(n) {
      utils::combn(covariates, n,
        FUN = function(row) {
          paste0(
            outcome, " ~ ", exposure,
            " + ", paste0(row, collapse = " + "),
            paste0(" + ", "strata(", match_ID, ")")
          )
        }
      )
    }))

    ## run models
    models <- lapply(formulas, function(frml) {
      survival::clogit(as.formula(frml), data = data, method = "exact", na.action = "na.exclude")
    })
    names(models) <- formulas # rename models using formula

    ## model selection paramaters
    model_selection <- AICcmodavg::aictab(cand.set = models)
  }
  return(model_selection)
}
