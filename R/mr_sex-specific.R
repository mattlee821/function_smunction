#' @title harmonise_data_sex: perform harmonise_data() using sex specific data
#' @return a data frame
#' @param data_exposure exposure data frame
#' @param data_outcome outcome data frame
#' @param harmonise_action action for harmonising alleles (1,2,3)
#' @note
#' you must have a single column called 'sex' in your exposure and outcome data

harmonise_data_sex <- function(data_exposure, data_outcome, harmonise_action){
  levels_sex <- unique(c(data_exposure$sex, data_outcome$sex))
  list_harmonise <- list()
  for (i in levels_sex) {
    a <- data_exposure %>% filter(sex == i)
    b <- data_outcome %>% filter(sex == i)
    data_harmonise <- harmonise_data(a, b, action = harmonise_action)
    list_harmonise[[i]] <- data_harmonise
  }
  data_harmonise <- bind_rows(list_harmonise)
  rm(levels_sex, list_harmonise, i, a, b)
  return(data_harmonise)
}

#' @title mr_sex: perform mr() using sex specific data
#' @return a data frame
#' @param data_harmonise hamronised data frame
#' @param methods a list of methds supplied to "method_list"
#' @note
#' you must have a single column called 'sex' in your harmonised data frame

mr_sex <- function(data_harmonise, methods){
  levels_sex <- unique(data_harmonise$sex)
  list_mr <- list()
  for (i in levels_sex) {
    a <- data_harmonise %>% filter(sex == i)
    data_mr <- mr(a, method_list = methods)
    list_mr[[i]] <- data_mr
  }
  data_mr <- bind_rows(list_mr)
  rm(levels_sex, list_mr, i, a)
  return(data_mr)
}

#' @title mr_singlesnp_sex: perform mr_singlesnp() using sex specific data
#' @return a data frame
#' @param data_harmonise harmonised data frame
#' @note
#' you must have a single column called 'sex' in your harmonise data frame
mr_singlesnp_sex <- function(data_harmonise){
  levels_sex <- unique(data_harmonise$sex)
  list_sensitivity <- list()
  for (i in levels_sex) {
    a <- data_harmonise %>% filter(sex == i)
    data_sensitivity <- mr_singlesnp(a)
    list_sensitivity[[i]] <- data_sensitivity
  }
  data_sensitivity <- bind_rows(list_sensitivity)
  return(data_sensitivity)
}

#' @title mr_heterogeneity_sex: perform mr_heterogeneity() using sex specific data
#' @return a data frame
#' @param data_harmonise harmonised data frame
#' @param methods_heterogeneity a list of methods to supply to method_list
#' @note
#' you must have a single column called 'sex' in your harmonise data frame
mr_heterogeneity_sex <- function(data_harmonise, methods_heterogeneity){
  levels_sex <- unique(data_harmonise$sex)
  list_sensitivity <- list()
  for (i in levels_sex) {
    a <- data_harmonise %>% filter(sex == i)
    data_sensitivity <- mr_heterogeneity(a, method_list = methods_heterogeneity)
    list_sensitivity[[i]] <- data_sensitivity
  }
  data_sensitivity <- bind_rows(list_sensitivity)
  return(data_sensitivity)
}

#' @title mr_pleiotropy_test_sex: perform mr_pleiotropy_test() using sex specific data
#' @return a data frame
#' @param data_harmonise harmonised data frame
#' @note
#' you must have a single column called 'sex' in your harmonise data frame
mr_pleiotropy_test_sex <- function(data_harmonise){
  levels_sex <- unique(data_harmonise$sex)
  list_sensitivity <- list()
  for (i in levels_sex) {
    a <- data_harmonise %>% filter(sex == i)
    data_sensitivity <- mr_pleiotropy_test(a)
    list_sensitivity[[i]] <- data_sensitivity
  }
  data_sensitivity <- bind_rows(list_sensitivity)
  return(data_sensitivity)
}

#' @title mr_leaveoneout_sex: perform mr_leaveoneout() using sex specific data
#' @return a data frame
#' @param data_harmonise harmonised data frame
#' @note
#' you must have a single column called 'sex' in your harmonise data frame
mr_leaveoneout_sex <- function(data_harmonise){
  levels_sex <- unique(data_harmonise$sex)
  list_sensitivity <- list()
  for (i in levels_sex) {
    a <- data_harmonise %>% filter(sex == i)
    data_sensitivity <- mr_leaveoneout(a)
    list_sensitivity[[i]] <- data_sensitivity
  }
  data_sensitivity <- bind_rows(list_sensitivity)
  return(data_sensitivity)
}
