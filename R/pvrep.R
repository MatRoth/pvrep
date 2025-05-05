#' Analysis with replicate weights and/or plausible values using the `survey` and `mitools` packages
#'
#' Currently tailored towards the PIAAC Cycle 2 dataset. Might be expanded in future.
#' The aim of the function is to make it simpler to write analysis code that works
#' in the replicate weights/plausible value framework. The function abstracts the loops
#' over replicate weights/plausible values and pools the results in an appropriate manner.
#'
#' @param cur_country  Character vector of length one. Must be the same as in the dataset. Defaults to NULL if dataset does not need to be subset.
#' @param country_var_name  Character vector of length one. Name of column containing country names. Defaults to NULL if dataset does not need to be subset.
#' @param cur_variables Character vector. Any column name in the dataset that is not a plausible value used in the analysis function.
#' @param cur_pv  List of character vectors or one character vector. Column names of the current plausible values.
#' @param cur_pv_name Character vector. The names the plausible values will be referenced by in the analysis function.
#' @param cur_func Analysis function. See Details for information on how to define the analysis function.
#' @param main_weight Character of length one. Name of the column with main weights.
#' @param rep_weights Character vector of the names of the replicate weights column names in the dataset.
#' @param rep_method Character vector of length one. "Fay" or "JK2" according to data
#' @param rho Numeric of length one. Fay`s method adjustment factor. Supplied by data provider.
#' @param dat Data with all variables, plausible values and replicate weights used in the analysis function.
#'
#' @details
#' The analysis function provided to `pvrep()` performs the analysis taking into account
#' replicate weights or replicate weights & plausible values.
#'
#' If only replicate weights are required (i.e. the analysis does not involve plausible values)
#' the function needs the following arguments in this order: dat, current_weight, current_variables.
#' * `dat` is the dataframe provided to pvrep.
#' * `current_weight` is the weight used for the analysis.
#' * `current_variables` are any other variable names (i.e. provided as a character vector) used in the analysis.
#'
#' This can be useful for performing loops over different combination of variables)
#' If plausible values are needed in addition to replicate weights the function needs the
#' following argument in this order: dat, cur_weights, cur_pv_name, cur_variables.
#' * `dat` is the dataframe provided to pvrep.
#' * `current_weight` is a variable containing the current weight used for the analysis.
#' * `cur_pv_name` is a variable containing the current plausible value(s) used in the analysis.
#' * `current_variables` are any other variable names (i.e. provided as a character vector) used in the analysis.
#'
#' @returns Dataframe with number of rows equal to computed statistics and the following columns.
#' * parameter: Name of the parameter.
#' * est: Point estimate of the parameter.
#' * se: Standard error of the parameter.
#' * df: Degrees of freedom according to Lumely (https://notstatschat.rbind.io/2019/06/08/design-degrees-of-freedom-brief-note/)
#' * ci_lwr: Lower bound of the 95% confidence interval.
#' * ci_upr: Upper bound of the 95% confidence interval.
#' * p_value: p value of a hypothesis test with H0 = 0.
#'
#' @author Matthias Roth
#'
#' @export
pvrep <- function(cur_country = NULL,
                  country_var_name = NULL,
                  cur_variables,
                  cur_pv = NULL,
                  cur_pv_name = NULL,
                  cur_func,
                  main_weight,
                  rep_weights,
                  rep_method,
                  rho = NULL,
                  dat){

  # Enforce list, if single vector of cur_pv is supplied
  if(!is.list(cur_pv)) cur_pv <- list(cur_pv)

  # Subset data if necessary
  if(is.null(cur_country) & is.null(country_var_name)){
    cur_analysis_data <- dat} else {
      cur_analysis_data <- dat[cur_country == dat[[country_var_name]],]
    }


  # Create survey object
  cur_svy_obj <- survey::svrepdesign(data = cur_analysis_data,
                                     repweights = cur_analysis_data[[rep_weights]],
                                     weights = cur_analysis_data[[main_weight]],
                                     variables = cur_variables,
                                     type = rep_method,
                                     rho = if(rep_method == "Fay") rho else NULL)

  # If PV are supplied
  if(!is.null(cur_pv) & !is.null(cur_pv_name)){
    # Outer loop -> Plausible values
    cur_pvs_formula_part <- map_chr(cur_pv,paste,collapse = "+")
    current_pv_formula <- map2(cur_pv_name,
                               cur_pvs_formula_part,
                               \(name_elem, form_elem) as.formula(glue::glue("{name_elem}~{form_elem}")))
    mi_results <- mitools::withPV(
      mapping = current_pv_formula,
      data = cur_analysis_data,
      rewrite = F,
      action = \(cur_data_pv){
        # Inner loop -> Replicate Weights
        rep_res <- survey::withReplicates(
          design = cur_svy_obj,
          theta = \(cur_weights, cur_data){
            cur_func(cur_data_pv,cur_weights, cur_pv_name, cur_variables) # <- observe function order
          })
        rep_res}) # Contains 10 PV estimates with SE from replicate weights
    #Pool pv results
    combined_para <- mitools::MIcombine(mi_results)
  }else{
    # If PV are not supplied
    combined_para <- survey::withReplicates(
      design = cur_svy_obj,
      theta = \(cur_weights, cur_data){
        cur_func(cur_analysis_data,cur_weights, cur_variables)}) # <- observe function order
  }

  results <- tibble::tibble(
    parameter = if(class(combined_para) == "svrepstat") names(combined_para) else names(combined_para$coefficients), # Parameter names
    est = if(class(combined_para) == "svrepstat") as.numeric(combined_para) else coef(combined_para), # Combines PV estimates
    se = sqrt(diag(vcov(combined_para))) |> # Combines SE estimates from replicate weights + PVs if avaliable
      as.numeric(),
    df = if(class(combined_para) == "svrepstat") qr(cur_svy_obj$repweights)$rank-1 else combined_para$df)

  # Add ci
  results["ci_lwr"] <- results$est-qt(0.975,results$df)*results$se
  results["ci_upr"] <- results$est+qt(0.975,results$df)*results$se
  results["p_value"]<- 2*pt(q = abs((results$est/results$se)),df = results$df,lower.tail = F)|> round(3)
  # Return tibble
  results
}
