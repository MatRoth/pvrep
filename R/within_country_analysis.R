#' Within country analysis with replicate weights and/or plausible values
#'
#' @param cur_country  Character vector of length one. Must be the same as in the dataset.
#' @param country_var_name  Character vector of length one. Name of column containing country names.
#' @param cur_variables Character vector. Any column name in the dataset that is not a plausible value used in the analysis function.
#' @param cur_pv  Character vector. Column names of the current plausible values.
#' @param cur_pv_name Character vector of length one. The name the plausible value will be referenced by in the analysis function.
#' @param cur_func Analysis function
#' @param main_weight Vector with main weights.
#' @param rep_weights Complete PIAAC dataframe
#' @param rep_method Character vector of length one. "Fay" or "JK2" according to data
#' @param rho # Numeric of length one. Supplied by data provider.
#' @param dat  # Dataframe with replicate weights.
#'
#' @returns Dataframe with number of rows equal to computed statistics and the follwing columns.
#'
#' @author Matthias Roth
#'
#' @export
within_country_analysis <- function(cur_country, # Character vector of length one. Must be the same as in the dataset.
                                    country_var_name, # Character vector of length one. Name of column containing country names.
                                    cur_variables, # Character vector. Any column name in the dataset that is not a plausible value used in the analysis function.
                                    cur_pv = NULL, # Character vector. Column names of the current plausible values.
                                    cur_pv_name = NULL, # Character vector of length one. The name the plausible value will be referenced by in the analysis function.
                                    cur_func, # Analysis function
                                    main_weight, # Vector with main weights.
                                    rep_weights, # Complete PIAAC dataframe
                                    rep_method, #FAY or JK2 according to data
                                    rho = NULL,
                                    dat ){ # Dataframe with replicate weights.

  # Subset data

  cur_analysis_data <- dat[cur_country == dat[[country_var_name]],]

  # Create survey object
  cur_svy_obj <- survey::svrepdesign(data = cur_analysis_data,
                                     repweights = rep_weights,
                                     weights = cur_analysis_data[[main_weight]],
                                     variables = cur_variables,
                                     type = rep_method,
                                     rho = if(rep_method == "Fay") rho else NULL)

  # If PV are supplied
  if(!is.null(cur_pv) & !is.null(cur_pv_name)){
    # Outer loop -> Plausible values
    cur_pvs_formula_part <- paste(cur_pv,collapse = "+")
    current_pv_formula <- as.formula(glue::glue("{cur_pv_name}~{cur_pvs_formula_part}"))
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
