#' Forward selection of traits or environmental variables using dc-CA.
#' 
#' @inheritParams FS.wrda
# @inheritParams dc_CAForwardSelectionEnv 
#' @param mod dc-CA model. 
#' @param select Character, the default \code{"traits"} selects
#' trait variables based on the environmental model of \code{mod}.
#' Set \code{select = "e"} to select environmental variables based 
#' on the trait model of \code{mod}.
#' @param consider character vector of names considered for addition, either 
#' trait names in \code{mod$data$dataTraits} or
#' environmental variable names in in \code{mod$data$dataEnv}.
#' The names may include transformations of predictor variables,
#' such as \code{log(.)},
#' if \code{consider} does not include factors or if \code{factor2categories=FALSE}.
#' If \code{consider} includes factors, such transformations
#' give in a error in the default setting (\code{factor2categories=TRUE}).
#' @param initial_model character specifying what should be inside 
#' \code{Condition()}. Default: \code{"1"} (nothing, the intercept only).
#' Examples: With selection of environmental varialbles, 
#' \code{"region"} for a within-region analysis or
#' \code{"A*B"} for a within analysis specified by the interaction 
#' of factors \code{A} and \code{B}, with \code{region, A, B} in the data.
#' Analogously, a trait or trait factor can be the initial model.
#' 
#' @details 
#' The selection is on the basis of the additional fit (inertia) of a variable 
#' given the variables already in the model.
#' 
#' The function \code{FS} does not implement the max test and may thus 
#' be liberal. It is recommended to test the final model (element
#' \code{model_ final} of the returned object) by \code{anova}.
#' 
#' @returns list with three elements: \code{final...} with selected variables,
#' \code{model_final}, and \code{process} with account of the selection process 
#' 
#' @seealso \code{\link{dc_CA}}
#' 
#' @example demo/dune_forward_slct.R
#' 
#' @export
#'
FS.dcca <- function(mod,  
                    ...,
                    select = c("traits", "env"),
                    consider,
                    permutations = 999,
                    n_axes = "all",
                    initial_model = "1",
                    factor2categories = TRUE,
                    test = TRUE, 
                    threshold_P = 0.10, 
                    PvalAdjustMethod = "holm",
                    max_step = 10,
                    verbose = FALSE) {
  select <- match.arg(select)
  FUNFS <- if (select == "traits") FStraits else if (select == "env") FSenv
  ret <- FUNFS(mod, 
               consider = consider,
               permutations = permutations,
               n_axes = n_axes,
               initial_model = initial_model,
               factor2categories = factor2categories,
               test = test, 
               threshold_P = threshold_P, 
               PvalAdjustMethod = PvalAdjustMethod,
               max_step = max_step,
               verbose = verbose)
  return(ret)
}

#' @noRd
#' @keywords internal
FStraits <- function(mod,
                     consider = NULL, 
                     permutations = 999,
                     n_axes = "all",
                     initial_model = "1",
                     factor2categories = TRUE,
                     test = TRUE, 
                     threshold_P = 0.10, 
                     PvalAdjustMethod = "holm",
                     max_step = 10,
                     verbose = FALSE) {
  # function for forward selection of predictor variables using wrda.
  if (is.null(mod$SNCs_orthonormal_env)) {
    SNCs_orthonormal_env <- fSNC_ortho(mod)
  } else {
    SNCs_orthonormal_env <- mod$SNCs_orthonormal_env
  }
  # perform a wrda of on the traits
  n <- nrow(SNCs_orthonormal_env)
  modf <- wrda(mod$formulaTraits, 
               response = SNCs_orthonormal_env * sqrt((n - 1) / n),
               data = mod$data$dataTraits,
               weights = mod$weights$columns,
               object4QR = mod$CCAonTraits)
  ret <- FS(modf, 
            consider = consider,
            permutations = permutations,
            n_axes = n_axes,
            initial_model = initial_model,
            factor2categories = factor2categories,
            test = test, 
            threshold_P = threshold_P, 
            PvalAdjustMethod = PvalAdjustMethod,
            max_step = max_step,
            verbose = verbose)
  model_final <- dc_CA(mod$formulaEnv, formulaTraits = ret$formula,
                       response = mod$data$Y, 
                       mod$data$dataEnv, ret$model_final$data,
                       divideBySiteTotals = FALSE,
                       verbose = FALSE)
  
  ret$model_final <- model_final
  return(ret)
}

#' @noRd
#' @keywords internal
FSenv <- function(mod,
                  consider = NULL, 
                  permutations = 999,
                  n_axes = "all",
                  initial_model = "1",
                  factor2categories = TRUE,
                  test = TRUE, 
                  threshold_P = 0.10, 
                  PvalAdjustMethod = "holm",
                  max_step = 10,
                  verbose = FALSE) {
  # function for forward selection of predictor variables using wrda.
  CWMs_orthonormal_traits <- mod$CWMs_orthonormal_traits
  # perform a wrda of on the traits
  modf <- wrda(mod$formulaEnv, 
               response = CWMs_orthonormal_traits,
               data = mod$data$dataEnv,
               weights = mod$weights$rows,
               object4QR = mod$RDAonEnv)
  ret <- FS(modf, 
            consider = consider,
            permutations = permutations,
            n_axes = n_axes,
            initial_model = initial_model,
            factor2categories = factor2categories,
            test = test, 
            threshold_P = threshold_P, 
            PvalAdjustMethod = PvalAdjustMethod,
            max_step = max_step,
            verbose = verbose)
  model_final <- dc_CA(formulaEnv = ret$formula, 
                       formulaTraits = mod$formulaTraits,
                       response = mod$data$Y, 
                       ret$model_final$data, 
                       mod$data$dataTraits,
                       divideBySiteTotals = FALSE,
                       verbose = FALSE)
  ret$model_final <- model_final
  return(ret)
}
