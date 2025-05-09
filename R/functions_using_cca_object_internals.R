#' @noRd
#' @keywords internal
f_inertia <- function(object,
                      env_explainL = TRUE) {
  # function f_inertia uses vegan 2.6-4 internal structure
  # object: a dccav object, results of dc_CA
  # value:  a matrix (currently with 1 column) with in first column the inertias
  conditionE <- NULL
  env_explain <- NA
  if (is.null(object$CCAonTraits)) {
    total <- object$tot.chi
    conditionT <- object$conditionT
    env_explain <- object$inertia["env_explain"]
    if (is.null(env_explain)) env_explain <- NA
  } else {
    total <- object$CCAonTraits$tot.chi
    conditionT <- object$CCAonTraits$pCCA$tot.chi
	if (env_explainL && !is.null(object$data$Y)) {
      env_explain <- cca0(object$formulaEnv, response = object$data$Y,
        data = object$data$dataEnv, traceonly = TRUE,
        cca_object = object$CCAonTraits,object4QR = object$RDAonEnv)
      names(env_explain) <- NULL
      conditionE <- if (env_explain[1] > 1e-10) env_explain[1]
      env_explain <- env_explain[2]
    }
  }
  if (is.na(env_explain)) env_explain <- NULL
  names(env_explain) <- NULL
  inertia <- cbind(c(total = total,
                     conditionT = conditionT,
					 conditionE = conditionE,
                     traits_explain = object$RDAonEnv$tot.chi,
                     env_explain = env_explain,
                     conditionTE = object$RDAonEnv$pCCA$tot.chi,
                     constraintsTE = object$RDAonEnv$CCA$tot.chi))
  colnames(inertia) <- "weighted variance"
  expla <- c("total inertia (= weighted variation)",
             "variation fitted by the trait condition",
			 "variation fitted by the environmental condition",
             "trait-constrained variation", 
             "environment-constrained variation",
             "trait-constrained variation explained by the condition in formulaEnv",
             "trait-constrained variation explained by the predictors in formulaEnv")
  names(expla) <- c("total", "conditionT", "conditionE", "traits_explain", 
                    "env_explain", "conditionTE", "constraintsTE")
  attr(inertia, which = "meaning") <- 
    matrix(expla[rownames(inertia)], ncol = 1, 
           dimnames = list(rownames(inertia), "meaning"))
  return(inertia)
}

#' @noRd
#' @keywords internal
get_QR <- function(object, model = "CCA"){
  # function get_QR uses vegan 2.6-4 internal structure
  # gets the qr decompostion of object
  # model = "CCA" or "pCCA"
  if (model == "CCA") {
    QR <- object$CCA$QR 
  } else if (model == "pCCA") {
    QR <- object$pCCA$QR
  } else {
    stop("model not supported.\n")
  }
  return(QR)
}
