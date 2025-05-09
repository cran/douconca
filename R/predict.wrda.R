#' @title Prediction from cca0 and wrda models
#'
#' @description
#' Prediction of response and lc scores from environment data
#' using \code{\link{cca0}} and \code{\link{wrda}} models.
#' 
#' @inheritParams scores.dcca
#'
#' @param object return value of \code{\link{cca0}} or \code{\link{wrda}}.
#' @param ...  Other arguments passed to the function (currently ignored).
#' @param type type of prediction, \code{c(
#' "response", "lc")} for 
#' response (expected abundance) and constrained scores for sites. 
#' @param newdata Data in which to look for variables with which to predict.
#' @param rank rank (number of axes to use). Default "full" for all axes 
#' (no rank-reduction).
#' @param weights list of weights of species and of sites in \code{newdata} 
#' when \code{type = "response"}, else ignored (default NULL
#' yielding equal species and site weights, both summing to 1). 
#' Example: weights = list(species = c(100, 1, 1), sites = c(1, 1, 1, 1)), in 
#' that order, with traits of three new species in newdata[[1]] and 
#' environmental values (and levels of factors) of four new sites in 
#' newdata[[2]]. Species weights are scaled to a sum of one.				  
#' 
#' @details
#' Variables that are in the model but not in \code{newdata} are set to their 
#' weighted means in the training data. Predictions are thus at the (weighted)
#' mean of the quantitative variables not included. Predictions with 
#' not-included factors are at the reference level (the first level of the 
#' factor).
#'
#' In a \code{\link{cca0}} model with \code{type = "response"},
#' many of the predicted values may be negative,
#' indicating expected absences (0) or small expected response values.
#' 
#' @returns a matrix with the predictions. The exact content of the matrix 
#' depends on the \code{type} of predictions that are being made.
#'
#' @example demo/dune_cca0.R
#' 
#' @export
predict.wrda <- function(object,
                         ...,
                         type = c("response",  "lc" ),
                         rank = "full",
                         newdata = NULL,
                         weights = NULL,
                         scaling = "symmetric") {
  type <- match.arg(type)
  object$formulaEnv <- object$formula
  if (rank == "full") {
    rank <- length(object$eigenvalues)
  }
  if (type == "response") {
    if (is.null(newdata)) {
      newdata <-  object$data
      weights <- object$weights
    } 
    if (is.null(weights[[2]])) weights$sites <- 
        rep(1 / nrow(newdata), nrow(newdata))
    
    if (length(weights[[2]]) != nrow(newdata)) {
      weights[[2]] <- rep(1 / nrow(newdata), nrow(newdata))
      warning("length of weights for sites does not match new environment data. ",
              "Site weights reset to equal weights.\n")
    }
  } else if (is.null(newdata)) {
    if (type %in% c("SNC", "lc")) {
      newdata <- if (inherits(object,"wrda"))object$data else object$data$dataEnv
    } 
  }
  ret <- switch(type,
                SNC = predict_env(object, newdata, rank),
                response = predict_response_wrda(object, newdata, rank, weights),
                lc = predict_lc(object, newdata, rank, scaling = scaling))
  return(ret)
}
