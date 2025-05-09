#' @title Performs (weighted) double constrained correspondence analysis (dc-CA)
#'
#' @description
#' Double constrained correspondence analysis (dc-CA) for analyzing
#' (multi-)trait (multi-)environment ecological data using library \code{vegan} 
#' and native R code. It has a \code{formula} interface which allows to assess,
#' for example, the importance of trait interactions in shaping ecological 
#' communities. The function \code{dc_CA} has an option to divide the abundance
#' data of a site by the site total, giving equal site weights. This division 
#' has the advantage that the multivariate analysis corresponds with an 
#' unweighted (multi-trait) community-level analysis, instead of being weighted
#' (Kleyer et al. 2012, ter Braak and van Rossum, 2025).
#'
#' @param formulaEnv two-sided or one-sided formula for the rows (samples) with 
#' row predictors in \code{dataEnv}. The left hand side of the formula is ignored if
#' it is specified in the \code{response} argument. Specify row covariates (if any) by
#' adding \code{+ Condition(covariate-formula)} to \code{formulaEnv} as 
#' in \code{\link[vegan]{rda}}. The \code{covariate-formula} should not contain 
#' a \code{~} (tilde). Default: \code{NULL} for \code{~.}, i.e. all variables 
#' in \code{dataEnv} are predictor variables.
#' @param formulaTraits formula or one-sided formula for the columns (species) 
#' with column predictors in \code{dataTraits}. When two-sided, the left hand 
#' side of the formula is not used. Specify column covariates (if any ) by 
#' adding \code{+ Condition(covariate-formula)} to \code{formulaTraits} as 
#' in \code{\link[vegan]{cca}}. The \code{covariate-formula} should not contain 
#' a \code{~} (tilde). Default: \code{NULL} for \code{~.}, i.e. all variables 
#' in \code{dataTraits} are predictor traits.
#' @param response matrix, data frame of the abundance data 
#' (dimension \emph{n} x \emph{m}) or list with community weighted means (CWMs)
#' from \code{\link{fCWM_SNC}}, NULL. If \code{NULL}, the response should be
#' at the left-hand side of \code{formulaEnv}. 
#' See Details for analyses starting from community
#' weighted means. Rownames of \code{response}, if any, are carried through.
#' @param dataEnv matrix or data frame of the row predictors, with rows 
#' corresponding to those in \code{response}. (dimension \emph{n} x \emph{p}).
#' @param dataTraits matrix or data frame of the column predictors, with rows 
#' corresponding to the columns in \code{response}.
#' (dimension \emph{m} x \emph{q}).
#' @param divideBySiteTotals logical; default \code{TRUE} for closing the 
#' data by dividing the rows in the \code{response} by their total.
#' However, the default is \code{FALSE}, when the species totals are 
#' proportional to \code{N2*(N-N2)} with  
#' \code{N2} the Hill numbers of order 2 of the species 
#' and \code{N} the number of sites,
#' as indicator that the response data have been pre-processed to
#' \code{N2}-based marginals using \code{\link{ipf2N2}}.												
#' @param dc_CA_object  optional object from an earlier run of this function. 
#' Useful if the same formula for the columns (\code{formulaTraits}), 
#' \code{dataTraits} and \code{response} are used with a new formula for the 
#' rows. If set, the data of the previous run is used and the result of its 
#' first step is taken for the new analysis and \code{env_explain} is set 
#' to \code{FALSE}. 
#' @param env_explain logical (default \code{TRUE})
#' for calculation of the inertia explained 
#' by the environmental variable (based on a CCA of abundance 
#' (with \code{divideBySiteTotals}, if true) on the environmental formula).
#' @param use_vegan_cca default \code{TRUE}.
#' @param verbose logical for printing a simple summary (default: TRUE)
#
#' @details
#' Empty (all zero) rows and columns in \code{response} are removed from the 
#' \code{response} and the corresponding rows from \code{dataEnv} and 
#' \code{dataTraits}. Subsequently, any columns with missing values are 
#' removed from  \code{dataEnv} and \code{dataTraits}. It gives an error 
#' ('name_of_variable' not found), if variables with missing entries are
#' specified in \code{formulaEnv} and \code{formulaTraits}.
#'
#' Computationally, dc-CA can be carried out by a single singular value 
#' decomposition (ter Braak et al. 2018), but it is here computed in two steps.
#' In the first step, the transpose of the \code{response} is regressed on to 
#' the traits (the column predictors) using \code{\link[vegan]{cca}} with 
#' \code{formulaTraits}. The column scores of this analysis (in scaling 1) are 
#' community weighted means (CWM) of the orthonormalized traits. These are then 
#' regressed on the environmental (row) predictors using \code{\link{wrda}} 
#' with \code{formulaEnv} or using \code{\link[vegan]{rda}}, if site weights 
#' are equal.
#'
#' A dc-CA can be carried out on, what statisticians call, the sufficient
#' statistics of the method. This is useful, when the abundance data are not 
#' available or could not be made public in a paper attempting reproducible 
#' research. In this case, \code{response} should be a list
#' with as first element community weighted means 
#' (e.g. \code{list(CWM = CWMs)}) with respect to the 
#' traits, and the trait data, and, optionally, further list elements, for functions
#' related to \code{dc_CA}. The minimum is a 
#' \code{list(CWM = CWMs, weight = list(columns = species_weights))} with CWM a matrix
#' or data.frame, but then \code{formulaEnv}, \code{formulaTraits}, 
#' \code{dataEnv}, \code{dataTraits} must be specified in the call to 
#' \code{dc_CA}. The function \code{\link{fCWM_SNC}} and its example
#' show how to set the
#' \code{response} for this and helps to create the \code{response} from 
#' abundance data in these non-standard applications of dc-CA. Species and site 
#' weights, if not set in \code{response$weights} can be set by a variable
#' \code{weight} in the data frames \code{dataTraits} and \code{dataEnv}, 
#' respectively, but formulas should then not be \code{~.}.
#'
#' The statistics and scores in the example \code{dune_dcCA.r}, have been 
#' checked against the results in Canoco 5.15 (ter Braak & Šmilauer, 2018).
#'
#' @returns
#' A list of \code{class} \code{dcca}; that is a list with elements
#' \describe{
#' \item{CCAonTraits}{a \code{\link[vegan]{cca.object}} from the 
#' \code{\link[vegan]{cca}} analysis of the transpose of the closed 
#' \code{response} using formula \code{formulaTraits}.}
#' \item{formulaTraits}{the argument \code{formulaTraits}. If the formula was
#' \code{~.}, it was changed to explicit trait names.}
#' \item{data}{a list of \code{Y}, \code{dataEnv} and \code{dataTraits}, 
#' after removing empty rows and columns in \code{response} and after closure if 
#' \code{divideBySiteTotals = TRUE} and with the corresponding rows in 
#' \code{dataEnv} and \code{dataTraits} removed.}
#' \item{weights}{a list of unit-sum weights of columns and rows. The names of 
#' the list are \code{c("columns", "rows")}, in that order.}
#' \item{Nobs}{number of sites (rows).}
#' \item{CWMs_orthonormal_traits}{Community weighted means w.r.t. 
#' orthonormalized traits.}
#' \item{RDAonEnv}{a \code{\link{wrda}} object or 
#' \code{\link[vegan]{cca.object}} from the 
#' \code{\link{wrda}} or, if with equal row weights, 
#' \code{\link[vegan]{rda}} analysis, respectively of the column scores of the
#' \code{cca}, which are the CWMs of orthonormalized traits, using formula 
#' \code{formulaEnv}.}
#' \item{formulaEnv}{the argument \code{formulaEnv}. If the formula was 
#' \code{~.}, it was changed to explicit environmental variable names.}
#' \item{eigenvalues}{the dc-CA eigenvalues (same as those of the 
#' \code{\link[vegan]{rda}} analysis).}
#' \item{c_traits_normed0}{mean, sd, VIF and (regression) coefficients of 
#' the traits that define the dc-CA axes in terms of the 
#' traits with t-ratios missing indicated by \code{NA}s for 'tval1'.}
#' \item{inertia}{a one-column matrix with, at most,
#' six inertias (weighted variances):
#' \itemize{
#' \item total: the total inertia.
#' \item conditionT: the inertia explained by the condition in 
#' \code{formulaTraits} if present (neglecting row constraints).
#' \item traits_explain: the trait-structured variation, \emph{i.e.}
#' the inertia explained by the traits (without constaints on the rows and 
#' conditional on the Condition in \code{formulaTraits}).
#' This is the  maximum that the row predictors could explain in dc-CA
#' (the sum of the last two items is thus less than this value).
#' \item env_explain: the environmentally structured variation, \emph{i.e.}
#' the inertia explained by the environment (without constraints on the 
#' columns but conditional on the Condition \code{formulaEnv}).
#' This is the maximum that the column predictors could explain in 
#' dc-CA (the item \code{constraintsTE} is thus less than this value).
#' The value is \code{NA}, if there is collinearity in the environmental data.
#' \item conditionTE: the trait-constrained variation explained by the 
#' condition in \code{formulaEnv}.
#' \item constraintsTE: the trait-constrained variation explained by the 
#' predictors (without the row covariates).
#' }
#' }
#' }
#' If \code{verbose} is \code{TRUE} (or after \code{out <- print(out)} is 
#' invoked) there are three more items.
#' \itemize{
#' \item \code{c_traits_normed}: mean, sd, VIF and (regression) coefficients of
#' the traits that define the dc-CA trait axes (composite traits), and their 
#' optimistic t-ratio.
#' \item \code{c_env_normed}: mean, sd, VIF and (regression) coefficients of 
#' the environmental variables that define the dc-CA axes in terms of the 
#' environmental variables (composite gradients), and their optimistic t-ratio.
#' \item \code{species_axes}: a list with four items
#'   \itemize{
#'   \item \code{species_scores}: a list with names 
#'   \code{c("species_scores_unconstrained", "lc_traits_scores")} with the
#'   matrix with species niche centroids along the dc-CA axes (composite 
#'   gradients) and the matrix with linear combinations of traits.
#'   \item \code{correlation}: a matrix with inter-set correlations of the 
#'   traits with their SNCs.
#'   \item \code{b_se}: a matrix with (unstandardized) regression coefficients 
#'   for traits and their optimistic standard errors.
#'   \item \code{R2_traits}: a vector with coefficient of determination (R2) 
#'   of the SNCs on to the traits. The square-root thereof could be called 
#'   the species-trait correlation in analogy with the species-environment 
#'   correlation in CCA.
#'  }
#'  \item \code{sites_axes}: a list with four items
#'    \itemize{
#'    \item \code{site_scores}: a list with names 
#'    \code{c("site_scores_unconstrained", "lc_env_scores")} with the matrix 
#'    with community weighted means (CWMs) along the dc-CA axes (composite 
#'    gradients) and the matrix with linear combinations of environmental 
#'    variables.
#'    \item \code{correlation}: a matrix with inter-set correlations of the 
#'    environmental variables with their CWMs.
#'    \item \code{b_se}: a matrix with (unstandardized) regression coefficients 
#'    for environmental variables and their optimistic standard errors.
#'    \item \code{R2_env}: a vector with coefficient of determination (R2) of 
#'    the CWMs on to the environmental variables. The square-root thereof 
#'    has been called the species-environmental correlation in CCA.
#'  }
#' }
#' All scores in the \code{dcca} object are in scaling \code{"sites"} (1): 
#' the scaling with \emph{Focus on Case distances} .
#'
#' @references
#' Kleyer, M., Dray, S., Bello, F., Lepš, J., Pakeman, R.J., Strauss, B., Thuiller,
#' W. & Lavorel, S. (2012) Assessing species and community functional responses to
#' environmental gradients: which multivariate methods?
#' Journal of Vegetation Science, 23, 805-821.
#' \doi{10.1111/j.1654-1103.2012.01402.x}
#'
#' ter Braak, CJF, Šmilauer P, and Dray S. (2018). Algorithms and biplots for
#' double constrained correspondence analysis.
#' Environmental and Ecological Statistics, 25(2), 171-197.
#' \doi{10.1007/s10651-017-0395-x}
#'
#' ter Braak C.J.F. and  P. Šmilauer  (2018). Canoco reference manual
#' and user's guide: software for ordination (version 5.1x).
#' Microcomputer Power, Ithaca, USA, 536 pp.
#' 
#' ter  Braak, C.J.F. and van Rossum, B. (2025).
#' Linking Multivariate Trait Variation to the Environment: 
#' Advantages of Double Constrained Correspondence Analysis 
#' with the R Package Douconca. Ecological Informatics, 88.
#' \doi{10.1016/j.ecoinf.2025.103143}
#' 
#' Oksanen, J., et al. (2024).
#' vegan: Community Ecology Package. R package version 2.6-6.1.
#' \url{https://CRAN.R-project.org/package=vegan}.
#'
#' @seealso \code{\link{plot.dcca}}, \code{\link{scores.dcca}}, 
#' \code{\link{print.dcca}} and \code{\link{anova.dcca}}
#' 
#' @example demo/dune_dcCA.R
#' @export
dc_CA <- function(formulaEnv = NULL, 
                  formulaTraits = NULL,
                  response = NULL, 
                  dataEnv = NULL,
                  dataTraits = NULL,
                  divideBySiteTotals = NULL,
                  dc_CA_object = NULL,
                  env_explain = TRUE,
                  use_vegan_cca = FALSE,
                  verbose = TRUE) {
  # response matrix or data frame, dataEnv and dataTraits data frames 
  # in which formualaE and formulaT are evaluated
  # dc_CA_object = result (value) of a previous run, 
  # can be used to save computing time for
  # runs that modify the formula for samples (step2: RDAonEnv) only
  # The step1 (CCAonTraits and the data and formulaTraits) are taken 
  # from dc_CA_object into the new result.
  # If set, formulaTraits, response, dataEnv, dataTraits are not used at all 
  # and have no efffect on the result
  call <- match.call()																					 
  if (is.null(response) && is.null(dc_CA_object)) {
    response <- get_response(formulaEnv)
  } else {
    if (is.list(response) && inherits(response[[1]], c("matrix", "data.frame"))) {
      # response is a list of CWMs_orthonormal_traits and a weights list
      if (any(is.na(response[[1]]))) {
        stop("The CWMs should not have missing entries.\n")
      }
      # create a sufficient dc_CA_object object
      dc_CA_object <- response
    }
  }
  if (is.null(dc_CA_object)) {
    out0 <- check_data_dc_CA(formulaEnv, formulaTraits, response, dataEnv, 
                             dataTraits, divideBySiteTotals, call)
    tY <- t(out0$data$Y)
    formulaTraits <- change_reponse(out0$formulaTraits, "tY", out0$data$dataTraits)
    environment(formulaTraits) <- environment()
    step1 <- if (use_vegan_cca) vegan::cca(formulaTraits, data = out0$data$dataTraits) else
      cca0(formulaTraits, response = tY, data = out0$data$dataTraits)
    if (is.null(vegan::eigenvals(step1, model = "constrained"))) {
      message("Trait constraints or CWMs are constant,", 
              "no dc-CA analysis possible.\n",
              "Inspect return value$CWM for issues.\n")
      ms <- f_wmean(out0$formulaTraits, tY = out0$data$Y, out0$data$dataTraits,
                    weights=out0$weights, name = "CWM", SNConly = TRUE)
      out0$CWM <- ms
      out0$CCAonTraits <- step1
      return(out0)
    }
    n <- out0$Nobs
    CWMs_orthonormal_traits <-
      scores(step1, display = "species", scaling = "species",
             choices = seq_len(rank_mod(step1))) * sqrt((n - 1) / n)
    rownames(CWMs_orthonormal_traits) <- 
      if (is.null(rownames(out0$data$Y))) paste0("Sam", seq_len((n))) else 
        rownames(out0$data$Y)
    if (rownames(CWMs_orthonormal_traits)[1] == "col1") {
      rownames(CWMs_orthonormal_traits) <- paste0("Sam", seq_len((n)))
    }
    out1 <- c(list(CCAonTraits = step1), 
              out0,
              list (CWMs_orthonormal_traits = CWMs_orthonormal_traits))
  } else {
    dc_CA_object$call <- call
    env_explain <- FALSE
    if (!is.null(formulaEnv)) {
      dc_CA_object$formulaEnv <- formulaEnv 
    } else if (is.null(dc_CA_object$formulaEnv)) {
      warning("formulaEnv set to ~. in dc_CA.\n")
      dc_CA_object$formulaEnv <- ~.
    }
    if (inherits(dc_CA_object, "dcca")) {
      out1 <- dc_CA_object[c("CCAonTraits", "formulaTraits", "formulaEnv",
                             "data", "call", "weights", "Nobs",
                             "CWMs_orthonormal_traits")]
    } else if (is.list(response) && 
               inherits(response[[1]], c("matrix", "data.frame"))) {
      out1 <- checkCWM2dc_CA(dc_CA_object, dataEnv, dataTraits, formulaTraits)
    } else{
      stop("The class of dc_CA_object should be dcca, whereas it is now:",
           class(dc_CA_object)[1], class(dc_CA_object)[2], ".\n")
    }
  }
  formulaEnv <- change_reponse(out1$formulaEnv, "out1$CWMs_orthonormal_traits", 
                               data = out1$data$dataEnv)
  environment(formulaEnv) <- environment()
  if (diff(range(out1$weights$rows)) < 1.0e-4 / out1$Nobs) {
    step2 <- vegan::rda(formulaEnv, data = out1$data$dataEnv)
    eigenvalues <- vegan::eigenvals(step2, model = "constrained")
  } else {
    step2 <- try(wrda(formulaEnv, response = out1$CWMs_orthonormal_traits, 
                      weights = out1$weights$rows, data = out1$data$dataEnv))
    eigenvalues <- if (!inherits(step2, "try-error")) step2$CCA$eig else NULL
  }
  if (is.null(eigenvalues)) {
    message("Environmental constraints or SNCs are constant,",
            "no dc-CA analysis possible.\n",
            "Inspect return value$SNC for issues.\n")
    mt <- f_wmean(formulaEnv = out1$formulaEnv, tY = t(out1$data$Y), out1$data$dataEnv,
                  weights = out1$weights, name = "SNC", SNConly = TRUE)
    out0$SNC <- mt
    out0$CCAonTraits <- step1
    out0$RDAonEnv <- step2
    return(out0)
  }
  if (length(eigenvalues)) {
    names(eigenvalues) <- paste0("dcCA", seq_along(eigenvalues))
  } else {
    warning(" no eigenvalue(s) obtained, no dc-CA axis found.\n")
  }
  out <- c(out1, list(RDAonEnv = step2,
                      eigenvalues =  eigenvalues))
  out$c_traits_normed0 <- try(f_canonical_coef_traits2(out))
  inertia <- try(f_inertia(out, env_explain))
  if (inherits(inertia, "try-error")) {
    warning("could not obtain inertias.\n") 
    print(inertia)
  }
  out$inertia <- inertia
  if (inherits(out$RDAonEnv, "rda")) {
    class(out) <- c("dccav", "dcca", "list") 
  } else {
    class(out) <- c("dcca", "list")
  }
  if (verbose) {
    out <- print(out)
  }
  return(out)
}

#' @noRd
#' @keywords internal
get_response <- function(formula) {
  trms <- terms(formula, specials = "Condition", keep.order = TRUE) 
  response_nam <- rownames(attr(trms, "factors"))[attr(trms, "response")]
  if (length(response_nam)) {
    response <- eval(parse(text = response_nam), envir = environment(formula))
  } else {
    stop("Specify the response.\n")
  }
  return(response)
}