#' @title Performs a weighted redundancy analysis
#'
#' @description
#' \code{wrda} is formula-based implementation of weighted redundancy analysis.
#' 
#' @inheritParams cca0
#'
#' @param response matrix or data frame of the abundance data (dimension 
#' \emph{n} x \emph{m}). Rownames of \code{response}, if any, are carried 
#' through. Can be \code{NULL} if \code{cca_object} is supplied or 
#' if the response is \code{formula} is two-sided.							  
#' @param cca_object a vegan-type cca-object of \emph{transposed} 
#' \code{response}, from which centred abundance values
#'  and row and column weights can be obtained.
#' @param weights row weights (a vector). If not specified unit weights are 
#' used.
#
#' @details
#' The algorithm is a modified version of published R-code for weighted 
#' redundancy analysis (ter Braak, 2022).
#'
#' Compared to  \code{\link[vegan]{rda}}, \code{wrda} does not have residual 
#' axes, \emph{i.e.} no SVD or PCA of the residuals is performed.
#'
#' @returns
#' All scores in the \code{wrda} object are in scaling \code{"sites"} (1): 
#' the scaling with \emph{Focus on Case distances}.
#'
#' @references
#' ter Braak C.J.F. and  P. Šmilauer  (2018). Canoco reference manual
#' and user's guide: software for ordination (version 5.1x).
#' Microcomputer Power, Ithaca, USA, 536 pp.
#'
#' Oksanen, J., et al. (2022)
#' vegan: Community Ecology Package. R package version 2.6-4.
#' \url{https://CRAN.R-project.org/package=vegan}.
#'
#' @seealso \code{\link{scores.wrda}}, \code{\link{anova.wrda}},
#' \code{\link{print.wrda}}
#' 
#' @example demo/dune_wrda.R
#' 
#' @export
wrda <- function(formula, 
                 response = NULL, 
                 data, 
                 weights = rep(1 / nrow(data), nrow(data)),
                 traceonly = FALSE,
                 cca_object = NULL,
                 object4QR = NULL) {
  call <- match.call()
  if (!is.null(cca_object)) {
    # check whether eY works with and without trait covariates
    Yw <- cca_object$Ybar
    transp <- f_transpose(Yw, data, response)
    if (is.na(transp)) {
      warning("cca_object not usable.\n")
      cca_object <- NULL
    }
  }
  if (is.null(cca_object)) {
    if (is.null(response)) response <- get_response(formula)														
    Yn <- as.matrix(response) 
    K <- rep(1 / ncol(Yn), ncol(Yn))
    Wn <- weights / sum(weights)
    sWn <- sqrt(Wn)
    # transform to the unweighted case
    Yw <- as.matrix(response) * sWn
    # center
    Yw <- unweighted_lm_Orthnorm(Yw, matrix(sWn)) * sqrt(nrow(Yn) / (nrow(Yn) - 1))
    total_variance <- sum(Yw ^ 2) #total_inertia
  } else {  # use cca_object
    if (transp) {
      Yw <- t(Yw) 
      weights <- rev(cca_object$weights)
    } else {
      weights <- cca_object$weights
    }
    K <- weights[[1]]
    Wn <- weights[[2]]
    sWn <-sqrt(Wn)
    sumY <- cca_object$sumY
    total_variance <- cca_object$tot.chi
  } # cca_object
  msqr <- msdvif(formula, data = data, weights = Wn, XZ = TRUE,
                 object4QR = object4QR)
  eY <- qr.resid(msqr$qrZ, Yw)
  Yfit_X <- qr.fitted(msqr$qrXZ, eY)
  ssY_gZ <- sum(eY ^ 2) 
  ssY_XgZ <- sum(Yfit_X ^ 2)
  if (traceonly) {
    inert <- c(Condition_inertia = total_variance - ssY_gZ, Fit_inertia = ssY_XgZ)
    attr(inert, which = "rank") <-
      c(Cond = msqr$qrZ$rank - 1, Fit = msqr$qrXZ$rank - msqr$qrZ$rank + 1)
    return(inert)
  }
  svd_Yfit_X <- SVDfull(Yfit_X)
  biplot <- NULL
  eig <- svd_Yfit_X$d ^ 2
  names(eig) <- paste0("wRDA", seq_along(eig))
  CCA <- with(svd_Yfit_X, list(eig = eig, poseig = eig, u = u, v = v, 
                               rank = rank, qrank = msqr$qrXZ$rank,
                               tot.chi = ssY_XgZ, QR = msqr$qrXZ, 
                               biplot = biplot, envcentre = NULL, 
                               centroids = NULL))
  rank_CA <- min(nrow(Yw)-1, ncol(Yw)) 
  eig_CA <- rep(NA, rank_CA)
  names(eig_CA) <- paste0("wPCA", seq_len(rank_CA))
  if (ncol(msqr$Zw) == 1) {
    pCCA <- NULL 
  } else {
    pCCA <- list(rank = min(ncol(Yw), msqr$qrZ$rank), 
                 tot.chi = total_variance - ssY_gZ,
                 QR = msqr$qrZ, envcentre = NULL)
  }
  if (length(svd_Yfit_X$d) == 1) {
    diagd <- matrix(svd_Yfit_X$d)
  } else {
    diagd <- diag(svd_Yfit_X$d)
  }
  # need to be orthogonalized w.r.t Z
  site_axes <- list(
    site_scores = list(
      site_scores_unconstrained = qr.resid(msqr$qrZ, eY %*% svd_Yfit_X$v) / sWn,
      lc_env_scores = (svd_Yfit_X$u %*% diagd) /sWn
    )
  )
  species_axes <- list(species_scores = 
                         list(species_scores_unconstrained = svd_Yfit_X$v))
  object <- list(call = call, method = "wrda", tot.chi = total_variance,
                 formula = formula, site_axes = site_axes, 
                 species_axes = species_axes, Nobs = nrow(Yw), eigenvalues = eig,
                 weights = list(columns = K, rows = Wn),
                 data = data, Ybar = Yw, pCCA = pCCA, CCA = CCA, 
                 CA = list(tot.chi = ssY_gZ - ssY_XgZ, rank = rank_CA, eig = eig_CA),
                 inertia = "weighted variance"
  )
  class(object) <- "wrda"
  return(object)
}

#' @noRd
#' @keywords internal
f_transpose <- function(Yw, 
                        data,
                        response) {
  transpose <- if (nrow(data) == ncol(Yw)) TRUE else 
    if (nrow(data) == nrow(Yw)) FALSE else NA
  if (nrow(Yw) == ncol(Yw) && !is.null(response)) {
    transpose <- if (rownames(Yw)[1] == colnames(response)[1]) TRUE else
      if (rownames(Yw)[1] == rownames(response)[1]) FALSE else NA
  }
  return(transpose)
}		