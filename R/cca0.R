#' @title Performs a canonical correspondence analysis
#'
#' @description
#' \code{cca0} is formula-based implementation of canonical correspondence 
#' analysis.
#'
#' @param formula one or two-sided formula for the rows (samples) with row 
#' predictors in \code{data}. The left hand side of the formula is ignored if
#' it is specified in the next argument (\code{response}). Specify row 
#' covariates (if any ) by adding \code{+ Condition(covariate-formula)} to 
#' \code{formula} as in \code{\link[vegan]{rda}}. The \code{covariate-formula}
#' should not contain a \code{~} (tilde). 
#' @param response matrix or data frame of the abundance data (dimension 
#' \emph{n} x \emph{m}). Rownames of \code{response}, if any, are carried 
#' through. BEWARE: all rows and columns should have positive sums!
#' Can be \code{NULL} if \code{cca_object} is supplied or 
#' if the response is \code{formula} is two-sided.
#' @param data matrix or data frame of the row predictors, with rows 
#' corresponding to those in \code{response} (dimension \emph{n} x \emph{p}).
#' @param cca_object a vegan-type cca-object of \emph{transposed} 
#' \code{response}, from which chisq_residuals and row and column weights can 
#' be obtained.
#' @param object4QR a vegan-type cca-object with weighted QR's for 
#' \code{formula}, i.e. \code{qr(Z)} and \code{qr(XZ)} obtainable via 
#' \code{get_QR(object4QR, model = "pCCA")} and
#' \code{get_QR(object4QR, model = "CCA")}, respectively. 
#' @param traceonly logical, default \code{FALSE}. If \code{TRUE}, only 
#' the explained variance of the predictors and the \code{Condition()}
#' are returned, \emph{i.e} without performing a singular value
#' decompostion.
#' 
#' @details
#' The algorithm is a wrda on the abundance data
#' after transformation to chi-square residuals. 
#' 
#' It is much slower than \code{\link[vegan]{cca}}. The only reason to use
#' it, is that \code{\link{anova.cca0}} does residualized predictor permutation.
#' It is unknown to the authors of \code{douconca} which method
#' \code{\link[vegan]{anova.cca}} implements. See \code{\link{anova.cca0}}.
#'
#' Compared to  \code{\link[vegan]{cca}}, \code{cca0} does not have residual 
#' axes, \emph{i.e.} no CA of the residuals is performed.
#'
#' @returns
#' All scores in the \code{cca0} object are in scaling \code{"sites"} (1): 
#' the scaling with \emph{Focus on Case distances}.
#' 
#' The returned object has class \code{c("cca0" "wrda")} so that
#' the methods \code{print}, \code{predict} and \code{scores}
#' can use the \code{wrda} variant.
#'
#' @references
#' ter Braak C.J.F. and  P. Šmilauer  (2018). Canoco reference manual
#' and user's guide: software for ordination (version 5.1x).
#' Microcomputer Power, Ithaca, USA, 536 pp.
#'
#' Oksanen, J., et al. (2022)
#' vegan: Community Ecology Package. R package version 2.6-8.
#' \url{https://CRAN.R-project.org/package=vegan}.
#'
#' @seealso \code{\link{scores.wrda}}, \code{\link{anova.cca0}},
#' \code{\link{print.wrda}} and \code{\link{predict.wrda}}
#' 
#' @example demo/dune_cca0.R
#' 
#' @export
cca0 <- function(formula, 
                 response = NULL, 
                 data, 
                 traceonly = FALSE,
                 cca_object = NULL,
                 object4QR = NULL) {
  call <- match.call()
  if (!is.null(cca_object)) {
    Yw <- if (inherits(cca_object, "cca0")) cca_object$Ybar else
      # check whether eY works with and without trait covariates 
      vegan::ordiYbar(cca_object, model = "CA")
    transp <- f_transpose(Yw, data, response)
    if (is.na(transp)) {
      warning("cca_object not usable.\n")
      cca_object <- NULL
    }
  }  
  if (is.null(cca_object)) {
    if (is.null(response)) response <- get_response(formula)
	if (is.null(formula) || is_rhs_dot(formula)) data <- sanitize_df(data)
    sumY <- sum(response)
    Yn <- as.matrix(response) / sumY
    R <- rowSums(Yn)
    sWn <- sqrt(R)
    K <- colSums(Yn)
    if (any(R == 0)) {
      stop("Some sites do not have species.\n", names(R)[which(R==0)], "\n")
    }
    if (any(K == 0)) {
      stop("Some species are absent in every site.\n", 
           names(K)[which(K==0)], "\n")
    }
    sWn <- sqrt(R)
    # transform to the unweighted case
    # Yw <-diag(1/sWn)%*% ((Yn - R%*% t(K))%*% diag(1/sqrt(K)))# much much slower!!
    # chi-square residuals
    Yw <- (Yn - R %*% t(K)) *
      matrix(1 / sqrt(K), nrow = length(R), ncol= length(K), byrow = TRUE) / sWn
    total_variance <- sum(Yw ^ 2) #total_inertia
  } else { # use cca_object
    if (transp) {
      Yw <- t(Yw) 
      weights <- if (inherits(cca_object, "cca0"))
        rev(cca_object$weights) else
          list(columns= cca_object$rowsum, rows = cca_object$colsum)
    } else {
      weights <- if (inherits(cca_object, "cca0"))
        cca_object$weights else
          list(columns = cca_object$colsum, rows = cca_object$rowsum)
    }
    K <- weights[[1]]
    R <- weights[[2]]
    sWn <-sqrt(R)
    sumY <- cca_object$sumY
    total_variance <- cca_object$tot.chi
  } # cca_object
  msqr <- msdvif(formula, data = data, weights = R, XZ = TRUE,
                 object4QR = object4QR)
  eY <- qr.resid(msqr$qrZ, Yw)
  Yfit_X <- qr.fitted(msqr$qrXZ, eY)
  ssY_gZ <- sum(eY ^ 2) 
  ssY_XgZ <- sum(Yfit_X ^ 2) 
  if (traceonly) {
    inert <- c(Condition_inertia = total_variance - ssY_gZ, Fit_inertia = ssY_XgZ)
    attr(inert, which = "rank") <-
      c(Cond = msqr$qrZ$rank - 1, Fit = msqr$qrXZ$rank - msqr$qrZ$rank+1)
    return(inert)		 
  }
  svd_Yfit_X <- SVDfull(Yfit_X)
  biplot <- NULL
  eig <- svd_Yfit_X$d ^ 2
  names(eig) <- paste0("CCA", seq_along(eig))
  CCA <- with(svd_Yfit_X, list(eig = eig, poseig = eig, u = u, v = v, 
                               rank = rank, qrank = msqr$qrXZ$rank,
                               tot.chi = ssY_XgZ, QR = msqr$qrXZ, 
                               biplot = biplot, envcentre = NULL, 
                               centroids = NULL))
  rank_CA <- min(dim(Yw)) - 1
  eig_CA <- rep(NA, rank_CA)
  names(eig_CA) <- paste("CA", seq_len(rank_CA), sep = "")
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
      lc_env_scores = (svd_Yfit_X$u %*% diagd) / sWn
    )
  )
  species_axes <- list(species_scores = 
                         list(species_scores_unconstrained = svd_Yfit_X$v / sqrt(K)))
  object <- list(call = call, 
                 method = "cca", 
                 tot.chi = total_variance,
                 formula = formula, 
                 site_axes = site_axes, 
                 species_axes = species_axes, 
                 Nobs = length(R), 
                 eigenvalues = eig,
                 weights = list(columns = K, rows = R),
                 sumY = sumY,
                 data = data,
                 Ybar = Yw,
                 pCCA = pCCA,
                 CCA = CCA,
                 CA = list(tot.chi = ssY_gZ - ssY_XgZ, 
                           rank = rank_CA, 
                           eig = eig_CA),
                 inertia = "scaled Chi-square"
  )
  class(object) <- c("cca0", "wrda")
  return(object)
}

