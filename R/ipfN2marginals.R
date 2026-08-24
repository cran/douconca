#' Iterative proportional fitting of an abundance table to Hill-N2 marginals
#' 
#' Function for pre-processing/transforming an abundance table
#' by iterative proportional fitting,  
#' so that the transformed table has marginals
#' proportional to \code{N2} or \code{N2(1-N2/N)}
#' with \code{N} the number of elements in the margin.
#' Hill-N2 is the effective number of species. It is of intrinsic interest in
#' weighted averaging (CWM and SNC) as their variance is approximately 
#' inversely proportional to N2 (ter Braak 2019), 
#' and therefore of interest in \code{\link{dc_CA}}.
#'  
#' @param Y abundance table (matrix or dataframe-like), ideally, 
#' with names for rows and columns. 
#' @param updateN2 logical, default \code{TRUE}.
#' If \code{FALSE} the marginal sums are proportional to 
#' the N2-marginals of the initial table, but the N2-marginals of the returned 
#' matrix may not be equal to their marginal sum.
#' If \code{updateN2 = TRUE} and \code{N2N_N2_species=TRUE} (the default), 
#' the column marginals are \code{N2(N-N2)/N} with \code{N} the number of sites.
#' The row sums are then proportional to, what we term, the effective number of
#' informative species.
#' If \code{N2N_N2_species = FALSE},								
#' the returned transformed table has N2 columns marginals, 
#' \emph{i.e.} \code{colSums(Y2) = N2species(Y2)} with \code{Y2} 
#' the return value of \code{ipf2N2}.
#' If converged, N2 row marginals are equal to the row sums, \emph{i.e.} 
#' \code{rowSums(Y2) = approx. const*N2sites(Y2)}and \code{const} a constant.
#' @param N2N_N2_species Set species marginal to
#'  the value of \code{N2(1-N2/N)} for each species.
#' Default \code{TRUE}. If \code{FALSE}, the marginal is set to the \code{N2}
#' value of each species.
#' @param N2N_N2_sites Default \code{FALSE} sets the marginal proportional
#' to the \code{N2} value of each site. 
#' If \code{TRUE}, the marginal is set to \code{N2(1-N2/m)}, with \code{m} the 
#' number of species.				 
#' @param max_iter maximum number of iterative proportional fitting (ipf) 
#' iterations.
#' If \code{max_iter = 0}, the columns are divided by their effective number or
#' informativeness (\code{N2} or \code{N2(1-N2/N)}, depending on the setting
#' of \code{N2N_N2_species}) without
#' further pre-processing and the row sums are then,
#' with \code{N2N_N2_species = TRUE}, sums of
#' informativeness instead of effective number of informative species.
#'
#' @return a matrix of the same order as the input \code{Y},
#' obtained after ipf to N2-marginals.
#' 
#' @details
#' Applying \code{ipf2N2} with \code{N2N_N2_species=FALSE} 
#' to an presence-absence data table returns the same table.
#' However, a species that occurs everywhere (or in most of the sites)
#' is not very informative. This is acknowledged with the default option
#' \code{N2N_N2_species=TRUE}. Then,
#' with \code{N2N_N2_species=TRUE}, species that occur
#' in more than halve the number of sites are down-weighted, so that
#' the row sum is no longer equal to the richness of the site (the number of species),
#' but proportional to the number of informative species. 
#' The returned matrix has the intended species marginal (column sums),
#' by construction of the algorithm, even without convergence. 
#' On convergence, it has the intended site marginal (row sums).
#' 
#' @references
#' ter Braak, C.J.F. (2019). 
#' New robust weighted averaging- and model-based methods
#' for assessing trait-environment relationships. 
#' Methods in Ecology and Evolution, 10 (11), 1962-1971. 
#' \doi{10.1111/2041-210X.13278}
#' 
#' ter Braak, C.J.F. (2026).
#' Fourth-corner latent variable models overstate confidence in
#' trait–environment relationships and what to use instead
#' Environmental and Ecological Statistics.
#'\doi{10.1007/s10651-025-00696-0}
#' 
#' @example demo/dune_ipf2N2.R
#' 
#' @export
ipf2N2 <- function(Y, 
                   max_iter = 10000, 
                   updateN2 = TRUE, 
                   N2N_N2_species = TRUE, 
                   N2N_N2_sites = FALSE) {
  Y <- as.matrix(Y)
  rownams <- rownames(Y)
  colnams <- colnames(Y)
  R <- rowSums(Y)
  K <- colSums(Y)
  if (any(R == 0)) {
    warning("Some sites do not have species.\n ",
            paste(names(R)[which(R == 0)], collapse = " "), "\n")
  }
  if (any(K == 0)) {
    warning("Some species are absent in every site.\n", 
            paste(names(K)[which(K == 0)],  collapse = " "), "\n")
  }
  K[K < .Machine$double.eps] <- .Machine$double.eps
  R[R < .Machine$double.eps] <- .Machine$double.eps
  N2spp <- N2spp0 <- fN2N_N2(Y, 2, N2N_N2 = N2N_N2_species)
  N2sites <- N2sites0 <- fN2N_N2(Y, 1, N2N_N2 = N2N_N2_sites)
  N <- nrow(Y)
  crit0  <- 1.0e10
  crit1 <- crit0 - 1
  iter <- 0
  ratio <- 1
  if(max_iter == 0){
    Y <- Y %*% diag(N2spp / K)
  } else {
    while (crit1 < crit0 && iter < max_iter) { 
      iter <- iter + 1 
      crit0 <- crit1
      Y0 <- Y
      Y <- Y * ((N2sites / R) * ratio)
      K <- colSums(Y)
      if (updateN2) N2spp <- fN2N_N2(Y, 2, N2N_N2 = N2N_N2_species)
      Y <- Y %*% diag(N2spp / K)
      R <- rowSums(Y)
      if (updateN2) N2sites <- fN2N_N2(Y, 1, N2N_N2 = N2N_N2_sites)
      ratio <- sum(N2spp)/sum(N2sites)
      mm <- range((ratio * N2sites+1) / (R + 1))
      crit1 <- max(c(1/mm[1], mm[2]))
    }
    Y <- Y0
    if (updateN2) {
      N2sites <- fN2N_N2(Y, 1, N2N_N2 = N2N_N2_sites)
      N2spp <- fN2N_N2(Y, 2, N2N_N2 = N2N_N2_species)
    }
    if (iter == max_iter && max_iter > 0) {
      warning(paste0("No convergence in ", max_iter, " iterations."),"\n")
    }
    if (any(Y < 0)){
      warning("some values in preprocessed Y negative")
      Y[Y < 0] <- 0 # tiny non-negative values should not occur
    }
  }
  R <- rowSums(Y)
  attr(Y, which = "N2species_original") <- N2spp0
  attr(Y, which = "N2sites_original") <- N2sites0
  attr(Y, which = "N2species") <- N2spp
  attr(Y, which = "N2sites") <- N2sites
  attr(Y, which = "R/N2") <- (R / N2sites) / (sum(R) / sum(N2sites))
  attr(Y, which = "iter") <- iter
  attr(Y, which = "crit") <- crit0
  colnames(Y) <- colnams
  rownames(Y) <- rownams
  if (fN2(R) / length(R) < 0.5) {
    message("Warning: unbalanced site totals in return value:",
            "N2 of row sums less than halve the number of rows.\n")
  }
  if (any(N2sites <= .Machine$double.eps)) {
    warning("After processing: Some sites do not have species.\n ",
            paste(rownames(Y)[which(N2sites <= .Machine$double.eps)], 
                  collapse = " "), "\n")
  }
  if (any(N2spp <= .Machine$double.eps)) {
    warning("After processing: Some species are absent in every site.\n ",
            paste(colnames(Y)[which(N2spp <= .Machine$double.eps)],  
                  collapse = " "), "\n")
  }
  return(Y)
}

#' Hill number of order 2: N2
#' 
#' Calculate Hill number N2.
#' 
#' @param x a numeric vector.
#' @return scalar: Hill number N2.
#' @example demo/dune_ipf2N2.R
#' @references
#' Hill,M.O. (1973).
#' Diversity and evenness: a unifying notation and its consequences.
#' Ecology, 54, 427-432.
#' \doi{10.2307/1934352}.
#' 
#' ter Braak, C.J.F. (2019).
#' New robust weighted averaging- and model-based methods 
#' for assessing trait-environment relationships.
#' Methods in Ecology and Evolution, 10 (11), 1962-1971.
#' \doi{10.1111/2041-210X.13278}.
#' @export
# @noRd
# @keywords internal
fN2 <- function(x) {
  sx <- sum(x)
  if (is.na(sx) || sx < 1.0e-6) {
    return(.Machine$double.eps)
  }
  x <- x / sx		 
  return(1 / sum(x * x))
}

# @noRd
# @keywords internal
fN2N_N2 <- function(Y,
                    margin,
                    N2N_N2 = TRUE) {
  N2 <- apply(X = Y, MARGIN = margin, FUN = fN2)
  margin1 <- if (margin == 1) 2 else 1
  if (N2N_N2) N2 <- N2 * pmax(1 - N2 / dim(Y)[margin1], .Machine$double.eps)
  return(N2)
}