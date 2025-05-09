#' Iterative proportional fitting of an abundance table to Hill-N2 marginals
#' 
#' Function for preprocessing/transforming an abundance table
#' by iterative proportional fitting,  
#' so that the transformed table has marginals \code{N2} or \code{N2(N-N2)}
#' with \code{N} the number of elements in the margin.
#' Hill-N2 is the effective number of species. It is of intrinsic interest in
#' weighted averaging (CWM and SNC) as their variance is approximately 
#' inversely proportional to N2 (ter Braak 2019), 
#' and therefore of interest in \code{\link{dc_CA}}.
#'  
#' @param Y abundance table (matrix or dataframe-like), ideally, 
#' with names for rows and columns. 
#' BEWARE: all rows and columns should have positive sums!
#' @param updateN2 logical, default \code{TRUE}.
#' If \code{FALSE} the marginal sums are proportional to 
#' the N2-marginals of the initial table, but the N2-marginals of the returned 
#' matrix may not be equal to their marginal sum.
#' If \code{updateN2 = TRUE} and \code{N2N_N2_species=TRUE} (the default), 
#' the column marginals are \code{N2(N-N2)/2} with \code{N} the number of sites. 
#' The row sums are then proportional to, what we term, the effective number of
#' informative species.
#' If \code{N2N_N2_species = FALSE},								
#' the returned transformed table has N2 columns marginals, 
#' \emph{i.e.} \code{colSums(Y2) = const*N2species(Y2)} with \code{Y2} 
#' the return value of \code{ipf2N2} and \code{const} a constant.
#' If converged, N2 row marginals are equal to the row sums, \emph{i.e.} 
#' \code{rowSums(Y2) = approx. N2sites(Y2)}.
#' @param N2N_N2_species Set marginals proportional to \code{N2(N-N2)} Default 
#' \code{TRUE}.
#' @param N2N_N2_sites Default \code{FALSE}. Do not change.
#' @param max_iter maximum number of iterative proportional fitting (ipf) 
#' iterations.
#' @param crit stopping criterion.
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
#' 
#' @references
#' ter Braak, C.J.F. (2019). 
#' New robust weighted averaging- and model-based methods
#' for assessing trait-environment relationships. 
#' Methods in Ecology and Evolution, 10 (11), 1962-1971. 
#' \doi{10.1111/2041-210X.13278}
#' 
#' @example demo/dune_ipf2N2.R
#' 
#' @export
ipf2N2 <- function(Y, 
                   max_iter = 100, 
                   updateN2 = TRUE, 
                   N2N_N2_species = TRUE, 
                   N2N_N2_sites = FALSE, 
                   crit = 1.0e-2) {
  if (N2N_N2_sites) {
    warning("Setting N2N_N2_sites to TRUE is not recommended.\n")
  }
  Y <- as.matrix(Y)
  rownams <- rownames(Y)
  colnams <- colnames(Y)
  R <- rowSums(Y)
  K <- colSums(Y)
  if (any(R == 0)) {
    stop("Some sites do not have species.\n", names(R)[which(R == 0)], "\n")
  }
  if (any(K == 0)) {
    stop("Some species are absent in every site.\n", 
         names(K)[which(K == 0)], "\n")
  }
  N2spp <- fN2N_N2(Y, 2, N2N_N2 = N2N_N2_species)
  N2sites <- fN2N_N2(Y, 1, N2N_N2 = N2N_N2_sites)
  N <- nrow(Y)
  tot <- sum(R)
  iter <- 0
  # the criterion is on the weights made relative (what counts in CA)
  while (iter == 0 || 
         (max(N2sites / sum(N2sites) - R / sum(R)) > crit / N && iter < max_iter)) { 
    iter <- iter + 1
    Y <- Y * (N2sites / R)
    K <- colSums(Y)
    if (updateN2) N2spp <- fN2N_N2(Y, 2, N2N_N2 = N2N_N2_species)
    Y <- Y %*% diag(N2spp / K)
    R <- rowSums(Y)
    if (updateN2) N2sites <- fN2N_N2(Y, 1, N2N_N2 = N2N_N2_sites)
  }
  if (iter == max_iter) {
    warning(paste0("No convergence in ", max_iter, " iterations."),"\n")
  }
  Y <- Y * (tot / sum(Y))
  attr(Y, which = "N2species") <- N2spp
  attr(Y, which = "N2sites") <- N2sites
  attr(Y, which = "R/N2") <- (R / N2sites) / (sum(R) / sum(N2sites))
  attr(Y, which = "iter") <- iter
  colnames(Y) <- colnams
  rownames(Y) <- rownams
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
  x <- x / sum(x)
  return(1 / sum(x * x))
}

# @noRd
# @keywords internal
fN2N_N2 <- function(Y,
                    margin,
                    N2N_N2 = TRUE) {
  N2 <- apply(X = Y, MARGIN = margin, FUN = fN2)
  # the reason for the 2 is to strech the number to N/2 when N2 = N/2
  margin1 <- if (margin == 1) 2 else 1
  if (N2N_N2) N2 <- 2 * N2 * (1 - N2 / dim(Y)[margin1])
  return(N2)
}