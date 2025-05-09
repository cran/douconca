#' @title Permutation Test for canonical correspondence analysis
#'
#' @description
#' \code{anova.cca0} performs residual predictor permutation for cca0,
#' which is robust against differences in the 
#' weights (ter Braak & te Beest, 2022). The arguments of the function are 
#' similar to those of \code{\link[vegan]{anova.cca}}, but more restricted.
#
#' @param object an object from \code{\link{dc_CA}}.
#' @param ... unused.
#' @param permutations a list of control values for the permutations as
#' returned by the function \code{\link[permute]{how}}, or the number of 
#' permutations required (default 999), or a permutation matrix where each row
#' gives the permuted indices.
#'
#' @param by character \code{"axis"} which sets the test statistic to the
#' first eigenvalue of the RDA model. Default: \code{NULL} which sets the test
#' statistic to the weighted variance fitted by the predictors
#' (=sum of all constrained eigenvalues). The default is quicker 
#' computationally as it avoids computation of an svd of permuted data sets.
#' @inheritParams anova_species
#' 
#' @details
#' The algorithm is based on published R-code for residual predictor 
#' permutation in canonical correspondence analysis (ter Braak & te Beest, 2022), 
#' but using QR-decomposition instead of ad-hoc least-squares functions.
#' 
#' Note that \code{anova.cca0} is much slower than \code{\link[vegan]{anova.cca}}. 
#' 
#' As \code{\link[vegan]{anova.cca}} is implemented in lower level language,
#' it is difficult to see what it implements.
#' In simulations (with and without \code{Condition()} ) 
#' \code{\link[vegan]{anova.cca}} gives results that are very similar
#' to residual predictor permutation (RPP),
#' but, with \code{by = "axis"}, it can be conservative 
#' and is then less powerful than RPP (as tested with vegan 2.6-8).
#' 
#' @returns
#' A list with two elements with names \code{table} and \code{eigenvalues}.
#' The \code{table} is as from \code{\link[vegan]{anova.cca}} and 
#' \code{eigenvalues} gives the CCA eigenvalues.
#' 
#' @references
#' ter Braak, C.J.F. & te Beest, D.E. (2022).
#' Testing environmental effects on taxonomic composition with canonical
#' correspondence analysis: alternative permutation tests are not equal. 
#' Environmental and Ecological Statistics. 29 (4), 849-868.
#'  \doi{10.1007/s10651-022-00545-4} 
#' 
#' @example demo/dune_cca0.R
#' 
#' @importFrom stats anova
#' @export
anova.cca0 <- function(object, 
                       ...,
                       permutations = 999, 
                       by = c("omnibus", "axis"),
                       n_axes = "all") {
  # permat a matrix of permutations. 
  # If set overrules permuations.
  by <- match.arg(by)
  N <- nrow(object$data) 
  if (inherits(permutations, c("numeric", "how", "matrix"))) {
    if (is.numeric(permutations) && !is.matrix(permutations)) {
      permutations <- permute::how(nperm = permutations[1])
    } else if (is.matrix(permutations) && ncol(permutations) != N) {
      stop("each row of permutations should have", N, "elements.\n")
    }
  } else {
    stop("argument permutations should be integer, matrix or ", 
         "specified by permute::how().\n")
  }
  # Perform a weighted RDAR(M^*~E): an RDA of M^* on the environmental variables
  # using row weights R.
  sWn <- sqrt(object$weights$rows)
  Yw <- object$Ybar
  msqr <- msdvif(object$formula, object$data, object$weights$rows, XZ = FALSE,
                 object4QR = object)
  Zw <- msqr$Zw
  Xw <- msqr$Xw
  Yw <- qr.resid(msqr$qrZ, Yw)
  dfpartial = msqr$qrZ$rank
  # residual predictor permutation.
  out_tes <- list()
  out_tes[[1]]  <- randperm_eX0sqrtw(Yw,Xw, Zw, sWn = sWn,
                                     permutations = permutations,
                                     by = by, n_axes  = n_axes,
                                     return = "all")
  if (by == "axis") {
    while (out_tes[[1]]$rank > length(out_tes)) {
      Zw <- cbind(Zw, out_tes[[length(out_tes)]]$EigVector1)
      out_tes[[length(out_tes) + 1]] <- 
        randperm_eX0sqrtw(Yw,Xw, Zw, sWn = sWn, 
                          permutations = permutations, by = by, return = "all")
    }
  }
  f_sites <- fanovatable(out_tes, Nobs = N, dfpartial = dfpartial, 
                         type = "cca", calltext = c(object$call),
                         n_axes = n_axes)
  result <- list(table = f_sites, eigenvalues = attr(f_sites, "eig"))
  return(result)
}
