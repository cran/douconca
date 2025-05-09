#' Print a summary of a wrda or cca0 object
#' 
#' @param x an object from \code{\link{wrda}} or \code{\link{cca0}}
#' @param ...  Other arguments passed to the function (currently ignored).
#' 
#' @returns No return value, results are printed to console.
#' 
#' @example demo/dune_wrda.r
#'
#' @export
print.wrda <- function(x, 
                       ...) {
  print_dcca(x)
}
