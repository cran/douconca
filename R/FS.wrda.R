#' Forward selection of predictor variables using wrda or cca0
#' 
#' @param mod initial wrda or cca0 model with at least on predictor variable, 
#' @param consider character vector of names in \code{mod$data}
#' to consider for addition.
#' 
#' @inheritParams anova.wrda
#' @param n_axes number of eigenvalues to select upon.
#' The sum of \code{n_axes} eigenvalues is taken as criterion.
#' Default \code{"full"} for selection without dimension reduction to 
#' \code{n_axes}.
#' If \code{n_axes =1}, selection is on the first eigenvalue 
#' for selection of variables that form an optimal one-dimensional model.
#' @param test logical; default: \code{TRUE}. 
#' @param threshold_P  significance level, after adjustment for testing 
#' multiplicity, for addition of a variable to the model.
#' @param initial_model character specifying what should be inside 
#' \code{Condition()}. Default: \code{"1"} (nothing, the intercept only).
#' Examples: \code{"region"} for a within-region analysis or
#' \code{"A*B"} for a within analysis specified by the interaction 
#' of factors \code{A} and \code{B}, with \code{region, A, B} 
#' in the data.
#' @param factor2categories logical, default \code{TRUE}, to convert
#' factors to their categories (set(s) of indicator values). 
#' If \code{FALSE}, the selection uses, the fit of a factor
#' divided by its number of categories minus 1.
#' @param PvalAdjustMethod method for correction for multiple testing 
#' in  \code{\link[stats]{p.adjust}}, default \code{"holm"}, which
#' is an improved version Bonferroni.
#' @param max_step maximal number of variables selected. 
#' @param verbose show progress, default: \code{TRUE}.
#' @details
#' 
#' The selection is on the basis of the additional fit (inertia) of a variable 
#' given the variables already in the model.
#' 
#' The names in \code{consider} may include 
#' transformations of predictor variables, such as \code{log(.)},
#' if \code{consider} does not include factors or if \code{factor2categories=FALSE}.
#' If \code{consider} does include factors, such transformations
#' give in a error in the default setting (\code{factor2categories=TRUE}).
#' 
#' 
#' @returns list with three elements: \code{final...} with selected variables
#'  and \code{model_final} and \code{process} with account of the selection process 
#'  If \code{is.numeric(n_axes)}, then the variance in the returned table is
#'  the sum of the n_axes eigenvalues of the current model
#'  (all variables so far included).
#' @example demo/dune_forward_slct.R
#' 
#' @export
FS.wrda <- function(mod,  
                    ...,
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
  if (is.null(consider)) {
    stop("Argument consider must contain variable names.\n")
  }	
  if (any(consider %in% names(mod$data))) {
    if (factor2categories && 
        length(intersect(names(mod$data)[sapply(mod$data, is.factor)], consider)) && 
        !all(consider %in% names(mod$data)))
      stop(paste0("Some variables to consider are not in the data.\n",
                  "Transformations of predictors are not allowed\n", 
                  "when considering factors",
                  " with the default setting of factor2categories.\n",
                  "The variables are:\n  ", paste(consider, collapse = ", "), 
                  "\nThe names in the data are:\n", 
                  paste(names(mod$data), collapse = ", "), "\n"))
  } else {
    stop(paste0("none of the variables to consider are found in the data.\n",
                "The variables are:\n  ", paste(consider, collapse = ", "), 
                "\nThe names in the data are:\n  ", 
                paste(names(mod$data), collapse = ","), "\n"))
  }
  # function for forward selection of predictor variables using wrda.
  FUN <- if (inherits(mod, "cca0")) cca0 else 
    if (inherits(mod, c("rda","wrda"))) wrda
  veganperm <- diff(range(mod$weights$rows)) < 1.0e-4 / mod$Nobs && 
    !is.numeric(n_axes)
  if (!length(consider)) consider <- names(mod$data)
  TF <- sapply(X = mod$data, FUN = is.factor)
  fcts <- intersect(names(mod$data)[TF], consider)
  if (factor2categories && length(fcts)) {
    id <- which(consider %in% fcts) 
    f <- as.formula(paste("~", paste(c(consider[-id], fcts), collapse = "+")))
    mydata <- as.data.frame(modelmatrixI(f, data = mod$data))
    consider <- colnames(mydata)
  } else {
    mydata <- mod$data
  }
  ## step 1
  considerk <- consider
  names(consider) <- consider
  considered <- initial_model
  fit_measures <- matrix(0, nrow = length(consider), ncol = 6)
  colnames(fit_measures) <- c("df", "variance", "cumdf", 
                              "cumvariance", "pval", "pvaladj")
  fit_measuresL <- list()
  attr(n_axes, which = "Names_initial_model") <- 
    colnames(model.matrix(make_formula(1, " ", initial_model),
                          data = mydata))
  traceonly <- !(test || is.numeric(n_axes) || !veganperm) 
  ### initial model
  if (initial_model == "1") {
    dfpartial <- 0
  } else {
    ini_model <- wrda(formula = make_formula("all", "1", initial_model),
                      data= mydata, traceonly = TRUE, cca_object = mod)
    dfpartial <- attr(ini_model, which = "rank")[2]
  }
  fit_measures0 <- fit_measures
  fit_measures0[,"cumdf"] <-dfpartial
  ### start selection
  stepk <- 0
  while(fit_measures[1, 6] < threshold_P && 
        length(considerk) && stepk < max_step) {
    stepk <- stepk + 1
    fit_measures <- matrix(nrow = length(considerk), ncol = 5)
    colnames(fit_measures) <- c("df", "variance", "cumdf", "cumvariance", "pval")
    rownames(fit_measures) <- considerk
    for (k in seq_along(considerk)) {
      out_FS <- try(FUN(formula = make_formula(n_axes, considered,
                                               considerk[k]),
                        data = mydata, traceonly = traceonly,
                        cca_object = mod))
      if (length(out_FS) == 2) {
        df <- attr(out_FS, which = "rank")[2]
        out_FS <- list(out_FS, eigenvalues = NA)
      } else if (!inherits(out_FS, "try-error")) {
        df <- out_FS$CCA$qrank
      }  else {
        df <- NA
      }
      if (inherits(out_FS, "try-error") || is.null(out_FS$eigenvalues)) {
        crit <- 0
        message("This error is taken care of; a variable may have zero variance")
      } else {
        crit <- selection_crit_wrda(out_FS, n_axes)
      }
      if (test && crit > 0) {
        if (is.numeric(n_axes)) {
          out_FS <- try(FUN(formula = make_formula("all", considered, considerk[k]),
                            data = mydata, traceonly = traceonly,
                            cca_object = mod ))
        }
        if (inherits(out_FS, "try-error") || is.null(out_FS$eigenvalues)) {
          pval <- NA 
        } else {
          if (test && out_FS$eigenvalues[1] > sqrt(.Machine$double.eps)) {
            if (veganperm) {
              ff <- change_reponse(make_formula("all", considered, 
                                                considerk[k]), "mod$Ybar")
              environment(ff) <- environment()
              out_FSv <- vegan::rda(ff, data = mydata)
              an <- anova(out_FSv, permutations = permutations)
            } else {
              an <- anova(out_FS, permutations = permutations, 
                          n_axes = n_axes)$table
            }
            pval <- an$`Pr(>F)`[1]
          } else {
            pval <- NA
          }
        }
      } else pval <- NA
      fit_measures[k, ] <- c(df, crit, df, crit, pval)
    }
    pvaladj <- if (test) stats::p.adjust(fit_measures[, "pval"], 
                                         method = PvalAdjustMethod) else 0
    fit_measures <- cbind(fit_measures, pvaladj)
    fit_measures <- setCrit(fit_measures, n_axes, fit_measures0)
    id <- order(fit_measures[,"variance"], decreasing = TRUE)
    fit_measures <- fit_measures[id, , drop = FALSE]
    fit_measures0 <- fit_measures
    if (verbose) print(round(fit_measures, 5))
    fit_measuresL[[stepk]] <- fit_measures
    best_trait <- rownames(fit_measures)[1]
    considered <- c(considered, consider[best_trait])
    if (verbose) {
      print("considered") 
      print(considered)
    }
    considerk <- rownames(fit_measures)[-1]
    fit_measuresL[[stepk]] <- fit_measures
  }
  names(fit_measuresL) <- considered[-1]
  # selected variables with last non-significant value
  tab <- t(sapply(X = fit_measuresL, FUN = function(x) {
    x[1, , drop = FALSE]
  }))
  colnames(tab) <- colnames(fit_measures)
  if (!test) tab <- tab[, -c(5, 6), drop = FALSE]
  if (!is.numeric(n_axes)) tab <- tab[, -c(3, 4), drop = FALSE]
  selected_vars <- if (nrow(tab) == 1 || !test) rownames(tab) else 
    rownames(tab)[-nrow(tab)]
  formula2 <- make_formula(n_axes = "all", initial_model, 
                           paste(selected_vars, collapse = "+"))
  mod2 <- FUN(formula = formula2, 
              data = mydata, 
              traceonly = FALSE,
              cca_object = mod)
  ret <- list(finalWithOneExtra = tab, model_final = mod2, formula = formula2,
              process = fit_measuresL)
  if (nrow(tab) == 1 && max_step > 1) {
    names(ret)[1] <- 
      if (length(consider) > 1) "final_Non_Significant" else "final" 
  }
  return(ret)
}

#' returns the selection_criterion.
#' 
#' @noRd
#' @keywords internal
selection_crit_wrda <- function(out_FS, 
                                n_axes = "full"){
  if (is.character(n_axes)) {
    crit <- if (names(out_FS)[1] == "call") out_FS$CCA$tot.chi else out_FS[[1]][2]
  } else {
    eig <- out_FS$eigenvalues
    n_eig <- length(eig)
    crit <- sum(eig[seq_len(min(n_eig, n_axes))])
  }
  crit <- if (crit < .Machine$double.eps) 0 else crit
  return(crit)
}

#' @param n_axes number of eigenvalues. The sum of \code{n_axes} 
#' eigenvalues is taken  as criterion. Default \code{"full"} for selection 
#' without n_axes reduction.
#' If \code{n_axes =1}, selection is on the first eigenvalue with 
#' all variables in the model.
#' 
#' @param considered variables currently in the model
#' @param considerkk variable_to_consider
#' 
#' @noRd
#' @keywords internal
make_formula <- function(n_axes, 
                         considered, 
                         considerkk) {
  # make formulas of all variables  if(is.numeric(n_axes )) else
  # a formula with Condition(variables currently in the model) + variable_to_consider.
  fconsidered <- paste(considered, collapse = "+")
  if (is.character(n_axes)) {
    formula_FS <-
      paste("~", considerkk, "+ Condition(", fconsidered , ")", collapse = " ")
  } else {
    formula_FS <- paste("~", fconsidered, " + ", considerkk , collapse = " ")
  }
  formula_FS <- as.formula(formula_FS)
  return(formula_FS)
}

#' Default forward selection function.
#' 
#' @param mod A fitted model.
#' @param ... Further arguments passed to other methods.
#' 
#' @return The results from applying forward selection on the fitted model. 
#' 
#' @export
FS <- function(mod, 
               ...) {
  UseMethod("FS")
}

#' @noRd
#' @keywords internal
setCrit <- function(fit_measures, 
                    n_axes,
                    fit_measures0) {
  if (is.numeric(n_axes)) { # adapt variance as it is cumulative with numeric n_axes
    fit_measures00 <- rep(fit_measures0[1,c(4)],each = nrow(fit_measures))
    fit_measures[, 2] <- 
      fit_measures[, 4] - matrix(fit_measures00, nrow(fit_measures))
  }
  # adapt df as it is cumulative in all versions 
  fit_measures00 <- rep(fit_measures0[1, 3], each = nrow(fit_measures))
  fit_measures[, 1] <- 
    fit_measures[, 3] - matrix(fit_measures00, nrow(fit_measures))
  fit_measures[, "variance"] <-  
    fit_measures[, "variance"] / fit_measures[, "df"]
  return(fit_measures)
}
