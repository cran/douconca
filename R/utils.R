# utility functions for f_trait_axes and f_env_axes in print.dcca and scores.dcca--------------

#' @noRd
#' @keywords internal
wcor <- function(X, 
                 Y = X,
                 w = rep(1, nrow(X))) {
  # weighted correlation between matrix X and Y
  w <- w / sum(w)
  Xstd <- standardize_w(X, w) * w
  Ystd <- standardize_w(Y, w)
  t(Xstd) %*% Ystd
}

#' @noRd
#' @keywords internal
mean_w <- function(X,
                   w = rep(1 / nrow(X), nrow(X))) {
  t(w / sum(w)) %*% X
}

#' @noRd
#' @keywords internal
center_w <- function(X,
                     w = rep(1 / nrow(X), nrow(X))) {
  X - matrix(t(w) %*% X, nrow = nrow(X), ncol= ncol(X), byrow = TRUE)
}

#' @noRd
#' @keywords internal
standardize_w <- function(X, 
                          w = rep(1 / nrow(X), nrow(X))) {
  w<- w / sum(w)						   
  Xc <- X - matrix(t(w) %*% X, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  std <- sqrt(colSums(Xc * Xc * w))
  Xstd <- Xc / matrix(std, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  return(Xstd)
}

#' @noRd
#' @keywords internal					 
mean_sd_w <- function(X,
                      w = rep(1 / nrow(X), nrow(X))){
  # NB requires w to be have sum 1
  mean_w <- t(w) %*% X
  Xc <- X - matrix(mean_w, nrow = nrow(X), ncol = ncol(X), byrow = TRUE) 
  sd_w <- sqrt(colSums(Xc * Xc * w))
  return(list(mean = mean_w, sd = sd_w))
}

#' @noRd
#' @keywords internal
msdvif <- function(formula = NULL, 
                   data, 
                   weights, 
                   XZ = FALSE, 
                   novif = FALSE, 
                   contrast = TRUE,
                   object4QR = NULL) {
  # calc mean variance and vif from for X given Z or XZ 
  # with qr of X|Z or of centered XZ
  # object4QR: cca_object with weighted QR of Z and XZ													  
  if (is.null(formula)) {
    f <- ~. 
  } else {
    f <- formula
  }
  ff <- get_Z_X_XZ_formula(f, data)
  if (XZ) {
    X <- model.matrix(ff$formula_XZ, data = data)[, -1, drop = FALSE]
  } else {
    if (contrast) {
      X <- model.matrix(ff$formula_X1, data = data)[, -1, drop = FALSE] 
    } else {
      X <- modelmatrixI(ff$formula_X1, data = data, XZ = FALSE)
    }
  }
  msd <- mean_sd_w(X, w = weights)
  if (novif) {
    return(list(msd = msd, ff_get = ff))
  }
  avg <- msd$mean
  sds <- msd$sd;
  sWn <- sqrt(weights)
  Zw <- model.matrix(ff$formula_Z, data) * sWn
  qrZ <- if (is.null(object4QR)) qr(Zw) else get_QR(object4QR, model = "pCCA")
  if (is.null(qrZ)) qrZ <- qr(matrix(sWn))
  if (XZ) {
    Xw <- qr.resid(qr(matrix(sWn)), X * sWn) 
  } else { 
    Xw <- qr.resid(qrZ, X * sWn)
  }
  qrX <- if (is.null(object4QR) || !XZ) qr(Xw) else get_QR(object4QR, model = "CCA")
  if (qrX$rank <0.5) message("A variable may have zero variance, giving the error:")
  diagXtX_inv <- diag(chol2inv(qrX$qr, size = qrX$rank))
  diagXtX_inv <- c(diagXtX_inv, rep(NA, length(sds)- length(diagXtX_inv)))
  VIF <- diagXtX_inv * sds ^ 2
  attr(qrX$qr, "scaled:center") <- c(msd$mean)
  attr(qrX$qr, "sd") <- c(msd$sd)
  meansdvif <- as.data.frame(t(rbind(avg, sds, VIF)))
  names(meansdvif) <- c("Avg", "SDS", "VIF")
  rr <- list(meansdvif = meansdvif, qrZ = qrZ, qrX = qrX, Zw = Zw, Xw = Xw)
  if (XZ) {
    names(rr)[c(3, 5)] <- c("qrXZ", "XZw")
  }
  return(rr)
}

#' @noRd
#' @keywords internal
modelmatrixI <- function(formula, 
                         data, 
                         XZ = TRUE){
  # model matrix with full identity contracts = full indicator coding
  ccontrasts <- function(x, data) {
    contrasts(data[[x]], contrasts = FALSE)
  }
  ff <- get_Z_X_XZ_formula(formula, data)
  if (XZ) {
    f <- ff$formula_XZ 
  } else {
    f <- ff$formula_X1
  }
  fcts <- c(ff$Condi_factor, ff$focal_factor)
  if (length(fcts)) {
    cntI <-  lapply(as.list(fcts), ccontrasts, data = data)
    names(cntI) <- fcts
    X <- model.matrix(f, contrasts.arg = cntI, data = data)
  } else {
    X <- model.matrix(ff$formula_XZ, data = data)
  }
  return(X[, -1, drop = FALSE])
}

#' @noRd
#' @keywords internal
is_rhs_dot <- function(f) {
   # For one-sided formulas (~ .), RHS is [[2]]; for two-sided (y ~ .), RHS is [[3]]
  rhs <- if (length(f) >= 3) f[[3]] else f[[2]]
  is.name(rhs) && identical(as.character(rhs), ".")
}

#' @noRd
#' @keywords internal
sanitize_df <- function(df) {
  names(df) <- make.names(names(df))
  for (nm in names(df)) {
    if (is.factor(df[[nm]])) {
      levels(df[[nm]]) <- make.names(levels(df[[nm]]))
    }
  }
  df
}

# from vegan 2.6-4 --------------------------------------------------------

#' @noRd
#' @keywords internal
centroids.cca <-  function(x, 
                           mf,
                           wt) {
  if (is.null(mf) || is.null(x)) {
    return(NULL)
  }
  facts <- sapply(mf, is.factor) | sapply(mf, is.character)
  if (!any(facts)) {
    return(NULL)
  }
  mf <- mf[, facts, drop = FALSE]
  ## Explicitly exclude NA as a level
  mf <- droplevels(mf, exclude = NA)
  if (missing(wt)) {
    wt <- rep(1, nrow(mf))
  }
  ind <- seq_len(nrow(mf))
  workhorse <- function(x, wt) {
    colSums(x * wt) / sum(wt)
  }
  ## As NA not a level, centroids only for non-NA levels of each factor
  tmp <- lapply(mf, function(fct) {
    tapply(ind, fct, function(i) workhorse(x[i, , drop = FALSE], wt[i]))
  })
  tmp <- lapply(tmp, function(z) sapply(z, rbind))
  pnam <- labels(tmp)
  out <- NULL
  if (ncol(x) == 1) {
    nm <- unlist(sapply(pnam,
                        function(nm) paste0(nm, names(tmp[[nm]]))),
                 use.names = FALSE)
    out <- matrix(unlist(tmp), nrow = 1, dimnames = list(NULL, nm))
  } else {
    for (i in seq_along(tmp)) {
      colnames(tmp[[i]]) <- paste0(pnam[i], colnames(tmp[[i]]))
      out <- cbind(out, tmp[[i]])
    }
  }
  out <- t(out)
  colnames(out) <- colnames(x)
  return(out)
}

