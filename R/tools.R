######################################################################
## These functions are minor modifications or directly copied from the
## glmnet package:
##        Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via Coordinate
#   Descent.
##        Journal of Statistical Software, 33(1), 1-22.
##        URL http://www.jstatsoft.org/v33/i01/.
## The reason they are copied here is because they are internal functions
## and hence are not exported into the global environment.


#'
#' get coefficients from and "FuncompCGL" object
#'
#' @description
#'          Computes the coefficients at the requested values for \code{lam} from a fitted
#'          \code{\link{FuncompCGL}} object.
#'
#' @param object fitted \code{\link{FuncompCGL}} object.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required. Default is the
#'          entire sequence used to create the model.
#'
#' @param \dots not used.
#'
#' @details
#'          \code{s} is the new vector at which predictions are requested. If \code{s} is not in the
#'          lambda sequence used for fitting the model, the \code{coef} function will use linear
#'          interpolation to make predictions. The new values are interpolated using a fraction of
#'          coefficients from both left and right \code{lam} indices.
#'
#' @return
#' The coefficients at the requested values for \code{s}.
#' @seealso \code{\link{FuncompCGL}}
#' @export
#'


coef.FuncompCGL <- function(object, s = NULL, ...) {
  beta <- object$beta
  if (!is.null(s)) {
    lam <- object$lam
    lamlist <- point.interp(lam, s)
    if(length(s) == 1)
    {
      beta = beta[, lamlist$left, drop=FALSE] *  (1 - lamlist$frac) +
        beta[, lamlist$right, drop=FALSE] * lamlist$frac
    } else {
      beta = beta[, lamlist$left, drop=FALSE] %*% diag(1 - lamlist$frac) +
        beta[, lamlist$right, drop=FALSE] %*% diag(lamlist$frac)
    }
    rownames(seq(s))
  }

  if(is.numeric(object$ref)) {
    # checked dimension after 'apply'
    beta <- apply(beta, 2, function(x, k, p1, ref)
                  c(x[ifelse(ref==1, 0, 1):((ref-1)*k)], -colSums(matrix(x[1:p1], byrow = TRUE, ncol = k)), x[((ref-1)*k+1):(length(x))]),
                  ref = object$ref, k = as.list(object$call)$k, p1 = ncol(object$Z) )
  }

  return(beta)
}


#'
#' get coefficients or make coefficient predictions from and "compCL" object
#'
#' @description
#'          Computes the coefficients at the requested values for \code{lam} from a fitted
#'          \code{\link{compCL}} object.
#'
#' @param object fitted \code{\link{compCL}} object.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required. Default is the
#'          entire sequence used to create the model.
#'
#' @param \dots not used.
#'
#' @details
#'          \code{s} is the new vector at which predictions are requested. If \code{s} is not in the
#'          lambda sequence used for fitting the model, the \code{coef} function will use linear
#'          interpolation to make predictions. The new values are interpolated using a fraction of
#'          coefficients from both left and right \code{lam} indices.
#'
#' @return
#' The coefficients at the requested values for \code{s}.
#' @seealso \code{\link{compCL}}
#' @export
#'
#'

# object <- compCL.object
# s <- 1
# s <- c(1, 0.5, 0.1)
# coef(compCL.object, s = 1)
# coef(compCL.object, s = c(1, 0.5, 0.1))
#coef(compCL.object)


# TO add baseline to coef.compCL
coef.compCL <- function(object, s = NULL, ...) {
  beta <- object$beta
  if (!is.null(s)) {
    lam <- object$lam
    lamlist <- point.interp(lam, s)
    if(length(s) == 1) {
      beta = beta[, lamlist$left, drop=FALSE] *  (1 - lamlist$frac) +
        beta[, lamlist$right, drop=FALSE] * lamlist$frac
    } else {
      beta = beta[, lamlist$left, drop=FALSE] %*% diag(1 - lamlist$frac) +
        beta[, lamlist$right, drop=FALSE] %*% diag(lamlist$frac)
    }
    colnames(beta) <- paste(seq(along = s))
  }
  return(beta)
}



#'
#' make predictions from a "FuncompCGL" object
#'
#' @description
#'  predicts fitted values and class labels from a fitted \code{\link{FuncompCGL}} object.
#'
#' @param object fitted \code{\link{FuncompCGL}} object.
#'
#' @param newx Data frame of new values for functional compositional data which predictions are to be made.
#' @param newZc Matrix of new values for time-invariate covariates which predictions are to be made.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required.
#'         Default is the entire sequence used to create the model.
#'
#' @param \dots Not used. Other arguments to predict.
#' @inheritParams FuncompCGL
#' @details
#'        \code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda
#'        sequence used for fitting the model, the \code{predict} function will use linear interpolation
#'        to make predictions. The new values are interpolated using a fraction of predicted values from
#'        both left and right \code{lam} indices.
#'
#' @return
#' prediction values at the requested values for \code{s}.
#' @seealso \code{\link{FuncompCGL}}
#' @examples
#'
#' df_beta = 5
#' p = 30
#' beta_C_true = matrix(0, nrow = p, ncol = df_beta)
#' beta_C_true[3, ] <- c(-0.8, -0.8 , 0.4 , 1 , 1)
#' beta_C_true[4, ] <- c(0.5, 0.5, -0.6  ,-0.6, -0.6)
#' beta_C_true[1, ] <- c(-0.5, -0.5, -0.5 , -1, -1)
#' beta_C_true[2, ] <- c(0.8, 0.8,  0.7,  0.6,  0.6)
#' Data <- Model(n = 50, p = p, m = 2, intercept = TRUE,
#'               SNR = 2, sigma = 2,
#'               rho_X = 0, rho_T = 0.5, df_beta = df_beta,
#'               n_T = 20, obs_spar = 0.8, theta.add = c(3,4,5),
#'               beta_C = as.vector(t(beta_C_true)))
#' TEST <- Model(n = 100, p = p, m = 2, intercept = TRUE,
#'               SNR = 2, sigma = 2,
#'               rho_X = 0, rho_T = 0.5, df_beta = df_beta,
#'               n_T = 20, obs_spar = 0.8, theta.add = c(3,4,5),
#'               beta_C = as.vector(t(beta_C_true)))
#' y <- Data$data$y
#' X <- Data$data$Comp
#' Zc <- Data$data$Zc
#' intercept <- Data$data$intercept
#'
#' k_use <- df_beta
#' m1 <- FuncompCGL(y = y, X = X , Zc = Zc, intercept = intercept,
#'                  k = k_use, basis_fun = "bs",
#'                  insert = "FALSE", method = "t",
#'                  dfmax = p, tol = 1e-6)
#' predict(m1, newx = TEST$data$Comp, newZc = TEST$data$Zc, intercept = intercept)

#' @export
#' @importFrom methods cbind2
#'


# predict.FuncompCGL <- function(object, newx, s = NULL, ...) {
#   beta <- object$beta

#   ## >>> Data pre-coding: integral <<<
#   #if(nrow(beta) -1 != ncol(newx) )

#   if(is.list(newx)) {
#     ## if newx is the raw data
#     ## newx = list(Comp = Comp, Zc = Zc)
#     ## composition is data.frame of longitudinal composition data and Zc is the matrix of time-invariant covariates
#     list_param = object$call) # list_param = as.list(object$call)
#   } else {
#     ## if new is the pre-coding data
#     ## newx is the cbind(Z, Zc)
#     ## Z is the matrix integral of compositional covariates
#   }
#   ## <<< Data pre-coding: integral >>>

#   ## >>> Take Extraction: lambda and beta <<<
#   if (!is.null(s)) {
#     lam <- object$lam
#     # linear interpolating s in lambda sequence
#     lamlist <- point.interp(lam, s)
#     if(length(s) == 1)
#     {
#       beta = beta[, lamlist$left, drop=FALSE] * (1 - lamlist$frac) +
#         beta[, lamlist$right, drop=FALSE] * lamlist$frac
#     } else {
#       beta = beta[, lamlist$left, drop=FALSE] %*% diag(1 - lamlist$frac) +
#         beta[, lamlist$right, drop=FALSE] %*% diag(lamlist$frac)
#     }
#     #dimnames(beta) <- list(vnames, paste(seq(along = s)))
#   }
#   ## <<< Take Extraction: lambda and beta >>>

#   ## >>> Make prediction: matrix mutiplication <<<
#   if (is.null(dim(newx))) newx = matrix(newx, nrow = 1)
#   fitting <- cbind2(newx, 1) %*% beta #as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
#   ## <<< Make prediction: matrix mutiplication >>>

#   return(fitting)
# }


predict.FuncompCGL <- function(object, newx, newZc = NULL, intercept = TRUE,
                               s = NULL, ref = NULL,
                               degree = 3, basis_fun, insert, method, interval, Trange,
                               T.name, ID.name, ...)
{




  # ## >>> Data pre-coding: integral <<<
  # if(is.list(newx)) {
  #   ## CASE I: newx is the raw data
  #   ## newx: longitudinal composition or categorical data is data.frame
  #   ## newZc: matrix of time-invariant covariates

  #   list_param = object$call # list_param = as.list(object$call)
  # } else {
  #   ## CASE II: newx is the coded integral data
  #   ## newx: matrix integral of compositional covariates
  #   ## newZc: matrix of time-invariant covariates
  # }
  # ## <<< Data pre-coding: integral >>>

  # ## >>> Take Extraction: lambda and beta <<<
  # if (!is.null(s)) {
  #   lam <- object$lam
  #   # linear interpolating s in lambda sequence
  #   lamlist <- point.interp(lam, s)
  #   if(length(s) == 1)
  #   {
  #     beta = beta[, lamlist$left, drop=FALSE] * (1 - lamlist$frac) +
  #       beta[, lamlist$right, drop=FALSE] * lamlist$frac
  #   } else {
  #     beta = beta[, lamlist$left, drop=FALSE] %*% diag(1 - lamlist$frac) +
  #       beta[, lamlist$right, drop=FALSE] %*% diag(lamlist$frac)
  #   }
  #   #dimnames(beta) <- list(vnames, paste(seq(along = s)))
  # }
  # ## <<< Take Extraction: lambda and beta >>>

  # ## >>> Make prediction: matrix mutiplication <<<
  # ------------------------------------------------------------
  # if (is.null(dim(newx))) newx = matrix(newx, nrow = 1)
  # ------------------------------------------------------------
  # fitting <- cbind2(newx, 1) %*% beta #as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
  # ## <<< Make prediction: matrix mutiplication >>>

  # return(fitting)
  newx = as.data.frame(newx)
  if(missing(T.name)) T.name = "TIME"
  if(missing(ID.name)) ID.name = "Subject_ID"
  if( !(T.name %in% names(newx)) || !(ID.name %in% names(newx))) stop("Please provide ID.name and T.name")
  # In case that X is not represented in numerical order of Subject_ID
  newx[, ID.name] <- factor(newx[, ID.name], levels = unique(newx[, ID.name]))
  X.names <- colnames(newx)[!colnames(newx) %in% c(T.name, ID.name)]
  newx[, X.names] = proc.comp(newx[, X.names])
  p = length(X.names)
  pc = ifelse(is.null(dim(newZc)), 0, ncol(newZc))
  beta = object$beta

  if(is.numeric(ref)) {
    if( !(ref %in% 1:p)) stop("Reference variable is out of range")
  }
  k = as.integer((nrow(beta) - 1 - pc) / (p - ifelse(is.numeric(ref), 1, 0)))
  #print(nrow(beta) - 1 - pc)
  #print(p - ifelse(is.numeric(ref), 1, 0))
  #print(k)
  #print(length(unique(newx[, ID.name])))
  list_param = object$call

  # ----------------------------------------------------------------------------

  # >>> integral domain <<<
  if(missing(Trange)) Trange = c(list_param$Trange[1], list_param$Trange[2])
  if(missing(interval)) interval = list_param$interval
  switch(interval,
         "Standard" = {
           #mapping time sequence on to [0,1]
           X[, T.name] = (X[, T.name] - min(Trange)) / diff(Trange) #(X[, T.name] - min(Trange[1])) / diff(Trange)
           interval = c(0, 1)
           Trange = c(0, 1)
         },
         "Original" = {
           interval = Trange
         })

  # <<< integral domain >>>


  # >>> integral interpolation: sseq <<<
  if(missing(insert)) insert = list_param$insert
  if(insert != "FALSE") {
    sseq <- round(seq(from = interval[1], to = interval[2],
                                            #by = length(sseq) * 2,
                                            by = min(diff(sseq))/10), #### digit = 10
                                        10)
  } else {
    sseq <- round(sort(unique(c(Trange, as.vector(X[, T.name])))), 10)
  }
  # <<< integral interpolation: sseq >>>

  # >>> integral basis <<<
  if(missing(basis_fun)) basis_fun = list_param$basis_fun
  nknots <- k - (degree + 1)
    if(nknots > 0) {
      knots <- ((1:nknots) / (nknots + 1)) * diff(interval) + interval[1]
    } else  knots <- NULL
  basis <- switch(basis_fun,
                  "bs" = bs(x = sseq, df = k, degree = degree,
                            Boundary.knots = interval, intercept = TRUE),
                  "fourier" = eval.basis(sseq,
                                         basisobj = create.fourier.basis(rangeval = interval, nbasis = k )),
                  "OBasis" = evaluate(OBasis(expand.knots(c(interval[1], knots, interval[2])),
                                             order = degree + 1),
                                      sseq))

  # <<< integral basis >>>


  # take log-ratio w.r.t. reference
  newx[, X.names] <- apply(newx[, X.names], 2, function(x, ref) x - ref,
                        ref = if(is.null(ref)) 0 else newx[, X.names[ref], drop = TRUE])
  # if ref is null, take the p+1-th component off (nothing);
  # ortherwise take the ref-th component off
  D <- split(newx[, c(T.name, X.names[-ifelse(is.null(ref), p + 1, ref)])], newx[, ID.name])
  n <- length(D)
  #print(length(D))
  p1 <- (p - length(ref)) * k
  #print(p1)
  # Z with (p - length(ref)) * k
  if(missing(method)) method = list_param$method
  newZ <- matrix(NA, nrow = n, ncol = p1)
  # method: integral method
  for(i in 1:n) newZ[i, ] <- ITG(D[[i]], basis, sseq, T.name, interval, insert, method)$integral
  # ----------------------------------------------------------------------------

  newX = cbind(newZ, newZc)
  #print(dim(newX))
  predict.linear(object, newX, s)


}


#'
#' make predictions from a "compCL" object
#'
#' @description
#'  predicts fitted values and class labels from a fitted \code{\link{compCL}} object.
#'
#' @param object fitted \code{\link{compCL}} object.
#'
#' @param newx,newzc newx is the matrix of new values for log composition and newzc is the matrix for confounding
#'                   matrix. Joint of newx and newzc is the design matrix at which predictions are to be made.
## at which predictions are to be made. Must be a matrix.
#'                   Default value for newzc is \code{NULL}.
#'
## @param  matrix for new values for Zc.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required.
#'          Default is the entire sequence used to create the model.
#'
#' @param \dots Not used. Other arguments to predict.
#'
#' @details
#'        \code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda
#'        sequence used for fitting the model, the \code{predict} function will use linear interpolation
#'        to make predictions. The new values are interpolated using a fraction of predicted values from
#'        both left and right \code{lam} indices.
#'
#' @return
#' prediction values at the requested values for \code{s}.
#' @seealso \code{\link{compCL}}
#' @export
#'

# object <- compCL.object
# Znew <- Z[1, ]
# Zcnew = NULL
# s <- 1
# s <- c(1, 0.5, 0.1)
# Znew <- Z[1:5,]
# predict(compCL.object, Znew = Z[1:5, ], s = 1)
# predict(compCL.object, Znew = Z[1:5, ])

# predict.compCL <- function(object, Znew, Zcnew = NULL,  s = NULL, ...) {
#   #print(object$beta)
#   beta <- object$beta
#   nvars <- dim(beta)[1] - 1

#   # >>> predictors <<<
#   ## Znew
#   if( is.null(dim(Znew)) ) Znew = matrix(Znew, nrow = 1)

#   if( any(abs(rowSums(Znew) - 1) > 1e-10) ) {
#     message("Z is transformed into compositional data by deviding rowSums")
#     Znew <- Znew / rowSums(Znew)
#   }
#   if(any(Znew == 0)) stop("There is zero entry in compositional data")
#   Znew <- log(Znew)

#   ## Zcnew
#   if( !is.null(Zcnew) && is.null(dim(Zcnew)) ) Zcnew = matrix(Zcnew, nrow = 1)
#   new <- cbind(Znew, Zcnew)
#   ## <<< predictors >>>

#   if( nvars != dim(new)[2] ) stop("numbers of variables in data is not consistant with estimated coefficients")


#   # >>> beta <<<
#   if (!is.null(s)) {
#     lam <- object$lam
#     lamlist <- point.interp(lam, s)
#     if(length(s) == 1) {
#       beta = beta[, lamlist$left, drop=FALSE] * (1 - lamlist$frac) +
#         beta[, lamlist$right, drop=FALSE] * lamlist$frac
#     } else {
#       beta = beta[, lamlist$left, drop=FALSE] %*% diag(1 - lamlist$frac) +
#         beta[, lamlist$right, drop=FALSE] %*% diag(lamlist$frac)
#     }
#     colnames(beta) <- paste(seq(along = s))
#   }
#   # <<< beta >>>

#   #colnames(beta) <- paste(seq(along = 1:dim(beta)[2]))
#   fitting <- cbind2(new, 1) %*% beta   #as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
#   return(fitting)
# }

predict.compCL <- function(object, newx, newzc = NULL,  s = NULL, ...) {
  newz <- proc.comp(newx)
  if( !is.null(newzc) && is.null(dim(newzc)) ) {
    #if only one new data point
    newzc = matrix(newzc, nrow = 1)
  }
  newX <- cbind(newz, newzc)
  predict.linear(object, newX, s)
}


#'
#' print a "FuncompCGL" object
#'
#' @description
#' Print the nonzero counts for composition varaibles at each lambda along the FuncompCGL path.
#'
#' @param x fitted \code{\link{FuncompCGL}} object.
#'
#' @param digits significant digits in printout.
#'
#' @param \dots Not used. Other arguments to predict.
#'
#' @return a two-column matrix, the first columns is the number of nonzero group counts and
#'         the second column is Lambda.
#' @seealso \code{\link{FuncompCGL}}
#' @export
#'

print.FuncompCGL <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  #print(cbind(Df = x$df, Lam = signif(x$lam, digits)))
  show <- cbind(Df = x$df, Lam = signif(x$lam, digits))
  colnames(show) <- c("DF", "Lam")
  rownames(show) <- seq(nrow(show))
  print(show)
}



#'
#' print a "compCL" object
#'
#' @description
#' Print the nonzero counts for composition varaibles at each lambda along the compCL path.
#'
#' @param x fitted \code{\link{compCL}} object.
#'
#' @param digits significant digits in printout.
#'
#' @param \dots Not used. Other arguments to predict.
#'
#' @return a two-column matrix, the first columns is the number of nonzero group counts and
#'         the second column is Lambda.
#' @seealso \code{\link{compCL}}
#' @export
#'

#print(compCL.object)

print.compCL <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  show <- cbind(Df = x$df, Lam = signif(x$lam, digits))
  colnames(show) <- c("Df", "Lam")
  rownames(show) <- paste0("L", seq(nrow(show)))
  print(show)
  #print(cbind(Df = object$df, Lam = signif(object$lam, digits)))
}


# print.compCL <- function(x, digits = max(3, getOption("digits") - 3), ...) {
#   cat("\nCall: ", deparse(x$call), "\n\n")
#   show <- cbind(Df = x$df, Lam = signif(x$lam, digits))
#   colnames(show) <- c("DF", "Lam")
#   rownames(show) <- seq(nrow(show))
#   print(show)
# }

#'
#' get coefficients from a cross-validation "cv.FuncompCGL" object.
#'
#' This function gets coefficients or makes coefficient predictions from a cross-validated
#' \code{FuncompCGL} model, using the stored \code{"FuncompCGL.fit"} object,
#' and the optimal value chosen for \code{lam} and\code{k} (degree of freedom of basis).
#' Questionable!!!
#'
#' @param object fitted \code{\link{cv.FuncompCGL}} object.
#'
#' @param trim logical, used the trimmed result or not.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required.
#'          Default is the value \code{s="lam.min"} stored on the CV \code{object}, it is the
#'          value of \code{lam} for optimal \code{k} such that gives minimum mean cross-validated error.
#'          If \code{s} is numeric, it is taken as the value(s) of \code{lam} to be used.
#'
#' @param k value of df of basis at which predictions are required. Default is \code{k = "k.min"},
#'          it is the value of \code{k} that gives minimum mean cross-validated error. If \code{k}
#'          is a scaler, it is taken as the value of \code{k} to be used and it must be stored in
#'          \code{\link{cv.FuncompCGL}} object.
#'
#' @param \dots not used.
#'
#' @return The coefficients at the requested values for \code{lam} and \code{k}.
#' @seealso \code{\link{cv.FuncompCGL}}
#' @export
#'


# s <- "lam.1se"
# object <- cv.m1
# trim = FALSE
coef.cv.FuncompCGL <- function(object, trim = FALSE, s = c("lam.min", "lam.1se"), k = NULL, ...) {
  trim <- ifelse(trim, "Ttrim", "Ftrim")

  if (is.numeric(s)) {
    lam_se <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    lam_se <- object[[trim]][[s]][1]
    k_se <- object[[trim]][[s]][2]
  } else stop("Invalid form for s")

  if(is.numeric(k)){
    if(! (k %in% as.integer(names(object$Funcomp.CGL.fit)) )) stop("K is out of range")
    k_se <- k
  }
  k_se <- as.character(k_se)
  cv.fit <- object$Funcomp.CGL.fit[[k_se]]

  #??????????????????????
  # pass arugment from inner function ??????????
  cv.fit$call <- as.list(cv.fit$call)
  cv.fit$call$k <- as.numeric(k_se)
  cv.fit$call$ref <-  as.numeric(cv.fit$ref)
  cv.fit$call
  #??????????????????????
  coef( cv.fit, s = lam_se) #coef.FuncompCGL

}

#'
#' get coefficients from a cross-validation "cv.compCL" object.
#'
#' This function gets coefficients or makes coefficient predictions from a cross-validated
#' \code{compCL} model, using the stored \code{"compCL.fit"} object,
#' and the optimal value chosen for \code{lam}.
#'
#' @param object fitted \code{\link{cv.compCL}} object.
#'
#' @param trim logical, used the trimmed result or not. Default is FASLE.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required.
#'         Default is the value \code{s="lam.min"} stored on the CV \code{object}, it is the
#'         optimal value of \code{lam} that gives minimum.
#'         Alternatively \code{s="lambda.min"} can be used,  it is the largest value of \code{lam}
#'         such that error is within 1 standard error of the minimum.
#'         If \code{s} is numeric, it is taken as the value(s) of \code{lam} to be used.
#'
#' @param \dots not used. Other arguments to predict.
#'
#' @return
#' The coefficients at the requested values for \code{lam}.
#' @seealso \code{\link{cv.compCL}}
#' @export
#'

# object <- cvm
# trim <- FALSE
# s <- 1
# coef(cvm, trim = FALSE, s = "lam.min")
#coef(cvm, trim = TRUE, s = "lam.min")

coef.cv.compCL <- function(object, trim = FALSE, s = c("lam.min", "lam.1se" ),...) {
  trim <- ifelse(trim, "Ttrim", "Ftrim")
  object_use <- object[[trim]]
  if (is.numeric(s)) {
    lam <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    # switch(s,
    #           "lam.min" = "lambda.min",
    #           "lam.1se" = "lambda.1se")
    lam <- object_use[[s]]
  } else {
    stop("Invalid form for s")
  }
  select <- list(beta = coef(object$compCL.fit, s = lam, ...), lam = lam)
  return(select)
}
#
#
# coef.cv.compCL <- function(object, trim = FALSE, s = "lam.min",...) {
#   trim <- ifelse(trim, "Ttrim", "Ftrim")
#
#   cv.fit <- object$compCGL.fit
#
#   if (is.numeric(s)) {
#     lam <- s
#   } else if (is.character(s)) {
#     lam <- object[[trim]]$lam.min
#   } else stop("Invalid form for s")
#
#   coef(cv.fit, s = lam)
#
# }


#' @title
#' make predictions from a "cv.FuncompCGL" object.
#' @description
#' This function makes prediction from a cross-validated \code{FuncompCGL} model,
#' using the stored \code{FuncompCGL.fit} object and the optimal value chose for \code{lambda}.
#'
#' @param Znew Data frame of new values for functional compositional data at which predictions are to be made.
#' @param Zcnew Matrix of new values for time-invariant covariates data at which predictions are to be made.
#' @param ... Other arguments are passed to \code{predict.FuncompCGL}
#' @inheritParams coef.cv.FuncompCGL
#' @seealso \code{\link{FuncompCGL}}
#' @examples
#' df_beta = 5
#' p = 30
#' beta_C_true = matrix(0, nrow = p, ncol = df_beta)
#' beta_C_true[3, ] <- c(-1, 0, 0, 0, -0.5)
#' beta_C_true[1, ] <- c(1, 0, 1 , 0, -0.5)
#' beta_C_true[2, ] <- c(0, 0,  -1,  0,  1)
#'
#' nfolds = 10
#' k_list <- c(4,5,6)
#' n_train = 100
#' n_test = 500
#'
#' Data <- Model(n = n_train, p = p, m = 0, intercept = TRUE,
#'               SNR = 2, sigma = 2,
#'               rho_X = 0, rho_T = 0.5,
#'               Corr_X = "CorrCS",
#'               df_beta = df_beta,
#'               n_T = 20, obs_spar = 1, theta.add = FALSE, #c(0,0,0),
#'               beta_C = as.vector(t(beta_C_true)))
#' y <- drop(Data$data$y)
#' n <- length(y)
#' X <- Data$data$Comp
#' Zc <- Data$data$Zc
#' intercept <- Data$data$intercept
#' m <- ifelse(is.null(Zc), 0, dim(Zc)[2])
#' arg_list <- as.list(Data$call)[-1]
#' arg_list$n <- n_test
#' Test <- do.call(Model, arg_list)
#' rule <- "lam.min"
#' # cv_cgl, Constrained group lasso
#' cv_cgl <-  cv.FuncompCGL(y = y, X = X, Zc = Zc, intercept = intercept,
#'                          W = rep(1, p), #W = function(x){ diag( 1 / apply(x, 2, sd) ) },
#'                          k = k_list, trim = 0,
#'                          nfolds = 10,
#'                          tol = 0, inner_eps = 1e-6, inner_maxiter = 1E3,
#'                          dfmax = 30, lambda.factor = 1e-3,
#'                          mu_ratio = 1, outer_eps = 1e-6,
#'                          keep = TRUE, Trange = c(0,1))
#' y_hat = predict(cv_cgl, Znew = Test$data$Comp, Zcnew = Test$data$Zc)
#'
#' @export

predict.cv.FuncompCGL <- function(object, Znew, Zcnew = NULL,
                                  s = c("lam.1se", "lam.min"), k = NULL,
                                  trim = FALSE, ...) {


  trim <- ifelse(trim, "Ttrim", "Ftrim")

  if (is.numeric(s)) {
    lam_se <- s
    if(is.null(k)) {
      stop("Need to provide df k")
    } else {
      if(is.numeric(k)){
        if(! (k %in% as.integer(names(object$Funcomp.CGL.fit)) )) stop("K is out of range")
        k_se <- k
  }
      k_se <- as.character(k_se)
    }
  } else if (is.character(s)) {
    s <- match.arg(s)
    lam_se <- object[[trim]][[s]][1]
    k_se <- object[[trim]][[s]][2]
  } else stop("Invalid form for s")

  predict(object$Funcomp.CGL.fit[[ which(names(object$Funcomp.CGL.fit) == k_se)]], Znew, Zcnew, s = lam_se, ...)

}



#' @title
#' make predictions from a \code{"cv.compCL"} object.
#' @description
#' This function makes prediction from a cross-validated \code{compCL} model,
#' using the stored \code{compCL.fit} object and the optimal value chose for \code{lambda}.
#'
#' @param object fitted \code{"cv.compCL"} model
#' @param Znew new compositional data
#' @param Zcnew new time-invariant covarites
#' @param s \code{"lam.min"} or \code{"lam.1se"}
#' @param trim \code{"FALSE"}
#' @param \dots not used. Other arguments to predict.
#'
# @inheritParams
#' @seealso \code{\link{cv.compCL}}
#' @examples
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c( beta, rep(0, times = p - length(beta)) )
#' Comp_data = comp_Model(n = n, p = p,
#'                        rho = 0.2, sigma = 0.5,
#'                        gamma  = 0.5, add.on = 1:5,
#'                        beta = beta, intercept = FALSE)
#' test_data = comp_Model(n = 100, p = p,
#'                        rho = 0.2, sigma = 0.5,
#'                        gamma  = 0.5, add.on = 1:5,
#'                        beta = beta, intercept = FALSE)
#' cvm <- cv.compCL(y = Comp_data$y,
#'                  Z = Comp_data$X.comp, Zc = Comp_data$Zc,
#'                  intercept = Comp_data$intercept,
#'                  lam = NULL, nfolds = 10, trim = 0.05, lambda.factor = 0.0001,
#'                  dfmax = p, mu_ratio = 1, outer_eps = 1e-10, inner_eps = 1e-8, inner_maxiter = 1e4)
#' y_hat = predict(cvm, Znew = test_data$X.comp, Zcnew = test_data$Zc)
#' plot(test_data$y, y_hat)
#' abline(a = 0, b = 1, col = "red")
#'
#' @export

predict.cv.compCL <- function(object, Znew, Zcnew = NULL, s = c("lam.min", "lam.1se" ), trim = FALSE, ...){
 if (is.numeric(s)) {
    lam <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    # switch(s,
    #           "lam.min" = "lambda.min",
    #           "lam.1se" = "lambda.1se")
    trim <- ifelse(trim, "Ttrim", "Ftrim")
    lam_use <- object[[trim]]
    lam <- lam_use[[s]]
  } else {
    stop("Invalid form for s")
  }
  predict(object$compCL.fit, Znew, Zcnew, s = lam, ...)
}



#' @title
#' plot coefficients from a \code{"FuncompCGL"} object
#'
#' @description
#' Produces a coefficient profile plot of the coefficient paths for a fitted
#' \code{"FuncompCGL"} object.
#'
#' @param x fitted \code{"FuncompCGL"} model.
#' @param p numbers of compositional variables.
#' @param k degree of freedom of basis used.
#' @param ylab What is on the Y-axis, L1-norm or L2-norm (default) of each group of coefficients.
#' @param xlab Either plot against \code{log(lambda)} (default) or \code{lambda}.
#' @param \dots Other graphical parameters to plot.
#'
#' @details
#' A plot is prodoced, and nothing is returned.
#' @seealso \code{\link{FuncompCGL}}
#' @export
#'


plot.FuncompCGL <- function(x, p, k,
                            ylab = c("L2", "L1"),
                            xlab = c("log", "lambda"),
                            ...) {
  obj <- x
  ylab = match.arg(ylab)
  xlab = match.arg(xlab)
  if( is.numeric(obj$ref) && !is.numeric(as.list(obj$call)$k) ) {
    obj$call <- as.list(obj$call)
    obj$call$k <- k
  }
  beta <- coef(obj) #obj$beta
  include_which <- sort(unique(unlist(apply(beta, 2, Nzero, p = p, k = k))))
  group <- matrix(1:(p*k), nrow = k)
  beta <- beta[group[, include_which], ]
  #cat("1\r\n")
  switch(ylab,
         "L1" = {
           yvalue = apply(beta, 2, function(x, k) {
             value = vector()
             for(l in 1:(length(x)/k)) value[l] <- sum(abs(x[(l*k-k+1):(l*k)]))
             return(value)
           }, k = k)
           ylab = "L1 norm"},
         "L2" = {
           yvalue = apply(beta, 2, function(x, k) {
             value = vector()
             for(l in 1:(length(x)/k)) value[l] <- sqrt(sum((x[(l*k-k+1):(l*k)])^2))
             return(value)
           }, k = k)
           ylab = "L2 norm"
         }
  )

  switch(xlab,
         "lambda" = {
           xlab = "Lambda"
           xvalue = drop(obj$lam)
         },
         "log" = {
           xlab = "Log(Lambda)"
           xvalue = log(drop(obj$lam))
         })


  dotlist=list(...) #list(...)
  type=dotlist$type
  if(is.null(type))
    matplot(xvalue,t(yvalue),lty=1,xlab=xlab,ylab=ylab,type="l", ylim = c(-0.1, max(yvalue)), ...
    )  else matplot(xvalue,t(yvalue),lty=1,xlab=xlab,ylab=ylab
                    , ylim = c(-0.1, max(yvalue)),...
    )

  #text(x = min(xvalue), y = yvalue[,dim(yvalue)[2]], labels = paste(include_which))
  # cvraw <- (drop(y) - cv.m1$fit.preval[, , 1])^2
  # N <- length(y) - apply(is.na(cv.m1$fit.preval[, , 1]), 2, sum)
  # cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  # cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))

  at.txt = apply(yvalue, 1, function(x) which(x>0)[1])
  at.index <- rbind(at.txt, sort(at.txt))
  #matplot(xvalue,t(yvalue),lty=1,xlab=xlab,ylab=ylab,type="l", ylim = c(-0.4, max(yvalue)))
  text(x = xvalue[at.txt[order(at.txt)[-((1:floor(length(at.txt)/2))*2)]]],
       y = -0.05, labels = paste(include_which[order(at.txt)][-((1:floor(length(at.txt)/2))*2)]),
       col = rep(c("black", "grey"), length.out = length(include_which) - length((1:floor(length(at.txt)/2))*2)) )
  text(x = xvalue[at.txt[order(at.txt)[((1:floor(length(at.txt)/2))*2)]]],
       y = -0.1, labels = paste(include_which[order(at.txt)][(1:floor(length(at.txt)/2))*2]),
       col = rep(c("black", "grey"), length.out = length((1:floor(length(at.txt)/2))*2)) )



}











#' @title
#' plot coefficients from a \code{"compCL"} object
#'
#' @description
#' Produces a coefficient profile plot of the coefficient paths for a fitted
#' \code{"compCL"} object.
#'
#' @param x fitted \code{"compCL"} model.
#' @param xlab What is on the X-axis. "\code{lambda} against \code{log(lambda)} (default)
#' or "\code{norm}" against \code{L1-norm} .
#' @param label if TRUE, lable the curve with variable sequence numbers.
#' @param \dots Other graphical parameters to plot.
#'
#'
#' @details
#' A plot is prodoced, and nothing is returned.
#' @seealso \code{\link{compCL}}
#' @export
#'

plot.compCL <- function(x, xlab=c("norm", "log"), label=FALSE, ...) {
  xlab <-  match.arg(xlab)
  beta <- x$beta
  lambda <- x$lam
  df <- x$df
  nr = nrow(beta)
  if(nr == 1) {
    index_which <- ifelse(any(abs(beta) > 0), 1, NULL)
  } else {
    # beta_index <- abs(beta) > 0
    # index_which <- seq(nr)
    # ones <- rep(1, ncol(beta))
    # nz = as.vector((beta %*% ones) > 0)
    # index_which = index_which[nz]
    index_which <- apply(beta, 1, function(x) ifelse(any(abs(x) > 0), TRUE, FALSE))
    index_which <- seq(nrow(beta))[index_which]
  }

  nwhich = length(index_which)

  if(nwhich == 0) {
    warning("No plot produced since all coefficients are zero")
    return()
  }

  beta = as.matrix(beta[index_which, , drop=FALSE])
  name_ylab = "Coefficients"
  switch(xlab,
         "norm" = {
           index_xlab = apply(abs(beta), 2, sum)
           name_xlab = "L1 Norm"
           approx.f = 1
         },
         "log" = {
           index_xlab = log(lambda)
           name_xlab = "log Lambda"
           approx.f = 0
         })
  dotlist = list(...)
  type = dotlist$type
  if(is.null(type)) {
    matplot(index_xlab, t(beta), lty=1, xlab = name_xlab, ylab = name_ylab, type = "l", ...)
  } else {
    matplot(index_xlab, t(beta), lty=1, xlab = name_xlab, ylab = name_ylab, ...)
  }

  atdf=pretty(index_xlab)
  ### compute df by interpolating to df at next smaller lambda
  ### thanks to R package glment and Yunyang Qian
  prettydf=approx(x=index_xlab,y=df,xout=atdf,rule=2,method="constant",f=approx.f)$y
  axis(3,at=atdf,labels=prettydf,tcl=NA)
  if(label){
    switch(xlab,
           "norm" = {
             xpos = max(index_xlab)
             pos = 4
           },
           "log" = {
             xpos = min(index_xlab)
             pos = 2
           })
    pos_xlab = rep(xpos, nwhich)
    pos_ylab = beta[, ncol(beta)]
    text(pos_xlab, pos_ylab, paste(index_which), cex = 0.8, pos = pos)
  }
}


























#' @title
#' plot the cross-validation curve produced by cv.FuncompCGL
#'
#' @description
#' Plots the cross-validation curve for different degree of freedom \code{k},
#' and upper and lower standard deviation curves,
#' as a function of the \code{lambda} values used.
#'
#' @param x fitted \code{"cv.FuncompCGL"}.
#' @param xlab Either plot against \code{log(lambda)} (default) or \code{lambda}.
#' @param trim Logic value, determining use trimmed result stored in \code{cv.FumcompCGL} or not. Default
#' @param k_list a vector, determining result of which \code{k} (degree freedom) to plot.
#'               If not provided, cross-validation curve for k that associate with \code{lambda.min} (blue)
#'               and \code{lambda.1se} (red) are plotted.
#' @param \dots Other graphical parameters to plot
#'
#' @details
#' A plot is prodoced, and nothing is returned.
#' @seealso \code{\link{cv.FuncompCGL}}
#' @export
#' @import graphics

plot.cv.FuncompCGL <- function(x, xlab = c("log", "lambda"),
                               trim = FALSE,
                               k_list,
                               ...) {
  #cat("1 \r\n")
  cvobj = x
  xlab <- match.arg(xlab)
  trim <- ifelse(trim, "Ttrim", "Ftrim")
  cvobj_use <- cvobj[[trim]]
  #xlab = "log"
  switch(xlab,
         "lambda" = {
           xlab = "Lambda"
           xvalue = drop(cvobj$lam)
         },
         "log" = {
           xlab = "Log(Lambda)"
           xvalue = log(drop(cvobj$lam))
         })
  #cat(xlab, "\r\n")
  k_listall <- as.numeric(names(cvobj$Funcomp.CGL.fit))
  rule_min <- cvobj_use$lam.min
  rule_1se <- cvobj_use$lam.1se
  if(missing(k_list)) k_list <- c(rule_min['df'], rule_1se['df'])
  if(is.numeric(k_list)) k_list <- k_list[k_list %in% k_listall]
  if(is.character(k_list)) k_list <- switch(k_list,
                                            "lam.1se" = rule_1se['df'],
                                            "lam.min" = rule_min['df'])

  #cat(k_list, "\r\n")
  N_list <- apply(cvobj_use$cvm[k_list- k_listall[1] + 1, ], 1, function(x) length(x[!is.na(x)]))
  plot.args = list(x = xvalue[seq(max(N_list))], xlab = xlab,
                   y = cvobj_use$cvm[k_list[which.max(N_list)] - k_listall[1] + 1, seq(max(N_list))],
                   #y = cvobj_use$cvm[k_list[1], ],
                   ylab = "MEAN-Squared Error", ##
                   type = "n",
                   ylim = range(cvobj_use$cvup[k_list- k_listall[1] + 1, ], cvobj_use$cvlo[k_list- k_listall[1] + 1, ], na.rm = TRUE)
  )
  #cat(plot.args[[1]], "\r\n")
  new.args = list(...) #list(...)
  #cat("2\r\n")
  if(length(new.args)) plot.args[names(new.args)] = new.args
  do.call("plot", plot.args)


  for(l in 1:length(k_list)) {
    if(k_list[l] != rule_min['df'] && k_list[l] != rule_1se['df']) {
      # ll <- N_list[l]
      # error.bars(xvalue[seq(ll)],
      #            cvobj_use$cvup[k_list[l] - k_listall[1] + 1, seq(ll)],
      #            cvobj_use$cvlo[k_list[l] - k_listall[1] + 1, seq(ll)],
      #            width=0.01,col="lightgrey")
      #

      error.bars(xvalue,
                 cvobj_use$cvup[k_list[l] - k_listall[1] + 1, ],
                 cvobj_use$cvlo[k_list[l] - k_listall[1] + 1, ],
                 width=0.01,col="lightgrey")
    }
  }

 if(rule_1se['df'] %in% k_list){
    # ll <- N_list[which(k_list == rule_1se['df'])[1]]
    # error.bars(xvalue[seq(ll)],
    #            cvobj_use$cvup[rule_1se['df'] - k_listall[1] + 1, seq(ll)],
    #            cvobj_use$cvlo[rule_1se['df'] - k_listall[1] + 1, seq(ll)],
    #            width=0.01,col="blue")


    error.bars(xvalue,
               cvobj_use$cvup[rule_1se['df'] - k_listall[1] + 1, ],
               cvobj_use$cvlo[rule_1se['df'] - k_listall[1] + 1, ],
               width=0.01,col="blue")
  }

  if(rule_min['df'] %in% k_list) {
    # ll <- N_list[which(k_list == rule_min['df'])[1]]
    # error.bars(xvalue[seq(ll)],
    #            cvobj_use$cvup[rule_min['df'] - k_listall[1] + 1, seq(ll)],
    #            cvobj_use$cvlo[rule_min['df'] - k_listall[1] + 1, seq(ll)],
    #            width=0.01,col="red")


    error.bars(xvalue,
               cvobj_use$cvup[rule_min['df'] - k_listall[1] + 1, ],
               cvobj_use$cvlo[rule_min['df'] - k_listall[1] + 1, ],
               width=0.01,col="red")
  }





  for(l in 1:length(k_list)) {
    if(k_list[l] != rule_min['df'] && k_list[l] != rule_1se['df']) {
      points(xvalue,
             cvobj_use$cvm[k_list[l] - k_listall[1] + 1, ],
             pch=18,
             col="limegreen")
    }
  }

  if(rule_1se['df'] %in% k_list){
    points(xvalue,
           cvobj$Ftrim$cvm[rule_1se['df'] - k_listall[1] + 1, ],
           pch=20,
           col="blue")
    axis(side=3,at=xvalue,
         labels=paste(cvobj$Funcomp.CGL.fit[[as.character(rule_1se['df'])]]$df[1:length(xvalue)]),
         tick=FALSE,
         line=0,
         pos = plot.args$ylim[2],
         col.axis = "blue")
  }

  if(rule_min['df'] %in% k_list) {
    points(xvalue,
           cvobj$Ftrim$cvm[rule_min['df'] - k_listall[1] + 1, ],
           pch=20,
           col="red")
    axis(side=3,at=xvalue,
         labels=paste(cvobj$Funcomp.CGL.fit[[as.character(rule_min['df'])]]$df[1:length(xvalue)]),
         tick=FALSE,
         line=0,
         col.axis = "red")

  }


  abline(v = switch(xlab,
                    "Lambda" = rule_1se["lam"],
                    "Log(Lambda)" = log(rule_1se["lam"])),
         lty = 3, col = "blue"
  )

  abline(v = switch(xlab,
                    "Lambda" = rule_min["lam"],
                    "Log(Lambda)" = log(rule_min["lam"])),
         lty = 3, col = "red"
  )



}



#' @title
#' plot the cross-validation curve produced by "cv.compCL" object.
#'
#' @description
#' Plots the cross-validation curve, and upper and lower standard deviation
#' curves, as a function of the \code{lambda} values used.
#'
#' @param x fitted \code{"cv.compCL"} object.
#' @param xlab Either plot against \code{log(lambda)} (default) or \code{lambda} or \code{-log(lambda)}.
#' @param \dots Other graphical parameters to plot.
#'
#' @details
#' A plot is produced, and nothing is returned.
#' @seealso \code{\link{cv.compCL}}
#' @export
#'


plot.cv.compCL <- function(x, xlab = c("log", "-log", "lambda"), ...) {
  cvobj <- x
  xlab <- match.arg(xlab)
  #trim <- ifelse(trim, "Ttrim", "Ftrim")
  cvobj_use <- cvobj[["Ftrim"]]
  switch(xlab,
        "log" = {
           xlab = "Log(Lambda)"
           xvalue = log(drop(cvobj$lam))
         },
         "-log" = {
           xlab = "-Log(Lambda)"
           xvalue = -log(drop(cvobj$lam))
         },
         "lambda" = {
           xlab = "Lambda"
           xvalue = drop(cvobj$lam)
         })


  plot.args = list(x = xvalue,y = cvobj_use$cvm,
                   ylim = range(cvobj_use$cvupper,cvobj_use$cvlo),
                   xlab=xlab,
                   ylab="MEAN-Squared Error",
                   type="n")
  new.args=list(...)
  if(length(new.args)) plot.args[names(new.args)]=new.args
  do.call("plot",plot.args)
  error.bars(xvalue, cvobj_use$cvupper, cvobj_use$cvlo, width=0.01, col="darkgrey")
  points(xvalue,cvobj_use$cvm,pch=20,col="limegreen")
  axis(side=3,at=xvalue,labels=paste(cvobj$compCL.fit$df),tick=FALSE,line=0)
  abline(v = switch(xlab,
                   "Log(Lambda)" = log(cvobj_use$lam.1se),
                   "-Log(Lambda)" = -log(cvobj_use$lam.1se),
                   "Lambda" = cvobj_use$lam.1se),
        lty=3, col = "blue")
  abline(v=switch(xlab,
                 "Log(Lambda)" = log(cvobj_use$lam.min),
                 "-Log(Lambda)" = -log(cvobj_use$lam.min),
                 "Lambda" = cvobj_use$lam.min),
        lty=3, col = "red")

  #invisible()
}



#'
#' get coefficients from a GIC "GIC.FuncompCGL" object.
#'
#' This function gets coefficients or makes coefficient predictions from a GIC
#' \code{FuncompCGL} model, using the stored \code{"FuncompCGL.fit"} object,
#' and the optimal value chosen for \code{lam} and\code{k} (degree of freedom of basis).
#'
#'
#' @param object fitted \code{\link{GIC.FuncompCGL}} object.
#'
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required.
#'          Default is the value \code{s="lam.min"} stored on the \code{GIC.FuncompCGL} \code{object}, it is the
#'          value of \code{lam} for optimal \code{k} such that gives minimum mean GIC value.
#'          If \code{s} is numeric, it is taken as the value(s) of \code{lam} to be used.
#'
#' @param k value of df of basis at which predictions are required. Default is the value of \code{k}
#'          that gives minimum mean GIC value. If \code{k}
#'          is a scaler, it is taken as the value of \code{k} to be used and it must be stored in
#'          \code{\link{GIC.FuncompCGL}} object.
#'
#' @param \dots not used.
#'
#' @return The coefficients at the requested values for \code{lam} and \code{k}.
#'
#' @export
#'


# s <- "lam.1se"
# object <- cv.m1
# trim = FALSE
coef.GIC.FuncompCGL <- function(object, s ="lam.min", k = NULL, ...) {

  if (is.numeric(s)) {
    lam_se <- s
  } else if (is.character(s)) {
    lam_se <- object$lammin['lam']
    k_se <- object$lammin['df']
  } else stop("Invalid form for s")

  if(is.numeric(k)){
    if(! (k %in% as.integer(names(object$Funcomp.CGL.fit)) )) stop("K is out of range")
    k_se <- k
  }
  k_se <- as.character(k_se)
  cv.fit <- object$Funcomp.CGL.fit[[k_se]]

  #??????????????????????
  # pass arugment from inner function ??????????
  cv.fit$call <- as.list(cv.fit$call)
  cv.fit$call$k <- as.numeric(k_se)
  cv.fit$call$ref <-  as.numeric(cv.fit$ref)
  cv.fit$call
  #??????????????????????
  coef( cv.fit, s = lam_se) #coef.FuncompCGL

}


#'
#' get coefficients or make coefficient predictions from a "GIC.compCL" object.
#'
#' This function gets coefficients or makes coefficient predictions from a regulaized fitting
#' \code{compCL} model, using the stored \code{"compCL.fit"} object,
#' and the optimal value chosen for \code{lam}.
#'
#' @param object fitted \code{\link{GIC.compCL}} object.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required.
#'         Default is the value \code{s="lam.min"} stored on the CV \code{object}, it is the
#'         optimal value of \code{lam} that gives minimum.
#'         Alternatively \code{s="lambda.min"} can be used,  it is the largest value of \code{lam}
#'         such that error is within 1 standard error of the minimum.
#'         If \code{s} is numeric, it is taken as the value(s) of \code{lam} to be used.
#'
#' @param \dots not used. Other arguments to predict.
#'
#' @return
#' The coefficients at the requested values for \code{lam}.
#' @seealso \code{\link{GIC.compCL}}
#' @export
#'

coef.GIC.compCL <- function(object, s = "lam.min",...) {

  if (is.numeric(s)) {
    lam <- s
  } else if (s == "lam.min") {
    lam <- object[[s]]
  } else {
    stop("Invalid form for s")
  }
  select <- list(beta = coef(object$compCL.fit, s = lam, ...), lam = lam)
  return(select)
}

#' @title
#' make predictions from a \code{"GIC.FuncompCGL"} object.
#' @description
#' This function makes prediction from a GIC \code{FuncompCGL} model,
#' using the stored \code{Funcomp.CGL.fit} object and the optimal value chose for \code{lambda}.
#'
#' @param object fitted \code{"GIC.FuncompCGL"} model
#' @param Znew new time-variate compositional data, stored as data frame
#' @param Zcnew new time-invariant covarites, stored as matrix
#' @param s \code{"lam.min"} or user provided value
#' @param k value of df of basis at which predictions are required. Default is the value of \code{k}
#'          that gives minimum mean GIC value. If \code{k}
#'          is a scaler, it is taken as the value of \code{k} to be used and it must be stored in
#'          \code{\link{GIC.FuncompCGL}} object.
#' @param \dots Other arguments to predict.
#'
#' @seealso \code{\link{GIC.FuncompCGL}}
#' @examples
#' df_beta = 5
#' p = 30 #30
#' beta_C_true = matrix(0, nrow = p, ncol = df_beta)
#' beta_C_true[3, ] <- c(-1, 0, 0, 0, -0.5) #c(-0.5, 0, 0, 0, -0.5)
#' beta_C_true[1, ] <- c(1, 0, 1 , 0, -0.5) #c(0.5, 0, 1 , 0, -0.5)
#' beta_C_true[2, ] <- c(0, 0,  -1,  0,  1)
#'
#' nfolds = 10
#' k_list <- c(4,5,6)
#' n_train = 100
#' n_test = 500
#'
#' Data <- Model(n = n_train, p = p, m = 0, intercept = TRUE,
#'               SNR = 3, sigma = 3,
#'               rho_X = 0, rho_T = 0,
#'               Corr_X = "CorrCS", Corr_T = "CorrAR",
#'               df_beta = df_beta,
#'               n_T = 20, obs_spar = 1, theta.add = FALSE, #c(0,0,0),
#'               beta_C = as.vector(t(beta_C_true)))
#' y <- drop(Data$data$y)
#' n <- length(y)
#' X <- Data$data$Comp
#' Zc <- Data$data$Zc
#' intercept <- Data$data$intercept
#' m <- ifelse(is.null(Zc), 0, dim(Zc)[2]) #+ as.integer(intercept)
#' #m1 <- m + as.integer(intercept)
#' sseq <- Data$basis.info[,1]
#' beta_C.true <- matrix(Data$beta[1:(p*(df_beta))],
#'                       nrow = p, ncol = df_beta, byrow = TRUE)
#' beta_curve.true <- Data$basis.info[,-1] %*% t(beta_C.true)
#' Non_zero.true <- (1:p)[apply(beta_C.true, 1, function(x) max(abs(x)) > 0)]
#' foldid <- sample(rep(seq(nfolds), length = n))
#'
#' arg_list <- as.list(Data$call)[-1]
#' arg_list$n <- n_test
#' Test <- do.call(Model, arg_list)
#' y_test <- drop(Test$data$y)
#'
#'
#'
#' GIC_m1 <- GIC.FuncompCGL( y = y, X = X, Zc = Zc, ref = NULL,
#'                           inner_eps = 1e-8, outer_eps = 1e-8, tol = 1e-8,
#'                           k = k_list)
#' y_hat = predict(GIC_m1, Znew = Test$data$Comp, Zcnew = Test$data$Zc)
#'
#' @export
#'

predict.GIC.FuncompCGL <- function(object, Znew, Zcnew = NULL,
                                   s = "lam.min", k = NULL, ...) {



  if (is.numeric(s)) {
    lam_se <- s
    if(is.null(k)) {
      stop("Need to provide df k")
    } else {
      if(is.numeric(k)){
        if(! (k %in% as.integer(names(object$Funcomp.CGL.fit)) )) stop("K is out of range")
        k_se <- k
      }
      k_se <- as.character(k_se)
    }
  } else if (s == "lam.min") {
    lam_se <- object$lammin[1]
    k_se <- object$lammin[2]
  } else stop("Invalid form for s")

  predict(object$Funcomp.CGL.fit[[ which(names(object$Funcomp.CGL.fit) == k_se)]], Znew, Zcnew, s = lam_se, ...)

}






#' @title
#' make predictions from a \code{"GIC.compCL"} object.
#' @description
#' This function makes prediction from a GIC  \code{compCL} model,
#' using the stored \code{compCL.fit} object and the optimal value chose for \code{lambda}.
#'
#' @param object fitted \code{"GIC.compCL"} model
#' @param Znew new compositional data
#' @param Zcnew new time-invariant covarites
#' @param s \code{"lam.min"} or user provided value
#' @param \dots not used. Other arguments to predict.
#'
#' @inheritParams coef.GIC.FuncompCGL
#' @seealso \code{\link{GIC.compCL}}
#' @examples
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c(beta, rep(0, times = p - length(beta)))
#' Comp_data = comp_Model(n = n, p = p,
#'                             rho = 0.2, sigma = 0.5,
#'                             gamma  = 0.5, add.on = 1:5,
#'                             beta = beta, intercept = FALSE)
#' test_data = comp_Model(n = 100, p = p,
#'                             rho = 0.2, sigma = 0.5,
#'                             gamma  = 0.5, add.on = 1:5,
#'                             beta = beta, intercept = FALSE)
#' GICm <- GIC.compCL(y = Comp_data$y,
#'                    Z = Comp_data$X.comp, Zc = Comp_data$Zc,
#'                    intercept = Comp_data$intercept,
#'                    lam = NULL,lambda.factor = 0.0001,
#'                    dfmax = p, outer_eps = 1e-10, mu_ratio = 1)
#' y_hat = predict(GICm, Znew = test_data$X.comp, Zcnew = test_data$Zc)
#' plot(y_hat, test_data$y)
#' abline(a = 0, b = 1, col = "red")
#' @export


predict.GIC.compCL <- function(object, Znew, Zcnew = NULL, s = "lam.min", ...){
  if (is.numeric(s)) {
    lam <- s
  } else if (s == "lam.min") {
    lam <- object$lam.min
  } else {
    stop("Invalid form for s")
  }

  predict(object$compCL.fit, Znew, Zcnew, s = lam, ...)
}














# #' @title
# #' plot the GIC curve produced by "GIC.FuncompCGL" object.
# #'
# #' @description
# #' Plots the GIC curve, as a function of the \code{lambda} values used.
# #'
# #' @param x fitted "\code{GIC.FuncompCGL}" object.
# #' @param xlab Either plot against \code{log(lambda)} (default) or \code{lambda}.
# #' @param \dots Other graphical parameters to plot.
# #' @param alpha scaler to penalty model selection.
# #' @details
# #' $$GIC(\eqn{\lambda}) = log(MSE) + Selection * \code{alpha}$$, for normal error.
# #' If include linear constraints, like in log contract and constraints group lasso model,
# #' selection = None zero group - 1; otherwise with consideration of linear constraints,
# #' selection = None zero group.
# #'
# #' BIC:alpha = log(n)*df, GIC:alpha=log(log(n)) /n * log(max(p*df,n)) * df
# #' @export
# #'
#
#
# plot.GIC.FuncompCGL <- function(x, xlab = c("log", "lambda"), alpha, ...) {
#   cvobj <- x
#   xlab <- match.arg(xlab)
#
#   switch(xlab,
#          "lambda" = {
#            xlab = "Lambda"
#            xvalue = drop(cvobj$lam)
#          },
#          "log" = {
#            xlab = "Log(Lambda)"
#            xvalue = log(drop(cvobj$lam))
#          })
#   GIC_val <- matrix(NA, nrow = length(cvobj$fit), length(cvobj$lam) )
#   for(i in length(cvobj$fit)) {
#
#   }
#   plot.args = list(x = xvalue,y = cvobj_use$cvm,
#                    ylim = range(cvobj_use$cvupper,cvobj_use$cvlo),
#                    xlab=xlab,
#                    ylab="GIC",
#                    type="n")
#   new.args=list(...)
#   if(length(new.args)) plot.args[names(new.args)]=new.args
#   do.call("plot",plot.args)
#   error.bars(xvalue, cvobj_use$cvupper, cvobj_use$culo, width=0.01, col="darkgrey")
#   points(xvalue,cvobj_use$cvm,pch=20,col="red")
#   axis(side=3,at=xvalue,labels=paste(cvobj$compCL.fit$df),tick=FALSE,line=0)
#   abline(v=switch(xlab,
#                   "Lambda" = cvobj_use$lam.min,
#                   "Log(Lambda)" = log(cvobj_use$lam.min)
#   ), lty=3)
#   abline(v = switch(xlab,
#                     "Lambda" = cvobj_use$lam.1se,
#                     "Log(Lambda)" = log(cvobj_use$lam.1se)
#   ),lty=3)
#   invisible()
# }
#
#




#' @title
#' plot the GIC curve produced by "GIC.FuncompCGL" object.
#'
#' @description
#' Plots the CIC curve as a function of the \code{lambda} values used.
#'
#' @param x fitted \code{"GIC.FuncompCGL"} object.
#' @param xlab Either plot against \code{log(lambda)} (default) or \code{-log{lambda}} or \code{lambda}.
#' @param \dots Other graphical parameters to plot.
#' @param k_list a vector, determining result of which \code{k} (degree freedom) to plot.
#'               If not provided, GIC curve for k that associate with \code{lambda.min} (red)
#'               is plotted.
#' @details
#' A plot is produced, and nothing is returned.
#'
#' @export
#' @importFrom grDevices rainbow
#'


plot.GIC.FuncompCGL <- function(x, xlab = c("log", "-log", "lambda"),
                                k_list, ...) {
  ## >>> x-axis <<<
  xlab <- match.arg(xlab)
  switch(xlab,
         "lambda" = {
           xlab = "Lambda"
           xvalue = drop(x$lam)
         },
         "log" = {
           xlab = "Log(Lambda)"
           xvalue = log(drop(x$lam))}
         ,
          "-log" = {
            xlab = "-Log(Lambda)"
            xvalue = -log(drop(x$lam))
          })

  ## <<< x-axis >>>

  ## >>> K list <<<
  k_listall <- as.numeric(names(x$Funcomp.CGL.fit))
  k_opt <- x$lammin['df']
  if(missing(k_list)) k_list <- k_opt
  if(is.numeric(k_list)) {
    k_list <- k_list[k_list %in% k_listall]
    if(!(k_opt %in% k_list)) warning("Optimal df 'k' is not included in the provied 'k_list' ")
  }

  k_list <- c(k_list[k_opt == k_list], k_list[k_opt != k_list])
  ## <<< K list >>>


  N_list <- apply(x$GIC[k_list- k_listall[1] + 1, , drop = FALSE], 1, function(x) length(x[!is.na(x)]))


  ## >>> plot <<<
  plot.args = list(x = xvalue[seq(max(N_list))], xlab = xlab,
                   y = x$GIC[k_list[which.max(N_list)] - k_listall[1] + 1, seq(max(N_list))],
                   #y = cvobj_use$cvm[k_list[1], ],
                   ylab = "GIC", ##
                   type = "n",
                   ylim = range(x$GIC[k_list- k_listall[1] + 1, ], x$GIC[k_list- k_listall[1] + 1, ], na.rm = TRUE)
  )

  new.args = list(...) #list(...)

  if(length(new.args)) plot.args[names(new.args)] = new.args
  do.call("plot", plot.args)

  pch_list = c(0, seq(length(k_list))[-1])
  col_list = rainbow(length(k_list)+1)[-1]

  for(l in 1:length(k_list)) {
    if(k_list[l] != k_opt) {
      points(xvalue,
             x$GIC[k_list[l] - k_listall[1] + 1, ],
             pch = pch_list[l],
             col = col_list[l] #col="limegreen",
             )
    } else {
      points(xvalue,
             x$GIC[k_opt - k_listall[1] + 1, ],
             pch=1,
             col="red")
      axis(side=3,at=xvalue,
           labels=paste(x$Funcomp.CGL.fit[[as.character(k_opt)]]$df[1:length(xvalue)]),
           tick=FALSE,
           line=0,
           #pos = plot.args$ylim[2],
           col.axis = "red")
      abline(v = switch(xlab,
                        "Lambda" = x$lammin['lam'],
                        "Log(Lambda)" = log(x$lammin['lam']),
                        "-Log(lambda)" = -log(x$lammin['lam'])),
            lty = 3, col = "red")
    }
  }

  if(k_opt %in% k_list){
    if(length(k_list) > 1) {
      legend(switch(xlab,
                    "Log(Lambda)" = "topright",
                    "Lambda" = "topright",
                    "-Log(Lambda)" = "topleft"),
             inset = 0.03,
             legend = c(paste("k", k_opt, sep = "="),
                        paste("k", k_list[k_opt != k_list], sep="=")),
             col = c("red", col_list[k_opt != k_list]),
             pch = c(1, pch_list[k_opt != k_list]),
             cex = 0.8, box.lty = 0
      )
    } else {
      legend(switch(xlab,
                    "Log(Lambda)" = "topright",
                    "Lambda" = "topright",
                    "-Log(Lambda)" = "topleft"),
             inset = 0.03,
             legend = paste("k", k_opt, sep = "="),
             col = "red",
             pch = 1,
             cex = 0.8, box.lty = 0
      )
    }

  } else {
    legend(switch(xlab,
                  "Log(Lambda)" = "topright",
                  "Lambda" = "topright",
                  "-Log(Lambda)" = "topleft"),
           inset = 0.03,
           legend = paste("k", k_list, sep="="),
           col = col_list,
           pch = pch_list,
           cex = 0.8, box.lty = 0
    )

  }
  ## <<< plot >>>
  invisible()

}


#' @title
#' plot the GIC curve produced by "GIC.compCL" object.
#'
#' @description
#' Plots the CIC curve as a function of the \code{lambda} values used.
#'
#' @param x fitted \code{"GIC.compCL"} object.
#' @param xlab Either plot against \code{log(lambda)} (default) or \code{-log{lambda}} or \code{lambda}.
#' @param \dots Other graphical parameters to plot.
#'
#' @details
#' A plot is produced, and nothing is returned.
#' @seealso \code{\link{GIC.compCL}}
#' @export
#'


plot.GIC.compCL <- function(x, xlab = c("log", "-log", "lambda"), ...) {
  xlab <- match.arg(xlab)
  switch(xlab,
         "lambda" = {
           xlab = "Lambda"
           xvalue = drop(x$lam)
         },
         "log" = {
           xlab = "Log(Lambda)"
           xvalue = log(drop(x$lam))
         },
         "-log" = {
           xlab = "-Log(Lambda)"
           xvalue = -log(drop(x$lam))
         })


  plot.args = list(x = xvalue,y = x$GIC,
                   ylim = range(x$GIC),
                   xlab=xlab,
                   ylab="GIC",
                   type="n")
  new.args=list(...)
  if(length(new.args)) plot.args[names(new.args)]=new.args
  do.call("plot",plot.args)
  points(xvalue, x$GIC)
  axis(side=3,at=xvalue,labels=paste(x$compCL.fit$df),tick=FALSE,line=0)

  abline(v=switch(xlab,
                  "Lambda" = x$lam.min,
                  "Log(Lambda)" = log(x$lam.min),
                  "-Log(Lambda)" = -log(x$lam.min)
  ), lty=3, col = "red")

  invisible()
}
