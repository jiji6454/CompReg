#' @title
#' Fits regularization paths for longitudinal compositional data with group-lasso penalty.
#' @description
#' Fits regularization paths for longitudinal compositional data with group-lasso penalty at a sequence of regularization parameters lambda and fixed degree of freedom of basis.
#'
# @usage
# FuncompCGL <- function(y, X, Zc = NULL, intercept = TRUE,
#                        ref = NULL, k, degree = 3,
#                        basis_fun = c("bs", "OBasis", "fourier"),
#                        insert = c("FALSE", "X", "basis"),
#                        method = c("trapezoidal", "step"),
#                        interval = c("Original", "Standard"), Trange,
#                        T.name = "TIME", ID.name = "Subject_ID",
#                        W = rep(1,times = p - length(ref)),
#                        dfmax = p - length(ref),
#                        pfmax = min(dfmax * 1.5, p - length(ref)),
#                        lam = NULL, nlam = 100,
#                        lambda.factor = ifelse(n < p1, 0.05, 0.001),
#                        tol = 0, mu_ratio = 1.01,
#                        outer_maxiter = 1e+08, outer_eps = 1e-8,
#                        inner_maxiter = 1e+4, inner_eps = 1e-8)
#
#' @param y a vector of response variable.
#'
#' @param X a data frame or matrix.
#'       \itemize{
#'       \item If \code{dim(X)[1]} > \eqn{n}, \eqn{n} is the sample size,
#'             \code{X} should be a data frame of longitudinal compositinal predictors with number \eqn{p},
#'             including subject ID and time variable. Order of subject ID should be the same as that of \code{y}.
#'       \item If \code{dim(X)[1]}=\eqn{n}, \code{X} is considered as after taken integration, a \eqn{n*(p*k)} matrix.
#'       }
#'
#' @param T.name,ID.name characters specifying names of time varaible and Subject ID respectively in X,
#'                       only needed as X is data frame of longitudinal compositinal varaibles.
#'                       Default are \code{"TIME"} and \code{"Subject_ID"}.
#'
#' @param Zc A design matrix for control variables, could be missing. Default is NULL. No penalty is imposed.
#'
#' @param k a scaler, degree of freedom of basis.
#' @param ref reference variable. If \code{ref} is set to a scalar between \code{[1,p]}, log-contract method is applied with the variable
#'            \code{ref} as baseline. If \code{ref} = \code{NULL} (default value), constrained group lasso method is applied
#' @param degree degree of basis - default value is 3.
#'
#' @param basis_fun a function of basis. For now one of the following three types,
#'        \itemize{
#'        \item \code{bs} B-splines see \code{\link{bs}}.
#'        \item \code{OBasis} Orthoganal B-splies, see \code{\link{orthogonalsplinebasis}}.
#'        \item \code{fourier} Fourier basis, see \code{\link{fda}}
#'        }
#'        Default is \code{"bs"}.
#'
#' @param interval a character string sepcifying domain of integral
#'        \itemize{
#'          \item "Original" On original time scale, interval = range(Time).
#'          \item "Standard" Time points are mapped onto [0,1], interval = (0,1).
#'        }
#'        Default is \code{"Original"}
#'
#' @param insert way to interpolation. If \code{insert} = \code{"X"} or \code{"basis"}, dense time sequence is generated, equally space
#'               by \code{min(diff(sseq))/20)}, where \code{sseq} is sorted set of all observed time points.
#'               \itemize{
#'               \item \code{"FALSE"} no interpolation.
#'               \item \code{"X"} linear interpolation of compositional data.
#'               \item \code{"basis"} compositional data is considered as step function, imposing basis on un-observed time points for each subject.
#'               }
#'               Default is \code{"FALSE"}
#'
#' @param method method used to approximate integral.
#'               \itemize{
#'               \item \code{"trapezoidal"} Sum up area under trapezoidal formulated by values of function at two adjacent observed time points. See \code{\link{ITG_trap}}.
#'               \item \code{"step"} Sum up area under rectangle formulated by step function at observed time points. See \code{\link{ITG_step}}.
#'               }
#'               Default is \code{"trapezoidal"}
#'
#' @param W a vector in length of p (the total number of groups), matrix with dimension \code{p1*p1} or character specifying function
#'          used to calculate inverted weight matrix for each group.
#'          \itemize{
#'          \item If vector, works as penalty factor. Separate penalty weights can be applied to each group of beta'ss.
#'                to allow differential shrinkage. Can be 0 for some groups, which implies no shrinkage, and results in that group
#'                always being included in the model.
#'          \item If matrix, a block diagonal matrix. Diagonal elements are inverted weights matrics for each group.
#'          \item if character, user should provide the function for inverted weights matrics.
#'          }
#'          Default value is rep(1, times = p).
#' @param Trange range of time points
#' @inheritParams cglasso
#'
#' @return An object with S3 calss \code{\link{FuncompCGL}}
#' \item{Z}{integral matrix for longitudinal compositinal predictors with dimension \eqn{n*(p*k)}.}
#' \item{lam}{the actual sequence of \code{lam} values used.}
#' \item{df}{the number of non-zero groups in estimated coefficients for \code{Z} at each value of \code{lam}}
#' \item{beta}{a matrix of coefficients for \code{cbind{Z, Zc, 1_n}}, with \code{nlam} rows.}
#' \item{dim}{dimension of coefficient matrix}
#' \item{call}{the call that produced this object.}
#'
#'
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
#' print(m1)
#' beta <- coef(m1, s = m1$lam[20])
#' #beta <- coef(m1)
#'
#' beta_C <- matrix(beta[1:(p*k_use)], nrow = p, byrow = TRUE)
#' colSums(beta_C)
#' Non.zero <- apply(beta_C, 1, function(x) ifelse(max(abs(x)) == 0, FALSE, TRUE))
#' Non.zero <- (1:p)[Non.zero]
#' Non.zero
#' plot(m1, ylab = "L2", p = p , k = k_use)
#'
#' @export
#'
#' @import stats
#' @importFrom splines bs
#' @importFrom fda eval.basis create.fourier.basis
#' @importFrom orthogonalsplinebasis evaluate OBasis expand.knots
#'
#'

FuncompCGL <- function(y, X, Zc = NULL, intercept = TRUE, ref = NULL,
                       k, degree = 3, basis_fun = c("bs", "OBasis", "fourier"),
                       insert = c("FALSE", "X", "basis"), method = c("trapezoidal", "step"),
                       interval = c("Original", "Standard"), Trange,
                       T.name = "TIME", ID.name = "Subject_ID",
                       W = rep(1,times = p - length(ref)),
                       dfmax = p - length(ref), pfmax = min(dfmax * 1.5, p - length(ref)),
                       lam = NULL, nlam = 100, lambda.factor = ifelse(n < p1, 0.05, 0.001),
                       tol = 0, mu_ratio = 1.01,
                       outer_maxiter = 1e+08, outer_eps = 1e-8,
                       inner_maxiter = 1e+4, inner_eps = 1e-8) {

  n <- length(y)
  this.call <- match.call()
  basis_fun <- match.arg(basis_fun)
  method <- match.arg(method)
  insert <- match.arg(insert)
  interval <- match.arg(interval)
  interval_choice <- interval
  #print(interval)

  if(!is.null(lam) && TRUE %in% (lam < 0)) stop("User provided lambda must be positive vector")
  if(dim(X)[1] == n) {
    # Case I: already take integral and ready to go
    ## 1.
    Z <-  as.matrix(X)
    p1 <- ncol(Z)
    p <- p1 / k + length(ref)
    if(length(ref)>0) mu_ratio = 0
    #p <- p1 / k
    Trange_choice <- Trange
    # ## 2.Subtracting baseline integral ??????
    # if(is.numeric(ref)) {
    #   if( !(ref %in% 1:p) ) stop("Reference variable is out of range")
    #   mu_ratio = 0
    #   group.index <- matrix(1:p1, nrow = k)
    #   Z_ref <- Z[, group.index[, ref]]
    #   Z <- Z[, -group.index[, ref]]
    #   for(j in 1:(p-1)) Z[, group.index[, j]] <- Z[, group.index[, j]] - Z_ref
    #   p1 <- ncol(Z)
    # }


  } else {
    # Case II: take integral
    X <- as.data.frame(X)
    X.names <- colnames(X)[!colnames(X) %in% c(T.name, ID.name)]
    if( any(X[, X.names] == 0) ) stop("There is entry with value 0")
    # transform X into percentage and take log
    X[, X.names] <- log(X[, X.names] / rowSums(X[, X.names]))
    p <- length(X.names)
    #if(is.numeric(ref) && !(ref %in% 1:p)) stop("Reference variable is out of range")
    if(is.numeric(ref)) {
      if( !(ref %in% 1:p) ) stop("Reference variable is out of range")
      mu_ratio = 0
    }
    # Time range provided in case that sample do not cover the entire range
    if(missing(Trange)) Trange <- range(X[, T.name])
    Trange_choice <- Trange
    # cat("Trange_choice is")
    # print(Trange_choice)
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



    # In case that X is not represented in numerical order of Subject_ID
    X[, ID.name] <- factor(X[, ID.name], levels = unique(X[, ID.name]))

    # generate discrete basis matrix
    sseq <- round(sort(unique(c(Trange, as.vector(X[, T.name])))), 10)
    if(insert != "FALSE") sseq <- round(seq(from = interval[1], to = interval[2],
                                            #by = length(sseq) * 2,
                                            by = min(diff(sseq))/10), #### digit = 10
                                        10)

    # generate knots equally
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



    # take log-ratio w.r.t. reference
    X[, X.names] <- apply(X[, X.names], 2, function(x, ref) x - ref,
                          ref = if(is.null(ref)) 0 else X[, X.names[ref], drop = TRUE])
    # if ref is null, take the p+1-th component off (nothing);
    # ortherwise take the ref-th component off
    D <- split(X[, c(T.name, X.names[-ifelse(is.null(ref), p + 1, ref)])], X[, ID.name])
    p1 <- (p - length(ref)) * k
    # Z with (p - length(ref)) * k
    Z <- matrix(NA, nrow = n, ncol = p1)
    for(i in 1:n) Z[i, ] <- ITG(D[[i]], basis, sseq, T.name, interval, insert, method)$integral

  }



  group.index <- matrix(1:p1, nrow = k)
  if( is.vector(W) ) {
    ###cat(" ,length(W) = ", length(W))
    if(length(W) != (p1 / k) ) stop("W shoulde be a vector of length p=", p1/ k)
  } else if( is.matrix(W) ) {
    if(any(dim(W) - c(p1, p1) != 0)) stop('Wrong dimensions of Weights matrix')
  } else if( is.function(W) ){
    W_fun <- W
    W <- diag(p1)
    for(i in 1:(p1/k)) W[group.index[, i], group.index[, i]] <- W_fun(Z[, group.index[, i]])
  }

  ###cat(" ,is.null(ref) = ", is.null(ref))
  ###cat(" ,dfmax = ", dfmax, "\r\n")
  output <- cglasso(y = y, Z = Z, Zc = Zc, k = k, W = W, intercept = intercept,
                    mu_ratio = mu_ratio,
                    lam = lam, nlam = nlam, lambda.factor = lambda.factor,
                    dfmax = dfmax, pfmax = pfmax,
                    tol = tol,
                    outer_maxiter = outer_maxiter, outer_eps = outer_eps,
                    inner_maxiter = inner_maxiter, inner_eps = inner_eps)
  output$Z <- Z
  output$W <- W
  #output$call <- this.call
  output$call <- as.list(this.call)
  output$call['basis_fun'] = basis_fun
  output$call['method'] = method
  output$call['insert'] = insert
  # output$call$interval0 = interval[1]
  # output$call$interval1 = interval[2]
  output$call['interval'] = interval_choice
  if(!missing(Trange_choice)) output$call$Trange = c(Trange_choice[1], Trange_choice[2])
  output$ref <- ref
  class(output) <- "FuncompCGL"
  # dim = ((p - length(ref)) * k  + p_c + 1 ) x nlam
  output$dim <- dim(output$beta)
  #if(!missing(Trange))print(Trange)
  #print(this.call)
  return(output)
}





# FuncompCGL <- function(y, X, Zc = NULL, intercept = TRUE, k,
#                        T.name = "TIME", ID.name = "Subject_ID",
#                        degree = 3, basis_fun = c("bs", "OBasis", "fourier"),
#                        insert = c("FALSE", "X", "basis"), method = c("trapezoidal", "step"),
#                        interval = c("Original", "Standard"), Trange,
#                        W = rep(1, times = p), #pf = rep(1, times = p),
#                        dfmax = p, pfmax = min(dfmax * 1.5, p),
#                        lam = NULL, nlam = 100, lambda.factor = ifelse(n < p1, 0.05, 0.001),
#                        tol = 0, mu_ratio = 1.01,
#                        outer_maxiter = 1e+08, outer_eps = 1e-8,
#                        inner_maxiter = 1e+4, inner_eps = 1e-8
#                        #, ...
#                        ) {
#
#   n <- length(y)
#   this.call <- match.call()
#
#
#   basis_fun <- match.arg(basis_fun)
#   method <- match.arg(method)
#   insert <- match.arg(insert)
#   interval <- match.arg(interval)
#
#   if(dim(X)[1] == n) {
#     # already take integral
#     Z <- X
#     p1 <- dim(X)[2]
#     p <- p1 / k
#   } else {
#     # take integral
#     X <- as.data.frame(X)
#     X.names <- colnames(X)[!colnames(X) %in% c(T.name, ID.name)]
#     p <- length(X.names)
#     if( any(X[, X.names] == 0) ) stop("There is entry with value 0")
#     p1 <- p * k
#     #if(missing(Time)) Time <- X[, T.name]
#     Time <- X[, T.name]
#     #if(missing(Subject_ID))
#     Subject_ID <- unique(X[, ID.name])
#     X[, ID.name] <- factor(X[, ID.name], levels = Subject_ID )
#
#     if(missing(Trange)) Trange <- range(Time)
#     if(interval == "Standard") {
#       interval <- c(0,1)
#       X[, T.name] <- (Time - min(Trange[1])) / diff(Trange)
#       #X[, T.name] <- (Time - min(Time)) / diff(range(Time))
#     } else {
#       interval <- Trange
#     }
#
#     sseq <- round(sort(unique(c(Trange, as.vector(X[, T.name])))), 10)
#     if(insert != FALSE) sseq <- round(seq(from = interval[1], to = interval[2], by = min(diff(sseq))/20), 10)
#
#
#
#     nknots <- k - (degree + 1)
#     if(nknots > 0) {
#       knots <- ((1:nknots) / (nknots + 1)) * diff(interval) + interval[1]
#     } else {
#       knots <- NULL
#     }
#
#     basis <- switch(basis_fun,
#
#                     "bs" = bs(x = sseq, df = k, degree = degree,
#                               Boundary.knots = interval, intercept = TRUE),
#
#                     "fourier" = eval.basis(sseq,
#                                            basisobj = create.fourier.basis(rangeval = interval, nbasis = k )),
#
#                     "OBasis" = evaluate(OBasis(expand.knots(c(interval[1], knots, interval[2])),
#                                                order = degree + 1),
#                                         sseq)
#
#     )
#
#
#
#
#
#     X[, X.names] <- log(X[, X.names]/ rowSums(X[, X.names]))
#     D <- split(X[, c(T.name, X.names)], X[, ID.name])
#
#     Z <- matrix(NA, nrow = n, ncol = p1)
#     for(i in 1:n) {
#       #cat(i)
#       d <- D[[i]]
#       Z[i, ] <- ITG(d, basis, sseq, T.name, interval, insert, method)$integral
#
#       #Z <- sapply(D, function(x, basis, sseq, T.name, interval, insert, method)
#       #           ITG(x, basis, sseq, T.name, interval, insert, method)
#       #          ,sseq = sseq, basis = basis, T.name = T.name, interval = interval, insert = insert, method = method)
#       #Z <- t(Z)
#
#     }
#
#
#
#   }
#
#
#   Z <- as.matrix(Z)
#   group.index <- matrix(1:p1, nrow = k)
#   if( is.vector(W) ) {
#     if(length(W) != p) stop("W shoulde be a vector of length p")
#   } else if( is.matrix(W) ) {
#     if(any(dim(W) - c(p1, p1) != 0)) stop('Wrong dimensions of Weights matrix')
#   } else if( is.function(W) ){
#     W_fun <- W
#     W <- diag(p1)
#     for(i in 1:p) W[group.index[, i], group.index[, i]] <- W_fun(Z[, group.index[, i]])
#   }
#
#   output <- cglasso(y = y, Z = Z, Zc = Zc, k = k, W = W, intercept = intercept,
#                     mu_ratio = mu_ratio,
#                     lam = lam, nlam = nlam, lambda.factor = lambda.factor,
#                     dfmax = dfmax, pfmax = pfmax,
#                     tol = tol,
#                     outer_maxiter = outer_maxiter, outer_eps = outer_eps,
#                     inner_maxiter = inner_maxiter, inner_eps = inner_eps
#                     #, ...
#                     )
#   output$Z <- Z
#   output$W <- W
#   output$call <- this.call
#   class(output) <- "FuncompCGL"
#   output$dim <- dim(output$beta)
#   return(output)
# }
