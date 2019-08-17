
#' Fits regularization paths for group-lasso penalized learning problems with linear constraints.
#'
#' Fits regularization paths for linear constraints group-lasso penalized learning problems at a sequence of regularization parameters lambda.
#'
#'
#' @param y a vector of response variable with length \eqn{n}.
#' @param Z a matrix, of dimension \eqn{n \times p1}{n*p1}, for predictors to be imposed with group lasso penalty.
#' @param k
########    ,p two scalers, \eqn{p} is the total number of groups that would be imposed with group-lasso penalty,
#'          \eqn{k} is the size for each group in \eqn{Z}. All groups are with the same size. Number of groups is
#'          \eqn{p = p1 / k }.
#'
#' @param Zc a design matrix for those group-lasso-penalty free predictors. Can be missing. Default value is NULL.
#'
#' @param A,b linear equality constraints \eqn{A\beta_p1 = b}, where \eqn{b} is a vector with length \eqn{k},
#'            and \eqn{A} is a \eqn{k \times p1}{k*p1} matrix. Default values, \eqn{b} is a vector of 0's and
#'            \code{A = kronecker(matrix(1, ncol = p), diag(k))}.

## @param W_fun a function calculating weight matrix for each penalizd group. If missing, each weight matix is
##              set as \eqn{I_k}.
#'
#' @param W a vector in length of p (the total number of groups), matrix with dimension \code{p1*p1} or character specifying function
#'          used to caculate weight matrix for each group.
#'          \itemize{
#'          \item If vector, works as penalty factor. Separate penalty weights can be applied to each group of beta'ss to allow differential shrinkage. Can be 0 for some groups, which implies no shrinkage, and results in that group always being included in the model.
#'          \item If matrix, a block diagonal matrix. Diagonal elements are inverted weights matrics for each group.
#'          \item if character, user should provide the function.
#'          }
#'          Default value is rep(1, times = p).
#' @param lambda.factor the factor for getting the minimal lambda in \code{lam} sequence, where
#'                      \code{min(lam)} = \code{lambda.factor} * \code{max(lam)}.
#'                      \code{max(lam)} is the smallest value of \code{lam} for which all penalized group are zero's.
#'                      The default depends on the relationship between \eqn{n}
#######                 (the number of rows in the matrix ofpredictors)
#'                      and \eqn{p1}
#'                      (the number of predictors to be penalized).
#'                      If \eqn{n >= p1}, the default is \code{0.001}, close to zero.
#'                      If \eqn{n < p1}, the default is \code{0.05}. A very small value of
#'                      \code{lambda.factor}
#'                      will lead to a saturated fit. It takes no effect if there is user-defined lambda sequence.
#'
#' @param nlam the length of \code{lam} sequence. Default is 100.
#'
#' @param lam a user supplied lambda sequence. Typically, by leaving this option unspecified users can have the
#'            program compute its own \code{lam} sequence based on \code{nlam} and \code{lambda.factor}.
#######       , with length \code{nlam} equally spaced on log-scale.
#'            Supplying a value of lambda overrides this.
#'            If \code{lam} is provided but a scaler and \code{nlam}\eqn{>}1,
#'            \code{lam} sequence is also created starting from \code{lam}.
#'            If a sequence of lambda is provided, it is better to supply a decreasing one,
#'            if not, the program will sort user-defined \code{lambda} sequence in decreasing order
#'            automatically.
#'
## @param pf (needed?? coulde be combined into Weigt matrices)
##           penalty factor, a vector in length group \code{p}. Separate penalty weights can be applied to each
##          group of beta'ss to allow different shrinkage.  Can be 0 for some groups, which implies no
##           shrinkage, and results in that group always being included in the model.
#'
#' @param dfmax limit the maximum number of groups in the model. Useful for very large \eqn{p},
#'              if a partial path is desired. Default is \eqn{p}.
#'
#' @param pfmax limit the maximum number of groups ever to be nonzero. For example once a group enters the
#'             model along the path, no matter how many times it exits or re-enters model through the path,
#'             it will be counted only once. Default is \code{min(dfmax*1.5, p)}.
#
#' @param u    \code{u} is the inital value for penalty parameter of augmented Lanrange method adopted in
#'                   outer loop - default value is 1.
#'
#' @param mu_ratio   \code{mu_ratio} is the increasing ratio for \code{u}. Default value is 1.01.
#'                   Inital values for scaled Lagrange multipliers are set as 0's.
#'                   If \code{mu_ratio} < 1,
########             \code{u} is set as 0 and \code{outer_maxiter} = 1.
#'                   there is no linear constraints included. Group lasso coefficients are estimated.
#'
#' @param tol tolerance for vectors betas to be considered as none zero's. For example, coefficient
#'                   \eqn{\beta_j} for group j, if \eqn{max(abs(\beta_j))} < \code{tol}, set \eqn{\beta_j} as 0's.
#'                   Default value is 0.
#'
#' @param outer_maxiter,outer_eps
#'                   \code{outer_maxiter} is the maximun munber of loops allowed for Augmented Lanrange method;
#'                   and \code{outer_eps} is the convergence termination tolerance.
#'
#' @param inner_maxiter,inner_eps
#'                   \code{inner_maxiter} is the maximun munber of loops allowed for blockwise-GMD;
#'                   and \code{inner_eps} is the convergence termination tolerance.
#'
#' @param intercept whether to include intercept. Default is TRUE.
#'
#' @return A list
#' \item{lam}{the actual sequence of lambda values used }
#' \item{df}{the number of nonzero groups in estimated coefficients for \code{Z} at each value of lambda}
#' \item{path}{a matrix of coefficients}
#'
#' @export
#'
#' @importFrom Rcpp evalCpp sourceCpp
#' @import RcppArmadillo
#' @useDynLib compReg, .registration = TRUE


cglasso <- function(y, Z, Zc = NULL, k,
                    #W_fun = NULL,
                    W = rep(1, times = p),
                    intercept = TRUE,
                    A =  kronecker(matrix(1, ncol = p), diag(k)), b = rep(0, times = k),
                    u = 1, mu_ratio = 1.01,
                    lam = NULL, nlam = 100,lambda.factor = ifelse(n < p1, 0.05, 0.001),
                    dfmax = p, pfmax = min(dfmax * 1.5, p),
                    #pf = rep(1, times = p),
                    tol = 0,
                    outer_maxiter = 3e+08, outer_eps = 1e-8,
                    inner_maxiter = 1e+6, inner_eps = 1e-8
                    #,estar = 1 + 1e-8
                    ) {


  y <- drop(y)
  Z <- as.matrix(Z)
  n <- length(y)
  p1 <- dim(Z)[2]
  p <- p1 / k
  group.index <- matrix(1:p1, nrow = k)



  inter <- as.integer(intercept)
  Zc <- as.matrix(cbind(Zc, rep(inter, n)))
  if(inter == 1) {
    Zc_proj <- solve(crossprod(Zc)) %*% t(Zc)
  } else {
    if(dim(Zc)[2] == 1) Zc_proj <- t(Zc) else Zc_proj <- rbind(solve(crossprod(Zc[, -ncol(Zc)])) %*% t(Zc[, -ncol(Zc)]), rep(inter, n))
  }
  beta_c <- as.vector(Zc_proj %*% y)
  beta.ini <- c(rep(0, times = p1), beta_c)

  if( is.vector(W) ){
    if(length(W) != p) stop("W should be a vector of length p")
    W_inver <- diag(p1)
    pf <- W
  } else if( is.matrix(W) ) {
    if(any(dim(W) - c(p1, p1) != 0)) stop('Wrong dimension of Weights matrix')
    pf <- rep(1, times = p)
    W_inver <- W
  }

  Z <- Z %*% W_inver
  A <- A %*% W_inver



  # if(is.null(lam)) {
  #   lam0 <- lam.ini(Z = Z, y = y - Zc %*% beta_c, ix = group.index[1, ], iy = group.index[k ,], pf = pf)
  #   lam <- exp(seq(from=log(lam0), to=log(lambda.factor * lam0),length=nlam))
  # } else if(length(lam) == 1) {
  #   lam <- exp(seq(from=log(lam), to=log(lambda.factor * lam),length=nlam))
  # }

  if(is.null(lam)) {
    lam0 <- lam.ini(Z = Z, y = y - Zc %*% beta_c, ix = group.index[1, ], iy = group.index[k ,], pf = pf)
    lam <- exp(seq(from=log(lam0), to=log(lambda.factor * lam0),length=nlam))
  } else {
    if(length(lam) == 1 && nlam > 1) {
    lam <- exp(seq(from=log(lam), to=log(lambda.factor * lam),length=nlam))
    } else {
    lam <- sort(lam, decreasing = TRUE)
    }
  }

  pfmax <- as.integer(pfmax)
  dfmax <- as.integer(dfmax)
  inner_maxiter <- as.integer(inner_maxiter)
  outer_maxiter <- as.integer(outer_maxiter)
  inner_eps <- as.double(inner_eps)
  outer_eps <- as.double(outer_eps)
  tol <- as.double(tol)
  u <- as.double(u)
  mu_ratio <- as.double(mu_ratio)
  #estar <- as.double(estar)

  output <- ALM_GMD(y = y, Z = Z, Zc = Zc, Zc_proj = Zc_proj, beta = beta.ini, lambda = lam,
                    pf = pf, dfmax = dfmax, pfmax = pfmax, A = A, b = b,
                    group_index = group.index, u_ini = u, mu_ratio = mu_ratio,
                    inner_eps = inner_eps, outer_eps = outer_eps,
                    inner_maxiter = inner_maxiter, outer_maxiter = outer_maxiter, tol = tol
                    #, estar = estar
                    )

  output$beta[1:p1, ] <-  W_inver %*% output$beta[1:p1, ]
  return(output)

}



#'
#' Fits regularization paths for lasso penalized learning problems with linear constraints.
#'
#' Fits regularization paths for linear constraints lasso penalized learning problems at a sequence of regularization parameters lambda.
#'
#'
#' @param y a vector of response variable with length n.
#'
#' @param Z a \eqn{n*p} matrix after taking log transformation on compositional data.
#'
#' @param Zc a design matrix of other covariates considered. Default is \code{NULL}.
#'
#' @param intercept Whether to include intercept in the model. Default is TRUE.
#'
#' @param pf penalty factor, a vector in length of p. Separate penalty weights can be applied to each coefficience
#'           \eqn{\beta} for composition variates to allow differential shrinkage. Can be 0 for some \eqn{\beta}'s,
#'           which implies no shrinkage, and results in that composition always being included in the model.
#'           Default value for each entry is the 1.
#'
#' @param A,b linear equality constraints \eqn{A\beta_p = b}, where \eqn{b} is a scaler,
#'            and \eqn{A} is a vector with length \code{p}. Default values, \eqn{b} is a vector of 0 and
#'            \code{A = rep(1, times = p)}.
#'
#' @param u    \code{u} is the inital value for penalty parameter of augmented Lanrange method adopted in
#'                   outer loop - default value is 1.
#' @param beta.ini inital value of beta
#' @inheritParams FuncompCGL
#'
#' @export
#'
#'


classo <- function(y, Z, Zc = NULL, intercept = TRUE,
                   pf = rep(1, times = p),
                   lam = NULL, nlam = 100,lambda.factor = ifelse(n < p, 0.05, 0.001),
                   dfmax = p, pfmax = min(dfmax * 1.5, p),
                   u = 1, mu_ratio = 1.01, tol = 1e-10,
                   outer_maxiter = 3e+08, outer_eps = 1e-8,
                   inner_maxiter = 1e+6, inner_eps = 1e-8,
                   A = rep(1, times = p), b = 0,
                   beta.ini

) {

  this.call <- match.call()
  y <- drop(y)
  Z <- as.matrix(Z)
  n <- length(y)
  p <- dim(Z)[2]
  if(length(A) != p) stop("Length of vector A in Ax = b is wrong")
  inter <- as.integer(intercept)
  Znames <- colnames(Z)
  if (is.null(Znames)) Znames <- paste0("Z", seq(p))

  Zc <- as.matrix(cbind(Zc, rep(inter, n)))
  if(inter == 1) {
    Zc_proj <- solve(crossprod(Zc)) %*% t(Zc)
  } else {
    if(dim(Zc)[2] == 1) Zc_proj <- t(Zc) else Zc_proj <- rbind(solve(crossprod(Zc[, -ncol(Zc)])) %*% t(Zc[, -ncol(Zc)]), rep(inter, n))
  }

  if(missing(beta.ini)) {
    beta_c <- as.vector(Zc_proj %*% y)
    beta.ini <- c(rep(0, times = p), beta_c)
  }

  m <- dim(Zc)[2]
  if(m > 1) {
    Zcnames <- colnames(Z)
    if (is.null(Zcnames)) Zcnames <- paste0("Zc", seq(m-1))
  } else {
    Zcnames <- NULL
  }

  if(is.null(lam)) {
    lam0 <- t(Z) %*% (y - Zc %*% beta_c) / n
    lam0 <- max(abs(lam0) / pf)
    lam <- exp(seq(from=log(lam0), to=log(lambda.factor * lam0),length=nlam))
  }

  pfmax <- as.integer(pfmax)
  dfmax <- as.integer(dfmax)
  inner_maxiter <- as.integer(inner_maxiter)
  outer_maxiter <- as.integer(outer_maxiter)
  inner_eps <- as.double(inner_eps)
  outer_eps <- as.double(outer_eps)
  tol <- as.double(tol)
  u <- as.double(u)
  mu_ratio <- as.double(mu_ratio)
  b <- as.double(b)
  #estar <- as.double(estar)

  output <- ALM_BG(y = y, Z = Z, Zc = Zc, Zc_proj = Zc_proj, beta = beta.ini,
                   lambda = lam, pf = pf, dfmax = dfmax, pfmax = pfmax,
                   inner_eps = inner_eps, outer_eps = outer_eps,
                   inner_maxiter = inner_maxiter, outer_maxiter = outer_maxiter,
                   u_ini = u, mu_ratio = mu_ratio, tol = tol,
                   A = A, b = b

  )

  output$call <- this.call
  class(output) <- "compCL"
  output$dim <- dim(output$beta)
  output$lam <- drop(output$lam)
  dimnames(output$beta) <- list(c(Znames, Zcnames, "Intercept"), paste0("L", seq(along = output$lam)) )
  return(output)

}


