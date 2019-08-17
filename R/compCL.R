#' @title
#' Fits regularization paths for compositional data with lasso penalty.
#'
#' @description
#' Fit compositional data via penalized \emph{log-contrast} model which is proposed in
#' \href{https://academic.oup.com/biomet/article/101/4/785/1775476}{Variable selection in regression with compositional covariates}.
#' The regularization paths is computed for the lasso penalty at a grid of values for the regularization parameter \code{lambda}.
#'
#'
#'
#'
# @usage
# compCL <- function(y, Z, Zc = NULL, intercept = TRUE, pf = rep(1, times = p),
#                      lam = NULL, nlam = 100, lambda.factor = ifelse(n < p, 0.05, 0.001),
#                      dfmax = p, pfmax = min(dfmax * 1.5, p),
#                      u = 1,mu_ratio = 1.01, tol = 1e-10,
#                      outer_maxiter = 1e+08, outer_eps = 1e-8,
#                      inner_maxiter = 1e+3, inner_eps = 1e-6)
#
#'
#'
#' @param y a vector of response variable with length n.
#'
#' @param Z a \eqn{n*p} matrix of compositional data or categorical data.
#'          If \code{Z} is categorical data, i.e., row sums of \code{Z} differ from 1, \code{Z} is transformed to compositional data
#'          by diving each row of it's corresponding sum.
#'          \code{Z} could NOT include entry of 0.
#'
#' @param Zc a \eqn{n*p_c} design matrix of other covariates considered except composition . Default is \code{NULL}.
#'
#' @param intercept Whether to include intercept in the model. Default is TRUE.
#'
#' @param pf penalty factor, a vector in length of p. Separate penalty weights can be applied to each coefficience
#'           \eqn{\beta} for composition variates to allow differential shrinkage. Can be 0 for some \eqn{\beta}'s,
#'           which implies no shrinkage, and results in that composition always being included in the model.
#'           Default value for each entry is the 1.
#'
#' @param u \code{u} is the inital value for penalty parameter of augmented Lanrange method.
#'                   Default value is 1.
#'
#' @param lambda.factor the factor for getting the minimal lambda in \code{lam} sequence, where
#'                      \code{min(lam)} = \code{lambda.factor} * \code{max(lam)}.
#'                      \code{max(lam)} is the smallest value of \code{lam} for which all penalized group are zero's.
#'                      The default depends on the relationship between \eqn{n} and \eqn{p}.
#'                      If \eqn{n >= p}, the default is \code{0.001}, close to zero.
#'                      If \eqn{n < p}, the default is \code{0.05}.
#'                      A very small value of \code{lambda.factor}
#'                      will lead to a saturated fit. It takes no effect if there is user-defined lambda sequence.
#'
#'
#' @param mu_ratio   \code{mu_ratio} is the increasing ratio for \code{u}. Default value is 1.01.
#'                   Inital values for scaled Lagrange multipliers are set as 0's.
#'                   If \code{mu_ratio} < 1, \code{u} is set as 0 and \code{outer_maxiter} = 1.
#'                   Then there is no linear constraints included. Lasso coefficients are estimated.
#'
#' @param tol tolerance for betas to be considered as none zero's. For example, coefficient
#'                   \eqn{\beta_j} for composition component j, if \eqn{abs(\beta_j)} < \code{tol}, set \eqn{\beta_j} as 0.
#'                   Default value is 0.
#'
#'
#' @inheritParams FuncompCGL
#'
#' @return An object with S3 calss \code{\link{compCL}.}
## \item{beta}{a matrix of coefficients for \code{cbind{Z, Zc, 1_n}}, with \code{nlam} rows.}
#' \item{beta}{a matrix of coefficients for \eqn{p+p_c+1} rows.
#'             If \code{intercept=FALSE}, then the last row of \code{beta} is set to 0's.}
#' \item{lam}{the actual sequence of \code{lam} values used.}
#' \item{df}{the number of non-zero \eqn{\beta}s in estimated coefficients for \code{Z} at each value of \code{lam}.}
#' \item{npass}{total iteration conducted in computation.}
#' \item{error}{error message for path at each each value of \code{lam}. If 0, no error occurs.}
#' \item{call}{the call that produced this object.}
#' \item{dim}{dimension of coefficient matrix.}
#'
#'
#' @examples
#'
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c( beta, rep(0, times = p - length(beta)) )
#' Comp_data = comp_Model(n = n, p = p, rho = 0.2, sigma = 0.5, gamma  = 0.5, add.on = 1:5,
#'                             beta = beta, intercept = FALSE)
#' m1 <- compCL(y = Comp_data$y, Z = Comp_data$X.comp, Zc = Comp_data$Zc,
#'              intercept = Comp_data$intercept,
#'              pf = rep(1, times = p),
#'              lam = NULL, nlam = 100,lambda.factor = ifelse(n < p, 0.05, 0.001),
#'              dfmax = 20, #pfmax = min(dfmax * 1.5, p),
#'              mu_ratio = 1, tol = 1e-10,
#'              outer_maxiter = 1e8, outer_eps = 1e-10,
#'              inner_maxiter = 1e3, inner_eps = 1e-6)
#'
#' print(m1)
#' coef(m1, s = m1$lam[50]) # extract coefficients at a single value of lambda
#' pm1 = predict(m1,newx=Comp_data$X.comp[1:10,],newz=Comp_data$Zc[1:10, ],s=m1$lam[45:55])
#' print(pm1)
#' plot(m1)
#' plot(pm1[, 1], Comp_data$y[1:10])
#' @export
#'
#'



# compCL <- function(y, Z, Zc = NULL, intercept = TRUE,
#                    pf = rep(1, times = p),
#                    lam = NULL, nlam = 100, lambda.factor = ifelse(n < p, 0.05, 0.001),
#                    dfmax = p, pfmax = min(dfmax * 1.5, p),
#                    u = 1,
#                    mu_ratio = 1.01, tol = 1e-10,
#                    outer_maxiter = 1e+08, outer_eps = 1e-8,
#                    inner_maxiter = 1e+3, inner_eps = 1e-6) {
#   #u <- 1
#   if(!is.null(lam) && TRUE %in% (lam < 0)) stop("User provided lambda must be positive vector")
#   this.call <- match.call()
#   y <- drop(y)
#   Z <- as.matrix(Z)

#   if( any(abs(rowSums(Z) - 1) > 1e-10) ) {
#     message("Z is transformed into compositional data by deviding rowSums")
#     Z <- Z / rowSums(Z)
#   }
#   if(any(Z == 0)) stop("There is zero entry in compositional data")
#   Z <- log(Z)

#   n <- length(y)
#   p <- dim(Z)[2]
#   inter <- as.integer(intercept)
#   Znames <- colnames(Z)
#   if (is.null(Znames)) Znames <- paste0("Z", seq(p))

#   Zc <- as.matrix(cbind(Zc, rep(inter, n)))
#   if(inter == 1) {
#     Zc_proj <- solve(crossprod(Zc)) %*% t(Zc)
#   } else {
#     if(dim(Zc)[2] == 1) Zc_proj <- t(Zc) else Zc_proj <- rbind(solve(crossprod(Zc[, -ncol(Zc)])) %*% t(Zc[, -ncol(Zc)]), rep(inter, n))
#   }
#   beta_c <- as.vector(Zc_proj %*% y)
#   beta.ini <- c(rep(0, times = p), beta_c)

#   m <- dim(Zc)[2]
#   if(m > 1) {
#     Zcnames <- colnames(Z)
#     if (is.null(Zcnames)) Zcnames <- paste0("Zc", seq(m-1))
#   } else {
#     Zcnames <- NULL
#   }

#   if(is.null(lam)) {
#     lam0 <- t(Z) %*% (y - Zc %*% beta_c) / n
#     lam0 <- max(abs(lam0) / pf)
#     lam <- exp(seq(from=log(lam0), to=log(lambda.factor * lam0),length=nlam))
#   }else{
#     if(length(lam) == 1 && nlam > 1) {
#       lam <- exp(seq(from=log(lam), to=log(lambda.factor * lam),length=nlam))
#     } else {
#       lam <- sort(lam, decreasing = TRUE)
#     }
#   }

#   pfmax <- as.integer(pfmax)
#   dfmax <- as.integer(dfmax)
#   inner_maxiter <- as.integer(inner_maxiter)
#   outer_maxiter <- as.integer(outer_maxiter)
#   inner_eps <- as.double(inner_eps)
#   outer_eps <- as.double(outer_eps)
#   tol <- as.double(tol)
#   u <- as.double(u)
#   mu_ratio <- as.double(mu_ratio)
#   #estar <- as.double(estar)

#   output <- ALM_B(y = y, Z = Z, Zc = Zc, Zc_proj = Zc_proj, beta = beta.ini,
#                   lambda = lam, pf = pf, dfmax = dfmax, pfmax = pfmax,
#                   inner_eps = inner_eps, outer_eps = outer_eps,
#                   inner_maxiter = inner_maxiter, outer_maxiter = outer_maxiter,
#                   u_ini = u, mu_ratio = mu_ratio, tol = tol

#   )

#   output$call <- this.call
#   class(output) <- "compCL"
#   output$dim <- dim(output$beta)
#   output$lam <- drop(output$lam)
#   output$Z_log <- Z

#   dimnames(output$beta) <- list(c(Znames, Zcnames, "Intercept"), paste0("L", seq(along = output$lam)) )
#   return(output)

# }



compCL <- function(y, Z, Zc = NULL, intercept = TRUE,
                   pf = rep(1, times = p),
                   lam = NULL, nlam = 100, lambda.factor = ifelse(n < p, 0.05, 0.001),
                   dfmax = p, pfmax = min(dfmax * 1.5, p),
                   u = 1,
                   mu_ratio = 1.01, tol = 1e-10,
                   outer_maxiter = 1e+08, outer_eps = 1e-8,
                   inner_maxiter = 1e+3, inner_eps = 1e-6) {
  #u <- 1
  if(!is.null(lam) && TRUE %in% (lam < 0)) stop("User provided lambda must be positive vector")

  this.call <- match.call()

  y <- drop(y)
  n <- length(y)

  # Z <- as.matrix(Z)
  # if( any(abs(rowSums(Z) - 1) > 1e-10) ) {
  #   message("Z is transformed into compositional data by deviding rowSums")
  #   Z <- Z / rowSums(Z)
  # }
  # if(any(Z == 0)) stop("There is zero entry in compositional data")
  # Z <- log(Z)
  Z <- proc.comp(Z)
  p <- ncol(Z) #dim(Z)[2]
  Znames <- colnames(Z)
  if (is.null(Znames)) Znames <- paste0("Z", seq(p))

  inter <- as.integer(intercept)

  # Zc <- as.matrix(cbind(Zc, rep(inter, n)))
  if(is.null(Zc)) Zc <- matrix(inter, nrow = n) else Zc <- cbind2(as.matrix(Zc), inter)

  if(inter == 1) {
    Zc_proj <- solve(crossprod(Zc)) %*% t(Zc)
  } else {
    if(dim(Zc)[2] == 1) Zc_proj <- t(Zc) else Zc_proj <- rbind(solve(crossprod(Zc[, -ncol(Zc)])) %*% t(Zc[, -ncol(Zc)]), rep(inter, n))
  }

  beta_c <- as.vector(Zc_proj %*% y)
  beta.ini <- c(rep(0, times = p), beta_c)

  m <- ncol(Zc) #dim(Zc)[2]
  if(m > 1) {
    Zcnames <- colnames(Zc)
    if (is.null(Zcnames)) Zcnames <- paste0("Zc", seq(m-1))
  } else {
    Zcnames <- NULL
  }

  if(is.null(lam)) {
    lam0 <- t(Z) %*% (y - Zc %*% beta_c) / n
    lam0 <- max(abs(lam0) / pf)
    lam <- exp(seq(from=log(lam0), to=log(lambda.factor*lam0), length=nlam))
  } else {
    if(length(lam) == 1 && nlam > 1) {
      lam <- exp(seq(from=log(lam), to=log(lambda.factor*lam), length=nlam))
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

  output <- ALM_B(y = y, Z = Z, Zc = Zc, Zc_proj = Zc_proj, beta = beta.ini,
                  lambda = lam, pf = pf, dfmax = dfmax, pfmax = pfmax,
                  inner_eps = inner_eps, outer_eps = outer_eps,
                  inner_maxiter = inner_maxiter, outer_maxiter = outer_maxiter,
                  u_ini = u, mu_ratio = mu_ratio, tol = tol

  )

  output$call <- this.call
  class(output) <- "compCL"
  output$dim <- dim(output$beta)
  output$lam <- drop(output$lam)
  output$Z_log <- Z
  dimnames(output$beta) <- list(c(Znames, Zcnames, "Intercept"), paste0("L", seq(along = output$lam)) )
  return(output)

}



