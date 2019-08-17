#'
#' GIC cirterion selection for compCL
#'
#' Calculate GIC for compCL, return value of \code{lam}.
#' The function follows Variable selection in regression with compositional covariates by
#' WEI LIN, PIXU SHI, RUI FENG AND HONGZHE LI
#'
#' @inheritParams compCL
#'
#' @param \dots other arguments that can be passed to compCL.
#'
#' @return an object of class \code{\link{GIC.compCL}} is returned.
#' \item{compCL.fit}{a fitted \code{\link{compCL}} object for the full data}
#' \item{lam}{the values of \code{lam} used in the fits}
#' \item{GIC}{a vector of GIC values for each \code{lam}}
#' \item{lam.min}{\code{lam} value such that minimize \eqn{GIC(\lambda)} }
#'
#' @details
#' \deqn{\textrm{GIC}(\lambda) = \log{\hat{\sigma}^2_\lambda}
#'                       + (s_\lambda - 1) \frac{\log{\log{n}}}{n} \log{max(p, n)} },
#' where \eqn{\hat{\sigma}^2_\lambda} is the MSE for fitted path.
#'
#'
#' \deqn{F_{n} =\frac{\phi^{n} - \psi^{n}}{\sqrt{5}}.}
#
# deqn ASCII example
#'
#\deqn{ \sigma = \sqrt{ \frac{Z}{n} \sum
#  \left[ \textstyle\frac{1}{2}\displaystyle
#    \left( \log \frac{H_i}{L_i} \right)^2  - (2\log 2-1)
#    \left( \log \frac{C_i}{O_i} \right)^2 \right] }
#}{sqrt(N/n * runSum(0.5 * log(OHLC[,2]/OHLC[,3])^2 -
#           (2*log(2)-1) * log(OHLC[,4]/OHLC[,1])^2, n))}
#'
#' @examples
#'
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c(beta, rep(0, times = p - length(beta)))
#' Comp_data = comp_Model(n = n, p = p,
#'                             rho = 0.2, sigma = 0.5,
#'                             gamma  = 0.5, add.on = 1:5,
#'                             beta = beta, intercept = FALSE)
#' Comp_data$Zc
#' GICm <- GIC.compCL(y = Comp_data$y,
#'                    Z = Comp_data$X.comp, Zc = Comp_data$Zc,
#'                    intercept = Comp_data$intercept,
#'                    lam = NULL,lambda.factor = 0.0001,
#'                    dfmax = p, outer_eps = 1e-10, mu_ratio = 1)
#' coef(GICm)
#' plot(GICm)
#' @export
#'


GIC.compCL <- function(y, Z, Zc = NULL, intercept = FALSE,
                       lam = NULL, ...) {
  digits = 5
  this.call <- match.call()

  y <- drop(y)
  n <- length(y)

  p <- ncol(Z)
  pc <- ifelse(is.null(Zc), 0, dim(Zc)[2])
  pc <- pc + as.integer(intercept)


  compCL.object <- compCL(y = y, Z = Z, Zc = Zc, intercept = intercept, lam = lam, ...)
  lam <- compCL.object$lam

  # >>> MSE <<<
  # pred_fun <- paste("predict", class(compCL.object), sep = "_")
  # pred_fun <- get(pred_fun)
  newx <- cbind(compCL.object$Z_log, Zc)
  #predmat <- suppressMessages(do.call(pred_fun, list(object, new_com, s = NULL)))
  predmat <- predict.linear(compCL.object, newx, s = NULL)
  cvraw <- (y - predmat)^2
  MSE <- colSums(cvraw) / n
  # <<< MSE >>>

  # >>> GIC curve <<<
  scaler <- log(log(n)) * log(max(p + pc, n)) / n
  S <- apply(abs(compCL.object$beta[1:p, ]) > 0, 2, sum) ## support
  GIC <- log(MSE) + scaler * ( ifelse(S>=2, S-1, 0) + pc)
  GIC <- round(GIC, digits = digits)
  # <<< GIC curve >>>

  # >>> lambda selection <<<
  GIC.min <- min(GIC[drop(compCL.object$df) > 0])
  idmin <- GIC <= GIC.min
  idmin[drop(compCL.object$df) < 2 ] <- FALSE
  lam.min <- max(lam[idmin])
  # <<< lambda selection >>>

  result <- list(compCL.fit = compCL.object,
                 lam = lam,
                 GIC = GIC,
                 lam.min = lam.min)

  class(result) <- "GIC.compCL"
  result$call <- this.call
  return(result)
}



#' @title
#' GIC cirterion selection for FuncompCGL
#'
#' @description
#' Calculate GIC for compCL, return value of \code{lam}.
#'
#### @usage
#### GIC.FuncompCGL <- function(y, X, Zc = NULL, lam = NULL,
####                            ref = NULL,
####                            W = rep(1,times = p - length(ref)),
####                            k = 4:10, nlam = 100, outer_maxiter = 1e+6,
####                            cut_off = c("Curve","Matrix", "NULL"),
####                            lower_tri = 0.01, ...)
####
#'
#'
#' @inheritParams FuncompCGL
#'
#' @param \dots other arguments that could be passed to FuncompCL.
#'
#' @return an object of class \code{\link{GIC.FuncompCGL}} is returned.
#' \item{Funcomp.CGL.fit}{a list, length of \code{k},
#'                        of fitted \code{\link{FuncompCGL}} object for the full data.
#'                        objects with S3 calss \code{\link{FuncompCGL}}}
#' \item{lam}{the values of \code{lam} used in the fits}
#' \item{MSE}{matrix of mean squared error with size \code{k} by \code{nlam} (the length of
#'            actually used \code{lambda} sequence, migth pre-stop by \code{dfmax} or
#'            \code{pfmax}). MSE is equivalent to likelihood under normal error model. \cr
#'            \strong{Could be edited for other linkage.}}
#' \item{Nzero}{a \code{k} by nlam matrix for Nzero group cut-off by \code{cut_off} and \code{lower_tri}}
#'
#' @examples
#'
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
#'
#' GIC_curve <- GIC_m1$GIC
#' k_opt <- GIC_m1$lammin['df']
#' beta_GIC <- coef(GIC_m1)
#' plot(GIC_m1, xlab = "log", k_list = k_list)
#' plot.args = list(x = seq(length(GIC_m1$lam)), #GIC_m1$lam, #log(GIC_m1$lam),
#'                  y = GIC_curve[1, ],
#'                  ylim = range(GIC_curve),
#'                  xlab= "lambda Index",#"lambda", #"log(lambda)",
#'                  ylab="GIC",
#'                  type="n")
#' #do.call("plot",plot.args)
#' # for(i in 1:length(k_list)) {
#' #
#' #   points(x = seq(length(GIC_m1$lam)), #GIC_m1$lam, #log(GIC_m1$lam),
#' #          y = GIC_curve[i, ], col = rainbow(length(k_list))[i])
#' #   text(length(GIC_m1$lam), #tail(log(GIC_m1$lam), 1),
#' #        GIC_curve[i, length(GIC_m1$lam)], labels=paste(k_list[i]),
#' #        cex= 1, pos= 4, col = rainbow(length(k_list))[i])
#' # }
#' # axis(3, at = pretty(seq(length(GIC_m1$lam))), labels = rev(pretty(GIC_m1$lam)))
#' # loc  = which(GIC_curve == min(GIC_curve), arr.ind = TRUE)
#'
#'
#'
#'
#' beta_C <- matrix(beta_GIC[1:(p*k_opt)], byrow = TRUE, nrow = p)
#' cat("colSums:", colSums(beta_C) , "\r\n")
#' #Non.zero <- which(abs(beta_C[,1]) > 0)
#' Non.zero <- apply(beta_C, 1, function(x) ifelse(max(abs(x)) >0, TRUE, FALSE))
#' Non.zero <- (1:p)[Non.zero]
#' cat("None zero groups:", Non.zero)
#' #vet(beta, p = p, k = k_opt)
#'
#' par(mfrow=c(1,4))
#' do.call("plot",plot.args)
#' for(i in 1:length(k_list)) {
#'   points(x = seq(length(GIC_m1$lam)), #log(GIC_m1$lam),
#'          y = GIC_curve[i, ], col = rainbow(length(k_list))[i], pch = seq(length(k_list))[i])
#'   text(length(GIC_m1$lam), #tail(log(GIC_m1$lam), 1),
#'        GIC_curve[i, length(GIC_m1$lam)], labels=paste(k_list[i]),
#'        cex= 1, pos= 4, col = rainbow(length(k_list))[i])
#' }
#' #axis(3, at = pretty(seq(length(GIC_m1$lam))), labels = rev(pretty(GIC_m1$lam)))
#'
#' matplot(sseq, beta_curve.true,
#'         ylab = "coeffcients curve", xlab = "TIME", #main = "TRUE",
#'         ylim = range(Data$beta[1:(p*df_beta)]),
#'         type = "l")
#' abline(a = 0, b = 0, col = "grey", lwd = 2)
#' title("TRUE", line = 0.5)
#' text(0, beta_curve.true[1, Non_zero.true], labels = paste(Non_zero.true))
#'
#' B <- splines::bs(Data$basis.info[,1], df = k_opt, intercept = TRUE)
#' beta_curve <- B %*% t(beta_C)
#' matplot(sseq, beta_curve,
#'         ylab = "coef", xlab = "TIME", #main = "ESTI",
#'         ylim = range(Data$beta[1:(p*df_beta)])#,
#'         #type = "l"
#' )
#' abline(a = 0, b = 0, col = "grey", lwd = 2)
#' title("Estimate", line = 0.5)
#' text(0, beta_curve[1, Non.zero], labels = paste(Non.zero))
#' text(tail(sseq, 1), beta_curve[dim(beta_curve)[1], Non.zero], labels = paste(Non.zero))
#' plot(apply(abs(beta_C),1,sum))
#' text(seq(length(GIC_m1$lam))[which(apply(abs(beta_C),1,sum) > 0)], #tail(log(GIC_m1$lam), 1),
#'      apply(abs(beta_C),1,sum)[which(apply(abs(beta_C),1,sum) > 0)],
#'      labels=paste(seq(length(GIC_m1$lam))[which(apply(abs(beta_C),1,sum) > 0)]),
#'      cex= 1, pos= 4)
#'
#' title(paste0("k=", k_opt), line = 0.5)
#' title(paste0("Method cgl"), outer=TRUE, line = -2)
#' par(mfrow=c(1,1))
#'
#' ##set a cutoff when you compute nonzeros
#' Non.zero <- apply(beta_C, 1, function(x)
#'   ifelse(sqrt(sum(x^2)) > sqrt(sum(beta_C^2))/100, TRUE, FALSE))
#' Non.zero <- (1:p)[Non.zero]
#' Non.zero
#' @seealso \code{\link{FuncompCGL}}
#' @export
#'
#'

GIC.FuncompCGL <- function(y, X, Zc = NULL, intercept = TRUE, ref = NULL,
                           lam = NULL, nlam = 100,
                           W = rep(1,times = p - length(ref)),
                           k = 4:10, outer_maxiter = 1e+6, mu_ratio=1.01,
                           ...) {
  y <- drop(y)
  n <- length(y)
  object <- as.list(seq(length(k)))
  this.call <- match.call()
  names(object) <- k
  if(!is.null(lam) || length(k) == 1) {

    # Case I
    if(dim(X)[1] == n ) p = dim(X)[2] / k else p <- dim(X)[2] - 2

    for(i in 1:length(k)){
      object[[i]] <- FuncompCGL(y = y, X = X, Zc = Zc, lam = lam,
                                W = W, ref = ref,
                                k = k[i],
                                nlam = nlam, outer_maxiter = outer_maxiter,
                                intercept = intercept, mu_ratio = mu_ratio,
                                ...)
    }
  } else {

    # Case II: find commom lambda first
    p <- ncol(X) - 2
    for(i in 1:length(k)){
      # Caculate integral Z and W matrix (if W is a functoin)
      object[[i]] <- FuncompCGL(y = y, X = X, Zc = Zc, nlam = 1,
                                W = W, ref = ref,
                                k = k[i], outer_maxiter = 0,
                                intercept = intercept, mu_ratio = mu_ratio,
                                ...)
    }

    lam0 <- max(sapply(object, "[[", "lam")) # Shared lam0 for different k
    dotlist = list(...)
    Trange_choice = dotlist$Trange
    # Solution path for each df k
    for(i in 1:length(k)) {
      #print(i)
      if(is.null(Trange_choice)){
      object[[i]] <- FuncompCGL(y = y, X = object[[i]]$Z, Zc = Zc,
                                W = object[[i]]$W, ref = ref,
                                lam = lam0, nlam = nlam,
                                outer_maxiter = outer_maxiter,
                                k = k[i], Trange = c(object[[i]]$call$Trange[1], object[[i]]$call$Trange[2]), ...)

      } else {

      object[[i]] <- FuncompCGL(y = y, X = object[[i]]$Z, Zc = Zc,
                                W = object[[i]]$W, ref = ref,
                                lam = lam0, nlam = nlam,
                                outer_maxiter = outer_maxiter,
                                k = k[i], ...)

      }
      object[[i]]$call <- as.list(object[[i]]$call)
      object[[i]]$call$k <- k[i]
    }
  }

  # lam_start = max of all lam_start
  # lam_end might diff beacuse of stopping criterion
  GIC.nlam <- sapply(object, function(x) length(drop(x$lam)))
  GIC.nlam_id <- which.min(GIC.nlam)
  GIC.lam <- object[[GIC.nlam_id]]$lam
  GIC.nlam <- GIC.nlam[GIC.nlam_id]


  ### >>> get GIC curves <<<
  GIC_curve <- matrix(NA, nrow = length(k), ncol = GIC.nlam)
  MSE <- GIC_curve
  # df of time-invariant covariates
  pc = ifelse(is.null(Zc), 0, ncol(Zc)) + as.integer(intercept)
  for(i in seq(length(k))) {
    df = k[i]
    scaler = log(max(p*df+ pc, n)) * log(log(n)) / n
    predmat <- predict.linear(object = object[[i]], newX = cbind(object[[i]]$Z, Zc))
    GIC_curve[i, ] <- apply(predmat[, 1:GIC.nlam] - y, 2, function(x) mean(x^2))
    GIC_curve[i, ] <- log(GIC_curve[i, ] / n)
    N_zero <- object[[i]]$df[1:GIC.nlam]


    if(mu_ratio != 0) {
      # clg method
      GIC_curve[i, ] = GIC_curve[i, ] +
        (ifelse(N_zero == 0, 0, N_zero - 1) * df + pc) * scaler
    } else {

      # if(is.null(ref)) {
      #   # naive method
      #   GIC_curve[i, ] = GIC_curve[i, ] + (N_zero * df + pc) * scaler

      # } else {
      #   # base method
      #   # df of beta is N_zero
      #   GIC_curve[i, ] = GIC_curve[i, ] +
      #                 (ifelse(N_zero == 0, 0, N_zero) * df + pc) * scaler
      # }
      GIC_curve[i, ] = GIC_curve[i, ] + (N_zero * df + pc) * scaler

    }
  }

  ## select the minimun of GIC_curve matrix
  lammin <- ggetmin(lam=GIC.lam, cvm = GIC_curve, k_list = k)
  lammin <- lammin$lam.min
  # result <- list(Funcomp.CGL.fit = object, lam = GIC.lam, MSE = MSE)
  result <- list(Funcomp.CGL.fit = object, lam = GIC.lam,
                 GIC = GIC_curve, lammin = lammin#, MSE = MSE
  )
  class(result) <- "GIC.FuncompCGL"
  result$call <- this.call
  return(result)
}





