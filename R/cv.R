

#' @title
#'
#' Cross-validation for FuncompCGL
#'
#' @description
#' Does nfolds cross-validation for FuncompCGL, return value of \code{lam}
#' and df \code{k}. The function is modified based on the \code{cv} function
#' from \code{glmnet} package
#'
#' @usage
#' cv.FuncompCGL(y, X, Zc = NULL, lam = NULL, ref = NULL,
#'               W = rep(1,times = p - length(ref)),
#'               k = 4:10,
#'               foldid, nfolds = 10, nlam = 100,
#'               trim = 0, outer_maxiter = 1e+06,
#'               keep = FALSE, ...)
#'
#'
#' @param y a vector of response variable.
#'
#' @param X a data frame of longitudinal compositinal predictors with number \eqn{p},
#'          subject ID and time variable. Order of subject ID should be the same as that of \code{y}.
#'          If df \code{k} is a scalar, X could be the matrix with dimension \code{n*(k*p - length(ref))}
#'          after integral.
#  In this case no selection on df(\eqn{k}).
#'
#' @param k a vector or scalar consists of df for basis - default is 4:10.
#'
#' @param nfolds number of folds - default is 10. Smallest value allowable is nfolds=3.
#'
#' @param foldid an optional vector of values between 1 and \code{nfolds} identifying what fold each
#'               observation is in. If supplied, \code{nfold} can be missing.
#'
#' @param trim a scaler specifying percentage to be trimmed off for prediction error - default is 0. \cr
#'             \strong{This feature could be deleted later.}
#'
# @param T.name,ID.name characters specifying names of time and Subject ID respectively.
#'
#' @param ... other arguments that can be passed to \code{\link{FuncompCGL}}.
#'
#' @param nlam the length of \code{lam} sequence - default is 100.
#'
#' @param outer_maxiter  maximun munber of loops allowed for Augmented Lanrange method.
#'
#' @param keep If \code{keep=TRUE}, a prevalidated array is returned containing fitted values for each observation,
#'            of lambda and k. This means these fits are computed with this observation and the rest of
#'            its fold omitted. The folid vector is also returned. Default is keep=FALSE
#' @inheritParams FuncompCGL
#'
#' @return a list
#' \item{Funcomp.CGL.fit}{a list, length of \code{k},
#'                        of fitted \code{\link{FuncompCGL}} object for the full data.
#'                        objects with S3 calss \code{\link{FuncompCGL}}}
#' \item{lam}{the values of \code{lam} used in the fits}
#' \item{Ftrim}{a list for cross-validation result with trim = 0.
#'                \itemize{
#'                \item \code{cvm} the mean cross-validated error without trimming - a matrix of  \code{length(k)*length(lam)}
#'                \item \code{cvsd} estimate of standard error of cvm without trimming- a matrix of  \code{length(k)*length(lam)}
#'                \item \code{cvup} a matrix of  \code{length(k)*length(lam)}
#'                \item \code{cvlo} a matrix of  \code{length(k)*length(lam)}
#'                \item \code{lam.min} The optimal value of \code{k} and \code{lam} that gives minimum cross validation error cvm without trimming
#'                \item \code{lam.1se} The optimal value of \code{k} and \code{lam} for "lam.1se" without trimming}
#'                }
#'
#' \item{Ttrim}{a list for cross-validation result with trimming (if \code{trim} provided).
#'              Same as these for Ftrim. \cr
#'              \strong{This feature could be deleted later.}}
#'
#' \item{fit.preval, foldid}{fitting matrix and foldid. Only keept when \code{keep=TRUE}.}
#'
#'
#' @examples
#'
#' df_beta = 5
#' p = 30
#' beta_C_true = matrix(0, nrow = p, ncol = df_beta)
#' beta_C_true[3, ] <- c(-1, 0, 0, 0, -0.5)
#' beta_C_true[1, ] <- c(1, 0, 1 , 0, -0.5)
#' beta_C_true[2, ] <- c(0, 0,  -1,  0,  1)
#'
#' nfolds = 10
#' k_list <- c(4,5)
#' n_train = 100
#' n_test = 500
#'
#' Data <- Model(n = 100, p = p, m = 0, intercept = TRUE,
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
#' m <- ifelse(is.null(Zc), 0, dim(Zc)[2]) #+ as.integer(intercept)
#' m1 <- m + as.integer(intercept)
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
#' rule <- "lam.min"
#' # cv_cgl, Constrained group lasso
#' cv_cgl <-  cv.FuncompCGL(y = y, X = X, Zc = Zc, intercept = intercept,
#'                          W = rep(1, p), #W = function(x){ diag( 1 / apply(x, 2, sd) ) },
#'                          k = k_list, trim = 0,
#'                          foldid = foldid, #nfolds = 10,
#'                          tol = 0, inner_eps = 1e-6, inner_maxiter = 1E3,
#'                          dfmax = 30, lambda.factor = 1e-3,
#'                          mu_ratio = 1, outer_eps = 1e-6,
#'                          keep = TRUE, Trange = c(0,1)
#'                          #lam = c(exp(seq(log(0.01),log(0.00001),length=100)),0)
#'                          )
#' plot(cv_cgl,k_list = k_list)
#' cv_cgl$Ftrim[c("lam.min", "lam.1se")]
#' beta <-  coef(cv_cgl, trim = FALSE, s = rule)
#' k_opt <- (length(drop(beta)) - (m + 1)) / p
#' #cv_cgl$Funcomp.CGL.fit[[as.character(k_opt)]]
#'
#'
#' par(mfrow=c(1,4))
#' plot(cv_cgl,k_list = k_list)
#' plot(cv_cgl$Funcomp.CGL.fit[[as.character(k_opt)]],
#'      p = p, k = k_opt, ylab = "L2")
#' plot(cv_cgl$Ftrim$cvm[k_opt - k_list[1] + 1, ])
#' title(paste0("k=", k_opt), line = 0.5)
#' # apply(cv_cgl$Funcomp.CGL.fit[[k_opt - 3]]$beta, 2,
#' #       function(x, p, k) {
#' #         #which(abs(x[seq(1, (p-1)*k+1, by = k)])>0),
#' #         (1:p)[apply(matrix(x[1:(p*k)], byrow = TRUE, nrow = p), 1,
#' #                     function(x) max(abs(x)) > 0)]
#' #       },p = p , k = k_opt)
#' if(k_opt == df_beta) {
#'   plot(Data$beta, col = "red", pch = 19,
#'        ylim = range(c(range(Data$beta), range(beta)))) #range(Data$beta))
#'   abline(v= seq(from = 0, to = (p*df_beta), by = df_beta ))
#'   abline(h = 0)
#'   points(beta)
#'   if(m1 > 0) points(p*df_beta + 1:m1, tail(Data$beta, m1),
#'                     col = "blue", pch = 19)
#' } else {
#'   plot(beta, ylim = range(c(range(Data$beta), range(beta))) )
#'   abline(v= seq(from = 0, to = (p*k_opt), by = k_opt ))
#'   abline(h = 0, col = "red")
#'   if(m1 > 0) points(p*k_opt + 1:m1, tail(Data$beta, m1),
#'                     col = "blue", pch = 19)
#' }
#' title(paste0("Method cgl"), outer=TRUE, line = -2)
#' par(mfrow=c(1,1))
#'
#' beta_C <- matrix(beta[1:(p*k_opt)], byrow = TRUE, nrow = p)
#' cat("colSums:", colSums(beta_C))
#' #Non.zero <- which(abs(beta_C[,1]) > 0)
#' Non.zero <- apply(beta_C, 1, function(x) ifelse(max(abs(x)) >0, TRUE, FALSE))
#' Non.zero <- (1:p)[Non.zero]
#' cat("None zero groups:", Non.zero)
#' #vet(beta, p = p, k = k_opt)
#'
#' par(mfrow=c(1,4))
#' plot(cv_cgl) #plot(cv_cgl,k_list = k_list)
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
#' title(paste0("k=", k_opt), line = 0.5)
#' title(paste0("Method cgl"), outer=TRUE, line = -2)
#' par(mfrow=c(1,1))
#' ##set a cutoff when you compute nonzeros
#' Non.zero <- apply(beta_C, 1, function(x)
#'                  ifelse(sqrt(sum(x^2)) > sqrt(sum(beta_C^2))/100, TRUE, FALSE))
#' Non.zero <- (1:p)[Non.zero]
#' Non.zero
#' Non.zero <- apply(beta_C, 1, function(x)
#'                  ifelse(sum(x^2) > sum(beta_C^2)/100, TRUE, FALSE))
#' Non.zero <- (1:p)[Non.zero]
#' Non.zero
#' ## cut by curve
#' Curve_L2 <- colSums(beta_curve^2)
#' Curve_L2 <- Curve_L2 - colSums(beta_curve[c(1, nrow(Data$basis.info)), ]^2) / 2
#' Curve_L2 <- Curve_L2 * (Data$basis.info[2,1] - Data$basis.info[1,1])
#' plot(Curve_L2)
#' which(sqrt(Curve_L2) > sqrt(sum(Curve_L2)) / 100)
#'
#' cgl <- list()
#' # MSE <- crossprod(y -  cbind2(cbind(cv_cgl$Funcomp.CGL.fit[[k_opt - k_list[1]+ 1]]$Z,
#' # Zc), 1) %*% beta) / length(y)
#' X_train <- cbind2(cbind(cv_cgl$Funcomp.CGL.fit[[k_opt - k_list[1]+ 1]]$Z, Zc), 1)
#' MSE <- sum((y -  X_train %*% beta)^2) / length(y)
#' # R_sqr <- 1 - crossprod(y -  cbind2(cbind(cv_cgl$Funcomp.CGL.fit[[k_opt - k_list[1]+ 1]]$Z,
#' # Zc), 1) %*% beta) / crossprod(y -  mean(y))
#' R_sqr <- sum((y -  X_train %*% beta)^2)
#' R_sqr <- 1 - R_sqr / crossprod(y -  mean(y))
#'
#' obj <- FuncompCGL(y = y_test, X = Test$data$Comp, k = k_opt, nlam = 1, outer_maxiter = 0)
#' # PE <- sum((Test$data$y - cbind2(cbind(obj$Z, Test$data$Zc), 1) %*% beta)^2)
#' # / length(drop(Test$data$y))
#' X_test <- cbind2(cbind(obj$Z, Test$data$Zc), 1)
#' PE <- sum((y_test - X_test %*% beta)^2) / length(y_test)
#' cgl$pred_error <- c(MSE = MSE, PE = PE, Rsqr_train = R_sqr)
#'
#' cgl$Non.zero_cut <- Non.zero
#' cgl <- c(cgl,
#'              ERROR_fun(beta_fit = beta, beta_true = Data$beta,
#'                        basis_fit = B, basis_true = Data$basis.info[,-1],
#'                        sseq = Data$basis.info[, 1],
#'                        m = m, p = p, Nzero_group = length(Non_zero.true), tol = 0),
#'              k = k_opt)
#' cgl$coef <- list(beta_C = beta_C, beta_c = tail(beta, m1))
#'
#' \dontrun{
#' # naive model
#' # set mu_raio = 0 to identifying without linear constraints,
#' # no outer_loop for Lagrange augmented multiplier
#' # mu_ratio = 0
#' cv_naive <-  cv.FuncompCGL(y = y, X = X, Zc = Zc, intercept = intercept,
#'                            W = rep(1, p), #W = function(x){ diag( 1 / apply(x, 2, sd) ) },
#'                            k = k_list, nfolds = 10, trim = 0,
#'                            tol = 0, inner_eps = 1e-6, inner_maxiter = 1E3,
#'                            dfmax = 30, lambda.factor = 1e-3,
#'                            mu_ratio = 0, outer_eps = 1e-6,
#'                            keep = FALSE, Trange = c(0,1)
#'                            #lam = c(exp(seq(log(0.01),log(0.00001),length=100)),0)
#' )
#'
#' plot(cv_naive,k_list = k_list)
#' cv_naive$Ftrim[c("lam.min", "lam.1se")]
#' beta <-  coef(cv_naive, trim = FALSE, s = rule)
#' k_opt <- (length(drop(beta)) - (m + 1)) / p
#' #cv_naive$Funcomp.CGL.fit[[as.character(k_opt)]]
#'
#'
#' par(mfrow=c(1,4))
#' plot(cv_naive,k_list = k_list)
#' plot(cv_naive$Funcomp.CGL.fit[[as.character(k_opt)]],
#'      p = p, k = k_opt, ylab = "L2")
#' plot(cv_naive$Ftrim$cvm[k_opt - k_list[1] + 1, ])
#' title(paste0("k=", k_opt), line = 0.5)
#' # apply(cv_naive$Funcomp.CGL.fit[[k_opt - 3]]$beta, 2,
#' #       function(x, p, k) {
#' #         #which(abs(x[seq(1, (p-1)*k+1, by = k)])>0),
#' #         (1:p)[apply(matrix(x[1:(p*k)], byrow = TRUE, nrow = p), 1,
#' #                     function(x) max(abs(x)) > 0)]
#' #       },p = p , k = k_opt)
#' if(k_opt == df_beta) {
#'   plot(Data$beta, col = "red", pch = 19,
#'        ylim = range(c(range(Data$beta), range(beta)))) #range(Data$beta))
#'   abline(v= seq(from = 0, to = (p*df_beta), by = df_beta ))
#'   abline(h = 0)
#'   points(beta)
#'   if(m1 > 0) points(p*df_beta + 1:m1, tail(Data$beta, m1),
#'                     col = "blue", pch = 19)
#' } else {
#'   plot(beta, ylim = range(c(range(Data$beta), range(beta))) )
#'   abline(v= seq(from = 0, to = (p*k_opt), by = k_opt ))
#'   abline(h = 0, col = "red")
#'   if(m1 > 0) points(p*k_opt + 1:m1, tail(Data$beta, m1),
#'                     col = "blue", pch = 19)
#' }
#' title(paste0("Method naive"), outer=TRUE, line = -2)
#' par(mfrow=c(1,1))
#'
#' beta_C <- matrix(beta[1:(p*k_opt)], byrow = TRUE, nrow = p)
#' cat("colSums:", colSums(beta_C))
#' #Non.zero <- which(abs(beta_C[,1]) > 0)
#' Non.zero <- apply(beta_C, 1, function(x) ifelse(max(abs(x)) >0, TRUE, FALSE))
#' Non.zero <- (1:p)[Non.zero]
#' cat("None zero groups:", Non.zero)
#' #vet(beta, p = p, k = k_opt)
#'
#' par(mfrow=c(1,4))
#' plot(cv_naive) #plot(cv_naive,k_list = k_list)
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
#' title(paste0("k=", k_opt), line = 0.5)
#' title(paste0("Method naive"), outer=TRUE, line = -2)
#' par(mfrow=c(1,1))
#' ##set a cutoff when you compute nonzeros
#' Non.zero <- apply(beta_C, 1, function(x)
#'                   ifelse(sqrt(sum(x^2)) > sqrt(sum(beta_C^2))/100, TRUE, FALSE))
#' Non.zero <- (1:p)[Non.zero]
#' Non.zero
#'
#'
#' naive <- list()
#' # MSE <- crossprod(y -  cbind2(cbind(cv_naive$Funcomp.CGL.fit[[k_opt - k_list[1]+ 1]]$Z,
#' # Zc), 1) %*% beta)
#' X_train <- cbind2(cbind(cv_naive$Funcomp.CGL.fit[[k_opt - k_list[1]+ 1]]$Z, Zc), 1)
#' MSE <- sum((y -  X_train %*% beta)^2)/ length(y)
#' # R_sqr <- 1 - crossprod(y -  cbind2(cbind(cv_naive$Funcomp.CGL.fit[[k_opt - k_list[1]+ 1]]$Z,
#' # Zc), 1) %*% beta) / crossprod(y -  mean(y))
#' R_sqr <- sum((y -  X_train %*% beta)^2)
#' R_sqr <- 1 - R_sqr / crossprod(y -  mean(y))
#'
#' obj <- FuncompCGL(y = Test$data$y, X = Test$data$Comp, k = k_opt, nlam = 1, outer_maxiter = 0)
#' # PE <- sum((Test$data$y - cbind2(cbind(obj$Z, Test$data$Zc), 1) %*% beta)^2)
#' # / length(drop(Test$data$y))
#' X_test <- cbind2(cbind(obj$Z, Test$data$Zc), 1)
#' PE <- sum((y_test - X_test %*% beta)^2) / length(y_test)
#' naive$pred_error <- c(MSE = MSE, PE = PE, Rsqr_train = R_sqr)
#' naive$Non.zero_cut <- Non.zero
#' naive <- c(naive,
#'            ERROR_fun(beta_fit = beta, beta_true = Data$beta,
#'                      basis_fit = B, basis_true = Data$basis.info[,-1],
#'                      sseq = Data$basis.info[, 1],
#'                      m = m, p = p, Nzero_group = length(Non_zero.true), tol = 0),
#'            k = k_opt)
#' naive$coef <- list(beta_C = beta_C, beta_c = tail(beta, m1))
#'
#'
#'
#' # log contract model
#' # set reference variable and once ref is set to a scalar in range,
#' # mu_ratio is set to 0 automatically
#' # ref = sample(1:p, 1)
#' cv_base <- cv.FuncompCGL(y = y, X = X, Zc = Zc, intercept = intercept, ref = sample(1:p, 1),
#'                          W = rep(1, p - 1), #W = function(x){ diag( 1 / apply(x, 2, sd) ) },
#'                          k = k_list, nfolds = 10, trim = 0,
#'                          tol = 0, inner_eps = 1e-6, inner_maxiter = 1E3,
#'                          dfmax = 30, lambda.factor = 1e-3,
#'                          mu_ratio = 0, outer_eps = 1e-6,
#'                          keep = FALSE, Trange = c(0,1)
#'                          #,lam = c(exp(seq(log(0.01),log(0.00001),length=100)),0)
#' )
#'
#' plot(cv_base,k_list = k_list)
#' cv_base$Ftrim[c("lam.min", "lam.1se")]
#' beta <-  coef(cv_base, trim = FALSE, s = rule)
#' k_opt <- (length(drop(beta)) - (m + 1)) / p
#' #cv_base$Funcomp.CGL.fit[[as.character(k_opt)]]
#'
#'
#' par(mfrow=c(1,4))
#' plot(cv_base,k_list = k_list)
#' plot(cv_base$Funcomp.CGL.fit[[as.character(k_opt)]],
#'      p = p, k = k_opt, ylab = "L2")
#' plot(cv_base$Ftrim$cvm[k_opt - k_list[1] + 1, ])
#' title(paste0("k=", k_opt), line = 0.5)
#' # apply(cv_base$Funcomp.CGL.fit[[k_opt - 3]]$beta, 2,
#' #       function(x, p, k) {
#' #         #which(abs(x[seq(1, (p-1)*k+1, by = k)])>0),
#' #         (1:p)[apply(matrix(x[1:(p*k)], byrow = TRUE, nrow = p), 1,
#' #                     function(x) max(abs(x)) > 0)]
#' #       },p = p , k = k_opt)
#' if(k_opt == df_beta) {
#'   plot(Data$beta, col = "red", pch = 19,
#'        ylim = range(c(range(Data$beta), range(beta)))) #range(Data$beta))
#'   abline(v= seq(from = 0, to = (p*df_beta), by = df_beta ))
#'   abline(h = 0)
#'   points(beta)
#'   if(m1 > 0) points(p*df_beta + 1:m1, tail(Data$beta, m1),
#'                     col = "blue", pch = 19)
#' } else {
#'   plot(beta, ylim = range(c(range(Data$beta), range(beta))) )
#'   abline(v= seq(from = 0, to = (p*k_opt), by = k_opt ))
#'   abline(h = 0, col = "red")
#'   if(m1 > 0) points(p*k_opt + 1:m1, tail(Data$beta, m1),
#'                     col = "blue", pch = 19)
#' }
#' title(paste0("Method baseline=", cv_base$Funcomp.CGL.fit[[1]]$ref), outer=TRUE, line = -2)
#' par(mfrow=c(1,1))
#'
#' beta_C <- matrix(beta[1:(p*k_opt)], byrow = TRUE, nrow = p)
#' cat("colSums:", colSums(beta_C))
#' #Non.zero <- which(abs(beta_C[,1]) > 0)
#' Non.zero <- apply(beta_C, 1, function(x) ifelse(max(abs(x)) >0, TRUE, FALSE))
#' Non.zero <- (1:p)[Non.zero]
#' cat("None zero groups:", Non.zero)
#' #vet(beta, p = p, k = k_opt)
#'
#' par(mfrow=c(1,4))
#' plot(cv_base) #plot(cv_base,k_list = k_list)
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
#' title(paste0("k=", k_opt), line = 0.5)
#' title(paste0("Method baseline=", cv_base$Funcomp.CGL.fit[[1]]$ref), outer=TRUE, line = -2)
#' par(mfrow=c(1,1))
#' ##set a cutoff when you compute nonzeros
#' Non.zero <- apply(beta_C, 1, function(x)
#'                  ifelse(sqrt(sum(x^2)) > sqrt(sum(beta_C^2))/100, TRUE, FALSE))
#' Non.zero <- (1:p)[Non.zero]
#' Non.zero
#'
#'
#' base <- list()
#' # MSE <- crossprod(y -  cbind2(cbind(cv_cgl$Funcomp.CGL.fit[[k_opt - k_list[1]+ 1]]$Z,
#' # Zc), 1) %*% beta) / length(y)
#' X_train <- cbind2(cbind(cv_cgl$Funcomp.CGL.fit[[k_opt - k_list[1]+ 1]]$Z, Zc), 1)
#' MSE <- crossprod(y -  X_train %*% beta) / length(y)
#' # R_sqr <- 1 - crossprod(y -  cbind2(cbind(cv_cgl$Funcomp.CGL.fit[[k_opt - k_list[1]+ 1]]$Z,
#' # Zc), 1) %*% beta) / crossprod(y -  mean(y))
#' R_sqr <- crossprod(y -  X_train %*% beta)
#' R_sqr <- 1 - R_sqr / crossprod(y -  mean(y))
#'
#' obj <- FuncompCGL(y = Test$data$y, X = Test$data$Comp, k = k_opt, nlam = 1, outer_maxiter = 0)
#' # PE <- sum((Test$data$y - cbind2(cbind(obj$Z, Test$data$Zc), 1) %*% beta)^2)
#' # / length(drop(Test$data$y))
#' X_test <- cbind2(cbind(obj$Z, Test$data$Zc), 1)
#' PE <- sum((y_test - X_test %*% beta)^2) / length(y_test)
#' base$pred_error <- c(MSE = MSE, PE = PE, Rsqr_train = R_sqr)
#'
#' base$Non.zero_cut <- Non.zero
#' base <- c(base,
#'           ERROR_fun(beta_fit = beta, beta_true = Data$beta,
#'                     basis_fit = B, basis_true = Data$basis.info[,-1],
#'                     sseq = Data$basis.info[, 1],
#'                     m = m, p = p, Nzero_group = length(Non_zero.true), tol = 0),
#'           k = k_opt)
#' base$coef <- list(beta_C = beta_C, beta_c = tail(beta, m1))
#'
#' }
#'
#' @export
#'


cv.FuncompCGL <- function(y, X, Zc = NULL, lam = NULL, ref = NULL,
                          W = rep(1,times = p - length(ref)),
                          k = 4:10,
                          foldid, nfolds = 10, nlam = 100,
                          trim = 0,outer_maxiter = 1e+6,
                          keep = FALSE, ...) {
  y <- drop(y)
  n <- length(y)

  if (missing(foldid)) foldid <- sample(rep(seq(nfolds), length = n)) else nfolds <- length(unique(foldid)) #max(foldid)
  if (nfolds < 3) stop("nfolds must be bigger than 3; nfolds=10 recommended")


  object <- as.list(seq(length(k))) # list of data for different k
  names(object) <- k
  if(!is.null(lam) || length(k) == 1) {

    # Case I
    if(dim(X)[1] == n ) p = dim(X)[2] / k else p <- dim(X)[2] - 2

    for(i in 1:length(k)){
      ###cat("1", k[i], "\r\n")
      object[[i]] <- FuncompCGL(y = y, X = X, Zc = Zc, lam = lam,
                                W = W, ref = ref,
                                k = k[i],
                                nlam = nlam, outer_maxiter = outer_maxiter,
                                ...)
    }
  } else {
    # Case II: find commom lambda first
    p <- ncol(X) - 2
    ###cat(p)
    ###cat(W)
    ###cat("length(W)", length(W), "\r\n")
    ###cat("is.missing(ref)", is.null(ref), "\r\n")

    for(i in 1:length(k)){
      ###cat("B", k[i], "\n")
      # Caculate integral Z and W matrix (if W is a functoin)

      object[[i]] <- FuncompCGL(y = y, X = X, Zc = Zc, nlam = 1,
                                W = W, ref = ref,
                                k = k[i], outer_maxiter = 0,
                                ...)
    }

    lam0 <- max(sapply(object, "[[", "lam")) # Shared lam0 for different k
    # Solution path for each df k
    dotlist = list(...)
    Trange_choice = dotlist$Trange
    for(i in 1:length(k)) {
      if(is.null(Trange_choice))
      {
      object[[i]] <- FuncompCGL(y = y, X = object[[i]]$Z, Zc = Zc,
                                W = object[[i]]$W, ref = ref,
                                lam = lam0, nlam = nlam,
                                outer_maxiter = outer_maxiter,
                                k = k[i], Trange = c(object[[i]]$Trange[1], object[[i]]$Trange[2]), ...)
      }
      else{
      object[[i]] <- FuncompCGL(y = y, X = object[[i]]$Z, Zc = Zc,
                                W = object[[i]]$W, ref = ref,
                                lam = lam0, nlam = nlam,
                                outer_maxiter = outer_maxiter,
                                k = k[i], ...)
      }
    }
  }



  cv.nlam <- sapply(object, function(x) length(drop(x$lam))) # different stopping lambda sequence by dfmax/pfmax
  cv.nlam_id <- which.min(cv.nlam)
  cv.lam <- object[[cv.nlam_id]]$lam # shared lambda sequence for different k
  cv.nlam <- cv.nlam[cv.nlam_id]


  cvm <- matrix(NA, nrow = length(k), ncol = max(cv.nlam))
  rownames(cvm) <- paste0("df=", k)
  colnames(cvm) <- seq(cv.nlam)
  cvsd <- cvm

  # deleting later
  if(trim > 0) {
    cvm.trim <- cvm
    cvsd.trim <- cvm
  }

  if(keep) {
    preval <- array(NA, dim = c(n, cv.nlam, length(k)))
  }

  # Cross-validatoin via different folders
  for(l in 1:length(k)) {

    outlist <- as.list(seq(nfolds))

    for(i in 1:nfolds) {
      ###cat("folds", i, "\r\n")
      which <- foldid == i
      y_train <- y[!which]
      Z <- object[[l]]$Z
      Z_train <- Z[!which, , drop = FALSE]
      Zc_train <- Zc[!which, , drop = FALSE]
      outlist[[i]] <- FuncompCGL(y = y_train, X = Z_train, Zc = Zc_train,
                                 k = k[l],
                                 lam = cv.lam, W = object[[l]]$W, ref = ref,
                                 outer_maxiter = outer_maxiter,
                                 ...)
    }


    cvstuff <- cv.test(outlist, y, X = cbind(Z, Zc), foldid, lam = cv.lam, trim = trim, keep = keep) #define diffrent cv.test for GLM
    cvm[l, ] <- cvstuff$cvm
    cvsd[l, ] <- cvstuff$cvsd
    if(keep) preval[, , l] <- cvstuff$fit.preval

    # delet later
    if(trim > 0) {
      cvm.trim[l, ] <- cvstuff$cvmtrim
      cvsd.trim[l, ] <- cvstuff$cvsdtrim
    }

  }


  # Select lambda
  Ftrim = list(cvm = cvm, cvsd = cvsd)
  Ftrim$cvup = Ftrim$cvm + Ftrim$cvsd
  Ftrim$cvlo = Ftrim$cvm - Ftrim$cvsd
  lammin <- ggetmin(lam = cv.lam, cvm = Ftrim$cvm, cvsd = Ftrim$cvsd, k_list = k)
  Ftrim <- c(Ftrim, lammin)

  # Delect later
  if(trim > 0) {
    Ttrim <- list(cvm = cvm.trim, cvsd = cvsd.trim)
    Ttrim$cvup = Ttrim$cvm + Ttrim$cvsd
    Ttrim$cvlo = Ttrim$cvm - Ttrim$cvsd
    lammin <- ggetmin(lam = cv.lam, cvm = Ttrim$cvm, cvsd = Ttrim$cvsd, k_list = k)
    Ttrim <- c(Ttrim, lammin)
  } else {
    Ttrim <- NULL
  }

  result <- list(Funcomp.CGL.fit = object,
                 lam = cv.lam,
                 Ftrim = Ftrim,
                 Ttrim = Ttrim)
  if(keep) result <- c(result, list(fit.preval = preval, foldid = foldid))
  class(result) <- "cv.FuncompCGL"
  return(result)
}



# cv.FuncompCGL <- function(y, X, Zc = NULL, lam = NULL, W = rep(1, p),
#                           T.name = "TIME", ID.name = "Subject_ID",
#                           k = 4:10,
#                           foldid, nfolds = 10, nlam = 100,
#                           trim = 0.05,
#                           outer_maxiter = 1e+6, keep = FALSE,
#                           ...) {
#   y <- drop(y)
#   n <- length(y)
#   #p <- length(X.names)
#
#   if (missing(foldid)) foldid <- sample(rep(seq(nfolds), length = n)) else nfolds <- max(foldid)
#   if (nfolds < 3) stop("nfolds must be bigger than 3; nfolds=10 recommended")
#
#
#   object <- as.list(seq(length(k)))
#   names(object) <- k
#   if(!is.null(lam) || length(k) == 1) {
#
#     if(dim(X)[1] == n ) p = dim(X)[2] / k else p <- dim(X)[2] - 2
#
#     for(i in 1:length(k)){
#       object[[i]] <- FuncompCGL(y = y, X = X, Zc = Zc, lam = lam, W = W,
#                                 k = k[i],
#                                 T.name = T.name, ID.name = ID.name,
#                                 nlam = nlam, outer_maxiter = outer_maxiter,
#                                 ...)
#     }
#   } else {
#     p <- dim(X)[2] - 2
#     for(i in 1:length(k)){
#
#       object[[i]] <- FuncompCGL(y = y, X = X, Zc = Zc, nlam = 1, W = W,
#                                 k = k[i], outer_maxiter = 0,
#                                 T.name = T.name, ID.name = ID.name, ...)
#     }
#
#     lam0 <- max(sapply(object, "[[", "lam"))
#
#     for(i in 1:length(k)){
#       object[[i]] <- FuncompCGL(y = y, X = object[[i]]$Z, Zc = Zc, lam = lam0, W = object[[i]]$W,
#                                 T.name = T.name, ID.name = ID.name, nlam = nlam, outer_maxiter = outer_maxiter,
#                                 k = k[i], ...)
#     }
#   }
#
#
#
#   cv.nlam <- sapply(object, function(x) length(drop(x$lam)))
#   cv.nlam_id <- which.min(cv.nlam)
#   cv.lam <- object[[cv.nlam_id]]$lam
#   cv.nlam <- cv.nlam[cv.nlam_id]
#
#
#   cvm <- matrix(NA, nrow = length(k), ncol = max(cv.nlam))
#   rownames(cvm) <- paste0("df=", k)
#   colnames(cvm) <- seq(cv.nlam)
#   cvsd <- cvm
#   if(trim > 0) {
#     cvm.trim <- cvm
#     cvsd.trim <- cvm
#   }
#
#   if(keep) {
#     preval <- array(NA, dim = c(n, cv.nlam, length(k)))
#   }
#
#   for(l in 1:length(k)) {
#
#     outlist <- as.list(seq(nfolds))
#
#     for(i in 1:nfolds) {
#       #cat("folds", i, "\r\n")
#       which <- foldid == i
#       y_train <- y[!which]
#       Z <- object[[l]]$Z
#       Z_train <- Z[!which, , drop = FALSE]
#       Zc_train <- Zc[!which, , drop = FALSE]
#       outlist[[i]] <- FuncompCGL(y = y_train, X = Z_train, Zc = Zc_train, k = k[l], lam = cv.lam, W = object[[l]]$W,
#                                  T.name = T.name, ID.name = ID.name,outer_maxiter = outer_maxiter,
#                                  ...)
#     }
#
#     #X <- cbind(Z, Zc)
#     cvstuff <- cv.test(outlist, y, X = cbind(Z, Zc), foldid, lam = cv.lam, trim = trim, keep = keep)
#     cvm[l, ] <- cvstuff$cvm
#     cvsd[l, ] <- cvstuff$cvsd
#     if(keep) preval[, , l] <- cvstuff$fit.preval
#     if(trim > 0) {
#       cvm.trim[l, ] <- cvstuff$cvmtrim
#       cvsd.trim[l, ] <- cvstuff$cvsdtrim
#     }
#
#   }
#
#   # mincvm_crossk <- apply(cvm, 1, min, na.rm = TRUE)
#   # k_select <- which.min(mincvm_crossk)
#   # lam_select <- which.min(cvm[k_select, ])
#   # k_select <- k[k_select]
#   # lam_select <- cv.lam[lam_select]
#   # Ftrim = list(cvm = drop(cvm), cvsd = drop(cvsd)
#   #              #, lam.min = lam_select, k = k_select
#   #              )
#
#   Ftrim = list(cvm = cvm, cvsd = cvsd
#                #, lam.min = lam_select, k = k_select
#   )
#   Ftrim$cvup = Ftrim$cvm + Ftrim$cvsd
#   Ftrim$cvlo = Ftrim$cvm - Ftrim$cvsd
#   lammin <- ggetmin(lam = cv.lam, cvm = Ftrim$cvm, cvsd = Ftrim$cvsd, k_list = k)
#   Ftrim <- c(Ftrim, lammin)
#   #selection <- as.list(1:2)
#   #names(selection) <- c("lam.min", "lam.1sd")
#   #selection$lam.min <- aslist()
#
#   if(trim > 0) {
#     # mincvm_crossk <- apply(cvm.trim, 1, min, na.rm = TRUE)
#     # k_select <- which.min(mincvm_crossk)
#     # lam_select <- which.min(cvm.trim[k_select, ])
#     # k_select <- k[k_select]
#     # lam_select <- cv.lam[lam_select]
#
#     # Ttrim <- list(cvm = drop(cvm.trim), cvsd = drop(cvsd.trim)
#     #               #, lam.min = lam_select, k = k_select
#     #               )
#     #
#     Ttrim <- list(cvm = cvm.trim, cvsd = cvsd.trim
#                   #, lam.min = lam_select, k = k_select
#     )
#
#     Ttrim$cvup = Ttrim$cvm + Ttrim$cvsd
#     Ttrim$cvlo = Ttrim$cvm - Ttrim$cvsd
#     lammin <- ggetmin(lam = cv.lam, cvm = Ttrim$cvm, cvsd = Ttrim$cvsd, k_list = k)
#     Ttrim <- c(Ttrim, lammin)
#   } else {
#     Ttrim <- NULL
#   }
#
#   result <- list(Funcomp.CGL.fit = object,
#                  lam = cv.lam,
#                  Ftrim = Ftrim,
#                  Ttrim = Ttrim
#                  #,foldid = foldid
#                  )
#   if(keep) result <- c(result, list(fit.preval = preval, foldid = foldid))
#   class(result) <- "cv.FuncompCGL"
#   return(result)
# }




#' @title
#' Cross-validation for compCL
#'
#' @description
#' Does nfolds cross-validation for compCL, return value of \code{lam}.
#' The function is modified based on the \code{cv} function from \code{glmnet} package
#'
# @usage
# cv.compCL <- function(y, Z, Zc = NULL, intercept = FALSE,
#                       lam = NULL, nfolds = 10, foldid,
#                       trim = 0.1, ...)
#'
#' @inheritParams compCL
#'
#' @param nfolds number of folds - default is 10. Smallest value allowable is nfolds=3.
#'
#' @param foldid an optional vector of values between 1 and \code{nfolds} identifying what fold each
#'               observation is in. If supplied, \code{nfold} can be missing.
#'
#' @param trim a scaler specifying percentage to be trimmed off for prediction error - default is 0.
#' @param keep If \code{keep=TRUE}, a prevalidated array is returned containing fitted values for each observation,
#'            of lambda. This means these fits are computed with this observation and the rest of
#'            its fold omitted. Default is keep=FALSE
#'
#' @param ... other arguments that can be passed to compCL.
#'
#' @return an object of class \code{\link{cv.compCL}} is returned.
#' \item{compCL.fit}{a fitted \code{\link{compCL}} object for the full data}
#' \item{lam}{the values of \code{lam} used in the fits}
#' \item{Ftrim}{a list of cross-validation result without trimming.
#'                \itemize{
#'                \item \code{cvm} the mean cross-validated error without trimming -  a vector of  \code{length(lam)}
#'                \item \code{cvsd} estimate of standard error of cvm without trimming- a vector of  \code{llength(lam)}
#'                \item \code{cvupper} upper curve = \code{cvm+cvsd}.
#'                \item \code{cvlo} lower curve = \code{cvm-cvsd}.
#'                \item \code{lam.min} The optimal value of \code{lam} that gives minimum cross validation error \code{cvm}
#'                \item \code{lam.1se} The largest value of lam such that error is within 1 standard error of the minimum \code{cvm}
#'                }
#'             }
#'
#' \item{Ttrim}{a list of cross-validation result with \code{trim*100\%}, if provided, of tails trimmed off for cross validation error.
#'                \itemize{
#'                \item \code{cvm} the mean cross-validated error with with \code{trim*100\%} trimmed - a vector of  \code{length(lam)}
#'                \item \code{cvsd} estimate of standard error of cvm with \code{trim*100\%} trimmed - a vector of  \code{length(lam)}
#'                \item \code{cvupper} upper curve = \code{cvm+cvsd}.
#'                \item \code{cvlo} lower curve = \code{cvm-cvsd}.
#'                \item \code{lam.min} The optimal value of \code{lam} that gives minimum cross validation error cvm with \code{trim*100\%} trimmed
#'                \item \code{lam.1se} The largest value of lam such that error is within 1 standard error of the minimum \code{cvm} after \code{trim*100\%} trimmed.
#'                }
#'             }
#'
#' \item{foldid}{the values of \code{folidi} used in fits.}
#'
#'
#' @examples
#'
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c( beta, rep(0, times = p - length(beta)) )
#' Comp_data = comp_Model(n = n, p = p,
#'                             rho = 0.2, sigma = 0.5,
#'                             gamma  = 0.5, add.on = 1:5,
#'                             beta = beta, intercept = FALSE)
#' Comp_data$Zc
#' cvm <- cv.compCL(y = Comp_data$y,
#'                  Z = Comp_data$X.comp, Zc = Comp_data$Zc,
#'                  intercept = Comp_data$intercept,
#'                  lam = NULL, nfolds = 10, trim = 0.05, lambda.factor = 0.0001,
#'                  dfmax = p, mu_ratio = 1, outer_eps = 1e-10, inner_eps = 1e-8, inner_maxiter = 1e4)
#'
#' plot(cvm)
#' coef(cvm, s = "lam.min")
#' cvm$compCL.fit
#' which(abs(coef(cvm, s = "lam.min")$beta[1:p]) > 0)
#' which(abs(coef(cvm, s= "lam.1se")$beta[1:p]) > 0)
#'
#' @export
#'


# cv.compCL <- function(y, Z, Zc = NULL, intercept = FALSE,
#                       lam = NULL,
#                       nfolds = 10, foldid,
#                       trim = 0.1, ...) {
#   y <- drop(y)
#   n <- length(y)
#   if (missing(foldid)) foldid <- sample(rep(seq(nfolds), length = n)) else nfolds <- max(foldid)
#   if (nfolds < 3) stop("nfolds must be bigger than 3; nfolds=10 recommended")

#   compCL.object <- compCL(y = y, Z = Z, Zc = Zc, intercept = intercept,
#                           lam = lam, ...)
#   cv.lam <- compCL.object$lam

#   outlist <- as.list(seq(nfolds))

#   ###Now fit the nfold models and store them
#   for (i in seq(nfolds)) {
#     which <- foldid == i
#     outlist[[i]] <- suppressMessages(compCL(y = y[!which], Z = Z[!which, , drop = FALSE], Zc = Zc[!which, , drop = FALSE], intercept = intercept,
#                            lam = cv.lam,  ...))
#   }

#   #cvstuff <- cv.test2(outlist, y, Z = compCL.object$Z_log, Zc = Zc, foldid, lam = cv.lam, trim = trim)
#   cvstuff <- suppressMessages(cv.test2(outlist, y, Z = Z, Zc = Zc, foldid, lam = cv.lam, trim = trim))
#   cvm <- drop(cvstuff$cvm)
#   cvsd <- drop(cvstuff$cvsd)
#   if(trim > 0) {
#     cvm.trim <- drop(cvstuff$cvmtrim)
#     cvsd.trim <- drop(cvstuff$cvsdtrim)
#   }


#   Ftrim = list(cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, culo = cvm - cvsd)
#   lam.min <- getmin(lam = cv.lam, cvm = cvm, cvsd = cvsd)
#   Ftrim <- c(Ftrim, lam.min)

#   if(trim > 0) {
#     Ttrim <- list(cvm = cvm.trim, cvsd = cvsd.trim, cvupper = cvm.trim + cvsd.trim,
#                   cvlo = cvm.trim-cvsd.trim)
#     lam.min <- getmin(lam = cv.lam, cvm = cvm.trim, cvsd = cvsd.trim)
#     Ttrim <- c(Ttrim, lam.min)
#   } else {
#     Ttrim <- NULL
#   }

#   result <- list(compCL.fit = compCL.object,
#                  lam = cv.lam,
#                  Ftrim = Ftrim,
#                  Ttrim = Ttrim,
#                  foldid = foldid)


#   class(result) <- "cv.compCL"
#   return(result)
# }


cv.compCL <- function(y, Z, Zc = NULL, intercept = FALSE,
                       lam = NULL,
                       nfolds = 10, foldid,
                       trim = 0.1, keep = FALSE,...) {
  y <- drop(y)
  n <- length(y)
  if (missing(foldid)) foldid <- sample(rep(seq(nfolds), length = n)) else nfolds <- max(foldid)
  if (nfolds < 3) stop("nfolds must be bigger than 3; nfolds=10 recommended")

  compCL.object <- compCL(y = y, Z = Z, Zc = Zc, intercept = intercept,
                          lam = lam, ...)
  cv.lam <- compCL.object$lam

  outlist <- as.list(seq(nfolds))

  ###Now fit the nfold models and store them
  for (i in seq(nfolds)) {
    which <- foldid == i
    outlist[[i]] <- suppressMessages((compCL(y = y[!which], Z = Z[!which, , drop = FALSE], Zc = Zc[!which, , drop = FALSE], intercept = intercept,
                           lam = cv.lam,  ...)))
  }

  newx <- cbind(compCL.object$Z_log, Zc)
  cvstuff <- cv.test(outlist, y, X = newx, foldid, lam = cv.lam, trim = trim, keep = keep)
  #cvstuff <- cv.test(outlist, y, newx, foldid, lam = cv.lam, trim = trim)
  cvm <- drop(cvstuff$cvm)
  cvsd <- drop(cvstuff$cvsd)
  if(trim > 0) {
    cvm.trim <- drop(cvstuff$cvmtrim)
    cvsd.trim <- drop(cvstuff$cvsdtrim)
  }


  Ftrim = list(cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, cvlo = cvm - cvsd)
  lam.min <- getmin(lam = cv.lam, cvm = cvm, cvsd = cvsd)
  Ftrim <- c(Ftrim, lam.min)

  if(trim > 0) {
    Ttrim <- list(cvm = cvm.trim, cvsd = cvsd.trim, cvupper = cvm.trim + cvsd.trim,
                  cvlo = cvm.trim-cvsd.trim)
    lam.min <- getmin(lam = cv.lam, cvm = cvm.trim, cvsd = cvsd.trim)
    Ttrim <- c(Ttrim, lam.min)
  } else {
    Ttrim <- NULL
  }

  result <- list(compCL.fit = compCL.object,
                 lam = cv.lam,
                 Ftrim = Ftrim,
                 Ttrim = Ttrim,
                 foldid = foldid)


  class(result) <- "cv.compCL"
  return(result)
}




