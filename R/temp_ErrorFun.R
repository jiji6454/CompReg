
#' error function
#'
#' calculate classification, prediction and accuracy errors for fitted
#' \code{FuncompCGL} model in simulation.
#'
#'
#' @param beta_fit fitted coefficients vector
#' @param beta_true true coefficientse vector
#' @param basis_fit basis matrix with selected df \code{k} and \code{basis_fun}
#' @param basis_true basis matrix with true df \code{k} and \code{basis_fun}
#' @param sseq sequence of time points used to generate basis matrix
#' @param m number of time-invariant variables
#' @param p number of compositional varialbes
#' @param Nzero_group number of true None zero group
#' @param tol tolerance
#'
#' @return a list
#' \item{Non.zero}{None zero rows in coefficient matrix \code{beta_C}. }
#' \item{class_error}{classification errors, including FP, FN, FPR, FNR, etc.}
#' \item{coef_error}{\itemize{
#'                            \item \code{L2.cont} L2 norm for time-invariant coefficients error
#'                            \item \code{L2.comp, L2_inf.comp, L2_L1.comp} if selected df \code{k} is
#'                                  the same as true degree freedom of basis, coefficients matrix error is
#'                                  calculated}
#'                  }
#' \item{curve_error}{norms for curve error by function norms}
#' \item{group_norm}{L2 norm for rows in coefficients matrix \code{beta_C}}
#'
#'
#' @export
#'
#'

ERROR_fun <- function(beta_fit, beta_true,
                      basis_fit, basis_true, sseq,
                      m, p, Nzero_group, tol = 0) {

  error.list <- list()

  pos_group <- 1:Nzero_group
  neg_group <- (1:p)[-pos_group]
  df_fit <- (length(beta_fit) - 1 - m) / p
  df_true <- (length(beta_true) -1 -m) / p
  #cat(df_fit, "\r\n")
  Non.zero_select <- Nzero(beta = beta_fit, p = p, k = df_fit, tol = tol)
  error.list$Non.zero <- Non.zero_select
  error.list$class_error <- Class_error(pos_select = Non.zero_select,
                                        pos_group = pos_group, neg_group = neg_group)

  error.list$coef_error <- Coef_error(coef_true = beta_true, coef_esti = beta_fit,
                                      p = p, k_fit = df_fit, k_true = df_true)


  error.list$curve_error <- Curve_error(coef_true = beta_true, coef_esti = beta_fit, p = p,
                                        k_fit = df_fit, k_true = df_true,
                                        basis_true = basis_true, basis_fit = basis_fit, sseq = sseq)


  error.list$group_norm <- apply(vet(beta = beta_fit, p = p, k = df_fit)$C, 1, function(x, X) sqrt(sum(x^2)))
  return(error.list)
}






#
# curve error norm
#
# @description
# calculate function norm for curves, L1, L2, L1_1, L1_inf,
# L2_1, L2_inf norms. If information for true beta is not provides,
# then calculate these norms for estimated coefficient beta.
#
# @param coef_true true beta vector
# @param coef_esti estimated beta vector
# @param k_fit selected \code{df} of basis
# @param k_true true \code{df} for beta
# @inheritParams ERROR_fun
#
# @export
#
Curve_error <- function(coef_true, coef_esti, p,
                        k_fit, k_true, basis_true, basis_fit, sseq) {

  coef_esti.comp <- vet(beta = coef_esti, p = p, k = k_fit)$C
  curve_esti <- coef_esti.comp %*% t(basis_fit)

  if(missing(coef_true) || missing(basis_true)){
    curve_true <- matrix(0, nrow = nrow(curve_esti), ncol = ncol(curve_esti))
  } else {
    coef_true.comp <- vet(beta = coef_true, p = p, k = k_true)$C
    curve_true <- coef_true.comp %*% t(basis_true)
  }

  curve_diff <- abs(curve_esti - curve_true) ## p by length(sseq)

  ns <- length(sseq)
  time_diff <- sseq[2] - sseq[1]
  extra_sum <- rowSums(curve_diff[, c(1, ns)]^2) * time_diff / 2## crossprod(curve_diff[, c(1, ns)]) * time_diff/2
  add_sum <- apply(curve_diff, 1, function(x) sum(x^2)) * time_diff
  ITG <- add_sum - extra_sum
  L2 <- sqrt(ITG)
  L2_L1 <- sum(L2)
  L2_inf <- max(L2)


  extra_sum <- rowSums(curve_diff[, c(1, ns)]) * time_diff / 2## crossprod(curve_diff[, c(1, ns)]) * time_diff/2
  add_sum <- rowSums(curve_diff) * time_diff
  ITG <- add_sum - extra_sum
  L1 <- ITG
  L1_L1 <- sum(L1)
  L1_inf <- max(L1)

  L1_each.max <- max(curve_diff)
  if(missing(coef_true) || missing(basis_true)){
    error <- c(L1 = L1, L2 = L2,
               L1_each.max = L1_each.max,
               L1_inf = L1_inf, L1_L1 = L1_L1,
               L2_inf = L2_inf, L2_L1 = L2_L1 )
  } else {
    error <- c(L1_each.max = L1_each.max,
               L1_inf = L1_inf, L1_L1 = L1_L1,
               L2_inf = L2_inf, L2_L1 = L2_L1 )
  }

  return(error)
}



##classification error
Class_error <- function(pos_select, pos_group, neg_group) {
  pos <- length(pos_group)
  neg <- length(neg_group)
  p <- pos + neg
  FP_select <- pos_select[which( pos_select %in% neg_group )]
  FN_select <- pos_group[which( !(pos_group %in% pos_select) )]
  TP_select <- pos_select[which( pos_select %in% pos_group )]
  TN_select <- neg_group[which( !(neg_group %in% pos_select) )]

  FP <- length(FP_select)
  FN <- length(FN_select)
  TP <- length(TP_select)
  TN <- length(TN_select)

  FPR <- FP / neg # false positive rate
  FNR <- FN / pos # false negative rate
  Sensitivity <- TP / pos # true positive rate
  Specificity <- TN / neg # true negative rate
  PPV <- TP / (TP + FP)   # precision/positive predictive value
  NPV <- TN/ (TN + FN)    # negative predictive value
  FDR <- 1 - PPV          # false discovery rate
  ACC <- (TP + TN) / p
  F1_score <- (2 * TP) / (2 * TP + FP + FN)
  MCC <- (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) # Matthews correlation coefficient
  Informedness <- Sensitivity + Specificity - 1 # Youden's J statistic
  Markedness <- PPV + NPV - 1

  error <- c(FP = FP, FN = FN, TP = TP, TN = TN,
             FPR = FPR, FNR = FNR, Sensi = Sensitivity, Speci = Specificity,
             PPV = PPV, NPV = NPV, FDR = FDR, ACC = ACC, F1_score = F1_score,
             MCC = MCC, Informedness = Informedness, Markedness = Markedness)

  return(error)
}

##Coefficients estimation Norm error
Coef_error <- function(coef_true, coef_esti, p, k_fit, k_true) {
  coef_esti.comp <- vet(beta = coef_esti, p = p, k = k_fit)$C
  coef_esti.cont <- vet(beta = coef_esti, p = p, k = k_fit)$b
  coef_true.comp <- vet(beta = coef_true, p = p, k = k_true)$C
  coef_true.cont <- vet(beta = coef_true, p = p, k = k_true)$b
  #cat("1")

  error <- vector()
  error <- c(L2.cont = sqrt(sum( (coef_esti.cont - coef_true.cont)^2 )) )
  if(k_fit == k_true) {
    L2_diff.beta <- sqrt(sum( (coef_true - coef_esti)^2 ))
    L2_diff.comp <- apply(coef_esti.comp - coef_true.comp, 1,
                          function(x) sqrt(sum(x^2)) )
    L2_diff_inf.comp <- max(L2_diff.comp)
    L2_diff_L1.comp <- sum(L2_diff.comp)
    L2_diff.comp <- sqrt(sum( (as.vector(coef_esti.comp - coef_true.comp))^2 ))
    error <- c(error, L2.beta = L2_diff.beta, L2.comp = L2_diff.comp,
               L2_inf.comp = L2_diff_inf.comp, L2_L1.comp = L2_diff_L1.comp)
  }
  #cat('2')
  return(error)
}




#' @title
#' Basis matrix generation
#' @description
#' generate the basis matrix for polynomial spline, fourier, and orthogonalsplinebasis B-splines.
#'
#' @param sseq the predictor variable.
#' @param df degree of freedom.
#' @param degree degree of piecewise polynomial - default is 3. Used in \code{bs} and \code{OBasis}.
#' @param type type of basis
#'
#' @return basis matrix of size \code{length(sseq)} by \code{df}
#'
#' @details
#' For type \code{bs} and \code{OBasis}, \code{intercept} is set to \code{TRUE}, minimum value for
#' \code{df} is 4.
#' For type \code{fourier}, \code{df} should be a even number. If not, \code{df} is automatically
#' set to \code{df} - 1.
#'
#' @export

Genbasis <- function(sseq, df, degree, type = c("bs", "fourier", "OBasis")) {

  type <- match.arg(type)
  interval <- range(sseq)

  nknots <- df - (degree + 1)
  if(nknots > 0) {
    knots <- ((1:nknots) / (nknots + 1)) * diff(interval) + interval[1]
  } else {
    knots <- NULL
  }

  basis <- switch(type,
                  "bs" = bs(x = sseq, df =  df, degree = degree,
                            Boundary.knots = interval, intercept = TRUE),

                  "fourier" = eval.basis(sseq,
                                         basisobj = create.fourier.basis(rangeval = interval, nbasis = df)
                  ),

                  "OBasis" = evaluate(OBasis(expand.knots(c(interval[1], knots, interval[2])),
                                             order = degree + 1),
                                      sseq)
  )
  return(basis)
}

vet <- function(beta, p, k) {
  p1 <- p * k
  coef <- matrix(beta[1:p1], byrow = TRUE, nrow = p)
  result<- list(C = coef)
  coef <- beta[(p1+1):length(beta)]
  result$b <- coef
  return(result)
}

Nzero <- function(beta, p, k, tol = 0) {
  coef <- vet(beta, p = p, k = k)$C
  group <- apply(coef, 1, function(x) ifelse(max(abs(x)) > tol , TRUE, FALSE))
  group <- (1:p)[group]
  return(group)
}



# Nzero2 <- function(beta, p, k, tol) {
#   coef <- vet(beta, p = p, k = k)$C
#   group <- apply(coef, 1, function(x, X) sqrt(sum(x^2)))
#   group <- ifelse(group > sqrt(sum(coef^2)) /100, TRUE, FALSE)
#   group <- (1:p)[group]
#   return(group)
# }
















#
# W_fun <- function(x){
#   x <- apply(x, 2, sd)
#   x <- 1/x
#   x <- diag(x)
#   return(x)
# }




# Modelplain <- function(n,p,df_beta = 5, beta_C_matrix, sigma, ns = 100, SNR, obs_spar) {
#
#   if(missing(beta_C_matrix)){
#
#     #cat(1)
#     beta_C = matrix(0, nrow = p, ncol = df_beta)
#     # beta_C[1, ] <- rep(0.7, times = df_beta)
#     # beta_C[2, ] <- rep(0.5, times = df_beta)
#     # beta_C[3, ] <- rep(-0.8, times = df_beta)
#     # beta_C[4, ] <- rep(-0.6, times = df_beta)
#     beta_C[3, ] <- c(-0.8, -0.8 , 0.4 , 1 , 1)
#     beta_C[4, ] <- c(0.5, 0.5, -0.6  ,-0.6, -0.6)
#     beta_C[1, ] <- c(-0.5, -0.5, -0.5 , -1, -1)
#     beta_C[2, ] <- c(0.8, 0.8,  0.7,  0.6,  0.6)
#     #cat("Column Sums of beta equal 0's: ", all.equal(drop(colSums(beta_C)), rep(0, times = df_beta)), "\r\n")
#     Nzero_group = 4
#     beta_C_matrix = beta_C
#     beta_C = as.vector(t(beta_C))
#   } else {
#     cat(2)
#     Nzero_group = which( abs(beta_C_matrix[, 1]) > 0)
#     beta_C <- as.vector(t(beta_C_matrix))
#   }
#
#   XX <- array(NA, dim = c(n,p, ns))
#   for(s in 1:100){
#     X <- matrix(rnorm(n*p, mean = 0, sd = sigma), nrow = n)
#     #X[X==0] <- 0.5
#     X <- exp(X)
#     X <- X / rowSums(X)
#     XX[, , s] <- X
#   }
#
#   sseq <- round(seq(from = 0, to = 1, length.out = ns ), 5)
#   BS <- splines::bs(sseq, df = df_beta, intercept = TRUE)
#   beta <- BS %*% t(beta_C_matrix)
#   D <- plyr::alply(XX, .margins = 1,
#                    function(x, sseq) data.frame(t(rbind(TIME = sseq, x))),
#                    sseq = sseq )
#
#   Z_ITG <- sapply(D, function(x, sseq, beta.basis, insert, method){
#     x[, -1] <- log(x[, -1])
#     return(ITG(X = x, basis = beta.basis, sseq = sseq, insert = insert, method = method)$integral)
#   } ,sseq = sseq, beta.basis = BS, insert = "FALSE", method = "trapezoidal")
#
#   Z_ITG <- t(Z_ITG)
#
#   Z_t.full <- plyr::ldply(D[1:n], data.frame, .id = "Subject_ID")
#
#   Y.tru <- Z_ITG %*% as.vector(t(beta_C_matrix))
#   error <- rnorm(n, 0, 1)
#   sigma_e <- sd(Y.tru) / (sd(error) * SNR)
#   error <- sigma_e*error
#   Y.obs <- Y.tru + error
#
#   Z_t.obs <- lapply(D, function(x, obs_spar) {
#     n <- dim(x)[1]
#     #lambda <- obs_spar * n
#     #n.obs <- rpois(1, lambda)
#     #T.obs <- sample(n, size = n.obs)
#     T.obs <- replicate(n, sample(c(0,1), 1, prob = c(1 - obs_spar, obs_spar)))
#     T.obs <- which(T.obs == 1)
#     #T.obs = sort(sample(seq(n), round(obs_spar*n)))
#     x.obs <- x[T.obs, ]
#     return(x.obs)
#   }, obs_spar = obs_spar)
#
#   Z_t.obs <- plyr::ldply(Z_t.obs, data.frame, .id = "Subject_ID")
#
#
#   data <- list(y = Y.obs, Comp = Z_t.obs
#                #, Zc = Z_control, intercept = intercept
#   )
#   beta <- c(beta_C, 0)#, beta_c, ifelse(intercept, beta0, 0))
#   data.raw <- list(Z_t.full = Z_t.full, Z_ITG = Z_ITG,
#                    Y.tru = Y.tru)
#   basis.info <- cbind(sseq, BS)
#
#   output <- list(data = data, beta = beta, basis.info = basis.info, data.raw = data.raw#,
#                  #parameter = parameter
#   )
#   #### Output ####
#
#   return(output)
#
#
# }
#
#
#
# Modelplain2 <- function(n,p,df_beta = 5, beta_C_matrix, sigma, ns = 100, SNR, obs_spar) {
#
#   if(missing(beta_C_matrix)){
#
#     #cat(1)
#     beta_C = matrix(0, nrow = p, ncol = df_beta)
#     # beta_C[1, ] <- rep(0.7, times = df_beta)
#     # beta_C[2, ] <- rep(0.5, times = df_beta)
#     # beta_C[3, ] <- rep(-0.8, times = df_beta)
#     # beta_C[4, ] <- rep(-0.6, times = df_beta)
#     beta_C[3, ] <- c(-0.8, -0.8 , 0.4 , 1 , 1)
#     beta_C[4, ] <- c(0.5, 0.5, -0.6  ,-0.6, -0.6)
#     beta_C[1, ] <- c(-0.5, -0.5, -0.5 , -1, -1)
#     beta_C[2, ] <- c(0.8, 0.8,  0.7,  0.6,  0.6)
#     beta_C = beta_C * 10
#     #cat("Column Sums of beta equal 0's: ", all.equal(drop(colSums(beta_C)), rep(0, times = df_beta)), "\r\n")
#     Nzero_group = 4
#     beta_C_matrix = beta_C
#     beta_C = as.vector(t(beta_C))
#   } else {
#     cat(2)
#     Nzero_group = which( abs(beta_C_matrix[, 1]) > 0)
#     beta_C <- as.vector(t(beta_C_matrix))
#   }
#
#   XX <- array(NA, dim = c(n,p, ns))
#   for(s in 1:100){
#     X <- matrix(rnorm(n*p, mean = 0, sd = sigma), nrow = n)
#     #X[X==0] <- 0.5
#     X <- exp(X)
#     X <- X / rowSums(X)
#     XX[, , s] <- X
#   }
#
#   sseq <- round(seq(from = 0, to = 1, length.out = ns ), 5)
#   BS <- splines::bs(sseq, df = df_beta, intercept = TRUE)
#   beta <- BS %*% t(beta_C_matrix)
#   D <- plyr::alply(XX, .margins = 1,
#                    function(x, sseq) data.frame(t(rbind(TIME = sseq, x))),
#                    sseq = sseq )
#
#   Z_ITG <- sapply(D, function(x, sseq, beta.basis, insert, method){
#     x[, -1] <- log(x[, -1])
#     return(ITG(X = x, basis = beta.basis, sseq = sseq, insert = insert, method = method)$integral)
#   } ,sseq = sseq, beta.basis = BS, insert = "FALSE", method = "trapezoidal")
#
#   Z_ITG <- t(Z_ITG)
#
#   Z_t.full <- plyr::ldply(D[1:n], data.frame, .id = "Subject_ID")
#
#   Y.tru <- Z_ITG %*% as.vector(t(beta_C_matrix))
#   error <- rnorm(n, 0, 1)
#   sigma_e <- sd(Y.tru) / (sd(error) * SNR)
#   error <- sigma_e*error
#   Y.obs <- Y.tru + error
#
#   Z_t.obs <- lapply(D, function(x, obs_spar) {
#     n <- dim(x)[1]
#     #lambda <- obs_spar * n
#     #n.obs <- rpois(1, lambda)
#     #T.obs <- sample(n, size = n.obs)
#     T.obs <- replicate(n, sample(c(0,1), 1, prob = c(1 - obs_spar, obs_spar)))
#     T.obs <- which(T.obs == 1)
#     #T.obs = sort(sample(seq(n), round(obs_spar*n)))
#     x.obs <- x[T.obs, ]
#     return(x.obs)
#   }, obs_spar = obs_spar)
#
#   Z_t.obs <- plyr::ldply(Z_t.obs, data.frame, .id = "Subject_ID")
#
#
#   data <- list(y = Y.obs, Comp = Z_t.obs
#                #, Zc = Z_control, intercept = intercept
#   )
#   beta <- c(beta_C, 0)#, beta_c, ifelse(intercept, beta0, 0))
#   data.raw <- list(Z_t.full = Z_t.full, Z_ITG = Z_ITG,
#                    Y.tru = Y.tru)
#   basis.info <- cbind(sseq, BS)
#
#   output <- list(data = data, beta = beta, basis.info = basis.info, data.raw = data.raw#,
#                  #parameter = parameter
#   )
#   #### Output ####
#
#   return(output)
#
#
# }









