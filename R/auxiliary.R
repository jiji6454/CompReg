#### NOT export


# vet <- function(X, p, k) {
#   b <- X[-(1: (p * k))]
#   C <- matrix(X[1 : (p * k)], byrow = TRUE, nrow = p)
#   rownames(C) <- paste0("X", 1:p)
#   return(list(b = b, C = C ))
# }
#
#
# Nzero <- function(X, p, k) {
#   Non.zero <- apply(vet(X, p, k)$C, 1, function(x) {
#     test <- max(abs(x)) == 0 #all.equal(x, rep(0, k))
#     #ifelse(test %in% "TRUE", TRUE, FALSE)
#   })
#   Non.zero <- (1:p)[ !Non.zero]
#   return(Non.zero)
# }


### initial value for lam
### maximun lam is designed
### minmun lam is designed by lambda.factor of maximun lambda
lam.ini <- function(y, Z, ix, iy, pf) {
  lam_max <- 0
  n <- length(y)
  X <- t(Z) %*% y / n
  for(i in 1:length(ix)) {
    a <- sqrt(sum(X[ix[i]:iy[i]]^2))
    if(pf[i] > 0) lam_max <- max(lam_max / pf[i], a)
  }
  return(lam_max)
}


cv.test <- function(outlist, y, X, foldid, lam, trim = 0, keep = FALSE) {
  nlam <- length(lam)
  n <- length(y)
  predmat <- matrix(NA, n, nlam) # prediction matrix
  nfolds <- max(foldid)

  for (i in seq(nfolds)) {
    which <- foldid == i
    #pred <- X_test %*% outlist[[i]]$path
    #pred <- apply(outlist[[i]], 1 ,function(beta,X_test) X_test %*% beta, X_test=X[which, ])
    pred <- predict.linear(outlist[[i]], X[which, , drop = FALSE], s = NULL) #predict(outlist[[i]], X[which, , drop = FALSE]) #
    nlami <- length(outlist[[i]]$lam)
    predmat[which, seq(nlami)] <- pred }

  cvraw <- (y - predmat)^2
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))

  if(trim > 0) {
    cv.trim <- apply(cvraw, 2, function(x) {
      x <- x[!is.na(x)]
      x.boundary <- quantile(x, probs = c(trim / 2, 1 - trim / 2))
      x <- x[x < x.boundary[2]]
      x <- x[x >= x.boundary[1]]
      x.mean <- mean(x)
      x.sd <- sd(x)
      return(c(MEAN = x.mean, SD = x.sd))
    })
  } else {cv.trim <- NULL}

  output <- list()
  output$cvm <- cvm
  output$cvsd <- cvsd
  output$cvmtrim <- cv.trim[1, ]
  output$cvsdtrim <- cv.trim[2, ]
  if(keep) output$fit.preval <- predmat
  return(output)
}



proc.comp <- function(x) {
  # transform compositional data 
  
  if(any(x <= 0)) stop("Data must be positive")
  if( is.null(dim(x)) ) x <- matrix(x, nrow = 1) else x <- as.matrix(x)
  
  if( any(abs(rowSums(x) - 1) > 1e-10) ) {
    message("Data is transformed into compositional data by deviding rowSums")
    x <- x / rowSums(x)
  }
  x <- log(x)
  
  return(x)
}

lamtobeta <- function(lambda, beta, s) {
  lamlist <- point.interp(lambda, s)
  if(length(s) == 1) {
    beta = beta[, lamlist$left, drop=FALSE] * (1 - lamlist$frac) +
      beta[, lamlist$right, drop=FALSE] * lamlist$frac
  } else {
    beta = beta[, lamlist$left, drop=FALSE] %*% diag(1 - lamlist$frac) +
      beta[, lamlist$right, drop=FALSE] %*% diag(lamlist$frac)
  }
  colnames(beta) <- paste(seq(along = s))
  return(beta)
}


predict.linear <- function(object, newX, s = NULL) {
  beta <- object$beta
  #if(is.null(dim(newX))) newX = matrix(newX, nrow = 1)
  if( (nrow(beta) - 1) != ncol(newX) ) {
    stop("numbers of variables in data is not consistant with estimated coefficients")
  }
  if(!is.null(s)) {
    beta <- lamtobeta(lambda = object$lam, beta = beta, s = s)
  }
  fitting <- cbind2(newX, 1) %*% beta
  return(fitting)
}


# cv.test2 <- function(outlist, y, Z, Zc = NULL, foldid, lam, trim = 0) {
#   nlam <- length(lam)
#   n <- length(y)
#   predmat <- matrix(NA, n, nlam)
#   nfolds <- max(foldid)

#   for (i in seq(nfolds)) {
#     which <- foldid == i
#     pred <- predict(outlist[[i]], Znew = Z[which, , drop = FALSE], Zcnew = Zc[which, , drop = FALSE])
#     nlami <- length(outlist[[i]]$lam)
#     predmat[which, seq(nlami)] <- pred }

#   cvraw <- (y - predmat)^2
#   N <- length(y) - apply(is.na(predmat), 2, sum)
#   cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
#   cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))

#   if(trim > 0) {
#     cv.trim <- apply(cvraw, 2, function(x) {
#       x <- x[!is.na(x)]
#       x.boundary <- quantile(x, probs = c(trim / 2, 1 - trim / 2))
#       x <- x[x < x.boundary[2]]
#       x <- x[x >= x.boundary[1]]
#       x.mean <- mean(x)
#       x.sd <- sd(x)
#       return(c(MEAN = x.mean, SD = x.sd))
#     })
#   } else {cv.trim <- NULL}

#   output <- list()
#   output$cvm <- cvm
#   output$cvsd <- cvsd
#   output$cvmtrim <- cv.trim[1, ]
#   output$cvsdtrim <- cv.trim[2, ]
#   return(output)
# }

# cv.test2 <- function(outlist, y, x, foldid, lam, trim = 0) {
#   nlam <- length(lam)
#   n <- length(y)
#   predmat <- matrix(NA, n, nlam)
#   nfolds <- max(foldid)
  
  
#   # pred_fun <- paste("predict", class(outlist[[1]]), sep = "_")
#   # pred_fun <- get(pred_fun)
  
#   for (i in seq(nfolds)) {
#     which <- foldid == i
#     #pred <- do.call(pred_fun, list(outlist[[i]], x[which, , drop = FALSE], s = NULL)) #predict(outlist[[i]], Znew = Z[which, , drop = FALSE], Zcnew = Zc[which, , drop = FALSE])
#     pred <- predict.linear(outlist[[i]], x[which, , drop = FALSE], s = NULL)
#     nlami <- length(outlist[[i]]$lam)
#     predmat[which, seq(nlami)] <- pred }
  
#   cvraw <- (y - predmat)^2
#   N <- length(y) - apply(is.na(predmat), 2, sum)
#   cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
#   cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
  
#   if(trim > 0) {
#     cv.trim <- apply(cvraw, 2, function(x) {
#       x <- x[!is.na(x)]
#       x.boundary <- quantile(x, probs = c(trim / 2, 1 - trim / 2))
#       x <- x[x < x.boundary[2]]
#       x <- x[x >= x.boundary[1]]
#       x.mean <- mean(x)
#       x.sd <- sd(x)
#       return(c(MEAN = x.mean, SD = x.sd))
#     })
#   } else {cv.trim <- NULL}
  
#   output <- list()
#   output$cvm <- cvm
#   output$cvsd <- cvsd
#   output$cvmtrim <- cv.trim[1, ]
#   output$cvsdtrim <- cv.trim[2, ]
#   return(output)
# }


# GIC.test2 <- function(object, y, Z, Zc = NULL, intercept = intercept) {
#   n <- length(y)
#   p <- dim(Z)[2]
#   scaler <- log(log(n)) * log(max(p, n)) / n
#   fix <- as.integer(intercept) + ifelse(is.null(Zc), 0, dim(Zc)[2])
#   predmat <- predict(object, Znew = Z, Zcnew = Zc)
#   cvraw <- (y - predmat)^2
#   MSE <- apply(cvraw, 2, mean) / n
#   S <- apply(abs(object$beta[1:p, ]) > 0, 2, sum)
#   GIC <- log(MSE) + scaler * (ifelse(S>=2, S-1, 0) + fix) ######
#   return(GIC)
# }

# GIC.test2 <- function(object, y, Z, Zc = NULL, intercept = intercept) {
#   n <- length(y)
#   p <- dim(Z)[2]
#   fix <- as.integer(intercept) + ifelse(is.null(Zc), 0, dim(Zc)[2])
#   scaler <- log(log(n)) * log(max(p + fix, n)) / n
#   pred_fun <- paste("predict", class(object), sep = "_")
#   pred_fun <- get(pred_fun)
#   new_com <- cbind(Z, Zc)
#   predmat <- suppressMessages(do.call(pred_fun, list(object, new_com, s = NULL)))
#   cvraw <- (y - predmat)^2
#   #MSE <- apply(cvraw, 2, mean) / n
#   MSE <- colSums(cvraw) / n
#   S <- apply(abs(object$beta[1:p, ]) > 0, 2, sum)
#   GIC <- log(MSE) + scaler * (ifelse(S>=2, S-1, 0) + fix) ######
#   return(GIC)
# }


getmin <- function(lam, cvm, cvsd, digits = 5) {
  idx <- !is.na(cvm)
  cvm1 <- cvm[idx]
  cvsd1 <- cvsd[idx]
  cvm1 <- round(cvm1, digits = digits)
  cvmin <- min(cvm1, na.rm = TRUE)
  idmin <- cvm1 <= cvmin
  lam.min <- max(lam[idmin])
  idmin <- match(lam.min, lam)
  semin <- round((cvm1 + cvsd1), digits = digits)[idmin]
  idmin <- cvm1 <= semin
  lam.1se <- max(lam[idmin])
  output <- list(lam.min = lam.min, lam.1se = lam.1se)
  return(output)
}



ggetmin <- function(lam, cvm, cvsd, digits = 5, k_list) {
  #cvmin <- apply(cvm, 1, min, na.rm = TRUE)
  cvm1 <- round(cvm, digits = digits)
  cvmin <- min(cvm1, na.rm = TRUE)
  lam.min <- apply(cvm1 <= cvmin, 1, function(x, lam)
    ifelse(length(lam[x]) == 0, 0, max(lam[x], na.rm = TRUE)), lam=lam)
# sapply(apply(cvm1 <= cvmin, 1, function(x, lam) lam[x], lam=lam),
#                     function(x) ifelse(length(x) > 0, max(x, na.rm = TRUE), 0))
  lam.min_k <- which.max(lam.min)
  lam.min <- max(lam.min)
  output <- list(lam.min = c("lam" = lam.min, "df" = k_list[lam.min_k]))
  if(!missing(cvsd)) {
    idmin <- match(lam.min, lam)
    semin <- round((cvm + cvsd), digits = digits)[lam.min_k, idmin]
    lam.1se <- apply(cvm1 <= semin, 1, function(x, lam)
      ifelse(length(lam[x]) == 0, 0, max(lam[x], na.rm = TRUE)), lam=lam)
      # sapply(apply(cvm1 <= semin, 1, function(x, lam) lam[x], lam=lam),
      #                 function(x) ifelse(length(x) > 0, max(x, na.rm = TRUE), 0))
    lam.1se_k <- which.max(lam.1se)
    lam.1se <- max(lam.1se)
    output <- c(output, 
                list(lam.1se = c("lam" = lam.1se, "df" = k_list[lam.1se_k])))
  }

  # idmin <- match(lam.min, lam)
  # semin <- round((cvm + cvsd), digits = digits)[lam.min_k, idmin]
  # lam.1se <- apply(cvm1 <= semin, 1, function(x, lam)
  #   ifelse(length(lam[x]) == 0, 0, max(lam[x], na.rm = TRUE)), lam=lam)
  #   # sapply(apply(cvm1 <= semin, 1, function(x, lam) lam[x], lam=lam),
  #   #                 function(x) ifelse(length(x) > 0, max(x, na.rm = TRUE), 0))
  # lam.1se_k <- which.max(lam.1se)
  # lam.1se <- max(lam.1se)

  # output <- list(lam.min = c("lam" = lam.min, "df" = k_list[lam.min_k]),
  #                lam.1se = c("lam" = lam.1se, "df" = k_list[lam.1se_k]))
  return(output)
}


point.interp <- function(sseq_obs, sseq_full) {
  ### sfrac*right+(1-sfrac)*left
  if (length(sseq_obs) == 1) {
    nums <- length(sseq_full)
    left <- rep(1, nums)
    right <- left
    sfrac <- rep(1, nums)
  } else {
    sseq_full[sseq_full > max(sseq_obs)] <- max(sseq_obs)
    sseq_full[sseq_full < min(sseq_obs)] <- min(sseq_obs)
    k <- length(sseq_obs)
    #sfrac <- (sseq_full - sseq_obs[1])/(sseq_obs[k] - sseq_obs[1])
    #sseq_obs <- (sseq_obs - sseq_obs[1])/(sseq_obs[k] - sseq_obs[1])
    #coord <- approx(sseq_obs, seq(sseq_obs), sfrac, method = "linear")$y
    coord <- approx(sseq_obs, seq(sseq_obs), sseq_full, method = "linear")$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sseq_full - sseq_obs[left])/(sseq_obs[right] - sseq_obs[left])
    sfrac[left == right] <- 1
  }
  list(left = left, right = right, frac = sfrac)
}

#'
#' check KKT condition for solution for "FuncompCGL".
#' check KKT condition for solution for "FuncompCGL".
#' @param y response
#' @param Z column combination of compositional variables integration and control variables
#' @param p numbers of compositional variables
#' @param k degree of freedom for basis
#' @param A,b linear constraints
#' @param lam \code{lambda} used
#' @param beta coefficients vector
#' @param tol tolerance de
#'
#' @return
#' \item{summary}{logicl, KKT condition satisfied or not}
#' \item{feasible}{feasible condition}
#' \item{primary}{primary condition}
#' \item{KKT_mutiplier}{KKT_mutiplier}
#' @export

KKT_glasso <- function(y, Z, p, k,
                       A, b,
                       lam, beta,
                       #path,
                       tol = 1e-10){

  #if(missing(lam)) lam <- drop(path$lam)
  #if(missing(beta)) beta <- path$beta
  if(missing(A)) A <- kronecker(matrix(1, ncol = p), diag(k))
  if(missing(b)) b <- rep(0, times = k)
  n <- length(drop(y))
  group <- matrix(1:(p*k), nrow = k)
  #vet(beta, p = p, k = k)$C
  NZ <- Nzero(beta, p = p, k = k)
  residual <- (y - cbind2(Z, 1 )%*% beta) / n
  p1 <- p*k
  m <- dim(Z)[2] - p1
  tao_matrix <- matrix(NA, nrow = k, ncol = length(NZ))
  for(i in 1:length(NZ)) {
    tao_matrix[, i] <- -t(Z[, group[, NZ[i]]]) %*% residual + lam[NZ[i]] * beta[group[, NZ[i]]] / sqrt(sum(beta[group[, NZ[i]]]^2))
    tao_matrix[, i] <- solve(t(A[, group[, NZ[i]]])) %*% tao_matrix[, i]
  }
  tao_matrix <- -tao_matrix
  tao <- rowMeans(tao_matrix)
  tao_matrix_error <- matrix(apply(tao_matrix, 1, function(x) max(x) - min(x)), ncol = 1)
  KKT_multiplier <- all(apply(tao_matrix[,2:dim(tao_matrix)[2], drop = FALSE],
                              2, function(x, y) sum((x-y)^2), y = tao_matrix[, 1, drop = FALSE]) < tol)
  #sum((tao_matrix_error)^2) < tol #all(tao_matrix_error < 1e-5)
  KKT_matrix2 <- matrix(NA, ncol = 3, nrow = length(NZ))
  KKT_matrix2[, 1] <- NZ
  colnames(KKT_matrix2) <- c("group", "value", "condition")
  KKT_matrix2 <- as.data.frame(KKT_matrix2)

  for(i in 1:length(NZ)) {
    KKT_matrix2[i, 2] <- sum( (t(A[, group[, NZ[i]]]) %*% (-tao_matrix[, i] + tao))^2 ) #sqrt(sum( (-tao_matrix[, i] + tao)^2 ))  #tao_matrix[, 1]

  }
  KKT_matrix2[, 3] <- KKT_matrix2[, 2] < tol
  if(m > 0) {
    KKT_control <- -t(Z[, p1+(1:m)]) %*% residual
    KKT_control <- sum(KKT_control^2)
    KKT_matrix2 <- rbind(KKT_matrix2, c(p1+1, KKT_control,KKT_control < tol))
  }

  if((p - length(NZ)) > 0) {
    KKT_matrix <- matrix(NA, ncol = 3, nrow = p - length(NZ))
    KKT_matrix[, 1] <- (1:p)[-NZ]
    colnames(KKT_matrix) <- c("group", "value", "condition")
    KKT_matrix <- as.data.frame(KKT_matrix)
    for(i in 1:(p - length(NZ))){
      KKT_matrix[i, 2] <-  sqrt(sum((-t(Z[, group[, KKT_matrix[i, 1]]]) %*% residual + t(A[, group[, KKT_matrix[i, 1]]]) %*% tao )^2))
    }

    KKT_matrix[, 3] <- KKT_matrix[, 2] <= lam[(1:p)[-NZ]] #+ tol
    KKT_primary <- rbind(KKT_matrix2, KKT_matrix)
  } else {
    KKT_primary <- KKT_matrix2
  }



  KKT_primary$condition <- as.logical(KKT_primary$condition)
  #KKT_primary <- KKT_primary[order(KKT_primary$group), ]
  rownames(KKT_primary) <- KKT_primary[, 1]
  KKT_feasible <- matrix(NA, nrow = k, ncol = 3)
  KKT_feasible[, 1] <- 1:k
  colnames(KKT_feasible) <- c("Constraint", "value", "Condition")
  KKT_feasible <- as.data.frame(KKT_feasible)
  KKT_feasible[, 2] <- A %*% beta[1:p1]
  KKT_feasible[, 3] <- abs(KKT_feasible[, 2]) < tol

  KKT_check <- list( summary = all(c(KKT_primary$condition, KKT_feasible$Condition,  KKT_multiplier) ),
                     feasbile = KKT_feasible,
                     primary = KKT_primary,
                     KKT_mutiplier = tao_matrix)

  return(KKT_check)
}



KKT_lasso <- function(y, Z, p, #path,
                      lam, beta,
                      A, b, tol = 1e-6){

  #if(missing(lam)) lam <- drop(path$lam)
  #if(missing(beta)) beta <- path$beta
  if(missing(A)) A <- rep(1, times = p)
  if(missing(b)) b <- 0
  n <- length(drop(y))
  NZ <- which(abs(beta[1:p]) > 0)
  residual <- (y - cbind2(Z, 1 )%*% beta) / n
  m <- dim(Z)[2] - p

  if(length(NZ) > 0) {
    tao_matrix <- matrix(NA, nrow = 1, ncol = length(NZ))
    for(i in 1:length(NZ)) {
      tao_matrix[, i] <- -t(Z[, NZ[i]]) %*% residual + lam[NZ[i]] * ifelse(beta[NZ[i]] > 0, 1, -1)
      tao_matrix[, i] <- 1 / A[NZ[i]] * tao_matrix[, i]
    }
    tao_matrix <- -tao_matrix
    tao_matrix_error <- max(tao_matrix) - min(tao_matrix)

    tao <- rowMeans(tao_matrix)
    KKT_multiplier <- tao_matrix_error < 1e-5

    KKT_matrix2 <- matrix(NA, ncol = 3, nrow = length(NZ))
    KKT_matrix2[, 1] <- NZ
    colnames(KKT_matrix2) <- c("X", "value", "condition")
    KKT_matrix2 <- as.data.frame(KKT_matrix2)
    for(i in 1:length(NZ)) {
      KKT_matrix2[i, 2] <- abs(A[NZ[i]] * (-tao_matrix[, i] + tao))
    }
    KKT_matrix2[, 3] <- KKT_matrix2[, 2] < tol
  } else {
    KKT_matrix2 <- NULL
  }


  if(m > 0) {
    KKT_control <- -t(Z[, p+(1:m)]) %*% residual
    KKT_control <- sqrt(sum(KKT_control^2))
    KKT_matrix2 <- rbind(KKT_matrix2, c(p+1, KKT_control,KKT_control < tol))
  }

  if((p - length(NZ)) > 0) {
    KKT_matrix <- matrix(NA, ncol = 3, nrow = p - length(NZ))
    KKT_matrix[, 1] <- (1:p)[-NZ]
    colnames(KKT_matrix) <- c("X", "value", "condition")
    KKT_matrix <- as.data.frame(KKT_matrix)
    for(i in 1:(p - length(NZ))){
      KKT_matrix[i, 2] <-  abs(-t(Z[,  KKT_matrix[i, 1]]) %*% residual + A[KKT_matrix[i, 1]] * tao)
    }

    KKT_matrix[, 3] <- KKT_matrix[, 2] <= lam[(1:p)[-NZ]]
    KKT_primary <- rbind(KKT_matrix2, KKT_matrix)
  } else {
    KKT_primary <- KKT_matrix2
  }
  #KKT_primary$condition <- as.logical(KKT_primary$condition)
  #rownames(KKT_primary) <- KKT_primary[, 1]


  KKT_feasible <- matrix(NA, nrow = 1, ncol = 3)
  KKT_feasible[, 1] <- 1
  colnames(KKT_feasible) <- c("Constraint", "value", "Condition")
  KKT_feasible <- as.data.frame(KKT_feasible)
  KKT_feasible[, 2] <- A %*% beta[1:p]
  KKT_feasible[, 3] <- abs(KKT_feasible[, 2]) < tol

  KKT_check <- list( summary = all(c(KKT_primary$condition, KKT_feasible$Condition,  KKT_multiplier) ),
                     feasbile = KKT_feasible,
                     primary = KKT_primary,
                     KKT_mutiplier = tao_matrix)

  return(KKT_check)
}


error.bars <- function(x, upper, lower, width = 0.02, ...)
{
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}
