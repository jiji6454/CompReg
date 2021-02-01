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