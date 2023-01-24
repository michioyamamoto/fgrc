##☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★
##   Functional Generalized Reduced Clustering
##
##  ファイル名：FGRC.R
##  ファイル内容：
##  作成者：YAMAMOTO, Michio
##  作成日：2013年09月20日
##  最終更新日：2014年11月26日
##  コメント：grcパッケージを参考に再作成．
##           現状は経時測定データのみを想定している．
##           OptimGRC_Cへの例外処理を行った（141126）
##☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★

##### arguments #####
##X: 過去のデータにも適用可能なように，リスト，行列，配列すべてに対応可能とする
##   N.sub * N.time * N.varを想定している
##N.comp1: dimensionality of a subspace for cluster structure
##N.comp2: dimensionality of a subspace for irrelevant variables
##rho1, rho2: the values of tuning parameters (rho1 > rho2)
##N.random: the number of initial random starts for the overall ALS algorithm
##show.random.ite: if TRUE, the number of iteration is shown
##nstart: the number of initial random starts for the k-means algorithm for each implementation
##Basis: basis function used for basis function expansions (only "Bsp" is available)
##par.Bsp: the number of knots and the degree of B-spline basis function
##lambda: the value of smoothing parameter
##maxit: the maximum number of iteration of the overall ALS algorithm
##eps: if the difference between the loss functions in iterations becomes less than eps, the algorithm will be finised.
##mc.core: the number of cores which is used for the multicore package.


##### dependency #####
## source("BasisExpand.R") ## comment out when creating the package
## dyn.load("OptimGRC_C.so") ## comment out when creating the package
OptimGRC_C <- function (X, U, A, info) {
  .Call("OptimGRC_C",
        as.double(X),
        as.double(U),
        as.double(A),
        as.double(info), ##rho1, rho2も含む
        pacakge = "fgrc")
}


##### main function #####
## X <- X; N.comp1 <- 2; N.comp2 <- 0; N.clust <- 4; rho1 <- 1; rho2 <- 0; N.random <- 10; show.random.ite <- TRUE; nstart <- 100; Basis <- "Bsp"; par.Bsp <- c(8, 4); lambda <- 0; maxit=100; eps=1e-05; mc.cores=1
FGRC <- function (X, N.comp1, N.comp2, N.clust, rho1=1, rho2=0, N.random=1, show.random.ite=FALSE, nstart=100, Basis="Bsp", par.Bsp=NULL, lambda=0, maxit=100, eps=1e-05, mc.cores=1) {

  if(Basis=="Bsp" && is.null(par.Bsp))
    stop("par.Bsp should be specified if 'Bsp' is used as basis functions")
  if(nchar(system.file(package="fda")) == 0)
    stop("Package 'fda' is required to implement FGRC.")
  ## require(fda, quietly=TRUE)

  if (is.list(X)) {
    N.sub <- dim(X[[1]])[1]
    N.time <- dim(X[[1]])[2]
    N.var <- length(X)
    X.temp <- X
    X <- array(0, dim=c(N.sub, N.time, N.var))
    for (i in 1:N.var)
      X[,,i] <- X.temp[[i]]
  } else if (is.matrix(X)) {
    N.sub <- dim(X)[1]
    N.time <- dim(X)[2]
    N.var <- 1
    X.temp <- X
    X <- array(0, dim=c(N.sub, N.time, 1))
    X[,,1] <- X.temp
  } else if (is.array(X)) {
    N.sub <- dim(X)[1]
    N.time <- dim(X)[2]
    N.var <- dim(X)[3]
  }
  N.comp <- N.comp1 + N.comp2

  ## smoothing based on basis function with lambda
  ## it is assumed that the measurement point is the same for all objects
  temp <- BasisExpand(X, Basis, par.Bsp, lambda)
  X.smooth <- temp$X.smooth
  coef.arr <- temp$coef.arr
  Phi <- temp$Phi
  S <- temp$S
  N.basis <- temp$N.basis
  K <- temp$K
  H.mat <- temp$H.mat
  Hat <- temp$Hat
  temp2 <- eigen(H.mat, symmetric=TRUE)
  H.half <- temp2$vectors %*% sqrt(diag(temp2$values)) %*% t(temp2$vectors)
  H.half.inv <- solve(H.half)


  ## preparation for applying the GRC algorithm
  G.H.P <- NULL
  for (n.var in 1:N.var)
    G.H.P <- cbind(G.H.P, coef.arr[,, n.var] %*% H.half)
  N.var.basis <- N.var * N.basis

  ## applying the GRC algorithm
  if (mc.cores > 1) {
    if(nchar(system.file(package="parallel")) == 0)
      stop("Package 'parallel' is required to calculate using multicores.")
##    require("parallel", quietly=TRUE)

    ##----- multiple random starts -----
    temp.solution <- parallel::mclapply(1:N.random, FUN = function(n.random) {
      if (show.random.ite)
        if (n.random %% (N.random / 10) == 0)
          cat(paste(n.random, "; ", sep=""))

      set.seed(n.random)

      A.H.P <- qr.Q(qr(matrix(rnorm(N.var.basis * N.comp), N.var.basis, N.comp))) ## initial values
      conv <- 0
      U <- matrix(0, N.sub, N.clust)
      cluster <- sample(1:N.clust, N.sub, replace=TRUE)

      ## Optimization based on K-means+EVD or K-means+GP algorithm
      info <- c(N.sub, N.var.basis, N.comp1, N.comp2, N.clust, nstart, maxit, eps, rho1, rho2)
      solution <- try(OptimGRC_C(G.H.P, U, A.H.P, info), silent=FALSE) ##Lapack関数への例外処理

      ans <- list()
      if (class(solution) == "try-error") {
        ans$loss <- Inf
      } else {
        ans$A.H.P <- solution[[1]]
        ans$U <- solution[[2]]
        ans$cluster <- c(ans$U %*% c(1:N.clust))
        ans$n.ite <- solution[[3]][1]
        ans$loss <- solution[[3]][2]
      }

      invisible({rm(list=c("solution")); gc(); gc(); gc()})

      return(ans)
    }, mc.cores=mc.cores)

    ## 並列処理からの値の受け渡し
    best.lossfunc <- sapply(temp.solution,
                            FUN = function(x) x$loss
                            )
    cluster.all <- sapply(temp.solution,
                          FUN = function(x) x$cluster
                          )
    if (all(is.na(best.lossfunc)))
      stop("\nmywarning-->Could not find a feasible starting point...exiting\n",
           call. = FALSE)
    nb <- which(best.lossfunc == min(best.lossfunc, na.rm = TRUE))[1]
    solution <- temp.solution[[nb]]
    A.H.P <- solution$A.H.P
    cluster <- solution$cluster
    conv <- solution$conv
    val.lossfunc <- solution$loss
    n.ite <- solution$n.ite

    invisible({rm(list=c("temp.solution")); gc(); gc()})


  } else {
    ## for the case where the parallel package is not used
    ret.val.lossfunc <- Inf

    ##----- multiple random starts -----
    for (n.random in 1:N.random) {
      if (show.random.ite)
        if (n.random %% (N.random / 10) == 0)
          cat(paste(n.random, "; ", sep=""))
      set.seed(n.random)

      A.H.P.temp <- qr.Q(qr(matrix(rnorm(N.var.basis * N.comp), N.var.basis, N.comp))) ## initial values
      conv <- 0
      U.temp <- matrix(0, N.sub, N.clust)
      cluster.temp <- sample(1:N.clust, N.sub, replace=TRUE)
      U.temp[col(U.temp) == cluster.temp] <- 1

      ## Optimization based on K-means+EVD or K-means+GP algorithm
      info <- c(N.sub, N.var.basis, N.comp1, N.comp2, N.clust, nstart, maxit, eps, rho1, rho2)
      solution <- OptimGRC_C(G.H.P, U.temp, A.H.P.temp, info)

      ## comparison of loss function values
      if (solution[[3]][2] < ret.val.lossfunc) {
        ret.val.lossfunc <- solution[[3]][2]
        A.H.P <- solution[[1]]
        U <- solution[[2]]
        cluster <- c(U %*% c(1:N.clust))
        n.ite <- solution[[3]][1]
        val.lossfunc <- solution[[3]][2]
      }
    }
  } ## End of optimization using GRC algorithm


  ## Estimated weight functions and coefficient matrix
  A <- array(0, dim=c(N.basis, N.comp, N.var))
  V <- array(0, dim=c(N.time, N.comp, N.var))
  if (is.matrix(A.H.P)) { ## 全く収束しなかった場合の例外処理
    for (n.var in 1:N.var) {
      num <- (n.var - 1) * N.basis + 1
      A[,, n.var] <- H.half.inv %*% A.H.P[num: (num + N.basis - 1), ]
      V[,, n.var] <- Phi %*% A[,, n.var]
    }

    if (N.comp == N.comp1) {
      A.H.P1 <- A.H.P
      A.H.P2 <- NULL
    } else {
      A.H.P1 <- A.H.P[, 1:N.comp1]
      A.H.P2 <- A.H.P[, (N.comp1 + 1):N.comp]
    }
    F <- G.H.P %*% A.H.P1
    U <- matrix(0, N.sub, N.clust)
    U[col(U) == cluster] <- 1


    ## projected cluster centroid functions
    if (any(diag(t(U) %*% U) == 0)) {
      F.mean <- NULL
      X.proj.clust.mean <- NULL
    } else {
      F.mean <- solve(t(U) %*% U) %*% t(U) %*% F
      X.proj.clust.mean <- array(0, dim=c(N.time, N.clust, N.var))
      for (n.var in 1:N.var) {
        num <- (n.var - 1) * N.basis + 1
        X.proj.clust.mean[,, n.var] <- Phi %*% H.half.inv %*% A.H.P[num:(num + N.basis - 1),] %*% t(A.H.P[num:(num + N.basis - 1),]) %*% t(G.H.P[,num:(num + N.basis - 1)]) %*% U %*% solve(t(U) %*% U)
      }
    }
  } else {
    V <- NULL
    cluster <- NULL
    F <- NULL
    A <- NULL
    X.proj.clust.mean <- NULL
    val.lossfunc <- NULL
    n.ite <- NULL
  }

  if (show.random.ite)
    cat("\n")

  return(list("V"=V, "cluster"=cluster, "F"=F, "A"=A, "F.clust.mean"=F.mean, "X.clust.mean"=X.proj.clust.mean, "lossfunc"=val.lossfunc, "n.ite"=n.ite))
}
