##¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú
##          Two-step approach for FFKM
##
##  ¥Õ¥¡¥¤¥ëÌ¾¡§FFKMts.R
##  ¥Õ¥¡¥¤¥ëÆâÍÆ¡§
##  ºîÀ®¼Ô¡§YAMAMOTO, Michio
##  ºîÀ®Æü¡§2013Ç¯09·î20Æü
##  ºÇ½ª¹¹¿·Æü¡§2013Ç¯09·î20Æü
##  ¥³¥á¥ó¥È¡§
##¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú¡ù¡ú

##### arguments #####
##X: ²áµî¤Î¥Ç¡¼¥¿¤Ë¤âÅ¬ÍÑ²ÄÇ½¤Ê¤è¤¦¤Ë¡¤¥ê¥¹¥È¡¤¹ÔÎó¡¤ÇÛÎó¤¹¤Ù¤Æ¤ËÂÐ±þ²ÄÇ½¤È¤¹¤ë
##   N.sub * N.time * N.var¤òÁÛÄê¤·¤Æ¤¤¤ë
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
        as.double(info), ##rho1, rho2¤â´Þ¤à
        pacakge = "fgrc")
}


##### main function #####
## X <- X; N.comp <- 2; N.clust <- 4; N.random <- 10; show.random.ite <- TRUE; nstart <- 100; Basis <- "Bsp"; par.Bsp <- c(8, 4); lambda <- 0; maxit <- 100; eps <- 1e-05; cr.cont <- 0.80; mc.cores <- 1
FFKMts <- function (X, N.comp, N.clust, N.random=1, show.random.ite=FALSE, nstart=100, Basis="Bsp", par.Bsp=NULL, lambda=0, maxit=100, eps=1e-05, cr.cont=0.80, mc.cores=1) {
  N.comp1 <- N.comp
  N.comp2 <- 0
  rho1 <- 2
  rho2 <- 1

  if(Basis=="Bsp" && is.null(par.Bsp))
    stop("par.Bsp should be specified if 'Bsp' is used as basis functions")
  if(nchar(system.file(package="fda")) == 0)
    stop("Package 'fda' is required to implement FGRC.")
##  require(fda, quietly=TRUE)

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
    X[,, 1] <- X.temp
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

  ## avoiding rank deficient using the ordinary PCA
  ret.eigen <- eigen(t(G.H.P) %*% G.H.P)
  SumCont <- function (eigen.val) {
    N.comp <- length(eigen.val)
    ret <- numeric(N.comp)
    for (n.comp in 1:N.comp) {
      ret[n.comp] <- sum(eigen.val[1:n.comp])
    }
    sum.cont <- ret / sum(eigen.val)
    return(sum.cont)
  }
  sum.cont <- SumCont(ret.eigen$values)
  N.reduced <- which(sum.cont > cr.cont)[1]
  if (N.reduced < N.comp)
    N.reduced <- N.comp
  ret.svd <- svd(G.H.P)
  F.pca <- ret.svd$u[,1:N.reduced] %*% diag(ret.svd$d[1:N.reduced])

  ## ¤â¤·¥é¥ó¥¯½èÍý¸å¤ÎÀ®Ê¬¿ô¤¬µá¤á¤¿¤¤À®Ê¬¿ô°Ê²¼¤Î¾ì¹ç¡¤N.reduced <- N.comp¤È¤·¤Æ¡¤F.pca¤ËÂÐ¤¹¤ëÃ±¤Ê¤ëk-means¤È¤Ê¤ë

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

      A.H.P <- qr.Q(qr(matrix(rnorm(N.reduced * N.comp), N.reduced, N.comp))) ## initial values
      val.lossfunc <- Inf
      conv <- 0
      U <- matrix(0, N.sub, N.clust)

      ## optimization based on K-means+EVD or K-means+GP algorithm
      info <- c(N.sub, N.reduced, N.comp1, N.comp2, N.clust, nstart, maxit, eps, rho1, rho2)
      solution <- OptimGRC_C(F.pca, U, A.H.P, info)

      ans <- list()
      ans$A.H.P <- solution[[1]]
      ans$U <- solution[[2]]
      ans$cluster <- c(ans$U %*% c(1:N.clust))
      ans$n.ite <- solution[[3]][1]
      ans$loss <- solution[[3]][2]

      return(ans)
    }, mc.cores=mc.cores)

    ## ÊÂÎó½èÍý¤«¤é¤ÎÃÍ¤Î¼õ¤±ÅÏ¤·
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
      set.seed(n.random)

      A.H.P.temp <- qr.Q(qr(matrix(rnorm(N.reduced * N.comp), N.reduced, N.comp))) ## initial values
      val.lossfunc <- Inf
      conv <- 0
      U.temp <- matrix(0, N.sub, N.clust)
      cluster.temp <- sample(1:N.clust, N.sub, replace=TRUE)
      U.temp[col(U.temp) == cluster.temp] <- 1

      ## optimization based on K-means+EVD or K-means+GP algorithm
      info <- c(N.sub, N.reduced, N.comp1, N.comp2, N.clust, nstart, maxit, eps, rho1, rho2)
      solution <- OptimGRC_C(G.H.P, U.temp, A.H.P.temp, info)

      ## comparison of loss functio values
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


  ## estimated weight functions
  B.H.P <- ret.svd$v[, 1:N.reduced]
  Z <- F.pca %*% A.H.P
  Y <- F.pca %*% t(B.H.P)
  ret.svd2 <- svd(t(Y) %*% Z)
  A.H.P.post <- ret.svd2$u %*% t(ret.svd2$v)
  A.post <- array(0, dim=c(N.basis, N.comp, N.var))
  V <- array(0, dim=c(N.time, N.comp, N.var))
  for (n.var in 1:N.var) {
    num <- (n.var - 1) * N.basis + 1
    A.post[,, n.var] <- H.half.inv %*% A.H.P.post[num: (num + N.basis - 1), ]
    V[,, n.var] <- Phi %*% A.post[,, n.var]
  }

  if (N.comp == N.comp1) {
    A.H.P1 <- A.H.P
    A.H.P2 <- NULL
    A.H.P1.post <- A.H.P.post
    A.H.P2.post <- NULL
  } else {
    A.H.P1 <- A.H.P[, 1:N.comp1]
    A.H.P2 <- A.H.P[, (N.comp1 + 1):N.comp]
    A.H.P1.post <- A.H.P.post[, 1:N.comp1]
    A.H.P2.post <- A.H.P.post[, (N.comp1 + 1):N.comp]
  }
  F <- F.pca %*% A.H.P1
  F.post <- G.H.P %*% A.H.P1.post
  U <- matrix(0, N.sub, N.clust)
  U[col(U) == cluster] <- 1


  ## projected cluster centroid functions
  F.mean <- solve(t(U) %*% U) %*% t(U) %*% F
  X.proj.clust.mean <- array(0, dim=c(N.time, N.clust, N.var))
  for (n.var in 1:N.var) {
    num <- (n.var - 1) * N.basis + 1
    X.proj.clust.mean[,, n.var] <- Phi %*% H.half.inv %*% A.H.P.post[num:(num + N.basis - 1),] %*% t(A.H.P.post[num:(num + N.basis - 1),]) %*% t(G.H.P[,num:(num + N.basis - 1)]) %*% U %*% solve(t(U) %*% U)
  }

  if (show.random.ite)
    cat("\n")

  return(list("V"=V, "cluster"=cluster, "F"=F, "A"=A.post, "F.post"=F.post, "F.clust.mean"=F.mean, "X.clust.mean"=X.proj.clust.mean, "lossfunc"=val.lossfunc, "n.ite"=n.ite))
}
