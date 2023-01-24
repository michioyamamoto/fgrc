##☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★
##     Warm starts for the values of rho1 and rho2
##
##  ファイル名：WarmStart.R
##  ファイル内容：
##  作成者：YAMAMOTO, Michio
##  作成日：2015年03月10日
##  最終更新日：2015年03月11日
##  コメント：L_C, L_D, N.clustはgivenとする
##☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★


##### arguments #####
##*X: N.sub*N.varのデータ行列 (matrix)
##*N.comp1: 成分数L_C
##*N.comp2: 成分数L_D
##*N.clust: クラスター数
##*rho1.vec: rho1の候補 (vector)
##*rho2.vec: rho2の候補 (vector)
##*by.rho1: rho1.vecの要素間の差
##*by.rho2: rho2.vecの要素間の差
##*show.ite: パラメータ番号を出力するかどうか
##*mc.cores: 最初の解の計算に利用するコア数


#### return ####
##*ret.list: 各パラメータの組 (rho1, rho2) におけるGRCの結果（list）
##*A.list: 各パラメータの組 (rho1, rho2) におけるA（list）
##*cluster.mat: 各パラメータの組 (rho1, rho2) におけるcluster（matrix）


##### dependency #####
if(nchar(system.file(package="grc")) == 0){
    library(devtools)
    install_github("michioyamamoto/grc")
}


##### main function #####
## N.comp1 <- 2; N.comp2 <- 0; N.clust <- 3; by.rho1 <- 0.5; by.rho2 <- 0.5; rho1.vec <- seq(-2, 2, by=by.rho1); rho2.vec <- seq(-3, 1, by=by.rho2); show.ite <- TRUE; mc.cores <- 8
warms <- function (X, N.comp1, N.comp2, N.clust, rho1.vec=NULL, rho2.vec=NULL, by.rho1=NULL, by.rho2=NULL, N.random=1, nstart=100, Basis="Bsp", par.Bsp=NULL, lambda=0, maxit=100, eps=1e-05, show.ite=FALSE, mc.cores=1) {
  if(Basis=="Bsp" && is.null(par.Bsp))
    stop("par.Bsp should be specified if 'Bsp' is used as basis functions")
  if(is.null(rho1.vec) || is.null(rho2.vec) || is.null(by.rho1) || is.null(by.rho2))
    stop("rho1.vec, rho2.vec, by.rho1, and by.rho2 should be specified")
  if(nchar(system.file(package="fda")) == 0)
    stop("Package 'fda' is required to implement 'FGRC'.")
##  require(fda, quietly=TRUE)
  if(nchar(system.file(package="grc")) == 0)
    stop("Package 'grc' is required to implement 'warms'.")

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

  ## parameters for warm starts
  X.orig <- X ## save the original dataset
  X <- G.H.P ## expanded dataset of matrix form
  par.set.temp <- expand.grid(N.comp1=N.comp1, N.comp2=N.comp2, N.clust=N.clust,
                              rho1=rho1.vec, rho2=rho2.vec)
  par.set <- subset(par.set.temp, rho1 > rho2)
  rownames(par.set) <- c(1:dim(par.set)[1])

  ## rho2の値ごとにリストに格納する
  N.rho2 <- length(rho2.vec)
  rho.list <- list()
  num.par.by.rho2 <- numeric(N.rho2)
  for (n.rho2 in 1:N.rho2) {
    par.subset <- subset(par.set, rho2 == rho2.vec[n.rho2])
    rho.list <- append(rho.list, list(par.subset))
    num.par.by.rho2[n.rho2] <- dim(par.subset)[1]
  }

  ## matrix for representing parameter patterns
  rho1.vec <- seq(0, 2, by=by.rho1)
  rho2.vec <- seq(-1, 1, by=by.rho2)
  rho.set.temp <- expand.grid(rho1=rho1.vec, rho2=rho2.vec)
  rho.set <- subset(rho.set.temp, rho1 > rho2)


  ## conduct FGRC for each parameter pattern
  A.arr <- array(0, dim=c(dim(par.set)[1], dim(X)[2], (N.comp1 + N.comp2)))
  ret.list <- NULL
  ret <- NULL
  ite.show <- 0

  if (show.ite)
    cat("N.par = ", dim(par.set)[1], "\n")
  for (n.list in 1:N.rho2) {
    for (n.par.by.rho2 in 1:num.par.by.rho2[n.list]) {
      ite.show <- ite.show + 1
      if (show.ite) {
        if (ite.show %% 100 == 0) {
          cat(ite.show, "\n")
        } else {
          cat(".")
        }
      }
      ret.old <- ret
      rho1.eval <- rho.list[[n.list]]$rho1[n.par.by.rho2]
      rho2.eval <- rho.list[[n.list]]$rho2[n.par.by.rho2]

      ## multiple initial values are given and executed only at the first time
      if (is.null(ret.old)) {
        ret <- list()
        class(ret) <- "try-error"
        while (class(ret) == "try-error")
          ret <- try(grc::GRC(X, N.comp1, N.comp2, N.clust, rho1.eval, rho2.eval, N.random, nstart, FALSE, maxit, eps, NULL, mc.cores), silent=FALSE)
        ret.list <- append(ret.list, list(ret))
      } else {
        ## 近い設定値でのsolutionを初期値とする
        if (n.list != 1 && n.par.by.rho2 == 1) {
          abs.diff1 <- abs(par.set$rho1 - rho1.eval) < by.rho1 / 10^7 ## 小数点の計算の誤差を考慮
          abs.diff2 <- abs(par.set$rho2 - (rho2.eval - by.rho2)) < by.rho2 / 10^7 ## 小数点の計算の誤差を考慮
          par.set.temp <- data.frame(par.set, abs.diff1=abs.diff1, abs.diff2=abs.diff2) ## 小数点の計算の誤差を考慮
          num.par.set <- as.numeric(rownames(subset(par.set.temp, abs.diff1==TRUE & abs.diff2==TRUE)))
          ret.old <- ret.list[[num.par.set]]
        } else {
          abs.diff1 <- abs(par.set$rho1 - rho1.eval) < by.rho1 / 10^7 ## 小数点の計算の誤差を考慮
          abs.diff2 <- abs(par.set$rho2 - rho2.eval) < by.rho2 / 10^7 ## 小数点の計算の誤差を考慮
          par.set.temp <- data.frame(par.set, abs.diff1=abs.diff1, abs.diff2=abs.diff2) ## 小数点の計算の誤差を考慮
          num.par.set <- as.numeric(rownames(subset(par.set.temp, abs.diff1==TRUE & abs.diff2==TRUE))) - 1
          ret.old <- ret.list[[num.par.set]]
        }

        ## past solutions are used as initial values
        class(ret) <- "try-error"
        while (class(ret) == "try-error")
          ret <- try(grc::GRC(X, N.comp1, N.comp2, N.clust, rho1.eval, rho2.eval, 100, 100, FALSE, 100, 1e-05, ret.old$A, 1), silent=FALSE)

        ret.list <- append(ret.list, list(ret))
      }
    }
  }

  ## obtain the values of A and cluster
  A.list <- lapply(ret.list, FUN=function(x) x$A)
  cluster.mat <- sapply(ret.list, FUN=function(x) x$cluster)

  if (show.ite)
    cat("\n")

  return(list("ret"=ret.list, "A"=A.list, "cluster"=cluster.mat, "rho.set"=rho.set))
}



