##☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★
##    Return results of the Basis Function Expansion
##
##  ファイル名：BasisExpand.R
##  ファイル内容：
##  作成者：YAMAMOTO, Michio
##  作成日：2013年04月06日
##  最終更新日：2023年01月日
##  コメント：  ペナルティも考慮しておくこと！！
##            Wavelet基底では2次元にのみ対応．
##            B-spline, Fourier基底は1次元にのみ対応．
##            新たに作り直した（121030）
##            時系列データのみを対象とした
##            とりあえずBスプライン基底関数のみを利用する
##            Myintegの積分範囲をc(1:N.time)に修正しておいた
##☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★

##### 引数 #####
##*X: 配列（N.sub*N.time*N.var）
##*Basis:利用する基底関数の種類（"Bsp":B-spline基底関数）
##*par.Bsp: B-spline基底関数に関する設定値（=c("節点の数", "基底の次数")）
##*lambda: 平滑化パラメータ（全ての次元，変数に対して一定としている）


#### 戻り値 ####
## list型変数 including
##  *X.smooth: data Xの基底関数展開したもの引数Xと同じサイズの配列．
##  *coef.mat: 基底関数展開した際の回帰係数．配列（N.sub*N.basis*N.var）．
##  *Phi: N.basis*N.varの基底関数行列．Matrix型．
##  *S: Spline Smoothing Operator S^2のSにあたるもの．Matrix型．
##      回帰スプラインの文脈では，Hをhat matrixとするとS = H^{1/2}と考えられる．
##  *B.basis: 基底関数展開における基底関数の数．
##  *K: B-spline基底関数の2階微分からなるGram matrix
##  *H.mat: 基底関数のGram matrix（lambdaも考慮されている）
##  *Hat: ハット行列


##### 依存関数 #####
Myinteg <- function(x, quad.point=NULL){
  ##If quad.point=NULL, integrate x on[0,1]
  n <- length(x)
  if(is.null(quad.point))
    quad.point <- seq(0, 1, len=n)
  h <- (quad.point[n] - quad.point[1]) / (n - 1)
  ret <- t(h * c(1/2, rep(1, len=n-2), 1/2)) %*% x

  return(ret)
}


##### 関数本体 #####
## X <- X; Basis <- "Bsp"; par.Bsp <- c(8, 4); lambda <- 0; N.fbasis <- NULL
BasisExpand <- function(X, Basis="Bsp", par.Bsp=NULL, lambda=0)
{
  if (!(is.array(X) && !is.matrix(X))) {
    stop("Error in BasisExpand: X should be an array form.")
  }

  N.sub <- dim(X)[1]
  N.time <- dim(X)[2]
  N.var <- dim(X)[3]

  X.smooth <- array(0, dim=c(N.sub, N.time, N.var))
  Phi <- NULL
  S <- NULL
  K <- NULL
  H.mat <- NULL ## with penalty term
  H <- NULL ## without penalty term

  if (Basis == "Bsp") {
    if (is.null(par.Bsp[1])) ## if N.knots is NULL
      N.knots <- N.time %/% 2 - 2
    if (is.null(par.Bsp[2])) ## if N.knots is NULL
      N.order <- 4

    N.knots <- par.Bsp[1]
    N.order <- par.Bsp[2]
    N.basis <- N.knots + N.order - 2
    coef.arr <- array(0, dim=c(N.sub, N.basis, N.var))

    ## for of B-spline basis functions
    Phi <- bsplineS(c(1:N.time), seq(1, N.time, len=N.knots), N.order)

    ## calculate S in Smoothing Spline Operator S^2 (= (Hat matrix)^{1/2})
    ## calculate K
    basisobj <- create.bspline.basis(c(1, N.time), N.basis, norder=N.order)
    K <- bsplinepen(basisobj)

    ## calculate H (here, H = (Hat matrix))
    Hat.temp <- Phi %*% solve(t(Phi) %*% Phi + lambda * K)
    Hat <- Hat.temp %*% t(Phi)

    ## calculate coef.arr, S, X.smooth
    for (n.var in 1:N.var) {
      ## coef for linear approximation using B-spline basis func with lambda
        coef.arr[,,n.var] <- X[,,n.var] %*% Hat.temp

        ## functions obtained by linear approximation using B-spline basis func
        X.smooth[,,n.var] <- coef.arr[,,n.var] %*% t(Phi)
    }

    temp.eig <- eigen(Hat, symmetric=TRUE)

    ## 固有値が0の部分が計算機上は負として表現される可能性があるため，固有値の絶対値をとっている
    S <- temp.eig$vectors %*% diag(sqrt(abs(temp.eig$values))) %*% t(temp.eig$vectors)

    ## calculate Gram matrix for basis functions
    H.mat <- array(0, dim=c(N.basis, N.basis, N.time))
    for (i in 1:N.basis) {
      for (j in i:N.basis) {
        H.mat[i,j,] <- H.mat[j,i,] <- Phi[,i] * Phi[,j]
      }
    }

    H.mat <- apply(H.mat, c(1,2), Myinteg, quad.point=c(1:N.time))
    H <- H.mat
    H.mat <- H.mat + lambda * K


  } else {
    stop("MY WARNING: Please check the argment 'Basis' in BasisExpand")
  }

  return(list("X.smooth"=X.smooth, "coef.arr"=coef.arr, "Phi"=Phi, "S"=S, "N.basis"=N.basis, "K"=K, "H.mat"=H.mat, "Hat"=Hat, "H"=H))
}


