##$B!y!z!y!z!y!z!y!z!y!z!y!z!y!z!y!z!y!z!y!z!y!z!y!z!y!z!y!z(B
##    Return results of the Basis Function Expansion
##
##  $B%U%!%$%kL>!'(BBasisExpand.R
##  $B%U%!%$%kFbMF!'(B
##  $B:n@.<T!'(BYAMAMOTO, Michio
##  $B:n@.F|!'(B2013$BG/(B04$B7n(B06$BF|(B
##  $B:G=*99?7F|!'(B2023$BG/(B01$B7nF|(B
##  $B%3%a%s%H!'(B  $B%Z%J%k%F%#$b9MN8$7$F$*$/$3$H!*!*(B
##            Wavelet$B4pDl$G$O(B2$B<!85$K$N$_BP1~!%(B
##            B-spline, Fourier$B4pDl$O(B1$B<!85$K$N$_BP1~!%(B
##            $B?7$?$K:n$jD>$7$?!J(B121030$B!K(B
##            $B;~7ONs%G!<%?$N$_$rBP>]$H$7$?(B
##            $B$H$j$"$($:(BB$B%9%W%i%$%s4pDl4X?t$N$_$rMxMQ$9$k(B
##            Myinteg$B$N@QJ,HO0O$r(Bc(1:N.time)$B$K=$@5$7$F$*$$$?(B
##$B!y!z!y!z!y!z!y!z!y!z!y!z!y!z!y!z!y!z!y!z!y!z!y!z!y!z!y!z(B

##### $B0z?t(B #####
##*X: $BG[Ns!J(BN.sub*N.time*N.var$B!K(B
##*Basis:$BMxMQ$9$k4pDl4X?t$N<oN`!J(B"Bsp":B-spline$B4pDl4X?t!K(B
##*par.Bsp: B-spline$B4pDl4X?t$K4X$9$k@_DjCM!J(B=c("$B@aE@$N?t(B", "$B4pDl$N<!?t(B")$B!K(B
##*lambda: $BJ?3j2=%Q%i%a!<%?!JA4$F$N<!85!$JQ?t$KBP$7$F0lDj$H$7$F$$$k!K(B


#### $BLa$jCM(B ####
## list$B7?JQ?t(B including
##  *X.smooth: data X$B$N4pDl4X?tE83+$7$?$b$N0z?t(BX$B$HF1$8%5%$%:$NG[Ns!%(B
##  *coef.mat: $B4pDl4X?tE83+$7$?:]$N2s5"78?t!%G[Ns!J(BN.sub*N.basis*N.var$B!K!%(B
##  *Phi: N.basis*N.var$B$N4pDl4X?t9TNs!%(BMatrix$B7?!%(B
##  *S: Spline Smoothing Operator S^2$B$N(BS$B$K$"$?$k$b$N!%(BMatrix$B7?!%(B
##      $B2s5"%9%W%i%$%s$NJ8L.$G$O!$(BH$B$r(Bhat matrix$B$H$9$k$H(BS = H^{1/2}$B$H9M$($i$l$k!%(B
##  *B.basis: $B4pDl4X?tE83+$K$*$1$k4pDl4X?t$N?t!%(B
##  *K: B-spline$B4pDl4X?t$N(B2$B3,HyJ,$+$i$J$k(BGram matrix
##  *H.mat: $B4pDl4X?t$N(BGram matrix$B!J(Blambda$B$b9MN8$5$l$F$$$k!K(B
##  *Hat: $B%O%C%H9TNs(B


##### $B0MB84X?t(B #####
Myinteg <- function(x, quad.point=NULL){
  ##If quad.point=NULL, integrate x on[0,1]
  n <- length(x)
  if(is.null(quad.point))
    quad.point <- seq(0, 1, len=n)
  h <- (quad.point[n] - quad.point[1]) / (n - 1)
  ret <- t(h * c(1/2, rep(1, len=n-2), 1/2)) %*% x

  return(ret)
}


##### $B4X?tK\BN(B #####
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

    ## $B8GM-CM$,(B0$B$NItJ,$,7W;;5!>e$OIi$H$7$FI=8=$5$l$k2DG=@-$,$"$k$?$a!$8GM-CM$N@dBPCM$r$H$C$F$$$k(B
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


