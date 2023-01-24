FPCK <- function (X, N.comp, N.clust, N.random=1, show.random.ite=FALSE, nstart=100, Basis="Bsp", par.Bsp=NULL, lambda=0, maxit=100, eps=1e-05, mc.cores=1)
{
  FGRC(X, N.comp, 0, N.clust, 1, 0, N.random, show.random.ite, nstart, Basis, par.Bsp, lambda, maxit, eps, mc.cores)
}
