\name{FGRC}
\alias{FGRC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Functional Generalized Reduced Clustering Analysis
}
\description{
	This function conducts the functional generalized reduced clustering analysis proposed
	by Yamamoto and Hwang (2017).
}
\usage{
FGRC(X, N.comp1, N.comp2, N.clust, rho1=1, rho2=0, N.random=1,
show.random.ite=FALSE, nstart=100, Basis="Bsp", par.Bsp=NULL, lambda=0,
maxit=100, eps=1e-05, mc.cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
		An object (matrix, array, or list) which includes multivariate
		functional data. If the data consist of a single variable, the
		object is allowed to be a matrix, where the row and the column
		correspond to modes of subject and time, respectively. When the data
		are multivariate, the object has to be a list or an array. See
		examples below.
	}
  \item{N.comp1}{
		The number of components that correspond to the subspace for a
		cluster structure.
	}
  \item{N.comp2}{
		The number of components taht correspond to the subspace for
		disturbing structure.
	}
  \item{N.clust}{
		The number of clusters.
	}
  \item{rho1}{
		A tuning parameter for the within-group variability in a subspace.
	}
  \item{rho2}{
		A tuning parameter for the between-group variability in a subspace.
	}
  \item{N.random}{
		The number of random sets of parmeters for initial random starts.
	}
  \item{show.random.ite}{
		If \code{"TRUE"}, the number of each iteration is shown.
	}
  \item{nstart}{
		The number of random starts for the k-means clustering in each
		optimization of U, which is a cluster assignment parameters.
	}
  \item{Basis}{
		A name of a basis function used. Only B-spline basis function, "Bsp", can
		be specified.
	}
  \item{par.Bsp}{
		Parameters values needed to use B-spline basis functions. Those can
		be specified as a vector of length 2, in which the first element
		represents the number of nodes and the second element represents the
		order of the B-spline basis function.
	}
  \item{lambda}{
		A value of a smoothing parameter.
	}
	\item{maxit}{
		A max number of iterations for ALS algorithm.
	}
	\item{eps}{
		A criterion for the convergence of the loss function.
	}
  \item{mc.cores}{
		If \code{"parallel"} package has been installed, \code{"FGRC"}
		adopts a multithread process for multiple initial random starts. If
		\code{"mc.cores"=1}, \code{"parallel"} package is not needed, and a
		single core process is conducted.
	}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
	\item{V}{Estimated weight functions (array of
		time*component*variable). For the i-th variable, \code{V[,,i]}
		consists of two components, \code{V[,,i]=(V_{C}|V_{D})}, that
		correspond to weight functions for cluster and disturbing
		structures, respectively.}
	\item{cluster}{A vector of cluster assignment parameters.}
	\item{F}{Estimated component score matrix which corresponds to a
		subspace for a cluster structure.}
	\item{A}{Estimated loading array of basis*component*variable. For the
		i-th variable, \code{A[,,i]} consists of two components,
		\code{A[,,i]=(A^{C}|A^{D})}, that correspond to coefficient matrices
		of basis function expansions of \code{V_{C}} and \code{V_{D}}, respectively.}
	\item{F.clust.mean}{Estimated component score matrix for cluster
		centroids.}
	\item{X.clust.mean}{Projected cluster centroid functions (time*cluster*variable).}
	\item{lossfunc}{The value of loss function.}
	\item{n.ite}{The number of iteration needed for convergence.}
}
\references{
	Yamamoto, M. and Hwang, H. (2017). Dimension-reduced clustering of
	functional data via subspace separation. Journal of Classification,
	34(2): 294--326. doi:10.1007/s00357-017-9232-z.
}
\author{
	Michio Yamamoto\cr
	\email{yamamoto.michio.hus@osaka-u.ac.jp}
}
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
##Random Functional Data (unmeaningful example)
##100 subjects with 50 sampling points (single-variate functional data)
##Consider three clusters in a two-dimensional subspace
X <- array(rnorm(100 * 50), dim=c(100, 50, 1))
res <- FGRC(X, 2, 0, 3, 1, 0, 10, FALSE, 100, "Bsp", c(8, 4), 0, 100,
1e-05, 1)

##bivariate functional data with 100 subjects and 50 sampling points
X1 <- matrix(rnorm(100 * 50), 100, 50)
X2 <- matrix(rnorm(100 * 50), 100, 50)
X <- list(X1, X2)
res <- FGRC(X, 2, 0, 3, 1, 0, 10, FALSE, 100, "Bsp", c(8, 4), 0, 100,
1e-05, 1)

##Another way for the use of the function (use an array form)
X <- array(rnorm(100 * 50), dim=c(100, 50, 2))
res <- FGRC(X, 2, 0, 3, 1, 0, 10, FALSE, 100, "Bsp", c(8, 4), 0, 100,
1e-05, 1)
}

