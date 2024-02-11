\name{FFKMts}
\alias{FFKMts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Functional Factorial K-means Clustering Analysis with two-step approach
}
\description{
	This function conducts the functional factorial k-means clustering
	analysis with two-step approach for rank-deficient data.
}
\usage{
FFKMts(X, N.comp, N.clust, N.random=1, show.random.ite=FALSE,
nstart=100, Basis="Bsp", par.Bsp=NULL, lambda=0, maxit=100, eps=1e-05, cr.cont=0.80, mc.cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
		An object (matrix, array, or list) which includes multivariate
		functional data. If the data consist of a single variable, the
		object is allowed to be a matrix, where the row and the column
		correspond to modes of subject and time, respectively. When the data
		are multivariate, the object has to be a list or an array. See also
		the examples in \code{\link{FGRC}}.
	}
  \item{N.comp}{
		The number of components.
	}
  \item{N.clust}{
		The number of clusters.
	}
  \item{N.random}{
		The number of random sets of parmeters for initial random starts.
	}
  \item{show.random.ite}{
		If \code{"TRUE"}, the number of each iteration is shown.
	}
  \item{nstart}{
		The number of random starts for the k-means clustering in each
		optimization of U.
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
		A smoothing parameter value.
	}
	\item{maxit}{
		A max number of iterations for ALS algorithm.
	}
	\item{eps}{
		A criterion for the convergence of the loss function.
	}
	\item{cr.cont}{
		The threshold to decide the dimensionality of F.pca, which is a
		reduced component score matrix. A value ranged from 0 to 1 should be
		specified.
	}
  \item{mc.cores}{
		If \code{"parallel"} package has been installed, \code{"FFKMts"}
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
		time*component*variable).}
	\item{cluster}{A vector of cluster assignment parameters.}
	\item{F}{Estimated component score matrix which corresponds to a
		subspace for a cluster structure by two-step approach.}
	\item{A}{Estimated loading array of basis*component*variable by
		two-step approach.}
	\item{F.post}{Estimated component score matrix which corresponds to
		the estimated loading array A.}
	\item{F.clust.mean}{Estimated component score matrix for cluster
		centroids of F.}
	H.half.inv %*% A.H.P[num:(num + N.basis - 1), ]
	\item{X.clust.mean}{Projected cluster centroid functions.}
	\item{lossfunc}{The value of loss function.}
	\item{n.ite}{The number of iteration needed for convergence.}
}
\references{
	Yamamoto, M. and Terada, Y. (2014). Functional factorial k-means
	analysis. Computational Statistics & Data Analysis, 79: 133--148.
}
\author{
	Michio Yamamoto\cr
	\email{yamamoto.michio.hus@osaka-u.ac.jp}
}
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{FGRC}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##Random Functional Data (unmeaningful example)
##100 subjects with 50 sampling point (single-variate functional data)
##Consider three clusters in a two-dimensional subspace
X <- array(rnorm(100 * 50), dim=c(100, 50, 1))
res <- FFKMts(X, 2, 2, 10, FALSE, 100, "Bsp", c(8, 4), 0, 100, 1e-05,
0.80, 1)
}

