\name{FPCK}
\alias{FPCK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Functional Principal Component K-means Analysis
}
\description{
	This function conducts the functional principal component k-means
	clustering (FPCK) analysis proposed by Yamamoto (2012).
}
\usage{
FPCK(X, N.comp, N.clust, N.random=1, show.random.ite=FALSE, nstart=100, Basis="Bsp", par.Bsp=NULL, lambda=0,
maxit=100, eps=1e-05, mc.cores=1)
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
		optimization of a cluster assignment parameter U.
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
  \item{mc.cores}{
		If \code{"parallel"} package has been installed, \code{"FFKM"}
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
		subspace for a cluster structure.}
	\item{A}{Estimated loading array of basis*component*variable.}
	\item{X.clust.mean}{Projected cluster centroid functions (time*cluster*variable).}
	\item{lossfunc}{The value of loss function.}
	\item{n.ite}{The number of iteration needed for convergence.}
}
\references{
	Yamamoto, M. (2012). Clustering of functional data in a
	low-dimensional subspace. Advances in Data Analysis and
	Classification, 6, 219-247.
}
\author{
	Michio Yamamoto\cr
	\email{michio.koko@gmail.com}
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
res <- FPCK(X, 2, 3, 10, FALSE, 100, "Bsp", c(8, 4), 0, 100, 1e-05, 1)
}

