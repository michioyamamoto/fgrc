\name{warms}
\alias{warms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Warm Starts in FGRC.
}
\description{
	This function conducts FGRC with warm staring values.
}
\usage{
warms(X, N.comp1, N.comp2, N.clust, rho1.vec=NULL, rho2.vec=NULL,
by.rho1=NULL, by.rho2=NULL, N.random=1,
nstart=100, Basis="Bsp", par.Bsp=NULL, lambda=0,
maxit=100, eps=1e-05, show.ite=FALSE, mc.cores=1)
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
		The number of components, which correspond to the subspace for a
		cluster structure.
	}
  \item{N.comp2}{
		The number of components, which correspond to the subspace for
		disturbing structure.
	}
  \item{N.clust}{
		The number of clusters.
	}
  \item{rho1.vec}{
		A tuning parameter vector for the within-group variability in a subspace.
	}
  \item{rho2.vec}{
		A tuning parameter vector for the between-group variability in a subspace.
	}
  \item{by.rho1}{
		A tuning parameter vector for the within-group variability in a subspace.
	}
  \item{by.rho2}{
		A tuning parameter vector for the within-group variability in a subspace.
	}
  \item{N.random}{
		The number of random sets of parmeters for initial random starts.
	}
  \item{show.ite}{
		If \code{"TRUE"}, the number of each iteration is shown.
	}
  \item{nstart}{
		The number of random starts for the k-means clustering in each
		optimization of U, which is a cluster assignment parameters.
	}
  \item{Basis}{
		A name of a basis function. Only B-spline basis function, "Bsp", is implemented.
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
		If \code{"parallel"} package has been installed, \code{"warms"}
		adopts a multithread process for multiple initial random starts. If
		\code{"mc.cores"=1}, \code{"parallel"} package is not needed, and a
		single core process is conducted.
	}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
	\item{ret}{
		A list of all results from generalized reduced clustering for
		coefficients for data matrices \code{G_{H}}. The order of the list
		is the number of a pair of rho's.
	}
	\item{A}{
		A list of all estimated coefficient matrices \code{A_{H}}. The order of the list
		is the number of a pair of rho's.
	}
	\item{cluster}{
		A matrix of all estimated cluster assignments (subject*a pair of rho's). The order of the column
		is the number of a pair of rho's.
	}
	\item{rho.set}{
		A matrix of values of rho1 and rho2. The row number corresponds to the
		number of a pair of rho's.
	}
}
\references{
	Yamamoto, M. and Hwang, H. (2017). Dimension-reduced clustering of
	functional data via subspace separation. Journal of Classification,
	34(2): 294--326. doi:10.1007/s00357-017-9232-z.
}
\author{
	Michio Yamamoto\cr
	\email{michio.koko@gmail.com}
}
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
## Random Functional Data (unmeaningful example)
## 100 subjects with 50 sampling points (single-variate functional data)
## Consider three clusters in a two-dimensional subspace
if(0){
X <- array(rnorm(100 * 50), dim=c(100, 50, 1))
by.rho1 <- 0.10
by.rho2 <- 0.10
rho1.vec <- seq(0, 2, by=by.rho1)
rho2.vec <- seq(-1, 1, by=by.rho2)
res <- warms(X, 2, 0, 3, rho1.vec, rho2.vec, by.rho1, by.rho2, 10, 100,
"Bsp", c(8,4), 0, 100, 1e-05, TRUE, 1)
}
}

