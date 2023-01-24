\name{BasisExpand}
\alias{BasisExpand}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	A function for a basis function expansion
}
\description{
	This function conducts the basis function expansion using the package 'fda'.
}
\usage{
BasisExpand(X, Basis="Bsp", par.Bsp=NULL, lambda=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
		A data array with the dimensionality, subjects*time*variable.
	}
  \item{Basis}{
		A name of a basis function. Only B-spline basis functions "Bsp" can
		be specified.
	}
  \item{par.Bsp}{
		Values of parameters needed to use B-spline basis functions.
	}
  \item{lambda}{
		A  smoothing parameter value.
	}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
	\item{A1}{Estimated loading matrix which corresponds to a subspace for
	a cluster structure.}
}
% \references{
% 	De Soete, G. and Carroll, J.D. (1994). K-means clustering in a
% 	low-dimensional Euclidean space. In Diday, E., Lechevallier, Y.,
% 	Schader, M., Bertrand, P., and Burtschy, B. (Eds.) New Approaches in
% 	Classification and Data Analysis (pp.212-219). Springer, Heidelberg.
% }
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
% \examples{
% }

