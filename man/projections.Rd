\encoding{latin1}
\name{projections}
\alias{project}
\alias{orth}
\alias{Corth}
\title{Projection and Orthogonalisation}
\description{
  Functions to project one matrix onto another, or to ortghogonalise it
  against the other.
}
\usage{
project(M, N)
orth(M, N)
Corth(M, N)
}
\arguments{
  \item{M}{matrix to be projected or orthogonalised}
  \item{N}{matrix to be projected onto or orthogonalised against}
}
\details{
  \code{project(M, N)} calculates the projection of \code{M} onto \code{N},
  i.e., \eqn{N (N^t N)^{-1} N^t M}{N (N' N)^(-1) N' M}.

  \code{orth(M, N)} orthogonalises \code{M} with respect to \code{N},
  i.e., it calculates the projection of \code{M} onto the orthogonal
  space of \code{N}: \eqn{M - N (N^t N)^{-1} N^t M}{M - N (N' N)^(-1) N' M}.

  \code{Corth(M, N)} calculates the coefficient matrix needed to
  orthogonalise future matrices, i.e.,
  \eqn{(N^t N)^{-1} N^t M}{(N' N)^(-1) N' M}.  Future
  matrices \code{m} and \code{n} can be orthogonalised with
  \code{m - n \%*\% Corth(M, N)}.
}
\value{
  A matrix.
}
\author{Bjørn-Helge Mevik}
\note{
  The functions need to be opitmised, both for speed and numerical
  accurracy.
}
\seealso{\code{\link{lspls}}, \code{\link{lsplsCv}},
  \code{\link{predict.lspls}}}
\examples{
##FIXME
}
\keyword{algebra}
