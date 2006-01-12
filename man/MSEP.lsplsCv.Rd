\encoding{latin1}
\name{MSEP.lsplsCv}
\alias{MSEP.lsplsCv}
\alias{RMSEP.lsplsCv}
\title{MSEP and RMSEP for LS-PLS}
\description{
  (Root) Mean Squared Error of Prediction ((R)MSEP) methods for LS-PLS
  cross-validations (\code{"lsplsCv"} objects).
}
\usage{
\method{MSEP}{lsplsCv}(object, \dots)
\method{RMSEP}{lsplsCv}(object, \dots)
}
\arguments{
  \item{object}{an \code{"lsplsCv"} object, typically the output from
    \code{lsplsCv}.}
  \item{\dots}{Further arguments.  Currently unused.}
}
\value{
  An array.  The first dimension corresponds to the responses (for
  single-response models, the length of this dimension is 1).  The rest
  of the dimensions correspond to the number of components from the
  PLS matrices.
}
\author{Bjørn-Helge Mevik}
\seealso{\code{\link{lsplsCv}}, \code{\link{plot.lsplsCv}}}
\examples{
}
\keyword{multivariate}
\keyword{regression}
