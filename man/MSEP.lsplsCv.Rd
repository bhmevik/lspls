\encoding{UTF-8}
\name{MSEP.lsplsCv}
\alias{MSEP.lsplsCv}
\alias{RMSEP.lsplsCv}
\alias{R2.lsplsCv}
\title{MSEP, RMSEP and R^2 for LS-PLS}
\description{
  (Root) Mean Squared Error of Prediction ((R)MSEP) and R^2 methods for LS-PLS
  cross-validations (\code{"lsplsCv"} objects).
}
\usage{
\method{MSEP}{lsplsCv}(object, scale = FALSE, \dots)
\method{RMSEP}{lsplsCv}(object, scale = FALSE, \dots)
\method{R2}{lsplsCv}(object, \dots)
}
\arguments{
  \item{object}{an \code{"lsplsCv"} object, typically the output from
    \code{lsplsCv}.}
  \item{scale}{logical.  Whether the responses and predicted values
    should be divided by the standard deviation of the response prior to
    calculating the measure.  This is most useful when comparing several
    responses.  Default is not to scale.  Note that this argument is
    ignored by the \code{R2} method, since \eqn{R^2} is independent of
    scale.}
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
\keyword{multivariate}
\keyword{regression}
