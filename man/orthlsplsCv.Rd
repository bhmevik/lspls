\encoding{latin1}
\name{orthlsplsCv}
\alias{orthlsplsCv}
\title{Low Level Cross-Validation Function}
\description{
  Low-level function to perform the cross-validation in \code{lsplsCv}.
}
\usage{
orthlsplsCv(Y, X, Z, ncomp, segments, \dots)
}
\arguments{
  \item{Y}{matrix.  Response matrix.}
  \item{X}{matrix.  The first predictor matrix (typically a design matrix).}
  \item{Z}{list.  List of predictor matrices.}
  \item{ncomp}{list.  The number of components to fit from each matrix.}
  \item{segments}{list.  The segments to use.}
  \item{\dots}{Further arguments.  Currently not used.}
}
\details{
  This function is not meant to be called directly by the user.  It
  performs cross-validation of ortogonalized LS-PLS-models without
  splitting of parallell matrices into common and unique components.
  See the references for details.
}
\value{
  An array of cross-validated predictions.  The first dimension
  corresponds to the observations, the second to the responses, and the
  rest to the number of components of the PLS models.
}
\references{
  Jørgensen, K., Segtnan, V. H., Thyholt, K., Næs, T. (2004)  A
  Comparison of Methods for Analysing Regression Models with Both
  Spectral and Designed Variables.
  \emph{Journal of Chemometrics}, \bold{18}(10), 451--464.

  Jørgensen, K., Mevik, B.-H., Næs, T.  Combining Designed Experiments
  with Several Blocks of Spectroscopic Data.
  (Submitted)

  Mevik, B.-H., Jørgensen, K., Måge, I., Næs, T.  LS-PLS: Combining
  Categorical Design Variables with Blocks of Spectroscopic
  Measurements.
  (Submitted)
}
\author{Bjørn-Helge Mevik}
\seealso{\code{\link{lspls}}, \code{\link{lsplsCv}},
  \code{\link{orthlspls.fit}}}
\keyword{multivariate}
\keyword{regression}
