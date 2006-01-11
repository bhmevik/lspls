\encoding{latin1}
\name{orthlspls.fit}
\alias{orthlspls.fit}
\title{Underlying LS-PLS Fit Function}
\description{
  Fits orthogonalized LS-PLS models.
}
\usage{
orthlspls.fit(Y, X, Z, ncomp)
}
\arguments{
  \item{Y}{matrix.  Response matrix.}
  \item{X}{matrix.  The first predictor matrix (typically a design matrix).}
  \item{Z}{list.  List of predictor matrices.}
  \item{ncomp}{list.  The number of components to fit from each matrix.}
}
\details{
  \code{orthlspls.fit} is not meant to be called by the user.  It is
  called by \code{\link{lspls}} to do the actual fitting.  See
  \code{\link{lspls}} for details about LS-PLS and \code{ncomp}.  Each
  element of the list \code{Z} should either be a matrix or a list of
  matrices.
}
\value{
  A list with components
  \item{coefficients}{matrix with the final prediction coefficients}
  \item{predictors}{matrix with variables and scores used in the final
    regression}
  \item{orthCoefs}{list of coefficient generating matrices, to be used
    when predicting new predictors.}
  \item{models}{list of fitted PLS models for the matrices}
  \item{ncomp}{list with the number of components used}
  \item{scores}{list of score matrices}
  \item{loadings}{list of loading matrices}
  \item{residuals}{matrix with fit residuals, one coloumn per response}
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
\note{
  The interface (arguments and return values) is likely to change in
  a future version.
}
\seealso{\code{\link{lspls}}}
\keyword{multivariate}
\keyword{regression}
