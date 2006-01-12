\encoding{latin1}
\name{lsplsCv}
\alias{lsplsCv}
\title{Cross-Validate LS-PLS Models}
\description{
  Calculate cross-validated predictions for LS-PLS models.
}
\usage{
lsplsCv(formula, ncomp, data, subset, na.action, segments = 10,
        segment.type = c("random", "consecutive", "interleaved"),
        length.seg, model = TRUE, \dots)
}
\arguments{
  \item{formula}{model formula.  See Details.}
  \item{ncomp}{list or vector of positive integers, giving the
    number of components to use for each PLS matrix.  See Details.}
  \item{data}{an optional data frame with the data to fit the model from.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when
    the data contain missing values.}
  \item{segments}{the number of segments to use, or a list with segments (see
          Details).}
  \item{segment.type}{the type of segments to use.  Ignored if
    \code{segments} is a list.}
  \item{length.seg}{Positive integer.  The length of the segments to
    use.  If specified, it overrides \code{segments} unless
    \code{segments} is a list.}
  \item{model}{logical.  If \code{TRUE}, the model frame is returned.}
  \item{\dots}{additional arguments, passed to the underlying
    cross-validation function.}
}
\details{
  The function performs a cross-validation, using the model and segments
  specified in the call.  It returns an object of class
  \code{"lsplsCv"}, which has a plot method (see
  \code{\link{plot.lsplsCv}}).  See \link{lspls-package} for typical
  usage and more about LS-PLS models.

  See \code{\link{lspls}} for details about specifying the model
  with \code{formula} and \code{ncomp}.  Note that \code{lsplsCv}
  cross-validates models with from 0 components to the numbers of
  components specified with \code{ncomp}.

  If \code{segments} is a list, the arguments \code{segment.type} and
  \code{length.seg} are ignored.  The elements of the list should be
  integer vectors specifying the indices of the segments.  See
  \code{\link[pls]{cvsegments}} for details.

  Otherwise, segments of type \code{segment.type} are generated.  How
  many segments to generate is selected by specifying the number of
  segments in \code{segments}, or giving the segment length in
  \code{length.seg}.  If both are specified, \code{segments} is
  ignored.
}
\value{
  An object of class \code{"lsplsCv"}, with components
  \item{pred}{the cross-validated predictions.  An array with one
    dimension for the observations, one for the responses, and one for
    each of the PLS matrices.}
  \item{segments}{the list of segments used in the cross-validation.}
  \item{na.action}{if observations with missing values were removed,
    \code{na.action} contains a vector with their indices.}
  \item{ncomp}{the list of number of components used in the model.}
  \item{call}{the function call.}
  \item{terms}{the model terms.}
  \item{model}{if \code{model = TRUE}, the model frame.}
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
  Currently, \code{lsplsCv} handles the formula and the data, and calls
  \code{\link{orthlsplsCv}} for the actual cross-validation.  The
  formula interface is experimental, and might change in future versions.  
}
\seealso{\code{\link{lspls}}, \code{\link{plot.lsplsCv}},
  \code{\link[pls]{cvsegments}}, \code{\link{orthlsplsCv}},
  \link{lspls-package}}
\examples{
##FIXME
}
\keyword{multivariate}
\keyword{regression}
