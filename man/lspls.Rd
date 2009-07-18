\encoding{latin1}
\name{lspls}
\alias{lspls}
\title{Fit LS-PLS Models}
\description{
  A function to fit LS-PLS (least squares--partial least squares) models.
}
\usage{
lspls(formula, ncomp, data, subset, na.action, model = TRUE, \dots)
}
\arguments{
  \item{formula}{model formula.  See Details.}
  \item{ncomp}{list or vector of positive integers, giving the
    number of components to use for each \sQuote{pls-matrix}.  See Details.}
  \item{data}{an optional data frame with the data to fit the model from.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when
    the data contain missing values.}
  \item{model}{logical.  If \code{TRUE}, the model frame is returned.}
  \item{\dots}{additional arguments, passed to the underlying PLSR fit
    function.}
}
\details{
  \code{lspls} fits LS-PLS models, in which matrices are added
  successively to the model.  The first matrix is fit
  with ordinary least squares (LS) regression.  The rest of the matrices
  are fit with partial least squares regression (PLSR), using the
  residuals from the preceeding model as response.  See
  \link{lspls-package} or the references for more details, and
  \link{lspls-package} for typical usage.

  The model formula is specified as \samp{\var{resp} ~ \var{term1} +
    \var{term2} + \dots}.  If \var{resp} is a matrix (with more than one
  coloumn), a multi-response model is fitted.  \var{term1} specifies the
  first matrix to be fitted, using LS.  Each of the remaining terms will
  be added sequentially in the order specified in the formula (from left
  to right).  Each term can either be a single matrix, which will be added
  by itself, or several matrices separated with \code{:}, e.g.,
  \code{Z:V:W}, which will be added simultaneously (these will be denoted
  \dfn{parallell} matrices).

  The first matrix, \var{term1}, is called the \dfn{LS matrix}, and
  the rest of the predictor matrices (whether parallell or not) are
  called \dfn{PLS matrices}.
  
  Note that an intercept is \emph{not} automatically added to the model.
  It should be included as a constant coloumn in the LS matrix, if
  desired.  (If no intercept is included, the PLS matrices should be
  centered.  This happens automatically if the LS matrix includes the
  intercept.)

  The number of components to use in each of the PLSR models is
  specified with the \code{ncomp} argument, which should be a
  list.  Each element of the list gives the number of components to
  use for the corresponding term in the formula.  If the term specifies
  parallell matrices (separated with \code{:}), the list element
  should be a vector with one integer for each matrix.  Otherwise, it
  should be a number.

  To simplify the specification of \code{ncomp}, the following
  conversions are made:  if \code{ncomp} is a vector, it will be
  converted to a list.  \code{ncomp} will also be recycled as neccessary to get
  one element for each term.  Finally, for a parallell term, the list
  element will be recycled as needed.  Thus, \code{ncomp = 4} will
  result in 4 components being fit for every PLS matrix.

  Currently, the function \code{lspls} itself handles the formula and
  the data, and calls the underlying fit function
  \code{\link{orthlspls.fit}} to do the actual fitting.  This implements
  the orthogonalized version of the LS-PLS algorithm, and without splitting
  of parallell matrices into common and unique components (see
  the references).  Extensions to non-orthogonalized algorithms, and
  splitting of parallell matrices are planned.
}
\value{
  An object of class \code{"lspls"}.  The object contains all components
  returned by the underlying fit function (currently
  \code{\link{orthlspls.fit}}).  In addition, it contains the following
  components:
  \item{fitted.values}{matrix with fitted values, one coloumn per response}
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
  The user interface (e.g. the model handling) is experimental, and
  might well change in later versions.

  The handling of \code{formula} (especially \code{:}) is non-standard.
  Note that the order of the terms is significant; terms are added
  from left to right.
}
\seealso{\link{lspls-package}, \code{\link{lsplsCv}},
  \code{\link{plot.lspls}}}
\examples{
##FIXME
}
\keyword{multivariate}
\keyword{regression}
