\encoding{latin1}
\name{lspls}
\alias{lspls}
\title{Fit LS-PLS Models}
\description{
  A function to fit LS-PLS (least squares-partial least squares) models.
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
  FIXME: Re-read.
  
  \code{lspls} fits LS-PLS models, in which matrices are added
  successively(FIXME:?) to the model.  The first matrix is fit
  with ordinary least squares (LS) regression.  The rest of the matrices
  are fit with partial least squares regression (PLSR), using the
  residuals from the preceeding model as response.  See the references
  for more details.

  The model is specified as \samp{\var{resp} ~ \var{term1} + \var{term2}
    + \dots}.  If \var{resp} is a matrix (with more than two coloumns),
  a multi-response model is fitted.  \var{term1} specifies the first
  matrix to be fitted, using LS.  Each of the remaining terms will be
  added sequentially in the order specified in the formula (from left to
  right).  Each term can either be a single matrix, which will be added
  by itself, or several matrices separated with \code{:}, e.g.,
  \code{Z:V:W}, which will be added simultaneously (this will be denoted
  \term{parallell} matrices).  Note that an intercept is \emph{not}
  automatically.  It should be included as a constant coloumn in
  \var{term1}, if desired.  (If no intercept is included, the rest of
  the matrices should be centered.  This happens automatically if
  \var{term1} includes the intercept.)

  The number of components to use in each of the PLSR models, is
  specified with the \code{ncomp} argument.  \code{ncomp} should be a
  list.  Each element of the list gives the number of components to
  use for the corresponding term in the formula.  If the term specifies
  parallell matrices (separated with \code{:}), the list element
  should be a vector with one integer for each matrix.  Otherwise, it
  should be a number.

  To simplify the specification of \code{ncomp}, the following
  conversions are made:  if \code{ncomp} is a vector, it will be
  converted to a list.  It will also be recycled as neccessary to get
  one element for each term.  Finally, for a parallell term, the list
  element will be recycled as needed.  Thus, \code{ncomp = 4} will
  result in 4 components being fit for every matrix.

  Currently, the function \code{lspls} itself handles the formula and
  the data, and calls the underlying fit function
  \code{\link{orthlspls.fit}} to do the actual fitting.  This implements
  the orthogonalized version of the LS-PLS algorithm, and without splitting
  of parallell matrices into common and unique components (see
  References).  Extensions to non-orthogonalized algorithms, and
  splitting of parallell matrices are planned.
}
\value{
  An object of class \code{"lspls"}.  The object contains all components
  returned by the underlying fit function.  In addition, it contains the
  following components:
  \item{na.action}{if observations with missing values were removed,
    \code{na.action} contains a vector with their indices.}
  \item{ncomp}{the list of number of components used in the model.}
  \item{call}{the function call.}
  \item{terms}{the model terms.}
  \item{model}{if \code{model = TRUE}, the model frame.}
}
\references{ FIXME: Jørgensen publ, Jørgensen subm, Mevik subm. }
\author{Bjørn-Helge Mevik}
\note{
  The user interface (e.g. the model handling) is experimental, and
  might well change in later versions.

  The handling of \code{model} (especially \code{:}) is non-standard.
  Note that the order of the terms is significant; matrices are added
  from left to right.
}
\seealso{\code{\link{lsplsCv}}, \code{\link{plot.lspls}}}
\examples{
##FIXME
}
\keyword{multivariate}
\keyword{regression}
