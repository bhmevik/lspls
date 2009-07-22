\encoding{latin1}
\name{plot.lsplsCv}
\alias{plot.lsplsCv}
\title{Plot Method for Cross-Validations}
\description{
  Plot method for \code{"lsplsCv"} objects.  It plots the
  cross-validated (R)MSEP or R^2 against the total number of components
  or the matrices included in the model.
}
\usage{
\method{plot}{lsplsCv}(x, which = c("RMSEP", "MSEP", "R2"), ncomp,
        separate = TRUE, scale = !isTRUE(separate), \dots)
}
\arguments{
  \item{x}{object of class \code{"lsplsCv"}.  Object to be plotted.
    Typically the output from \code{\link{lsplsCv}}.}
  \item{which}{character string.  Which measure to plot.}
  \item{ncomp}{list.  The number of components to use when plotting, for
    each PLS matrix in the model.  See Details.}
  \item{separate}{logical.  Whether separate plots should be generated
    for each response (default) or one plot with the sum of the measure
    for all responses.}
  \item{scale}{logical.  Whether the responses and predicted values
    should be divided by the standard deviation of the response prior to
    calculating the measure.  Default is to scale when producing a
    combined plot (\code{separate = FALSE}) and not to scale otherwise.}
  \item{\dots}{Further arguments, sent to the underlying plot function.}
}
\details{
  If \code{ncomp} is not specified, 
  the \code{plot} method generates a plot of the cross-validated (R)MSEP
  or R^2 values for all combinations of number of components.  The
  values are plotted against the total number of components.  Each point
  is labelled with the combination of number of components.  E.g., for
  a model with three PLS matrices, \samp{132} means one
  component from the first matrix, three from the second and two from
  the third. 
  Also, the lowest (R)MSEP or highest \eqn{R^2} values for each total
  number of components are joined by a line.
  
  If \code{ncomp} is specified, the \code{plot} method plots (R)MSEP
  or R^2 for models with the first matrix, with the two first matrices,
  etc.  \code{ncomp} should be specified as when running \code{lsplsCv},
  and is used for selecting the number of components for each PLS
  matrix.  For instance, if \code{mod <- lsplsCv(Y ~ X + Z + V:W, ...)},
  then \code{plot(mod, ncomp = list(2, c(1,3)))} would plot the RMSEPs
  for \code{Y ~ X}, \code{Y ~ X + Z} and \code{Y ~ X + Z + V:W}, using
  2, 1 and 3 components for \code{Z}, \code{V} and \code{W}, respectively.

  If \code{separate} is \code{TRUE}, a separate plot panel is produced
  for each response.  Otherwise the measure is added for all responses
  and shown in one plot.  If \code{scale} is \code{TRUE} (the default
  when producing a combined plot), the measures for each response are
  standardised by dividing the responses and predicted values by the
  standard deviation of the (corresponding) response prior to
  calculating the measure.  Note that \code{scale} is ignored when
  \code{which} is \code{"R2"} because \eqn{R^2} is independent of scale.)
}
\value{
  The function returns whatever the (last) underlying plot function
  returns.
}
\author{Bjørn-Helge Mevik}
\seealso{\code{\link{lsplsCv}}, \code{\link{lspls}}}
\examples{
##FIXME
}
\keyword{hplot}
