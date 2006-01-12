\encoding{latin1}
\name{plot.lsplsCv}
\alias{plot.lsplsCv}
\title{Plot Method for Cross-Validations}
\description{
  Plot method for \code{"lsplsCv"} objects.  It plots the
  cross-validated (R)MSEP against the total number of components.
}
\usage{
\method{plot}{lsplsCv}(x, which = c("RMSEP", "MSEP"), \dots)
}
\arguments{
  \item{x}{object of class \code{"lsplsCv"}.  Object to be plotted.
    Typically the output from \code{\link{lsplsCv}}.}
  \item{which}{character string.  Which measure to plot.}
  \item{\dots}{Further arguments, sent to the underlying plot function.}
}
\details{
  The \code{plot} method generates a plot of the cross-validated (R)MSEP
  values for all combinations of number of components.  The
  values are plotted against the total number of components.  Each point
  is labelled with the combination of number of components.  E.g., for
  a model with three PLS matrices, \samp{132} means one
  component from the first matrix, three from the second and two from
  the third.

  Also, the lowest (R)MSEP values for each total number of components
  are joined by a line.
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
