\encoding{latin1}
\name{plot.lspls}
\alias{plot.lspls}
\alias{scoreplot.lspls}
\alias{loadingplot.lspls}
\title{Plots of LS-PLS Models}
\description{
  Plot method for \code{"lspls"} objects.
}
\usage{
\method{plot}{lspls}(x, plottype = c("scores", "loadings"), \dots)
\method{scoreplot}{lspls}(object, \dots)
\method{loadingplot}{lspls}(object, \dots)
}
\arguments{
  \item{x, object}{Object of class \code{"lspls"}.  The model to be plotted.}
  \item{plottype}{character string.  What type of plot to generate.}
  \item{\dots}{Further arguments, passed on to underlying plot functions.}
}
\details{
  The \code{plot} method simply calls \code{scoreplot.lspls} or
  \code{loadingplot.lspls} depending on the \code{plottype} argument.

  \code{scoreplot.lspls} gives a series of score plots, one for each PLS
  model.  The user is asked to press Return between each plot.

  \code{loadingplot.lspls} shows a series of loading plots, one for each
  PLS model.  All plots are shown in the same plot window.
}
\value{
  The functions return whatever the (last) underlying plot function returns.
}
\author{Bjørn-Helge Mevik}
\seealso{\code{\link{lspls}}, \code{\link[pls]{scoreplot}},
  \code{\link[pls:scoreplot]{loadingplot}}, \code{\link{plot.lsplsCv}}}
\examples{
##FIXME
}
\keyword{hplot}
