\encoding{UTF-8}
\name{predict.lspls}
\alias{predict.lspls}
\title{Predict Method for LS-PLS Models}
\description{
  Predict method for \code{"lspls"} objects.  It predicts response
  values or scores from new data.
}
\usage{
\method{predict}{lspls}(object, newdata, type = c("response", "scores"),
        na.action = na.pass, \dots)
}
\arguments{
  \item{object}{object of class \code{"lspls"}.  The fitted model to
    predict with.}
  \item{newdata}{data frame.  The new data.}
  \item{type}{character.  Wether to predict responses or scores.}
  \item{na.action}{function determining what should be done with missing
    values in \code{newdata}.  The default is to predict \code{NA}.  See
    \code{\link{na.omit}} for alternatives.}
  \item{\dots}{further arguments.  Currently not used.}
}
%\details{
%}
\value{
  If \code{type = "response"}, a matrix with predicted response values
  is returned.  If \code{type = "scores"}, a matrix with predicted
  score values is returned.
}
\author{Bjørn-Helge Mevik}
\seealso{\code{\link{lspls}}}
\examples{
##FIXME
}
\keyword{multivariate}
\keyword{regression}
