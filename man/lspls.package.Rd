\encoding{latin1}
\name{lspls-package}
\alias{lspls-package}
\alias{lspls}
\docType{package}
\title{
LS-PLS Models
}
\description{
Functions to fit and validate LS-PLS models.  Includes predict and plot
methods.
}
\details{
\tabular{ll}{
Package: \tab lspls\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2005-12-02\cr
License: \tab GPL version 2\cr
}
LS-PLS models are written on the form
\eqn{Y = X\beta + S_1\gamma_1 + \dots + S_k\gamma_k + E},
where the terms \eqn{S_i} are one or more matrices \eqn{Z_{i,j}}
separated by a colon (:), e.g., \eqn{Z_{i,1}:Z_{i,2}:\dots:Z_{i,l_i}}.

The model is fitted from left to right.  First \eqn{Y} is fitted to
\eqn{X} using least squares (LS) regression and the residuals
calculated.
For each \eqn{i}, the the matrices \eqn{Z_{i,1}}, \dots, \eqn{Z_{i,k_i}}
are orthogonalised against the variables used in the regression sofar
(for \eqn{i = 1}, this means \eqn{X}).
The residuals from the LS regression are used as the response in PLS
regressions with the orthogonalised matrices as predictors (one PLS
regression for each matrix), and the desired number of PLS components
from each matrix are included among the LS prediction variables.
The LS regression is then updated and new residuals calculated.

The function to fit LS-PLS models is \code{lspls}.  A typical usage is
\code{mod <- lspls(y ~ X + Z + V1:V2 + W, ncomp = list(3, c(2,1), 2), data =
  mydata)}.  The first argument is the formula describing the model.
\code{X} is fit first, using LS.  Then PLS scores from \code{Z}
(orthogonalised) are added.  Then PLS scores from \code{V1} and
\code{V2} are added (simultaneously), and finally PLS scores from
\code{W}.  The next argument, \code{ncomp}, specifies the number of
components to use from each PLS: 3 \code{Z} score vectors, 2 \code{V1}
score vectors, 1 \code{V2} score vector and 2 \code{W} score vectors.
Finally, \code{mydata} should be a data frame with matrices \code{y},
\code{X}, \code{Z}, \code{V1}, \code{V2} and \code{W} (\code{y} can be a
vector).

Currently, score plots and loading plots of fitted models are
implemented.  \code{plot(mod, "scores")} gives score plots for each PLS
regression, and \code{plot(mod, "loadings")} gives loading plots.

There is a \code{predict} method to predict response or score values
from new data: \code{predict(mod, newdata = mynewdata)}.  (This predicts
response values.  Use \code{type = "scores"} to get scores.)

In order to determine the number of components to use from each matrix,
one can use cross-validation: \code{cvmod <- lsplsCv(y ~ X + Z + V1:V2 +
  W, ncomp = list(4, c(3,4), 3), segments = 12, data = mydata)}.  In
\code{lsplsCv}, \code{ncomp} gives the maximal number of components to
test.  The argument \code{segments} specifies the number of segments to
use.  One can specify the type of segments to use (random (default),
consequtive or interleaved) with the argument \code{segment.type}.
Alternatively, one can supply the segments explicitly with
\code{segments}.  See \code{\link{lsplsCv}} for details.

One can plot cross-validated RMSEP values with \code{plot(cvmod)}.
(Similarly, \code{plot(cvmod, "MSEP")} plots MSEP values.)  This makes
it easier to determine the optimal number of components for each PLS.
See \code{\link{plot.lsplsCv}} for details.  To calculate the RMSEP or
MSEP values explicitly, one can use the function \code{RMSEP} or
\code{MSEP}.
}
\author{
Bjørn-Helge Mevik

Maintainer: Bjørn-Helge Mevik <bhx5@mevik.net>
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
\keyword{ package }
\keyword{multivariate}
\keyword{regression}
\seealso{
  \code{\link{lspls}}, \code{\link{lsplsCv}}, \code{\link{plot.lspls}},
  \code{\link{plot.lsplsCv}}
}
\examples{
## FIXME
}
