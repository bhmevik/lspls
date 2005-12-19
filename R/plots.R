### plots.R:  Plot functions
### $Id$

###
### Plot method for lspls objects
###

plot.lspls <- function(x, plottype = c("scores", "loadings"), ...) {
    plottype <- match.arg(plottype)
    plotFunc <- switch(plottype,
                       scores = scoreplot.lspls,
                       loadings = loadingplot.lspls)
    plotFunc(x, ...)
}


###
### Scoreplot
###

scoreplot.lspls <- function(object, ...) {
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(ask = TRUE)
    for (i in seq(along = object$scores)) {
        if (is.matrix(object$scores[[i]])) {
            scoreplot(object$scores[[i]], comps = 1:object$ncomp[[i]], main = i, ...)
        } else {
            for (j in seq(along = object$scores[[i]])) {
                scoreplot(object$scores[[i]][[j]], comps = 1:object$ncomp[[i]][j], main = paste(i, j, sep = "."), ...)
            }
        }
    }
}


###
### Loadingplot
###

## Idea: Make `loadingplot' in pls generic, with methods for matrix,
## loadings(?), lspls and default (anything that has a 'loadings' method that
## gives a single matrix).

loadingplot.lspls <- function(object, ...) {
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(mfrow = n2mfrow(length(unlist(object$ncomp))))
    for (i in seq(along = object$loadings)) {
        if (is.matrix(object$loadings[[i]])) {
            loadingplot(object$loadings[[i]], comps = 1:object$ncomp[[i]], main = i, ...)
        } else {
            for (j in seq(along = object$loadings[[i]])) {
                loadingplot(object$loadings[[i]][[j]], comps = 1:object$ncomp[[i]][j], main = paste(i, j, sep = "."), ...)
            }
        }
    }
}
