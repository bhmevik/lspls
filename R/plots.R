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


###
### Plot method for lsplsCv objects:
###
## FIXME: Should maybe be a plot method for (R)MSEP objects...

plot.lsplsCv <- function(x, which = c("RMSEP", "MSEP"), ...) {
    which <- match.arg(which)
    val <- do.call(which, list(x))
    comps <- expand.grid(lapply(dimnames(val)[-1], as.numeric))
    ncomps <- rowSums(comps)
    ncombs <- nrow(comps)
    complabels <- apply(comps, 1, paste, collapse = "")
    mXlab <- "total number of components"
    mYlab <- which
    nResp <- dim(val)[1]
    if (nResp > 1) {
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        par(mfrow = n2mfrow(nResp), oma = c(1, 1, 0, 0) + 0.1,
            mar = c(3, 3, 3, 1) + 0.1)
        xlab <- ""
        ylab <- ""
    } else {
        xlab <- mXlab
        ylab <- mYlab
    }
    val <- aperm(val, c(2:length(dim(val)), 1)) # Make "resp" the last dimension
    for (i in 1:nResp) {
        cval <- c(val)[ncombs * (i - 1) + 1:ncombs]
        plot(ncomps, cval, type = "n", xlab = xlab, ylab = ylab, main = i, ...)
        text(ncomps, cval, labels = complabels)
        oncomps <- min(ncomps):max(ncomps)
        minval <- numeric(length(oncomps))
        for (i in seq(along = oncomps))
            minval[i] <- min(cval[ncomps == oncomps[i]])
        lines(oncomps, minval, lty = 2, col = 2)
    } ## for
    if (nResp > 1) {
        ## Add outer margin text:
        mtext(mXlab, side = 1, outer = TRUE)
        mtext(mYlab, side = 2, outer = TRUE)
    }
} ## function
