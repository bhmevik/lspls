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
    matnames <- strsplit(attr(object$terms, "term.labels")[-1], ":")
    for (i in seq(along = object$scores)) {
        if (is.matrix(object$scores[[i]])) {
            scoreplot(object$scores[[i]], comps = 1:object$ncomp[[i]],
                      main = matnames[[i]][1], ...)
        } else {
            for (j in seq(along = object$scores[[i]])) {
                scoreplot(object$scores[[i]][[j]],
                          comps = 1:object$ncomp[[i]][j],
                          main = matnames[[i]][j], ...)
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
    matnames <- strsplit(attr(object$terms, "term.labels")[-1], ":")
    for (i in seq(along = object$loadings)) {
        if (is.matrix(object$loadings[[i]])) {
            loadingplot(object$loadings[[i]], comps = 1:object$ncomp[[i]],
                        main = matnames[[i]][1], ...)
        } else {
            for (j in seq(along = object$loadings[[i]])) {
                loadingplot(object$loadings[[i]][[j]],
                            comps = 1:object$ncomp[[i]][j],
                            main = matnames[[i]][j], ...)
            }
        }
    }
}


###
### Plot method for lsplsCv objects:
###
## FIXME: Should maybe be a plot method for (R)MSEP objects...

plot.lsplsCv <- function(x, which = c("RMSEP", "MSEP", "R2"), separate = TRUE,
                         scale = !isTRUE(separate), ...) {
    which <- match.arg(which)
    val <- do.call(which, list(object = x, scale = scale))
    if (!isTRUE(separate)) {
        ## Aggregate over the responses, but keep a dummy dimension for
        ## the response (it simplifies the code below):
        dims <- c(1, dim(val)[-1])
        dns <- c(resp = "all responses", dimnames(val)[-1])
        if (which == "R2") {
            val <- array(colMeans(val), dim = dims, dimnames = dns)
        } else {
            ## FIXME: For RMSEP: take sqrt(colSums(MSEP))?
            val <- array(colSums(val), dim = dims, dimnames = dns)
        }
    }
    comps <- expand.grid(lapply(dimnames(val)[-1], as.numeric))
    ncomps <- rowSums(comps)
    ncombs <- nrow(comps)
    complabels <- apply(comps, 1, paste, collapse = "")
    mXlab <- "total number of components"
    mYlab <- if (isTRUE(scale)) paste(which, "(std. resp.)") else which
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
    respnames <-  dimnames(val)[[1]]
    val <- aperm(val, c(2:length(dim(val)), 1)) # Make "resp" the last dimension
    for (i in 1:nResp) {
        cval <- c(val)[ncombs * (i - 1) + 1:ncombs]
        plot(ncomps, cval, type = "n", xlab = xlab, ylab = ylab,
             main = respnames[i], ...)
        text(ncomps, cval, labels = complabels)
        oncomps <- min(ncomps):max(ncomps)
        bestval <- numeric(length(oncomps))
        for (i in seq(along = oncomps))
            bestval[i] <- if (which == "R2") max(cval[ncomps == oncomps[i]])
                          else min(cval[ncomps == oncomps[i]])
        lines(oncomps, bestval, lty = 2, col = 2)
    } ## for
    if (nResp > 1) {
        ## Add outer margin text:
        mtext(mXlab, side = 1, outer = TRUE)
        mtext(mYlab, side = 2, outer = TRUE)
    }
} ## function
