.First.lib <- function(libname, pkgname) {
    ## A small hack (MSEP should be made generic in pls):
    if (!exists("MSEP.default")) {
        MSEP.default <<- MSEP
        MSEP <<- function(object, ...) UseMethod("MSEP")
    }
    ## A small hack (RMSEP should be made generic in pls, or the mvrVal
    ## object should be changed to be a matrix):
    if (!exists("RMSEP.default")) {
        RMSEP.default <<- RMSEP
        RMSEP <<- function(object, ...) UseMethod("RMSEP")
    }
    ## Idea: Make `scoreplot' in pls generic, with methods for matrix,
    ## scores(?), lspls and default (anything that has a 'scores' method that
    ## gives a single matrix).
    ## Dirty hack:
    if (!exists("scoreplot.default")) {
        scoreplot.default <<- scoreplot
        scoreplot <<- function(object, ...) UseMethod("scoreplot")
    }
    ## Idea: Make `loadingplot' in pls generic, with methods for matrix,
    ## loadings(?), lspls and default (anything that has a 'loadings' method that
    ## gives a single matrix).
    ## Dirty hack:
    if (!exists("loadingplot.default")) {
        loadingplot.default <<- loadingplot
        loadingplot <<- function(object, ...) UseMethod("loadingplot")
    }
}
