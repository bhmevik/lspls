### MSEP.R: MSEP and RMSEP functions.
### $Id$

## A small hack (MSEP should be made generic in pls):
## This must go into .First.lib
if (!exists("MSEP.default")) {
    MSEP.default <- MSEP
    MSEP <- function(object, ...) UseMethod("MSEP")
}

## MSEP takes a CV-object, and calculates the MSEP
MSEP.lsplsCv <- function(object, ...) {
    if (is.null(object$mode))
        stop("`object' has no `model' component.  Recalculate with `model = TRUE'")
    colMeans((object$pred - c(model.response(model.frame(object))))^2)
}

## A small hack (RMSEP should be made generic in pls, or the mvrVal
## object should be changed to be a matrix):
## This must go into .First.lib
if (!exists("RMSEP.default")) {
    RMSEP.default <- RMSEP
    RMSEP <- function(object, ...) UseMethod("RMSEP")
}

## RMSEP is a wrapper around MSEP that returns its square root.
RMSEP.lsplsCv <- function(object, ...) sqrt(MSEP(object, ...))
