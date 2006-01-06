### MSEP.R: MSEP and RMSEP functions.
### $Id$


## MSEP takes a CV-object, and calculates the MSEP
MSEP.lsplsCv <- function(object, ...) {
    if (is.null(object$mode))
        stop("`object' has no `model' component.  Recalculate with `model = TRUE'")
    colMeans((object$pred - c(model.response(model.frame(object))))^2)
}


## RMSEP is a wrapper around MSEP that returns its square root.
RMSEP.lsplsCv <- function(object, ...) sqrt(MSEP(object, ...))
