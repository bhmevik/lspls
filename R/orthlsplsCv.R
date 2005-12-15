### orthlsplsCv.R: Cross-validation, orthogonalizing version
###
### $Id$

## The algorithm is based on recursion, after X has been handled.

orthlsplsCv <- function(Y, X, Z, A, method = getOption("pls.algorithm"),
                        segments = 10, segment.type = c("random",
                                       "consecutive", "interleaved"),
                        length.seg) {

    ## The recursive function:
    ## It uses the following variables from orthlsplsCv
    ## - Z: list of spectral matrices
    ## - A: list of #comps to use in the CV
    ## - segment: indices of the segment to be predicted
    ## - cvPreds: array of predictions; dim: c(nObs, nResp, unlist(A))
    ## - pls.fit: the pls fit function
    cvPredRest <- function(indices, prevCalib, prevPred, prevComps, prevRes) {
        ## indices is the indices of the remaining matrices
        ## prevCalib is a matrix with the X vars and scores used in the
        ##   previous calibrations
        ## prevPred is a matrix with the X vars and scores used in the
        ##   previous predictions
        ## prevComps is the numbers of components used in the previous
        ##   calibrations
        ## prevRes is the residuals from the previous calibrations

        ## The general idea is to handle the first matrix or list of
        ## matrices (Z[[indices[1]]]), and then recall itself on the rest of
        ## the matrices.

        ## The matrix/matrices to handle in this call:
        ind <- indices[1]
        M <- Z[[ind]]
        if (is.matrix(M)) {             # A single matrix
            ## Orthogonalise the calibration spectra wrt. prevVars
            Mcal <- M[-segment,]
            Mo <- orth(Mcal, prevCalib)

            ## Orthogonalise the prediction spectra
            Mpred <- M[segment,]
            Mpo <- Mpred - prevPred %*% Corth(Mcal, prevCalib)
            ## mal: Zorig[i,] - Xorig[i,] %*% Co(Xorig[-i,]) %*% Zorig[-i,]

            ## Estimate a model prevRes ~ orth. spectra + res
            plsM <- pls.fit(Mo, prevRes, A[[ind]])
            ## Save scores:
            calScores <- plsM$scores

            ## Predict new scores and response values
            predScores <- sweep(Mpo, 2, plsM$Xmeans) %*% plsM$projection
            ## FIXME: Only for orth.scores?
            predVals <- array(dim = c(nrow(predScores), dim(plsM$Yloadings)))
            for (a in 1:A[[ind]])
                predVals[,,a] <-
                    sweep(predScores[,1:a] %*% t(plsM$Yloadings[,1:a, drop=FALSE]),
                          2, plsM$Ymeans, "+")

            ## Add the predictions to the outer cvPreds variable
            ## Alt. 1:  Calculate the 1-index indices manuall, and use
            ## single indexing (probably quickest, but requires a loop).
            ## Alt. 2:  Use matrix indexing with an expanded grid.
            ##eg <- expand.grid(segment, 1:nResp, A[[ind]])
            ##indMat <- do.call("cbind", c(eg[1:2], as.list(prevComps), eg[3]))
            ##cvPreds[indMat] <- cvPreds[indMat] + predVals
            ## Alt. 3: Build and eval an expression which does what we want:
            ncomps <- length(prevComps)
            nrest <- length(dim(cvPreds)) - ncomps - 2
            dummy <- Quote(cvPreds[segment,])
            dummy[4 + seq(along = prevComps)] <- prevComps
            dummy[4 + ncomps + 1:nrest] <- dummy[rep(4, nrest)]
            eval(substitute(dummy <<- dummy + c(predVals), list(dummy = dummy)))

            ## Return if this is the last matrix/set of matrices
            if (length(indices) == 1) return()

            ## Calculate new residuals
            newResid <- - plsM$fitted.values + c(prevRes)

            ## To save space: drop the model object(s)
            rm(plsM)

            ## Recursively call ourself for each number of components in the
            ## present model
            for (i in seq(length = A[[ind]]))
                Recall(indices[-1], # Remove the index of the current matrix
                       cbind(prevCalib, calScores[,1:i]), # Add the scores we've used
                       cbind(prevPred, predScores[,1:i, drop=FALSE]), # Add the scores we've predicted
                       c(prevComps, i), # and the number of comps
                       newResid[,,i]) # update the residual

        } else {                        # List of parallell matrices
            Scal <- list()              # The current calibration scores
            Spred <- list()             # The current prediction scores
            for (j in seq(along = M)) {
                ## Orthogonalise the calibration spectra wrt. prevVars
                Mcal <- M[[j]][-segment,]
                Mo <- orth(Mcal, prevCalib)

                ## Orthogonalise the prediction spectra
                Mpred <- M[[j]][segment,]
                Mpo <- Mpred - prevPred %*% Corth(Mcal, prevCalib)

                ## Estimate a model prevRes ~ orth. spectra + res
                plsM <- pls.fit(Mo, prevRes, A[[ind]][j])
                ## Save scores:
                Scal[[j]] <- plsM$scores

                ## Predict new scores
                Spred[[j]] <- sweep(Mpo, 2, plsM$Xmeans) %*% plsM$projection
            }
            ## To save space: drop the model object
            rm(plsM)

            ## Loop over the different combinations of #comps:
            nComps <- expand.grid(lapply(as.list(A[[ind]]), seq))
            for (cind in 1:nrow(nComps)) {
                newComps <- nComps[cind,]
                comps <- c(prevComps, unlist(newComps))
                ## Predict new response values
                calScores <- mapply(function(B, b) B[,1:b], Scal, newComps)
                predScores <- mapply(function(B, b) B[,1:b], Spred, newComps)
                lsS <- lm.fit(calScores, prevRes) # FIXME: How about intercept?
                newResid <- lsS$residuals
                predVals <- predScores %*% lsS$coefficients

                ## Add the predictions to the outer cvPreds variable.  Build
                ## and eval an expression which does what we want:
                nc <- length(comps)
                nrest <- length(dim(cvPreds)) - nc - 2
                dummy <- Quote(cvPreds[segment,])
                dummy[4 + seq(along = comps)] <- comps
                dummy[4 + nc + 1:nrest] <- dummy[rep(4, nrest)]
                eval(substitute(dummy <<- dummy + c(predVals),
                                list(dummy = dummy)))

                if (length(indices) > 1) { # There are more matrices to fit
                    ## Recursively call ourself
                    Recall(indices[-1], # Remove the index of the current matrices
                           cbind(prevCalib, calScores), # Add the scores we've used
                           cbind(prevPred, predScores), # Add the scores we've predicted
                           c(comps), # and the number of comps
                           newResid) # use the new residual

                }
            } ## for
        } ## if
    } ## recursive function

    ## Setup:
    nObs <- nrow(X)
    nResp <- ncol(Y)
    ## cvPreds: the cross-validated predictions:
    cvPreds <- array(0, dim = c(nObs, nResp, unlist(A)))
    ## Build an unevaluated expression that will insert the predictions into
    ## cvPreds[segment,,...,] when evaluated:
    ndim <- length(dim(cvPreds))
    ## This creates an expression with the empty index argument repeated as
    ## many times as neccessary:
    dummy <- Quote(cvPreds[segment,])[c(1:3, rep(4, ndim - 1))]
    ## Substitute this in an assignment statement:
    addPredictions <- substitute(dummy <- predVals, list(dummy = dummy))

    ## FIXME: Maybe in outer function:
    method <- match.arg(method, c("oscorespls", "kernelpls", "simpls"))
    pls.fit <- switch(method, oscorespls = oscorespls.fit,
                      kernelpls = kernelpls.fit, simpls = simpls.fit)

    ## Set up segments:
    if (is.list(segments)) {
        if (is.null(attr(segments, "type")))
            attr(segments, "type") <- "user supplied"
    } else {
        if (missing(length.seg)) {
            segments <- cvsegments(nObs, k = segments, type = segment.type)
        } else {
            segments <- cvsegments(nObs, length.seg = length.seg,
                                   type = segment.type)
        }
    }

    ## The main cross-validation loop
    temp <- 0
    for (segment in segments) {
        cat(temp <- temp + 1, "")
        ## Handle X
        lsX <- lm.fit(X[-segment,, drop = FALSE], Y[-segment,, drop = FALSE])
        resid <- lsX$residuals
        predVals <- X[segment,, drop = FALSE] %*% lsX$coefficients

        ## Insert the predictions into the cvPred array:
        eval(addPredictions)

        ## Handle the rest of the matrices:
        cvPredRest(indices = 1:length(unlist(A)),
                   prevCalib = X[-segment,, drop = FALSE],
                   prevPred = X[segment,, drop = FALSE],
                   prevComps = c(),
                   prevRes = resid)
    }
    return(cvPreds)
} ## function
