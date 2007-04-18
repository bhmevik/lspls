###
### Orthogonal case, fitting of model
###

orig.orthlspls.fit <- function(Y, X, Z, ncomp) {
    ## Parametres:
    nObs <- nrow(X)
    totNumComps <- sum(unlist(ncomp))
    totNumCoefs <- ncol(X) + totNumComps
    if(totNumCoefs > nObs) stop("Too many variables/components selected.")

    ## Containers:
    B <- matrix(nrow = totNumCoefs, ncol = ncol(Y)) # Regr. coefs.
    V <- matrix(nrow = nObs, ncol = totNumCoefs) # Regr. variables
    models <- list()                    # plsr models
    orthCoefs <- list()                 # Orthogonalisation matrices
    ## These two are not strictly neccessary:
    S <- list()                         # Scores
    L <- list()                         # Loadings


    ## Choose PLS fit algorithm.  FIXME: Support other algs?
    pls.fit <- oscorespls.fit

    ## Start with X:
    lsX <- lm.fit(X, Y)
    nVar <- ncol(X)
    ## Extract
    V[,1:nVar] <- X
    B[1:nVar,] <- lsX$coefficients
    res <- lsX$residuals

    ## For testing:
    ##Balt <- B

    ## Walk through Z:
    for (i in seq(along = Z)) {
        ##cat("i =", i, "\n")
        M <- Z[[i]]
        if (is.matrix(M)) {             # Single matrix
            Mo <- orth(M, V[,1:nVar])   # Orth. M against all used variables
            orthCoefs[[i]] <- Corth(M, V[,1:nVar]) # For pred
            models[[i]] <- pls.fit(Mo, res, ncomp[[i]])# Could use Y
            V[,nVar + (1:ncomp[[i]])] <- S[[i]] <- models[[i]]$scores
            L[[i]] <- models[[i]]$loadings
            ## FIXME: Testing:
            ##lmZ <- lm.fit(models[[i]]$scores, res)
            ##Balt[nVar + (1:ncomp[[i]]),] <- lmZ$coefficients
            ## FIXME: This depends on the pls algorithm (at least the scaling):
            B[nVar + (1:ncomp[[i]]),] <- t(models[[i]]$Yloadings)
            res <- models[[i]]$residuals[,,ncomp[[i]]]
            nVar <- nVar + ncomp[[i]]
        } else {                        # Parallell matrices
            Vadd <- matrix(nrow = nObs, ncol = sum(ncomp[[i]])) # The variables to be added in the present step
            S[[i]] <- list()
            L[[i]] <- list()
            orthCoefs[[i]] <- list()
            models[[i]] <- list()
            if (length(M) == 2 && length(ncomp[[i]]) == 3) {
                ## Split information into common and unique components
                ## FIXME: Save the orthCoefs and models that are needed for
                ## prediction!
                nZ1u <- ncomp[[i]][[1]]
                nZ2u <- ncomp[[i]][[2]]
                nC <- ncomp[[i]][[3]]
                ## Step 2a: pls RvY ~ OvZ1:
                OvZ1 <- orth(M[[1]], V[,1:nVar])
                CvZ1 <- Corth(M[[1]], V[,1:nVar]) # For pred
                plsRvYOvZ1 <- pls.fit(OvZ1, res, nZ1u + nC)
                SZ1 <- plsRvYOvZ1$scores
                Rvz1 <- plsRvYOvZ1$residuals[,,nZ1u + nC]
                ## Step 2b: pls RvY ~ OvZ2:
                OvZ2 <- orth(M[[2]], V[,1:nVar])
                CvZ2 <- Corth(M[[2]], V[,1:nVar]) # For pred
                plsRvYOvZ2 <- pls.fit(OvZ2, res, nZ2u + nC)
                SZ2 <- plsRvYOvZ2$scores
                Rvz2 <- plsRvYOvZ2$residuals[,,nZ2u + nC]
                ## Step 3a: pls Rvz2Y ~ Oxz2Z1 (get unique scores):
                Ovz2Z1 <- orth(OvZ1, SZ2) # equiv: orth(Z1, cbind(V, SZ2))
                Cvz2Z1 <- Corth(OvZ1, SZ2) # For prediction
                plsRvz2YOvz2Z1 <- pls.fit(Ovz2Z1, Rvz2, nZ1u)
                SZ1u <- plsRvz2YOvz2Z1$scores
                Vadd[,1:nZ1u] <- S[[i]][[1]] <- SZ1u
                L[[i]][[1]] <- plsRvz2YOvz2Z1$loadings
                ## Step 3b: pls Rvz1Y ~ Oxz1Z2 (get unique scores):
                Ovz1Z2 <- orth(OvZ2, SZ1) # equiv: orth(Z2, cbind(V, SZ1))
                Cvz1Z2 <- Corth(OvZ2, SZ1) # For prediction
                plsRvz1YOvz1Z2 <- pls.fit(Ovz1Z2, Rvz1, nZ2u)
                SZ2u <- plsRvz1YOvz1Z2$scores
                Vadd[,nZ1u + 1:nZ2u] <- S[[i]][[2]] <- SZ2u
                L[[i]][[2]] <- plsRvz1YOvz1Z2$loadings
                ## Step 4: Calculate common components:
                SU <- cbind(SZ1u, SZ2u)
                SZ1c <- orth(SZ1, SU)
                CSZ1c <- Corth(SZ1, SU) # For prediction
                SZ2c <- orth(SZ2, SU)
                CSZ2c <- Corth(SZ2, SU) # For prediction
                ## Use cancor to get one set of scores:
                ## FIXME: Must handle matrices without full col. rank.
                ## See e.g. rank.condition(corpcor) or is.fullrank(limma)
                ccaC <- cancor(SZ1c, SZ2c)
                SC <- (SZ1c %*% ccaC$xcoef[,1:nC] +
                       SZ2c %*% ccaC$ycoef[,1:nC]) / 2
                Vadd[,nZ1u + nZ2u + 1:nC] <- S[[i]][[3]] <- SC
                L[[i]][[3]] <- PZ1c <- crossprod(OvZ1, SC)
                ## NB: crossprod(Ovz2Z1, SC) blir _ikke_ riktig, da Ovz2Z1
                ## er ortogonal til Z2 (egentlig SZ2), og det er ikke
                ## ønskelig for "felles" ladninger.
                L[[i]][[4]] <- PZ2c <- crossprod(OvZ2, SC)
            } else {
                ## Don't extract common/unique components
                added <- 0              # The number of comps. currently added
                for (j in seq(along = M)) {
                    ##cat("j =", j, "\n")
                    ## Walk through Z[[i]]
                    Mo <- orth(M[[j]], V[,1:nVar])
                    orthCoefs[[i]][[j]] <- Corth(M[[j]], V[,1:nVar]) # For pred
                    models[[i]][[j]] <- pls.fit(Mo, res, ncomp[[i]][[j]])# Could use Y
                    Vadd[,added + (1:ncomp[[i]][[j]])] <- S[[i]][[j]] <- models[[i]][[j]]$scores
                    L[[i]][[j]] <- models[[i]][[j]]$loadings
                    added <- added + ncomp[[i]][[j]]
                }
            } # if
            V[,nVar + (1:sum(ncomp[[i]]))] <- Vadd
            ## Only strictly neccessary for unique/common:
            lmZ <- lm.fit(Vadd, res)
            B[nVar + (1:sum(ncomp[[i]])),] <- lmZ$coefficients
            ##Balt[nVar + (1:sum(ncomp[[i]])),] <- lmZ$coefficients
            res <- lmZ$residuals
            nVar <- nVar + sum(ncomp[[i]])
        } # if
    } # for
    list(coefficients = B, predictors = V, orthCoefs = orthCoefs,
         models = models, ncomp = ncomp, scores = S, loadings = L, residuals = res)
} # function

## New version of common-components algorithm
orthlspls.fit <- function(Y, X, Z, ncomp) {
    ## Parametres:
    nObs <- nrow(X)
    ## Calculate the total number of components:
    tmp <- sapply(ncomp, function(x) switch(1+length(x), NULL, x, sum(x),
                                            sum(x), NULL, sum(x[1:3])))
    ## Check for valid lenghts in ncomp:
    if (any(sapply(tmp, is.null))) stop("Invalid ncomp")
    totNumComps <- sum(tmp)
    totNumCoefs <- ncol(X) + totNumComps
    if(totNumCoefs > nObs) stop("Too many variables/components selected.")

    ## Containers:
    B <- matrix(nrow = totNumCoefs, ncol = ncol(Y)) # Regr. coefs.
    V <- matrix(nrow = nObs, ncol = totNumCoefs) # Regr. variables
    models <- list()                    # plsr models
    orthCoefs <- list()                 # Orthogonalisation matrices
    ## These two are not strictly neccessary:
    S <- list()                         # Scores
    L <- list()                         # Loadings


    ## Choose PLS fit algorithm.  FIXME: Support other algs!
    pls.fit <- oscorespls.fit

    ## Start with X:
    lsX <- lm.fit(X, Y)
    nVar <- ncol(X)
    ## Extract
    V[,1:nVar] <- X
    B[1:nVar,] <- lsX$coefficients
    res <- lsX$residuals

    ## For testing:
    ##Balt <- B

    ## Walk through Z:
    for (i in seq(along = Z)) {
        ##cat("i =", i, "\n")
        M <- Z[[i]]
        if (is.matrix(M)) {             # Single matrix
            Mo <- orth(M, V[,1:nVar])   # Orth. M against all used variables
            orthCoefs[[i]] <- Corth(M, V[,1:nVar]) # For pred
            models[[i]] <- pls.fit(Mo, res, ncomp[[i]])# Could use Y
            V[,nVar + (1:ncomp[[i]])] <- S[[i]] <- models[[i]]$scores
            L[[i]] <- models[[i]]$loadings
            ## FIXME: Testing:
            ##lmZ <- lm.fit(models[[i]]$scores, res)
            ##Balt[nVar + (1:ncomp[[i]]),] <- lmZ$coefficients
            ## FIXME: This depends on the pls algorithm (at least the scaling):
            B[nVar + (1:ncomp[[i]]),] <- t(models[[i]]$Yloadings)
            res <- models[[i]]$residuals[,,ncomp[[i]]]
            nVar <- nVar + ncomp[[i]]
        } else {                        # Parallell matrices
            S[[i]] <- list()
            L[[i]] <- list()
            orthCoefs[[i]] <- list()
            models[[i]] <- list()
            if (length(M) == 2 && length(ncomp[[i]]) > 2) {
                ## Split information into common and unique components
                ## FIXME: Save the orthCoefs and models that are needed for
                ## prediction!
                nZ1u <- ncomp[[i]][1]
                nZ2u <- ncomp[[i]][2]
                nC <- ncomp[[i]][3]
                if (length(ncomp[[i]]) > 3) {
                    ## The user has specified how many Z1 and Z2 components
                    ## to use in the CCA.
                    nZ1c <- ncomp[[i]][4]
                    nZ2c <- ncomp[[i]][5]
                } else {
                    ## Use the same number of components as the number of
                    ## unique components
                    nZ1c <- nZ1u
                    nZ2c <- nZ2u
                }
                Vadd <- matrix(nrow = nObs, ncol = nC + nZ1u + nZ2u) # The variables to be added in the present step

                ## Calculate common components:
                ## pls RvY ~ OvZ1:
                OvZ1 <- orth(M[[1]], V[,1:nVar])
                ##CvZ1 <- Corth(M[[1]], V[,1:nVar]) # For pred
                plsRvYOvZ1 <- pls.fit(OvZ1, res, nZ1c)
                SZ1c <- plsRvYOvZ1$scores
                ##Rvz1 <- plsRvYOvZ1$residuals[,,nZ1c]
                ## pls RvY ~ OvZ2:
                OvZ2 <- orth(M[[2]], V[,1:nVar])
                ##CvZ2 <- Corth(M[[2]], V[,1:nVar]) # For pred
                plsRvYOvZ2 <- pls.fit(OvZ2, res, nZ2c)
                SZ2c <- plsRvYOvZ2$scores
                ##Rvz2 <- plsRvYOvZ2$residuals[,,nZ2c]
                ## Use cancor to get one set of scores:
                ## FIXME: Must handle matrices without full col. rank.
                ## See e.g. rank.condition(corpcor) or is.fullrank(limma)
                ccaC <- cancor(SZ1c, SZ2c)
                SC <- (SZ1c %*% ccaC$xcoef[,1:nC] +
                       SZ2c %*% ccaC$ycoef[,1:nC]) / 2
                Vadd[,nZ1u + nZ2u + 1:nC] <- S[[i]][[3]] <- SC
                L[[i]][[3]] <- PZ1c <- crossprod(OvZ1, SC)
                L[[i]][[4]] <- PZ2c <- crossprod(OvZ2, SC)

                ## Calculate unique components:
                ## pls RvY ~ OvcZ1 (equiv: RvcY ~ OvcZ1):
                OvcZ1 <- orth(OvZ1, SC) # equiv: orth(Z1, cbind(V, SC))
                CvcZ1 <- Corth(OvZ1, SC) # For prediction
                plsRvYOvcZ1 <- pls.fit(OvcZ1, res, nZ1u)
                SZ1u <- plsRvYOvcZ1$scores
                Vadd[,1:nZ1u] <- S[[i]][[1]] <- SZ1u
                L[[i]][[1]] <- plsRvYOvcZ1$loadings
                ## pls RvY ~ OvcZ2 (equiv: RvcY ~ OvcZ2):
                OvcZ2 <- orth(OvZ2, SC) # equiv: orth(Z2, cbind(V, SC))
                CvcZ2 <- Corth(OvZ2, SC) # For prediction
                plsRvYOvcZ2 <- pls.fit(OvcZ2, res, nZ2u)
                SZ2u <- plsRvYOvcZ2$scores
                Vadd[,nZ1u + 1:nZ2u] <- S[[i]][[2]] <- SZ2u
                L[[i]][[2]] <- plsRvYOvcZ2$loadings
                ## The number of components added to V in this step:
                added <- nZ1u + nZ2u + nC
            } else {
                ## Don't extract common/unique components
                Vadd <- matrix(nrow = nObs, ncol = sum(ncomp[[i]])) # The variables to be added in the present step
                added <- 0              # The number of comps. currently added
                for (j in seq(along = M)) {
                    ##cat("j =", j, "\n")
                    ## Walk through Z[[i]]
                    Mo <- orth(M[[j]], V[,1:nVar])
                    orthCoefs[[i]][[j]] <- Corth(M[[j]], V[,1:nVar]) # For pred
                    models[[i]][[j]] <- pls.fit(Mo, res, ncomp[[i]][j])# Could use Y
                    Vadd[,added + (1:ncomp[[i]][[j]])] <- S[[i]][[j]] <- models[[i]][[j]]$scores
                    L[[i]][[j]] <- models[[i]][[j]]$loadings
                    added <- added + ncomp[[i]][j]
                }
            } # if
            V[,nVar + 1:added] <- Vadd
            ## Only strictly neccessary for unique/common:
            lmZ <- lm.fit(Vadd, res)
            B[nVar + 1:added,] <- lmZ$coefficients
            ##Balt[nVar + 1:added,] <- lmZ$coefficients
            res <- lmZ$residuals
            nVar <- nVar + added
        } # if
    } # for
    list(coefficients = B, predictors = V, orthCoefs = orthCoefs,
         models = models, ncomp = ncomp, scores = S, loadings = L, residuals = res)
} # function

## Forenklingstanke: Gjoer om alle enkeltmatrisene i Z til lister med
## ett element.  Da blir algoritmene enklere (og burde ikke bli
## nevneverdig saktere).  Ved behov kan man teste paa length(M)
## (f.eks. ved beregning av B og ny res).

## FIXME:
## - Change name of 'models' component to 'plsmodels'?
## - Add fitted values?
