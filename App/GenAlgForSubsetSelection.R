GenAlgForSubsetSelection
function (P, Candidates, Test, ntoselect, npop = 100, nelite = 5, 
    keepbest = TRUE, tabu = TRUE, tabumemsize = 1, mutprob = 0.8, 
    mutintensity = 1, niterations = 500, minitbefstop = 200, 
    niterreg = 5, lambda = 1e-06, plotiters = FALSE, plottype = 1, 
    errorstat = "PEVMEAN2", C = NULL, mc.cores = 1, InitPop = NULL, 
    tolconv = 1e-07, Vg = NULL, Ve = NULL, Fedorov = FALSE) 
{
    if ((ncol(P) + 1) > ntoselect) {
        warning("The algorithm does not work well with p>ntrain, perhaps use a unsupervised dimension reduction on P.")
    }
    ridgereg <- function(y = CurrentPopFuncValues, x = designforCurrentPop, 
        lambda = lambda) {
        n <- nrow(x)
        p <- ncol(x)
        mindim <- min(p, n)
        rownames(x) <- NULL
        svdX <- svd(x, nu = mindim, nv = mindim)
        insvdnonzero <- 1:mindim
        diagvecforinv <- (svdX$d[insvdnonzero])/((svdX$d[insvdnonzero])^2 + 
            lambda)
        coef <- tcrossprod(svdX$v %*% diag(diagvecforinv), svdX$v) %*% 
            t(x) %*% (y - mean(y))
        return(list(coef = coef))
    }
    getsetsfedorov <- function(dummyx) {
        randint <- sample(1:10, 1)
        npc <- min(c(ntoselect - 10, ncol(P)))
        if (randint == 1) {
            optsol <- AlgDesign::optFederov(data = P[rownames(P) %in% 
                Candidates, 1:npc], nTrials = ntoselect, criterion = "D")
            outf <- Candidates[optsol$rows]
        }
        else if (randint == 2) {
            optsol <- AlgDesign::optFederov(data = P[rownames(P) %in% 
                Candidates, 1:npc], nTrials = ntoselect, criterion = "A")
            outf <- Candidates[optsol$rows]
        }
        else if (randint == 3) {
            optsol <- AlgDesign::optFederov(data = P[rownames(P) %in% 
                Candidates, 1:npc], nTrials = ntoselect, criterion = "I")
            outf <- Candidates[optsol$rows]
        }
        else if (randint == 4) {
            optsol <- AlgDesign::optFederov(data = as.data.frame(P[rownames(P) %in% 
                Candidates, 1:npc]), nTrials = ntoselect, criterion = "I", 
                space = as.data.frame(P[rownames(P) %in% Test, 
                  1:npc]), nRepeats = 1)
            outf <- Candidates[optsol$rows]
        }
        else {
            outf <- sample(Candidates, ntoselect)
        }
        return(outf)
    }
    if (is.null(InitPop)) {
        if (sum(errorstat %in% c("AOPT", "DOPT", 
            "CDMEAN", "PEVMEAN", "PEVMEAN2")) > 
            0) {
            if (Fedorov) {
                InitPop <- lapply(1:npop, function(q) {
                  return(getsetsfedorov(q))
                })
            }
            else {
                InitPop <- lapply(1:npop, function(q) {
                  return(sample(Candidates, ntoselect))
                })
            }
        }
    }
    mc.cores <- switch(Sys.info()[["sysname"]], Windows = {
        1
    }, Linux = {
        mc.cores
    }, Darwin = {
        mc.cores
    })
    nelite = nelite - 1
    if (tabu) {
        tabumemsize = tabumemsize
    }
    else {
        tabumemsize = 0
    }
    if (tabumemsize > 0) {
        memoryfortabu <- vector(mode = "list", length = tabumemsize)
    }
    else {
        memoryfortabu = NULL
    }
    completeton <- function(x) {
        if (length(x) < ntoselect) {
            x <- c(x, sample(setdiff(unique(unlist(InitPop)), 
                x), ntoselect - length(x)))
        }
        return(x)
    }
    if (is.null(InitPop)) {
        InitPop <- lapply(1:npop, function(x) {
            return(sample(Candidates, ntoselect))
        })
    }
    else {
        InitPop <- lapply(InitPop, completeton)
    }
    if (errorstat %in% c("CDMEANMM", "PEVMEANMM", 
        "GAUSSMEANMM", "userMM")) {
        K = solve(P)
    }
    linenames <- rownames(P)
    if (!is.null(errorstat)) {
        if (errorstat %in% c("CDMEANMM", "PEVMEANMM", 
            "GAUSSMEANMM", "userMM")) {
            InitPopFuncValues <- as.numeric(unlist(mclapply(InitPop, 
                FUN = function(x) {
                  do.call(errorstat, list(x, Test, P, K, lambda, 
                    C, Vg, Ve))
                }, mc.cores = mc.cores, mc.preschedule = T)))
        }
        else {
            InitPopFuncValues <- as.numeric(unlist(mclapply(InitPop, 
                FUN = function(x) {
                  do.call(errorstat, list(x, Test, P, lambda, 
                    C))
                }, mc.cores = mc.cores, mc.preschedule = T)))
        }
    }
    designforInitPop <- Reduce("rbind", lapply(InitPop, 
        function(x) {
            as.numeric(Candidates %in% x)
        }))
    if (length(InitPop) > 2) {
        sdofdesigns <- mean(apply(designforInitPop, 2, sd))
    }
    else {
        sdofdesigns = 0
    }
    betahat = ridgereg(y = InitPopFuncValues, x = designforInitPop, 
        lambda = lambda)$coef
    orderofInitPop <- order(InitPopFuncValues, decreasing = FALSE)
    ElitePop <- mclapply(orderofInitPop[1:nelite], FUN = function(x) {
        return(InitPop[[x]])
    }, mc.cores = mc.cores, mc.preschedule = T)
    ElitePopFuncValues <- InitPopFuncValues[orderofInitPop[1:nelite]]
    minvec <- c()
    sdvec <- c()
    designforElitePop <- Reduce("rbind", lapply(ElitePop, 
        function(x) {
            as.numeric(linenames %in% x)
        }))
    if (length(ElitePop) > 2) {
        sdofelitedesigns <- mean(apply(designforElitePop, 2, 
            sd))
    }
    else {
        sdofelitedesigns = 0
    }
    sdvec <- c(sdvec, sdofelitedesigns)
    minvec <- c(minvec, min(ElitePopFuncValues))
    for (iters in 1:niterations) {
        if (iters > minitbefstop) {
            maxmeans <- max(minvec[(length(minvec) - minitbefstop):length(minvec)])
            minmeans <- min(minvec[(length(minvec) - minitbefstop):length(minvec)])
            meansdiff <- maxmeans - minmeans
            if (meansdiff < tolconv) {
                (break)()
            }
        }
        if (iters > tabumemsize) {
            CurrentPop <- GenerateCrossesfromElites(Elites = ElitePop, 
                Candidates = Candidates, npop = npop, mutprob = mutprob, 
                mc.cores = mc.cores, mutintensity = mutintensity, 
                memoryfortabu = memoryfortabu)
        }
        else {
            CurrentPop <- GenerateCrossesfromElites(Elites = ElitePop, 
                Candidates = Candidates, npop = npop, mutprob = mutprob, 
                mc.cores = mc.cores, mutintensity = mutintensity, 
                memoryfortabu = NULL)
        }
        if (tabumemsize > 0) {
            memoryfortabu[[iters%%tabumemsize + 1]] <- CurrentPop
        }
        else {
            memoryfortabu = NULL
        }
        if (keepbest) {
            CurrentPop <- c(CurrentPop, ElitePop[1])
        }
        if (!is.null(errorstat)) {
            if (errorstat %in% c("CDMEANMM", "PEVMEANMM", 
                "GAUSSMEANMM", "userMM")) {
                CurrentPopFuncValues <- as.numeric(unlist(mclapply(CurrentPop, 
                  FUN = function(x) {
                    do.call(errorstat, list(x, Test, P, K, lambda, 
                      C, Vg, Ve))
                  }, mc.cores = mc.cores, mc.preschedule = T)))
            }
            else {
                CurrentPopFuncValues <- as.numeric(unlist(mclapply(CurrentPop, 
                  FUN = function(x) {
                    do.call(errorstat, list(x, Test, P, lambda, 
                      C))
                  }, mc.cores = mc.cores, mc.preschedule = T)))
            }
        }
        if (niterreg > iters) {
            designforCurrentPop <- Reduce("rbind", lapply(CurrentPop, 
                function(x) {
                  as.numeric(Candidates %in% x)
                }))
            if (length(CurrentPop) > 2) {
                sdofdesigns <- mean(apply(designforCurrentPop, 
                  2, sd))
            }
            else {
                sdofdesigns = 0
            }
            betahat = ridgereg(y = CurrentPopFuncValues, x = designforCurrentPop, 
                lambda = lambda)$coef
            orderofBetasforCurrentPop <- order(betahat, decreasing = FALSE)
            bestpredicted <- Candidates[orderofBetasforCurrentPop[1:ntoselect]]
            orderofCurrentPop <- order(CurrentPopFuncValues, 
                decreasing = FALSE)
            ElitePop <- lapply(orderofCurrentPop[1:nelite], FUN = function(x) {
                return(CurrentPop[[x]])
            })
            ElitePopFuncValues <- CurrentPopFuncValues[orderofCurrentPop[1:nelite]]
            designforElitePop <- Reduce("rbind", lapply(ElitePop, 
                function(x) {
                  as.numeric(linenames %in% x)
                }))
            if (length(ElitePop) > 2) {
                sdofelitedesigns <- mean(apply(designforElitePop, 
                  2, sd))
            }
            else {
                sdofelitedesigns = 0
            }
            sdvec <- c(sdvec, sdofelitedesigns)
            ElitePop[[nelite + 1]] <- bestpredicted
            ElitePopFuncValues[nelite + 1] <- mean(ElitePopFuncValues)
            minvec <- c(minvec, min(ElitePopFuncValues))
        }
        else {
            orderofCurrentPop <- order(CurrentPopFuncValues, 
                decreasing = FALSE)
            ElitePop <- lapply(orderofCurrentPop[1:nelite], FUN = function(x) {
                return(CurrentPop[[x]])
            })
            designforElitePop <- Reduce("rbind", lapply(ElitePop, 
                function(x) {
                  as.numeric(linenames %in% x)
                }))
            if (length(ElitePop) > 2) {
                sdofelitedesigns <- mean(apply(designforElitePop, 
                  2, sd))
            }
            else {
                sdofelitedesigns = 0
            }
            sdvec <- c(sdvec, sdofelitedesigns)
            ElitePopFuncValues <- CurrentPopFuncValues[orderofCurrentPop[1:(nelite + 
                1)]]
            minvec <- c(minvec, min(ElitePopFuncValues))
        }
        if (plotiters) {
            if (plottype == 1) {
                par(mfrow = c(1, 1))
                plot(minvec, type = "b")
                abline(h = min(minvec), col = "red")
                mtext(side = 3, line = -3, col = "blue", 
                  text = paste("\n     Minimum value of objective function is ", 
                    round(min(ElitePopFuncValues), 5), ".", 
                    sep = ""), adj = 0, outer = T)
            }
            else if (plottype == 2) {
                par(mfrow = c(2, 1))
                plot(minvec, type = "b")
                abline(h = min(minvec), col = "red")
                mtext(side = 3, line = -3, col = "blue", 
                  text = paste("\n     Minimum value of objective function is ", 
                    round(min(ElitePopFuncValues), 5), ".", 
                    sep = ""), adj = 0, outer = T)
                barplot(round(apply(designforElitePop, 2, mean), 
                  3))
                par(mfrow = c(1, 1))
            }
            else if (plottype == 3) {
                par(mfrow = c(3, 1))
                plot(minvec, type = "b")
                abline(h = min(minvec), col = "red")
                mtext(side = 3, line = -3, col = "blue", 
                  text = paste("\n    Divergence: ", round(sdofelitedesigns, 
                    5), ".\n     Minimum value of objective function is ", 
                    round(min(ElitePopFuncValues), 5), ".", 
                    sep = ""), adj = 0, outer = T)
                barplot(round(apply(designforElitePop, 2, mean), 
                  3))
                plot(sdvec, type = "b", col = (sdvec == 
                  0) + 1)
                par(mfrow = c(1, 1))
            }
            else {
            }
        }
    }
    ElitePop[[nelite + 2]] <- minvec
    for (i in 1:(nelite + 1)) {
        names(ElitePop)[i] <- paste("Solution with rank ", 
            i, sep = "")
    }
    names(ElitePop)[nelite + 2] <- "Best criterion values over iterarions"
    return(ElitePop)
}