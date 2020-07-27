#-----------------------------------------------------------------------------------------------
# Script: Functions for variation partitioning with a correction by means of MSR
# Project: Temporal scale on invertebrates
# Date: 24/06/2020
# Author: Castillo-Escrivà, A., Mesquita-Joanes, F. and Rueda, J.
# e-mail: acastilloescriva@gmail.com
#------------------------------------------------------------------------------------------------

# This code is based on Clappe et al., 2018. https://doi.org/10.1002/ecy.2376

# Variation partitioning with a correction by means of MSR for one (environment),
# two or three expanatory sets
vpmsr <- function(Y, X, ..., listw, method = "rda", nrepet.msr = 999, poly2 = TRUE,
                  perm.cca = 999){
  require(vegan)
  require(adespatial)
  require(Matrix)
  
  if(method == "rda"){
    Y <- decostand(Y, method = "hellinger")
  }
  if(method == "capscale"){
    Y <- vegdist(Y, method = "bray")
  }
  
  Xmsr <- msr(x = X, listwORorthobasis = listw, nrepet = nrepet.msr, simplify  = FALSE)
  
  if(isTRUE(poly2)){
    name <- names(X)
    ei <- 1:ncol(X)
    X <- as.matrix(sparse.model.matrix(as.formula(paste("~poly(X[,",ei,"],2)",
                                                        collapse = " + "))))[,-1]
    X <- as.data.frame(X)
    names(X) <- paste(rep(name, each = 2), rep(1:2, length(name)), sep= ".")
    
    Xmsr <- lapply(Xmsr, function(x) as.data.frame(as.matrix(
      sparse.model.matrix(as.formula(paste("~poly(x[,",ei,"],2)", collapse = " + "))))[,-1]))
    Xmsr <- lapply(Xmsr, setNames, names(X))
  }
  
  exp <- list(X, ...)
  if(length(exp) < 1 || length(exp) > 4){
    stop("needs 1 to three explanatory tables")
  }
  
  if(length(exp) == 1){
    vp <- vpmsr1(Y = Y, X = exp[[1]], Xmsr = Xmsr, listw = listw, method = method,
                 perm.cca = perm.cca)
  }
  if(length(exp) == 2){
    vp <- vpmsr2(Y = Y, X = exp[[1]], W = exp[[2]], Xmsr = Xmsr, listw = listw, method = method,
                 perm.cca = perm.cca)
  }
  if(length(exp) == 3){
    vp <- vpmsr3(Y = Y, X = exp[[1]], W = exp[[2]], Z = exp[[3]], Xmsr = Xmsr, listw = listw,
                 method = method, perm.cca = perm.cca)
  }
  return(vp)
}

#------------------------------------------------------------------------------------------------
# One explanatory set
vpmsr1 <- function(Y, X, Xmsr, listw, method, perm.cca){
  
  ## Standard Rsquared adjustment
  R2X <- RsquareAdj(do.call(method, list(Y~., X)), permutations = perm.cca)
  R2a.sta <- R2X$adj.r.squared
  names(R2a.sta) <- "a"
  
  ## R squared based on MSR
  # MSR of the environmental matrix
  
  R2dist <- function(Y, Xmsr, method){
    R2dist <- rep(NA, length(Xmsr))
    for(i in 1:length(R2dist)){
        R2dist[i] <- RsquareAdj(do.call(method, list(Y~., as.data.frame(Xmsr[[i]]))),
                                permutations = 0)$r.squared
        }
    return(R2dist)
    }
  
  # X whole part (e.g., environmental)
  R2X.msr <- R2dist(Y, Xmsr, method = method)
  R2aX.msr <- 1 - (1-R2X$r.squared) / (1-mean(R2X.msr))
  
  ## Standard and MSR-based table
  vp <- t(data.frame("standard" = R2a.sta, "msr" = R2aX.msr))
  return(vp)
}

#------------------------------------------------------------------------------------------------
# Two explanatory sets

vpmsr2 <- function(Y, X, W, Xmsr, listw, method, perm.cca){
  
  ## Standard variation partitioning
  if(method == "rda" || method == "capscale"){
    vp.sta <- varpart(Y, X, W)
  }
  if(method == "cca"){
    vp.sta <- varpart(Y, X, W, chisquare = TRUE, permutations = perm.cca)
  }
  
  ## Variation partitioning based on MSR
  # MSR of the environmental matrix
  
  R2dist <- function(Y, Xmsr, method){
    R2dist <- rep(NA, length(Xmsr))
    for(i in 1:length(R2dist)){
      R2dist[i] <- RsquareAdj(do.call(method, list(Y~., as.data.frame(Xmsr[[i]]))),
                                permutations = 0)$r.squared
      }
    return(R2dist)
  }
    
  
  # X whole part (e.g., environmental)
  R2X1 <- vp.sta$part$fract[1,2]
  R2X2 <- vp.sta$part$fract[2,2]
  R2aX2 <- vp.sta$part$fract[2,3]
  R2X1X2 <- vp.sta$part$fract[3,2]
  
  R2X1.msr <- R2dist(Y, Xmsr, method = method)
  R2aX1 <- 1 - (1-R2X1) / (1-mean(R2X1.msr))
  
  # X and W
  WXmsr <- lapply(Xmsr, cbind, W)
  R2X1X2.msr <- R2dist(Y, WXmsr, method = method)
  
  # Individual fractions
  a <- 1 - (1-(R2X1X2 - R2X2)) / (1-mean(R2X1X2.msr - R2X2))
  b <- R2aX1 -  a
  c <- R2aX2 - b
  d <- 1 - (a + b + c)
  
  vp.sta <- vp.sta$part$indfract[,3]
  vp.msr <- c(a, b, c, d)
  vp <- t(data.frame("standard" = vp.sta, "msr" = vp.msr))
  colnames(vp) <- letters[1:ncol(vp)]
  
  ## Standard and MSR-based table
  return(vp)
}


#------------------------------------------------------------------------------------------------
# Three explanatory sets

vpmsr3 <- function(Y, X, W, Z, Xmsr, listw, method, perm.cca){
  
  ## Standard variation partitioning
  if(method == "rda" || method == "capscale"){
    vp.sta <- suppressWarnings(varpart(Y, X, W, Z))
  }
  if(method == "cca"){
    vp.sta <- suppressWarnings(varpart(Y, X, W, Z, chisquare = TRUE, permutations = perm.cca))
  }
  
  ## Variation partitioning based on MSR
  # MSR of the environmental matrix
  R2dist <- function(Y, Xmsr, method){
    R2dist <- rep(NA, length(Xmsr))
    for(i in 1:length(R2dist)){
      R2dist[i] <- RsquareAdj(do.call(method, list(Y~., as.data.frame(Xmsr[[i]]))),
                                permutations = 0)$r.squared
      }
    return(R2dist)
    }
    
  # Partial partioning XW
  R2X1 <- vp.sta$part$fract[1,2]
  R2X2 <- vp.sta$part$fract[2,2]
  R2aX2 <- vp.sta$part$fract[2,3]
  R2X1X2 <- vp.sta$part$fract[4,2]
  
  R2X1.msr <- R2dist(Y, Xmsr, method = method)
  R2aX1 <- 1 - (1-R2X1) / (1-mean(R2X1.msr))
  
  WXmsr <- lapply(Xmsr, cbind, W)
  R2X1X2.msr <- R2dist(Y, WXmsr, method = method)
  
  a.W <- 1 - (1-(R2X1X2 - R2X2)) / (1-mean(R2X1X2.msr - R2X2))
  b.W <- R2aX1 -  a.W
  c.W <- R2aX2 - b.W
  R2aX1X2 <- a.W + b.W + c.W
  
  # Partial partioning XZ
  R2X3 <- vp.sta$part$fract[3,2]
  R2aX3 <- vp.sta$part$fract[3,3]
  R2X1X3 <- vp.sta$part$fract[5,2]
  
  ZXmsr <- lapply(Xmsr, cbind, Z)
  R2X1X3.msr <- R2dist(Y, ZXmsr, method = method)
  
  a.Z <- 1 - (1-(R2X1X3 - R2X3)) / (1-mean(R2X1X3.msr - R2X3))
  b.Z <- R2aX1 -  a.Z
  c.Z <- R2aX3 - b.Z
  R2aX1X3 <- a.Z + b.Z + c.Z
  
  # all the model
  R2X2X3 <- vp.sta$part$fract[6,2]
  R2aX2X3 <- vp.sta$part$fract[6,3]
  R2X1X2X3 <- vp.sta$part$fract[7,2]
  
  WZXmsr <- lapply(Xmsr, cbind, cbind(W,Z))
  R2X1X2X3.msr <- R2dist(Y, WZXmsr, method = method)
  
  a <- 1 - (1 - (R2X1X2X3 - R2X2X3)) / (1-mean(R2X1X2X3.msr - R2X2X3))
  all <- a + R2aX2X3
  b <- all - R2aX1X3
  c <- all - R2aX1X2 
  d <- all - R2aX3 - a - b
  e <- all - R2aX1 - b - c
  f <- all - R2X2 - a - c
  g <- all - a - b - c - d - e - f  
  h <- 1 - all  
  
  vp.sta <- vp.sta$part$indfract[,3]
  vp.msr <- c(a, b, c, d, e, f, g, h)
  vp <- t(data.frame("standard" = vp.sta, "msr" = vp.msr))
  colnames(vp) <- letters[1:ncol(vp)]
  
  ## Standard and MSR-based table
  return(vp)
}

