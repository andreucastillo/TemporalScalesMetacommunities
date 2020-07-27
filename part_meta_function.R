#-----------------------------------------------------------------------------------------------
# Script: Prepare data for basic metacommunity analysis (sampling data)
# Project: Temporal scale on invertebrates
# Date: 24/06/2020
# Author: Castillo-Escrivà, A., Mesquita-Joanes, F. and Rueda, J.
# e-mail: acastilloescriva@gmail.com
#------------------------------------------------------------------------------------------------

# Sampling simulations
sim.pre <- function(sim.dir, time.extent = 55, n.time = 10, n.site = 10, rm.start = TRUE){
  land.file <- list.files(path = sim.dir, pattern = "land", full.names = TRUE)
  land <- readRDS(land.file[1]) # load landscape
  
  # Sampling time
  sample <- land$land
  time.sel <- unique(sample$time)
  if(isTRUE(rm.start)){
    time.sel <- time.sel[-c(1:(land$start+land$burnin))]
  }
  time.sel <- round(seq(time.sel[1], time.sel[1] + time.extent, length.out = n.time))
  sample <- sample[sample$time %in% time.sel, ]
  if(isTRUE(rm.start)){
    sample$time <- sample$time - land$start - land$burnin
    time.sel <- time.sel - land$start - land$burnin
  }
  
  # sampling space
  site.sel <- sample(1:length(unique(sample$site)), n.site, replace = FALSE)
  site.sel <- sort(site.sel)
  sample <- sample[sample$site %in% unique(sample$site)[site.sel],]
  
  # sp matrix
  sp.files <- list.files(path = sim.dir, pattern = "sp", full.names = TRUE)
  sp.sel <- sp.files[time.sel]
  sp.in <- read.table(sp.sel[1], sep = ",", header = TRUE)
  sp.in <- sp.in[site.sel,]
  sp.out <- data.frame(matrix(0, nrow(sample), ncol = ncol(sp.in)))
  sp.out[sample$time == time.sel[1],] <- sp.in
  
  for(i in 2:length(time.sel)){
    sp.in <- read.table(sp.sel[i], sep = ",", header = TRUE)
    sp.in <- sp.in[site.sel, ,drop = FALSE]
    sp.out[sample$time == time.sel[i],] <- sp.in
  }
  colnames(sp.out) <- sprintf("sp%02d", 1:ncol(sp.out))
  rownames(sp.out) <- rownames(sample)
  
  env <- sample[,5, drop = FALSE]
  
  return(list("sample" = sample, "sp" = sp.out, "env" = env))
}

#----------------------------------------------------------------------------------------------
# subsampling empirical data
emp.pre <- function(sample, sp, env, trait, land, group, start, last, lag){
  
  # select land
  sample <- sample[sample$land == land, ]
  sp <- sp[rownames(sp) %in% rownames(sample), ]
  env <- env[rownames(env) %in% rownames(sample), ]
  env <- env[,!is.na(colSums(env))]
  
  # select species
  sp <- sp[ ,colnames(sp) %in% rownames(trait[trait$group == group, ])]
  sp[is.na(sp)] <- 0 # replace NA into 0
  sp <- sp[,colSums(sp) != 0]
  sp <- sp[rowSums(sp) != 0, ] # remove samples without records
  sample <- sample[rownames(sample) %in% rownames(sp),]
  env <- env[rownames(env) %in% rownames(sp),]
  
  # subsampling
  time.sel <- as.vector(unique(sample$time)[seq(start, last, lag)])
  sample <- sample[sample$time %in% time.sel, ]
  sp <- sp[rownames(sp) %in% rownames(sample),]
  env <- env[rownames(env) %in% rownames(sample),]
  
  # list with the data
  res <- list("sample" = sample, "sp" = sp, "env" = env)
  return(res)
}

#-----------------------------------------------------------------------------------------------
# forward selection

fw <- function(Y, X, mem.t = TRUE, space, time, method = "rda", perm.cca = 499, poly2 = FALSE){
  require(vegan)
  require(Matrix)
  if(method == "rda"){
    Y <- decostand(Y, method = "hellinger")
  }
  if(method == "capscale"){
    Y <- vegdist(Y, method = "bray")
  }
  
  name <- names(X)
  if(isTRUE(poly2)){
    ei <- 1:ncol(X)
    X.poly <- as.matrix(sparse.model.matrix(as.formula(paste("~poly(X[,",ei,"],2)",
                                                             collapse = " + "))))[,-1]
    X.poly <- as.data.frame(X.poly)
    names(X.poly) <- paste(rep(name, each = 2), rep(1:2, length(name)), sep= ".")
  } else {
    X.poly <- X
  }
  
  mod0 <- do.call(method, list(Y~1, X.poly))
  mod1 <- do.call(method, list(Y~., X.poly))
  fw <- ordiR2step(mod0, mod1, trace = FALSE, R2permutations = perm.cca, permutations = 2999) 
  X.sel <- X.poly[,colnames(X.poly) %in% rownames(scores(fw, display = "bp")), drop = FALSE]
  fw$anova$`Pr(>F)` <- p.adjust (fw$anova$`Pr(>F)`, method = "bonferroni", n = ncol(X.poly))
  adj <- fw$anova$`Pr(>F)` < 0.05
  X.sel <- X.sel[,adj[-length(adj)], drop = FALSE]
  if(isTRUE(poly2)){
    colnames(X.sel) <- sapply(strsplit(colnames(X.sel),"[.]"), `[`, 1)
    X.sel <- X[, colnames(X) %in% colnames(X.sel), drop = FALSE]
  }
  
  if(isTRUE(mem.t)){
    if(ncol(X.sel) > 0){
    p.value <- test.st(mem = X.sel, space = space, time = time)
    mem.res <- X.sel[,p.value$space > 0.05 && p.value$time > 0.05, drop = FALSE]
    mem.spa <- X.sel[,p.value$space < 0.05, drop = FALSE]
    mem.tim <- X.sel[,p.value$time < 0.05, drop = FALSE]
    
    mem.spa <- data.frame(mem.spa, mem.res)
    mem.tim <- data.frame(mem.tim, mem.res)
    if(ncol(mem.spa) > 0){
      colnames(mem.spa) <- paste("spa", 1:ncol(mem.spa), sep = "")
    }
    if(ncol(mem.tim) > 0){
      colnames(mem.tim) <- paste("tim", 1:ncol(mem.tim), sep = "")
    }
    X.sel <- list("mem.spa" = mem.spa, "mem.tim" = mem.tim)
    } else {
      X.sel <- list("mem.spa" = X.sel, "mem.tim" = X.sel)
    }
  }
  return(X.sel)
}

#-----------------------------------------------------------------------------------------------
# Complete metacommunity test routine - based on variation partitioning

test.meta <- function(sp, sample, env, method = "rda", nrepet.msr = 999, perm.cca = 999,
                      poly2 = TRUE, env.pca = TRUE){
  require(vegan)
  require(adespatial)
  require(spdep)
  
  sp <- sp[,colSums(sp)!=0, drop = FALSE]
  sp <- sp[rowSums(sp)!=0,, drop = FALSE]
  sample <- sample[rownames(sample) %in% rownames(sp), ]
  env <- env[rownames(env) %in% rownames(sp), , drop  = FALSE]
  
  if(isTRUE(env.pca)){
    pca <- rda(env, scale = TRUE)
    env <- as.data.frame(scores(pca, choices = 1:length(eigenvals(pca)), display = "sites"))
  }
  
  links <- xyz2links(names = sample[,4:3], coords = sample[,1:3])
  mem <- as.data.frame(scores.listw(links$listw, MEM.autocor = "positive"))
  mtest <- moran.randtest(mem, links$listw, nrepet = 99)
  mem <- mem[ ,mtest$pvalue < 0.05] # Select mems with significant positive eigenvalues
  
  # Forward selection
  env.sel <- fw(Y = sp, X = env, mem.t = FALSE, method = method, perm.cca = perm.cca, 
                poly2 = poly2)
  mem.sel <- fw(Y = sp, X = mem, mem.t = TRUE, space = sample$site, time = sample$time,
                method = method, perm.cca = perm.cca, poly2 = FALSE)
  
  n.var <- rep(0,4)
  n.var[1] <- ncol(env.sel)
  mem.cor <- round(cor(mem.sel$mem.spa, mem.sel$mem.tim))
  n.var[2] <- ncol(mem.sel$mem.spa[,rowSums(mem.cor) == 0, drop = FALSE])
  n.var[3] <- ncol(mem.sel$mem.tim[,colSums(mem.cor) == 0, drop = FALSE])
  n.var[4] <- ncol(mem.sel$mem.spa[,rowSums(mem.cor) != 0, drop = FALSE])
  names(n.var) <- c("env", "spa", "tim", "spa-tim")
  
  # Variation partitioning
  if(ncol(env.sel) > 0){
    if(ncol(mem.sel$mem.spa) > 0 && ncol(mem.sel$mem.tim) > 0){
      vp <- vpmsr(sp, env.sel, mem.sel$mem.spa, mem.sel$mem.tim, perm.cca = perm.cca,
                  listw = links$listw, method = method, nrepet.msr = nrepet.msr, poly2 = poly2)
      vp <- vp[2,]
    }
    if(ncol(mem.sel$mem.spa) == 0 && ncol(mem.sel$mem.tim) == 0){
      vp <- vpmsr(sp, env.sel, listw = links$listw, method = method, nrepet.msr = nrepet.msr,
                  poly2 = poly2, perm.cca = perm.cca)
      vp <- c(vp[2,], rep(0, 6), 1-vp[2,])
    }
    if(ncol(mem.sel$mem.spa) > 0 && ncol(mem.sel$mem.tim) == 0){
      vp <- vpmsr(sp, env.sel, mem.sel$mem.spa, listw = links$listw, method = method,
                  nrepet.msr = nrepet.msr, poly2 = poly2, perm.cca = perm.cca)
      vp <- c(vp[2,1], vp[2,3], 0, vp[2,2], rep(0, 3), vp[2,4])
    }
    if(ncol(mem.sel$mem.spa) == 0 && ncol(mem.sel$mem.tim) > 0){
      vp <- vpmsr(sp, env.sel, mem.sel$mem.tim, listw = links$listw, method = method,
                  nrepet.msr = nrepet.msr, poly2 = poly2, perm.cca = perm.cca)
      vp <- c(vp[2,1], 0, vp[2,3], rep(0, 2), vp[2,2], 0, vp[2,4])
    }  
  } else{
    if(ncol(mem.sel$mem.spa) > 0 && ncol(mem.sel$mem.tim) > 0){
      if(method == "rda" || method == "capscale"){
        if(method == "rda"){
          sp <- decostand(sp, method = "hellinger")
        }
        if(method == "capscale"){
          sp <- vegdist(sp, method = "bray")
        }
        vp <- suppressWarnings(varpart(sp, mem.sel$mem.spa, mem.sel$mem.tim))
      }
      if(method == "cca"){
        vp <- suppressWarnings(varpart(sp, mem.sel$mem.spa, mem.sel$mem.tim, chisquare = TRUE))
      }
      vp <- vp$part$indfract[,3]
      vp <- c(0, vp[1], vp[3], 0, vp[2], rep(0,2), vp[4])
    }
    if(ncol(mem.sel$mem.spa) == 0 && ncol(mem.sel$mem.tim) == 0){
      vp <- c(rep(0,7), 1)
    }
    if(ncol(mem.sel$mem.spa) > 0 && ncol(mem.sel$mem.tim) == 0){
      if(method == "rda"){
        sp <- decostand(sp, method = "hellinger")
      }
      if(method == "capscale"){
        sp <- vegdist(sp, method = "bray")
      }
      vp <- RsquareAdj(do.call(method, list(sp~., mem.sel$mem.spa)),
                       permutations = perm.cca)$adj.r.squared
      vp <- c(0, vp, rep(0,5), 1-vp)
    }
    if(ncol(mem.sel$mem.spa) == 0 && ncol(mem.sel$mem.tim) > 0){
      if(method == "rda"){
        sp <- decostand(sp, method = "hellinger")
      }
      if(method == "capscale"){
        sp <- vegdist(sp, method = "bray")
      }
      vp <- RsquareAdj(do.call(method, list(sp~., mem.sel$mem.tim)),
                       permutations = perm.cca)$adj.r.squared
      vp <- c(rep(0,2), vp, rep(0,4), 1-vp)
    }  
  }
  names(vp) <- letters[1:8]
  return(list("vp" = vp, "n.var" = n.var))
}

