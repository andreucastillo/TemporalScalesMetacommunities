#-----------------------------------------------------------------------------------------------
# Script: Functions for metacommunity simulations
# Project: Temporal scale on invertebrates
# Date: 24/06/2020
# Author: Castillo-Escrivà, A., Mesquita-Joanes, F. and Rueda, J.
# e-mail: acastilloescriva@gmail.com
#------------------------------------------------------------------------------------------------

# make landscape
make.land <- function(n.site = 10, space.side = 100, time.extent = 2000, start = 200,
                      burnin = 800, env.scale = 500){
  # env.type could be neu = neutral, spa = inly spatial change, tim = only temporal change,
  # st = spatiotemporal change.
  
  require(som.nn)
  require(RandomFields)
  require(vegan)
  t.total <- time.extent + start + burnin + 1
  
  # space (xy)
  coord <- replicate(2, runif(n.site, 0, space.side))
  
  # distances
  dist <- as.matrix(dist.torus(coord))
  
  # Environmental variable
  model <- RMexp(var = 0.05, scale = env.scale)
  
  repeat{
    sim <- RFsimulate(model = model, x = coord[,1]*10, y = coord[,2]*10, T = (start+1):t.total)
    env <- decostand(sim$variable1, "range")
    if((max(env[1:(n.site*100)]) - min(env[1:(n.site*100)])) > 0.6) {break}
  }
  
  ecum <- ecdf(env)
  env.cum <- ecum(env)
  
  # landscape
  land <- data.frame("x" = rep(coord[,1], t.total), "y" = rep(coord[,2], t.total))
  land$time <- rep(1:t.total, each = n.site)
  land$site <- rep(sprintf("s%02d", 1:n.site), t.total)
  land$env <- c(rep(env.cum[1:n.site], start), env.cum)
  
  plot(land$env[land$site == unique(land$site)[1]], type = "l", ylim = c(0,1), xlab = "time",
       ylab = "env")
  for(i in 2:n.site){
    lines(land$env[land$site == unique(land$site)[i]], col = i)
  }
  
  return(list("space.side" = space.side, "time.extent" = time.extent,
              "land" = land, "n.si" = n.site, "n.time" = t.total, "dist" = dist,
              "start" = start, "burnin" = burnin))
}

#------------------------------------------------------------------------------------
# make species regional pool
make.rp <- function(meta.name, sr = 10, fec = 10, intra = 1, inter = 0.95, sigma = 10,
                    a = 0.02, w = 10^-5, m = 0.0001){
  mu <- runif(sr, 0, 1)
  niche <- data.frame(mu, sigma)
  rp <- list("name" = as.character(meta.name), "sr" = sr, "fec" = fec, "intra" = intra*0.05,
             "inter" = inter * 0.05, "niche" = niche, "a" = a, "w" = w, "m" = m) 
  return(rp)
}

#------------------------------------------------------------------------------------------------
# Niche function
lambda <- function(env, mu = 0.5, sig = 10){
  lam <- exp((-(env - mu) ^ 2) / (2 * sig ^ 2))
  return(round(lam, 2))
}

#-------------------------------------------------------------------------------------
# Make sp matrices to preallocate
make.meta <- function(n.si, sr, n.time){
  sp <- matrix(0, n.si, sr)
  disp <- matrix(0, n.si, n.si)
  l <- matrix(0, n.si * n.time, sr)
  meta <- list("sp" = sp, "disp" = disp, "l" = l, "em" = sp)
  return(meta)
}

#-------------------------------------------------------------------------------------
# Dispersal function
# Dispersal function
wf <- function(w, dist) {
  wf <- as.vector(exp(- w * dist ^ 2))
  wf <- wf / sum(wf)
  return(wf)
}

dispersal <- function(em, disp){
  em <- disp %*% em
  return(em)
}

#-----------------------------------------------------------------------------------
# Immigration
immi <- function(sp, m){
  sp <- sp + rpois(length(sp), m) 
  return(sp)
}

#------------------------------------------------------------------------------------
# proj: competitive model
proj <- function(sp, fec, l, intra, inter){
  for(i in 1:nrow(sp)){
    sp[i,] <- sapply(fec * sp[i,] * l[i,], function(x) rpois(1,x)) /
      (1 + intra * sp[i,] + inter *(sum(sp[i,]) - sp[i,]) )
  }
  return(sp)
}

#------------------------------------------------------------------------------------------------
# Metacommunity simulation function
simul <- function(land, rp, meta, file.name){
  
  # species environmental performance (lambda)
  for(s in 1:rp$sr){
    meta$l[,s] <- lambda(land$land$env, rp$niche$mu[s], rp$niche$sigma[s])
  }
  
  # dispersal matrix
  meta$disp <- t(apply(land$dist, 1, function(x) wf(dist = x, rp$w)))
  
  # starting matrix
  meta$sp <- matrix(0, land$n.si, rp$sr)
  start.input <- seq(1, 100, 10)
  save <- (land$start + land$burnin + 1):land$n.time
  
  # run simulation
  for (t in 1:land$n.time){
    print(t)
    if(t %in% start.input){
      meta$sp <- meta$sp + rpois(length(meta$sp), 1)
    }
    meta$em <- apply(meta$sp, 1:2, function(x) rpois(1, x * rp$a)) 
    meta$sp <- meta$sp - meta$em
    meta$em[meta$sp<0] <-  meta$em[meta$sp<0] + meta$sp[meta$sp<0]
    meta$sp[meta$sp<0] <- meta$em[meta$sp<0] - meta$em[meta$sp<0]
    if(sum(meta$em) > 0){
      meta$em <- dispersal(em = meta$em, disp = meta$disp)
    }
    meta$sp <- meta$sp + meta$em
    meta$sp <- immi(meta$sp, rp$m)
    meta$sp <- proj(sp = meta$sp, fec = rp$fec, l = meta$l[land$land$time == t,],
                    intra = rp$intra, inter = rp$inter)
    meta$sp[meta$sp<0.1] <- 0
    if(t %in% save){
      write.table(meta$sp, paste(file.name, "_", sprintf("%06d", t), ".csv", sep = ""),
                  sep = ",", row.names = FALSE, col.names = TRUE)  
    }
  }
}

#-----------------------------------------------------------------------------------------------
# Simulation routine
metasim <- function(meta.name, n.rep, sigma, w, inter, ...){
  land <- make.land(...)
  rp <- make.rp(meta.name = meta.name, sigma = sigma, w = w, inter = inter)
  meta <- make.meta(land$n.si, rp$sr, land$n.time)
  
  sim.dir <- paste("sim", meta.name, sep = "_")
  dir.create(sim.dir)
  sim.wd <- paste(wd, "/", sim.dir, sep  = "")
  
  saveRDS(rp, file = paste(sim.wd, "/rp_", meta.name,".rds", sep = ""))
  sink(paste(sim.wd, "/rp_", meta.name,".txt", sep = ""))
  print(rp)
  sink()
  
  rep <- 1:n.rep
  for (r in rep){
    dir.create(paste(sim.wd, "/", sim.dir, "_", sprintf("%02d", rep[r]), sep = ""))
  }
  dir.l <- list.dirs(sim.wd)[-1]
  file.l <- dir.l
  for(r in 1:length(file.l)){
    file.l[r] <-
      paste(dir.l[r], "/sp_", meta.name, "_", sprintf("%02d", rep[r]), sep = "")
  }
  file.name <- file.l[1]
  
  for(r in 1:n.rep){
    file.name <- file.l[r]
    saveRDS(land, file = paste(dir.l[r], "/", "land_", meta.name, "_",
                               sprintf("%02d", rep[r]), ".rds", sep = ""))
    simul(land = land, rp = rp, meta = meta, file.name = file.name)
    if(r != n.rep){
      land <- make.land(...)  
    }
  }
}
