#-----------------------------------------------------------------------------------------------
# Script: Sampling and analyzing simulated metacommunities routine
# Project: Temporal scale on invertebrates
# Date: 29/06/2020
# Author: Castillo-Escrivà, A., Mesquita-Joanes, F. and Rueda, J.
# e-mail: acastilloescriva@gmail.com
#------------------------------------------------------------------------------------------------
rm(list = ls())
setwd("~/Escriptori") # here set your working directory
library(vegan)
library(adespatial)
source("part_meta_function.R")
source("links3D_function.R")
source("vpmsr_function.R")

#------------------------------------------------------------------------------------------------
# test only one simulation
#setwd("~/Escriptori/sim_test/") # here set your working directory
#dir <- "sim_test_01"
# Sampling metacommunities
#meta <- sim.pre(sim.dir = dir, time.extent = 20, n.time = 10, n.site = 10, rm.start = TRUE)
#colSums(meta$sp)

# Variation partitioning
#vp <- test.meta(sp = meta$sp, sample = meta$sample, env = meta$env, method = "rda",
#                nrepet.msr = 99, poly2 = FALSE, perm.cca = 99, env.pca = FALSE)
#round(vp$vp,2)

#-------------------------------------------------------------------------------------------------
# Variation partitioning in all the simulations at different temporal extent
# List simulations
setwd("~/Escriptori/simulations_w/") # here set your working directory
dir <- list.dirs(path = ".", recursive = FALSE)
dir <- unlist(lapply(dir, function(x) list.dirs(x)[-1]))
# Time extent
time.extent <- c(20, 40, 80, 160, 320, 640, 1280, 2000)

# Sampling metacommunities and analyzing
method <- "rda"
poly2 <- TRUE
file.name <- "vp1_sim_rda2.csv"

vp <- as.data.frame(matrix(0, nrow = length(dir)*length(time.extent), ncol = 17))
colnames(vp) <- c("name", "sigma", "disp", "time.extent", "rich", "n.env", "n.spa", "n.tim",
                  "n.st", letters[1:8]) 
vp$name <- rep(dir, each = length(time.extent))
vp$sigma <- rep(c(10, 0.5), each = length(dir)*length(time.extent)/2)
vp$disp <- rep(rep(c(10^-1, 10^-3, 10^-5), each = length(dir)*length(time.extent)/2/3), 2)
vp$time.extent <- rep(time.extent, length(dir))

meta <- sim.pre(sim.dir = vp$name[1], time.extent = vp$time.extent[1], n.time = 10,
                n.site = 10, rm.start = TRUE)
sp <- meta$sp
test <- test.meta(sp = meta$sp, sample = meta$sample, env = meta$env, method = method,
                  nrepet.msr = 99, perm.cca = 499, poly2 = poly2, env.pca = FALSE)

for(i in 1:nrow(vp)){
  # Sampling metacommunities
  print(i)
  meta <- sim.pre(sim.dir = vp$name[i], time.extent = vp$time.extent[i],
                  n.time = 10, n.site = 10, rm.start = TRUE)
  
  # total species in the metacommunity
  sp <- meta$sp[, colSums(meta$sp) !=0]
  vp$rich[i] <- ncol(sp)
  
  # Variation partitioning
  test <- test.meta(sp = meta$sp, sample = meta$sample, env = meta$env, method = method,
                    nrepet.msr = 99, perm.cca = 499, poly2 = poly2, env.pca = FALSE)
  vp[i, 6:9] <- test$n.var
  vp[i,10:17] <- test$vp
  print(vp[i, 10:17])
}

# check the filename
write.table(vp, file = file.name, sep = ",", row.names = TRUE, col.names = NA)
# check the table

rm(list = ls())

###############################################################################################
###############################################################################################
# Variation partiotning in empirical data
setwd("~/Escriptori")
source("part_meta_function.R")
source("links3D_function.R")
source("vpmsr_function.R")

# load data
sample <- read.table("https://docs.google.com/spreadsheets/d/e/2PACX-1vQuvpAMC9LIyPbXhQl-tE3JioMxErRgjqhTHIDcAbQliGnYn3Ta-FfC9NGwW8RmQJegQ9vZn3M_ceFl/pub?gid=2085174889&single=true&output=csv",
                     header = TRUE, row.names = 1, sep = ",")
env <- read.table("https://docs.google.com/spreadsheets/d/e/2PACX-1vQuvpAMC9LIyPbXhQl-tE3JioMxErRgjqhTHIDcAbQliGnYn3Ta-FfC9NGwW8RmQJegQ9vZn3M_ceFl/pub?gid=313905911&single=true&output=csv",
                  header = TRUE, row.names = 1, sep = ",")
sp <- read.table("https://docs.google.com/spreadsheets/d/e/2PACX-1vQuvpAMC9LIyPbXhQl-tE3JioMxErRgjqhTHIDcAbQliGnYn3Ta-FfC9NGwW8RmQJegQ9vZn3M_ceFl/pub?gid=1793471961&single=true&output=csv",
                 header = TRUE, row.names = 1, sep = ",")
trait <- read.table("https://docs.google.com/spreadsheets/d/e/2PACX-1vQuvpAMC9LIyPbXhQl-tE3JioMxErRgjqhTHIDcAbQliGnYn3Ta-FfC9NGwW8RmQJegQ9vZn3M_ceFl/pub?gid=206857959&single=true&output=csv",
                    header = TRUE, row.names = 1, sep = ",")

#-------------------------------------------------------------------------------------------------
# Subsampling empirical data - test
#meta <- emp.pre(sample = sample, sp = sp, env = env, trait = trait, land = "raco",
#                group = "dytiscidae", start = 1, last = 6, lag = 1)
#nrow(meta$sp) # number of samples
#mean(dist(scale(meta$env)))
#s <- matrix(0, nrow = nlevels(meta$sample$site), ncol = ncol(meta$env)) 
#for (i in 1:ncol(meta$env)){
#  s[,i] <- tapply(meta$env[,i], meta$sample$site, mean)
#}
#s <- s[!is.na(rowSums(s)), ]
#mean(dist(scale(s)))
#t <- matrix(0, nrow = nlevels(meta$sample$time), ncol = ncol(meta$env)) 
#for (i in 1:ncol(meta$env)){
#  t[,i] <- tapply(meta$env[,i], meta$sample$time, mean)
#}
#t <- t[!is.na(rowSums(t)), ]
#mean(dist(scale(t)))
#
#vp <- test.meta.emp(sp = meta$sp, sample = meta$sample, env = meta$env, method = "rda",
#                    nrepet.msr = 99, perm.cca = 99, poly2 = TRUE)
#vp

#--------------------------------------------------------------------------------------------------
# Subsampling empirical data
# table with the tests
method = "capscale"
poly2 = FALSE
file.name <- "vp_emp_dbrda.csv"

vp <- as.data.frame(matrix(NA, nrow = 36, ncol = 21))
colnames(vp) <- c("land", "group", "start", "last", "lag", "rich", "n.sample", "n.site",
                  "n.time", "n.env", "n.spa", "n.tim", "n.st",letters[1:8])
vp[,1] <- c(rep("magre", 3*4), rep("mallades", 3*4), rep("raco", 3*4))
vp[,2] <- rep(c("ostracoda", "odonata", "dytiscidae",
                "odonata", "dytiscidae", "diptera",
                "ostracoda", "dytiscidae", "diptera"), each = 4)
vp[,3] <- c(rep(c(1, 1, 6, 1), 3), rep(c(1, 1, 7, 1), 6))
vp[,4] <- c(rep(c(11, 6, 11, 11), 3), rep(c(12, 6, 12, 12), 6))
vp[,5] <- rep(c(1, 1, 1, 2), 9)

meta <- emp.pre(sample = sample, sp = sp, env = env, trait = trait, land = vp$land[1],
                group = vp$group[1], start = vp$start[1], last = vp$last[1], lag = vp$lag[1])

# variation partitioning
for(i in 1:nrow(vp)){
  print(i)
  print(paste(vp$land[i], vp$group[i], sep = "_"))
  meta <- emp.pre(sample = sample, sp = sp, env = env, trait = trait, land = vp$land[i],
                  group = vp$group[i], start = vp$start[i], last = vp$last[i], lag = vp$lag[i])
  
  vp[i,6] <- ncol(meta$sp)
  vp[i,7] <- nrow(meta$sample)
  vp[i,8] <- length(unique(meta$sample$site))
  vp[i,9] <- length(unique(meta$sample$time))
  
  if(nrow(meta$sample) > 12){
    test <- test.meta(sp = meta$sp, sample = meta$sample, env = meta$env,
                                  method = method, nrepet.msr = 99, perm.cca = 999, poly2 = poly2,
                                  env.pca = TRUE)  
    vp[i, 10:13] <- test$n.var
    vp[i,14:21] <- test$vp
  }
  print(vp[i, 14:21])
}

# check the filename
write.table(vp, file = file.name, sep = ",", row.names = TRUE, col.names = NA)
# check the table

rm(list = ls())
