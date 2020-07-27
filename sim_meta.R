#-----------------------------------------------------------------------------------------------
# Script: Metacommunity simulation routine
# Project: Temporal scale on invertebrates
# Date: 26/06/2020
# Author: Castillo-Escrivà, A., Mesquita-Joanes, F. and Rueda, J.
# e-mail: acastilloescriva@gmail.com
#------------------------------------------------------------------------------------------------
rm(list = ls())
wd <- "~/Escriptori" # change this to your working directory
setwd(wd)
source("sim_meta_function.R") # load additional functions

#-----------------------------------------------------------------------------------------------
# test - run one simulation
#meta.name <- "test"
#sigma <- 0.5
#inter <- 0.95
#a <- 0.4
#n.rep <- 1
#metasim(meta.name = meta.name, n.rep = n.rep, sigma = sigma, a = a, inter = inter)

#-----------------------------------------------------------------------------------------------
# run all the models
meta.name <- c("1ne_1lo", "1ne_2me", "1ne_3hi", "2ni_1lo", "2ni_2me", "2ni_3hi")
sigma <- rep(c(10, 0.5), each = 3)
inter <- rep(c(0.95, 0.55), each = 3)
w <- rep(c(10^-1, 10^-3, 10^-5), 2)

n.rep <- 10 # number of simulation replicates

wd <- "~/Escriptori/simulations_w" # change this to your working directory
dir.create(wd)
setwd(wd)
param <- data.frame(meta.name, sigma, w, inter)
for(i in 1:nrow(param)){
  metasim(meta.name = param$meta.name[i], n.rep = n.rep, sigma = param$sigma[i],
          w = param$w[i], inter = param$inter[i])
}
