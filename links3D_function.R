#-----------------------------------------------------------------------------------------------
# Script: Functions for spatio-temporal links
# Project: Temporal scale on invertebrates
# Date: 24/06/2020
# Author: Castillo-Escrivà, A., Mesquita-Joanes, F. and Rueda, J.
# e-mail: acastilloescriva@gmail.com
#------------------------------------------------------------------------------------------------

# create links in space (x,y) and time (z)
xyz2links <- function(names, coords){
  require(spdep)
  names <- as.data.frame(lapply(names, factor))
  
  coords <- as.data.frame(lapply(coords, as.numeric))
  
  n.site <- nlevels(names[,1])
  n.date <- nlevels(names[,2])
  
  # spatial links
  coord.list <- split(coords[,1:2], names[,2])
  nb.spa <- lapply(coord.list, function(x) graph2nb(gabrielneigh(as.matrix(x)), sym = TRUE))
  nb.spa <- nb.spa[which(sapply(nb.spa, function(x) sum(card(x))) != 0)]
  nb.spa <- lapply(nb.spa, function(x) which(unlist(nb2mat(x)) != 0, arr.ind = TRUE))
  pos <- sapply(1:length(coord.list), function(x) as.numeric(rownames(coord.list[[x]])[1]))
  pos <- pos[sapply(coord.list, nrow) > 1]
  nb.spa <- lapply(1:length(nb.spa), function(x) nb.spa[[x]] + pos[x] - 1)
  nb.spa <- do.call("rbind", nb.spa)
  
  # temporal links
  tim.list <- split(1:nrow(names), names[,1])
  nb.tim <- lapply(tim.list, function(x) cbind(x[-length(x)], x[-1]))
  nb.tim <- do.call("rbind", nb.tim)
  nb.tim <- rbind(nb.tim, nb.tim[,c(2,1)])
  
  # spatiotemporal links
  nb.st <- rbind(nb.spa, nb.tim)
  mat.st <- array(0, rep(nrow(names), 2))
  mat.st[nb.st] <- 1
  listw.st <- mat2listw(mat.st, style = "W")
  
  nb <- list("st" = nb.st, "spa" = nb.spa, "tim" = nb.tim)
  return(list("listw" = listw.st, "nb" = nb))
}

#-----------------------------------------------------------------------------------------------
# plot links in 3D

plot.st <- function(nb.st, coords, col = c("black", "red", "blue"), ...){
  require(plot3D)
  if(class(coords[,3]) != "numeric"){
    coords[,3] <- as.Date(coords[,3])  
  }
  coords <- as.data.frame(lapply(coords, as.numeric))
  scatter3D(x = coords[,1], y = coords[,2], z = coords[,3], col = col[1], zlab = "time", ...)
  for(i in 2:3){
    segments3D(x0 = coords[nb.st[[i]][,1],1], y0 = coords[nb.st[[i]][,1],2],
               z0 = coords[nb.st[[i]][,1],3], x1 = coords[nb.st[[i]][,2],1],
               y1 = coords[nb.st[[i]][,2],2], z1 = coords[nb.st[[i]][,2],3],
               col = col[i], add = TRUE, ...)
  }
}

#-----------------------------------------------------------------------------------------------
# plot a variable (like a MEM) through space and time

plot.st.var <- function(var, coords, col = c("black", "red", "blue"), ...){
  par(mar = rep(2,4))
  require(plot3D)
  if(class(coords[,3]) != "numeric"){
    coords[,3] <- as.Date(coords[,3])  
  }
  coords <- as.data.frame(lapply(coords, as.numeric))
  scatter3D(x = coords[,1], y = coords[,2], z = coords[,3], colvar = var,
            zlab = "time", pch = 19, ...)
}

#---------------------------------------------------------------------------------------------
# test spatiotemporal mems (two-way anova)

test.st <- function(mem, space, time){
  space <- as.factor(space)
  time <- as.factor(time)
  p.value <- matrix(0, nrow = ncol(mem), ncol = 3)
  for(i in 1:nrow(p.value)){
    p.value[i,] <- summary(aov(mem[,i] ~ space * time))[[1]][1:3,5]
  }
  colnames(p.value) <- c("space", "time", "interaction")
  rownames(p.value) <- colnames(mem)
  return(as.data.frame(p.value))
}

