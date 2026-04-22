# 4/15/24 smoothing
# 4/22/24 change average
# 9/17/24 average -> (floor(100x))/100
# 10/09/24 standardize dist and dem_30m
# 10/23/24 trees in 1985 <- 1985, 2020 <- 2020
# 8/12/25 data set create for paper
# 8/14/25 smooth by vectors

library(raster)
library(fields)
library(RColorBrewer)
library(rasterVis)
library(mvtnorm)
library(gridExtra)
library(fBasics)
library(coda)
library(tidyverse)
library(egg)

# number of grids N by N
N <- 300
# Length of grid
LoG <- 30
# number of years 1986-2020
TIME <- 35
# size in binomial
SIZE <- 100

setwd("/home/yawenguan/jiaqichen/Region88/paper/real data/data/N300")
dat.original <- read.csv("dat_8520_small_area_300_by_300_original.csv")

dat <- dat.original[, c("year", "s1", "s2", "Longitude", "Latitude", "Stream", "DEM_30m", "trees")]
index <- rep(1:(N*N), TIME+1)

dat <- cbind(index, dat)

# smoothing by average consecutive 3 years
trees.smth <- c()
prestyr <- dat[which(dat$year==1985), "trees"]
nextyr <- dat[which(dat$year==1986), "trees"]
smooth.trees.prest <- (prestyr + nextyr)/2
trees.smth <- c(trees.smth, round(smooth.trees.prest, digits = 0))

for (i in 1:(TIME-1)) {
  prest <- 1985+i
  prestyr <- dat[which(dat$year==prest), "trees"]
  prevyr <- dat[which(dat$year==(prest-1)), "trees"]
  nextyr <- dat[which(dat$year==(prest+1)), "trees"]
  
  smooth.trees.prest <- (prevyr + prestyr + nextyr)/3
  trees.smth <- c(trees.smth, round(smooth.trees.prest, digits = 0))
}

prestyr <- dat[which(dat$year==2020), "trees"]
prevyr <- dat[which(dat$year==2019), "trees"]
smooth.trees.prest <- (prevyr + prestyr)/2
trees.smth <- c(trees.smth, round(smooth.trees.prest, digits = 0))

dat$original.trees <- dat$trees
dat$trees <- trees.smth/SIZE

write.csv(dat, "dat_8520_small_area_300_by_300_smoothed_integered.csv")

# add a covariate as distance to the nearest stream
dat <- read.csv("dat_8520_small_area_300_by_300_smoothed_integered.csv")
data <- dat %>% filter(year == 1985)
coord <- cbind(data["s1"], data["s2"])
stream.cell <- data %>% filter(Stream == 1)
coord.stream <- cbind(stream.cell["s1"], stream.cell["s2"])
knnx <- get.knnx(coord.stream, coord, k = 1)
dist.str0 <- LoG*knnx$nn.dist
direc.h0 <- rep(0, dim(data)[1])
direc.v0 <- rep(0, dim(data)[1])
for (i in 1:dim(data)[1]) {
  str <- knnx$nn.index[i]
  direc.h0[i] <- LoG*(coord[i,"s1"]-coord.stream[str,"s1"])
  direc.v0[i] <- LoG*(coord[i,"s2"]-coord.stream[str,"s2"])
}

distance.str <- rep(dist.str0, TIME+1)
direc.to.str.h <- rep(direc.h0, TIME+1)
direc.to.str.v <- rep(direc.v0, TIME+1)

dat <- cbind(dat, distance.str, direc.to.str.h, direc.to.str.v)

# standardize dist, dem_30, Stream
dat$original.dist <- dat$distance.str
dat$original.DEM_30m <- dat$DEM_30m

dat_85 <- dat %>% filter(year == 1985)
dist_min <- min(dat_85$distance.str)
dist_max <- max(dat_85$distance.str)
dem30m_min <- min(dat_85$DEM_30m)
dem30m_max <- max(dat_85$DEM_30m)

dist85.std <- (dat_85$distance.str-dist_min)/(dist_max-dist_min)
dem30m.std <- (dat_85$DEM_30m-dem30m_min)/(dem30m_max-dem30m_min)

dat$distance.str <- rep(dist85.std, TIME+1)
dat$DEM_30m <- rep(dem30m.std, TIME+1)
dat$Stream <- rep(dat_85$Stream, TIME+1)

write.csv(dat, "dat_8520_small_area_300_by_300_smoothed_integered_standardized.csv")

# graph
dat85 <- dat %>% filter(year == 1985)
png(filename = "observation in 1985 smoothed 300 by 300.png")
ggplot(data = dat85, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 1985") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat90 <- dat %>% filter(year == 1990)
png(filename = "observation in 1990 smoothed 300 by 300.png")
ggplot(data = dat90, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 1990") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat95 <- dat %>% filter(year == 1995)
png(filename = "observation in 1995 smoothed 300 by 300.png")
ggplot(data = dat95, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 1995") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat00 <- dat %>% filter(year == 2000)
png(filename = "observation in 2000 smoothed 300 by 300.png")
ggplot(data = dat00, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 2000") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat05 <- dat %>% filter(year == 2005)
png(filename = "observation in 2005 smoothed 300 by 300.png")
ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 2005") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat10 <- dat %>% filter(year == 2010)
png(filename = "observation in 2010 smoothed 300 by 300.png")
ggplot(data = dat10, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 2010") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat15 <- dat %>% filter(year == 2015)
png(filename = "observation in 2015 smoothed 300 by 300.png")
ggplot(data = dat15, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 2015") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

dat20 <- dat %>% filter(year == 2020)
png(filename = "observation in 2020 smoothed 300 by 300.png")
ggplot(data = dat20, aes(s1*LoG, s2*LoG, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4") +
  labs(fill = "tree\ncoverage") +
  ggtitle("Observed Tree Coverage in 2020") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename="stream.png")
ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = Stream)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  ggtitle("Stream") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  coord_equal() 
dev.off()

png(filename="distance to stream.png")
ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = original.dist)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  labs(fill = "distance\nto\nstream") +
  ggtitle("Distance to Stream") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal() 
dev.off()

png(filename="elevation.png")
ggplot(data = dat05, aes(s1*LoG, s2*LoG, fill = original.DEM_30m)) +
  geom_tile() +
  ggtitle("Elevation") +
  scale_fill_gradientn(colours = terrain.colors(7)) +
  labs(fill = "Elevation\n(DEM_30m)") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal()
dev.off()

# aggregate to 100 by 100
setwd("/home/yawenguan/jiaqichen/Region88/paper/real data/data/N100")
N.100 <- 100
LoG.100 <- 90

Lmatrix <- matrix(0, nrow = N.100, ncol = N)
for (k in 1:N.100) {
  Lmatrix[k, 3*k-2] <- 1/3
  Lmatrix[k, 3*k-1] <- 1/3
  Lmatrix[k, 3*k] <- 1/3
}
Rmatrix <- t(Lmatrix)

# Stream
stream.matrix.85 <- matrix(dat85$Stream, nrow = N, ncol = N, byrow = F)
stream.matrix100.85 <- Lmatrix %*% stream.matrix.85 %*% Rmatrix
stream.vec100.85 <- as.vector(stream.matrix100.85)
stream.vec100.85 <- ceiling(stream.vec100.85)

# dem30m
dem30m.matrix.85 <- matrix(dat85$original.DEM_30m, nrow = N, ncol = N, byrow = F)
dem30m.matrix100.85 <- Lmatrix %*% dem30m.matrix.85 %*% Rmatrix
dem30m.vec100.85 <- as.vector(dem30m.matrix100.85)

# trees
trees <- c()
for (j in 1:(TIME+1)) {
  tree.vec300 <- dat[which(dat$year==1984+j), "original.trees"]
  tree.matrix <- matrix(tree.vec300, nrow = N, ncol = N, byrow = F)
  tree.matrix100 <- Lmatrix %*% tree.matrix %*% Rmatrix
  trees <- c(trees, as.vector(tree.matrix100))
}

s1.100 <- rep(rep(c(N.100:1), each = N.100), TIME+1)
s2.100 <- rep(rep(c(1:N.100), N.100), TIME+1)
dem30m <- rep(dem30m.vec100.85, TIME+1)
stream <- rep(stream.vec100.85, TIME+1)
year <- rep(c(1985:2020), each = N.100*N.100)
dat8520.N100 <- data.frame(trees = trees, DEM_30m = dem30m,
                      Stream = stream, year = year, s1 = s1.100, s2 = s2.100)

dat.N100 <- dat8520.N100

# convert coordinates to index, trees to percentage
dat.N100 <- dat.N100 %>% mutate(trees = round(trees, digits = 0)/100)

# add a covariate as distance to the nearest stream
# data <- dat %>% filter(year==1985)
# coord <- cbind(data["s1"], data["s2"])
# stream.cell <- data %>% filter(Stream == 1)
# coord.stream <- cbind(stream.cell["s1"], stream.cell["s2"])
# knnx <- get.knnx(coord.stream, coord, k = 1)
# dist.stream <- (N/100)*knnx$nn.dist
# direc.h <- rep(0, dim(data)[1])
# direc.v <- rep(0, dim(data)[1])
# for (i in 1:dim(data)[1]) {
#   str <- knnx$nn.index[i]
#   direc.h[i] <- (N/100)*(coord[i,"s1"]-coord.stream[str,"s1"])
#   direc.v[i] <- (N/100)*(coord[i,"s2"]-coord.stream[str,"s2"])
# }
# 
# distance.str <- rep(dist.stream, TIME+1)
# direc.to.str.h <- rep(direc.h, TIME+1)
# direc.to.str.v <- rep(direc.v, TIME+1)

# dat <- cbind(dat, distance.str, direc.to.str.h, direc.to.str.v)

index.N100 <- rep(c(1:(N.100 * N.100)), TIME+1)
dat.N100 <- cbind(index.N100, dat.N100)

write.csv(dat.N100, "dat_8520_small_area_100_by_100_smoothed_integered.csv")

pic_limit <- c(0,1)

data2015 <- dat.N100 %>% filter(year == 2015) %>% mutate(elevation = DEM_30m)
png(filename = "observation in 2015 with average aggreagation 100 by 100.png")
ggplot(data = data2015, aes(s1*LoG.100, s2*LoG.100, fill = trees)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="green4", limit = pic_limit) +
  ggtitle("Observation in 2015 with Average Aggregation 100 by 100") +
  labs(fill = "tree\ncoverage") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal()
dev.off()

png(filename = "stream with average aggregation 100 by 100.png")
ggplot(data = data2015, aes(s1*LoG.100, s2*LoG.100, fill = Stream)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  ggtitle("Stream with Average Aggregation 100 by 100") +
  theme_bw() +
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  coord_equal()
dev.off()

png(filename = "elevation with average aggregation 100 by 100.png")
ggplot(data = data2015, aes(s1*LoG.100, s2*LoG.100, fill = elevation)) +
  geom_tile() +
  ggtitle("Elevation with Average Aggregation 100 by 100") +
  labs(fill = "elevation\n(DEM_30m)") +
  theme_bw() +
  scale_fill_gradientn(colours = terrain.colors(7)) + 
  theme(plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal()
dev.off()

# dat <- read.csv("dat_8520_small_area_100_by_100_smoothed_integered_standardized.csv")

# dat <- read.csv("dat_8520_small_area_100_by_100_smoothed.csv")

# change average
# avg8520 <- dat %>% filter(trees != 0) %>% group_by(year) %>% summarise(Mean = mean(trees, na.rm = TRUE))
# plot(avg8520, type = "b", main = "Average trees without 0")
# avg.sl.8520 <- dat %>% group_by(year) %>% summarise(Mean.slope = mean(Slope_30m, na.rm = TRUE))
# plot(avg.sl.8520, type = "b", main = "Average Slope_30m")
# avg.tmean.8520 <- dat %>% group_by(year) %>% summarise(Mean.tmean = mean(tmean, na.rm = TRUE))
# plot(avg.tmean.8520, type = "b", main = "Average tmean")
# avg.dist.8520 <- dat %>% group_by(year) %>% summarise(Mean.dist = mean(distance.str, na.rm = TRUE))
# plot(avg.dist.8520, type = "b", main = "Average dist")
# avg.ppt.8520 <- dat %>% group_by(year) %>% summarise(Mean.ppt = mean(ppt, na.rm = TRUE))
# plot(avg.ppt.8520, type = "b", main = "Average ppt")

# avg.diff9196 <- (avg8520$Mean[12]-avg8520$Mean[7])/5
# avg.diff92 <- avg8520$Mean[7]+avg.diff9196-avg8520$Mean[8]
# avg.diff93 <- avg8520$Mean[7]+avg.diff9196*2-avg8520$Mean[9]
# avg.diff94 <- avg8520$Mean[7]+avg.diff9196*3-avg8520$Mean[10]
# avg.diff95 <- avg8520$Mean[7]+avg.diff9196*4-avg8520$Mean[11]
# 
# avg.diff0004 <- (avg8520$Mean[20]-avg8520$Mean[16])/4
# avg.diff01 <- avg8520$Mean[16]+avg.diff0004-avg8520$Mean[17]
# avg.diff02 <- avg8520$Mean[16]+avg.diff0004*2-avg8520$Mean[18]
# avg.diff03 <- avg8520$Mean[16]+avg.diff0004*3-avg8520$Mean[19]
# 
# avg.diff0813 <- (avg8520$Mean[29]-avg8520$Mean[24])/5
# avg.diff09 <- avg8520$Mean[24]+avg.diff0813-avg8520$Mean[25]
# avg.diff10 <- avg8520$Mean[24]+avg.diff0813*2-avg8520$Mean[26]
# avg.diff11 <- avg8520$Mean[24]+avg.diff0813*3-avg8520$Mean[27]
# avg.diff12 <- avg8520$Mean[24]+avg.diff0813*4-avg8520$Mean[28]

# dat <- dat %>% mutate(trees.smth = trees)
# for (i in 1:nrow(dat)) {
#   tree.value <- dat[i,"trees"]
#   if(dat[i,"Stream"] == 0 & tree.value !=0){
#     if(dat[i,"year"] == 1992){
#       dat[i,"trees"] <- ifelse(dat[i,"trees"]+avg.diff92 > 1, 1, dat[i,"trees"]+avg.diff92)
#     } else if(dat[i,"year"] == 1993){
#       dat[i,"trees"] <- ifelse(dat[i,"trees"]+avg.diff93 > 1, 1, dat[i,"trees"]+avg.diff93)
#     } else if(dat[i,"year"] == 1994){
#       dat[i,"trees"] <- ifelse(dat[i,"trees"]+avg.diff94 > 1, 1, dat[i,"trees"]+avg.diff94)
#     } else if(dat[i,"year"] == 1995){
#       dat[i,"trees"] <- ifelse(dat[i,"trees"]+avg.diff95 > 1, 1, dat[i,"trees"]+avg.diff95)
#     } else if(dat[i,"year"] == 2001){
#       dat[i,"trees"] <- ifelse(dat[i,"trees"]+avg.diff01 > 1, 1, dat[i,"trees"]+avg.diff01)
#     } else if(dat[i,"year"] == 2002){
#       dat[i,"trees"] <- ifelse(dat[i,"trees"]+avg.diff02 > 1, 1, dat[i,"trees"]+avg.diff02)
#     } else if(dat[i,"year"] == 2003){
#       dat[i,"trees"] <- ifelse(dat[i,"trees"]+avg.diff03 > 1, 1, dat[i,"trees"]+avg.diff03)
#     } else if(dat[i,"year"] == 2009){
#       dat[i,"trees"] <- ifelse(dat[i,"trees"]+avg.diff09 > 1, 1, dat[i,"trees"]+avg.diff09)
#     } else if(dat[i,"year"] == 2010){
#       dat[i,"trees"] <- ifelse(dat[i,"trees"]+avg.diff10 > 1, 1, dat[i,"trees"]+avg.diff10)
#     } else if(dat[i,"year"] == 2011){
#       dat[i,"trees"] <- ifelse(dat[i,"trees"]+avg.diff11 > 1, 1, dat[i,"trees"]+avg.diff11)
#     } else if(dat[i,"year"] == 2012){
#       dat[i,"trees"] <- ifelse(dat[i,"trees"]+avg.diff12 > 1, 1, dat[i,"trees"]+avg.diff12)
#     }
#   } else{
#     dat[i,"trees"] <- tree.value
#   }
# }
# 
# avg8520 <- dat %>% filter(trees != 0) %>% group_by(year) %>% summarise(Mean = mean(trees, na.rm = TRUE))
# plot(avg8520, type = "b", main = "Average trees without 0")
# write.csv(dat, "dat_8520_small_area_100_by_100_cleaned.csv")

# result display
# dat85 <- dat %>% filter(year == 1985)
# ggplot(data = dat85, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 1985") +
#   coord_equal() 
# 
# dat86 <- dat %>% filter(year == 1986)
# ggplot(data = dat86, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 1986") +
#   coord_equal()
# 
# dat87 <- dat %>% filter(year == 1987)
# ggplot(data = dat87, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 1987") +
#   coord_equal()
# 
# dat88 <- dat %>% filter(year == 1988)
# ggplot(data = dat88, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 1988") +
#   coord_equal() 
# 
# dat89 <- dat %>% filter(year == 1989)
# ggplot(data = dat89, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 1989") +
#   coord_equal() 
# 
# dat90 <- dat %>% filter(year == 1990)
# ggplot(data = dat90, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 1990") +
#   coord_equal() 
# 
# 
# dat91 <- dat %>% filter(year == 1991)
# ggplot(data = dat91, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 1991") +
#   coord_equal()
# 
# dat92 <- dat %>% filter(year == 1992)
# ggplot(data = dat92, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 1992") +
#   coord_equal()
# 
# dat93 <- dat %>% filter(year == 1993)
# ggplot(data = dat93, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 1993") +
#   coord_equal()
# 
# dat94 <- dat %>% filter(year == 1994)
# ggplot(data = dat94, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 1994") +
#   coord_equal()
# 
# dat95 <- dat %>% filter(year == 1995)
# ggplot(data = dat95, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 1995") +
#   coord_equal()
# 
# dat96 <- dat %>% filter(year == 1996)
# ggplot(data = dat96, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 1996") +
#   coord_equal()
# 
# dat97 <- dat %>% filter(year == 1997)
# ggplot(data = dat97, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 1997") +
#   coord_equal()
# 
# dat98 <- dat %>% filter(year == 1998)
# ggplot(data = dat98, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 1998") +
#   coord_equal()
# 
# dat99 <- dat %>% filter(year == 1999)
# ggplot(data = dat99, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 1999") +
#   coord_equal()
# 
# dat00 <- dat %>% filter(year == 2000)
# ggplot(data = dat00, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2000") +
#   coord_equal()
# 
# dat01 <- dat %>% filter(year == 2001)
# ggplot(data = dat01, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2001") +
#   coord_equal()
# 
# dat02 <- dat %>% filter(year == 2002)
# ggplot(data = dat02, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2002") +
#   coord_equal()
# 
# dat03 <- dat %>% filter(year == 2003)
# ggplot(data = dat03, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2003") +
#   coord_equal()
# 
# dat04 <- dat %>% filter(year == 2004)
# ggplot(data = dat04, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2004") +
#   coord_equal()
# 
# dat05 <- dat %>% filter(year == 2005)
# ggplot(data = dat05, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2005") +
#   coord_equal()
# 
# dat06 <- dat %>% filter(year == 2006)
# ggplot(data = dat06, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2006") +
#   coord_equal()
# 
# dat07 <- dat %>% filter(year == 2007)
# ggplot(data = dat07, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2007") +
#   coord_equal()
# 
# dat08 <- dat %>% filter(year == 2008)
# ggplot(data = dat08, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2008") +
#   coord_equal()
# 
# dat09 <- dat %>% filter(year == 2009)
# ggplot(data = dat09, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2009") +
#   coord_equal()
# 
# dat10 <- dat %>% filter(year == 2010)
# ggplot(data = dat10, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2010") +
#   coord_equal()
# 
# dat11 <- dat %>% filter(year == 2011)
# ggplot(data = dat11, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2011") +
#   coord_equal()
# 
# dat12 <- dat %>% filter(year == 2012)
# ggplot(data = dat12, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2012") +
#   coord_equal()
# 
# dat13 <- dat %>% filter(year == 2013)
# ggplot(data = dat13, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2013") +
#   coord_equal()
# 
# dat14 <- dat %>% filter(year == 2014)
# ggplot(data = dat14, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2014") +
#   coord_equal()
# 
# dat15 <- dat %>% filter(year == 2015)
# ggplot(data = dat15, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2015") +
#   coord_equal()
# 
# dat16 <- dat %>% filter(year == 2016)
# ggplot(data = dat16, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2016") +
#   coord_equal()
# 
# dat17 <- dat %>% filter(year == 2017)
# ggplot(data = dat17, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2017") +
#   coord_equal()
# 
# dat18 <- dat %>% filter(year == 2018)
# ggplot(data = dat18, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2018") +
#   coord_equal()
# 
# dat19 <- dat %>% filter(year == 2019)
# ggplot(data = dat19, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2019") +
#   coord_equal()
# 
# dat20 <- dat %>% filter(year == 2020)
# ggplot(data = dat20, aes(s1, s2, fill = trees)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="blue", limit = c(0,1)) +
#   labs(title = "observation in 2020") +
#   coord_equal()