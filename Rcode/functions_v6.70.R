# 12/30/23 remove edge in neighbor matrix
# 1/1/24 apply Hefley's neighborhood function
# 1/3/24 my harmonic function
# 1/9/24 if else fixed
# 1/13/24 numerical solution with homogenization
# 1/14/24 numerical solution of Fokker-Planck equation without homogenization
# 2/6/24 exponential growth term
# 2/9/24 logistic growth
# 2/28/24 H matrix move outside each time loop
# 3/5/24 zero-inflated beta likelihood and random sampling
# 3/5/24 zero-one-inflated beta likelihood and random sampling
# 3/26/24 optimization delete advection, death, exp growth, zero.rm
# 7/11/24 zero-inflated binomial
# 7/19/24 Rcpp
# 7/24/24 data.frame -> matrix
# 7/26/24 H2 matrix -> vec
# 8/5/24 loglikehood -> mapply?
# 8/5/24 n,epsilon,neighborhood function -> outside PDEsolution function
# 8/7/24 PDE solution function modification (sparse, ordinary, Rcpp)
# 8/12/24 neighborhood function fixed, ln, rn, bn, tn
# 8/12/24 Rcpp matrix, computation
# 8/18/24 stream is not edge, diffu.rate[edge/stream] not = 0 and related modification
# 8/18/24 modify mymeanfun(), myrstdevidefun(), myrstdevide2fun()
# 8/18/24 apply mymeanfun instead of mean
# 1/27/25 numerical solution of PDE without homogenization for simulation
# 1/28/25 neighborhood function update, ln->rn, rn->ln, to align with the longitude of the dataset, right hand side large longitude
# 2/3/25 homog-PDE solution function -> used homogenized stream vector, value converts to 0 at each step
# 2/3/25 nonhomog-PDE solution -> stream converts to 0 at each step
# 2/12/25 mynumericalPDEfun() error fixed
# 3/10/25 add fokker-planck exponential growth 
# 3/10/25 add PDE homogenization solution with exponential growth
# 3/14/25 adjust fokker-planck exponential growth 
# 3/14/25 adjust PDE homogenization solution with exponential growth
# 3/30/25 homogenization of mu0 in PDE output
# 6/20/25 keep zero-inflated binomial
# 6/20/25 keep one numericalPDE logistic growth with homogenization function and one numericalPDE logistic growth function
# 6/20/25 add length scale variable
# 6/21/25 c0 -> average of mu*D

library(raster, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
library(fields, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
library(RColorBrewer, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
library(rasterVis, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
library(mvtnorm, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
library(gridExtra, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
library(fBasics, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
# library(coda, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
library(dplyr, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
library(egg, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
# library(Rcpp, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
library(Rcpp, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
library(RcppArmadillo, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
library(Matrix, lib.loc = "/home/yawenguan/jiaqichen/test100/library")

sourceCpp("/home/yawenguan/jiaqichen/test100/Rcode/PDE_fun.cpp")
sourceCpp("/home/yawenguan/jiaqichen/test100/Rcode/H_matrix.cpp")

# the zero-inflated binomial likelihood
myZIbinomlikelihood <- function(x, gamma, p, s, log = TRUE){
  output <- NA
  if(sum(x<0|x>s)>0){
    likelihd <- 0
  } else{
    likelihd <- I(x==0)*(gamma + (1-gamma)*(1-p)^s) +
      I(x>0)*((1-gamma)*dbinom(x, size = s, prob = p, log = FALSE))
  }
  
  if(log==TRUE){
    loglikelihd <- log(likelihd)
    output <- loglikelihd
  } else{
    output <- likelihd
  }
  return(output)
}

loglikhd_apply_fun <- function(edge, resp, prob, gamma.v, size){
  llikhd <- ifelse(edge==0, 
                   myZIbinomlikelihood(x = resp * 100, 
                                       gamma = gamma.v,
                                       p = prob, 
                                       s = size,
                                       log = TRUE),
                   0)
  return(llikhd)
}

# the zero-inflated binomial sampler
myrZIbinom <- function(n, gamma, p, s){
  output <- rep(0, n)
  for (i in 1:n) {
    r <- runif(1, 0, 1)
    if(r <= gamma){
      output[i] <- 0
    } else{
      output[i] <- rbinom(n = 1, size = s, prob = p)
    }
  }
  
  return(output)
}

# the neighborhood information
neighborhood <- function(raster){
  nn <- matrix(rep(0, length(raster[ ])), nrow = length(raster[]), ncol = 4)
  for(i in 1:dim(nn)[1]){
    loc <- adjacent(raster,i)[ ,2]
    ln <- loc[which((loc-dim(raster)[1])==i)]
    rn <- loc[which((loc+dim(raster)[1])==i)]
    bn <- loc[which((loc+1)==i)]
    tn <- loc[which((loc-1)==i)]
    nn[i,1] <- if(length(ln)>0){ln} else{0}
    nn[i,2] <- if(length(rn)>0){rn} else{0}
    nn[i,3] <- if(length(bn)>0){bn} else{0}
    nn[i,4] <- if(length(tn)>0){tn} else{0}
  }
  nn
}

# H matrix
propagator.1st <- function(neighbor, diffu.rate, 
                           growth.rate, growth.rate2, dx, dy, dt){
  H <- matrix(0, nrow = dim(neighbor)[1], ncol = dim(neighbor)[1])
  for (i in 1:dim(neighbor)[1]) {
    if(sum(neighbor[i,]>0)==4)
    {
      H[i,i] <- 1-2*diffu.rate[i]*dt/dx^2-2*diffu.rate[i]*dt/dy^2+
        growth.rate[i]*dt
      H[i,neighbor[i,1]] <- diffu.rate[i]*dt/dx^2
      H[i,neighbor[i,2]] <- diffu.rate[i]*dt/dx^2
      H[i,neighbor[i,3]] <- diffu.rate[i]*dt/dy^2
      H[i,neighbor[i,4]] <- diffu.rate[i]*dt/dy^2
    }
  }
  return(H)
}

propagator.2nd <- function(neighbor, diffu.rate, 
                           growth.rate, growth.rate2, dx, dy, dt){
  H2 <- matrix(0, nrow = dim(neighbor)[1], ncol = dim(neighbor)[1])
  for (i in 1:dim(neighbor)[1]) {
    if(sum(neighbor[i,]>0)==4){
      H2[i,i] <- -growth.rate2[i]*dt
    }
  }
  return(H2)
}

propagator.2nd.vec <- function(neighbor, diffu.rate, 
                               growth.rate, growth.rate2, dx, dy, dt){
  H2.vec <- rep(0, dim(neighbor)[1])
  for (i in 1:dim(neighbor)[1]) {
    if(sum(neighbor[i,]>0)==4){
      H2.vec[i] <- -growth.rate2[i]*dt
    }
  }
  return(H2.vec)
}

# H2 matrix for exp growth
propagator.2nd.exp.growth <- function(neighbor, diffu.rate, 
                                      growth.rate, growth.rate2, dx, dy, dt){
  H2 <- matrix(0, nrow = dim(neighbor)[1], ncol = dim(neighbor)[1])
  for (i in 1:dim(neighbor)[1]) {
    if(sum(neighbor[i,]>0)==4){
      H2[i,i] <- 0
    }
  }
  return(H2)
}

# H matrix fokker planck
propagator.FP.1st <- function(neighbor, diffu.rate, 
                              growth.rate, dx, dy, dt){
  H <- matrix(0, nrow = dim(neighbor)[1], ncol = dim(neighbor)[1])
  for (i in 1:dim(neighbor)[1]) {
    if(sum(neighbor[i,]>0)==4)
    {
      H[i,i] <- 1-2*diffu.rate[i]*dt/dx^2-2*diffu.rate[i]*dt/dy^2+
        growth.rate[i]*dt
      H[i,neighbor[i,1]] <- diffu.rate[neighbor[i,1]]*dt/dx^2
      H[i,neighbor[i,2]] <- diffu.rate[neighbor[i,2]]*dt/dx^2
      H[i,neighbor[i,3]] <- diffu.rate[neighbor[i,3]]*dt/dy^2
      H[i,neighbor[i,4]] <- diffu.rate[neighbor[i,4]]*dt/dy^2
    }
  }
  return(H)
}

propagator.FP.2nd <- function(neighbor, diffu.rate, 
                              growth.rate, dx, dy, dt){
  H2 <- matrix(0, nrow = dim(neighbor)[1], ncol = dim(neighbor)[1])
  for (i in 1:dim(neighbor)[1]) {
    if(sum(neighbor[i,]>0)==4){
      H2[i,i] <- -growth.rate[i]*dt
    }
  }
  return(H2)
}

propagator.FP.2nd.exp.growth <- function(neighbor, diffu.rate, 
                                         growth.rate, dx, dy, dt){
  H2 <- matrix(0, nrow = dim(neighbor)[1], ncol = dim(neighbor)[1])
  for (i in 1:dim(neighbor)[1]) {
    if(sum(neighbor[i,]>0)==4){
      H2[i,i] <- 0
    }
  }
  return(H2)
}

# my harmonic mean function (only remove 0 grid)
myharmonicmeanfun <- function(x, na.rm){
  result <- NA
  n <- length(x)
  if(sum(x==0) == n){
    result <- 0
  } else{
    vec <- x[x!=0]
    result <- 1/mean(1/vec)
  }
  return(result)
}

# my mean function (only remove Inf grid)
mymeanfun <- function(x, na.rm){
  result <- NA
  n <- length(x)
  # if(sum(x==Inf | x==-Inf) == n){
  if(sum(x==Inf) + sum(x==-Inf) == n){
    result <- 0
  } else{
    vec <- x[x!=Inf & x!=-Inf]
    result <- mean(vec)
  }
  return(result)
}

# my division function on raster
myrstdividefun <- function(x, y, na.rm){
  n1 <- dim(x)[1]
  n2 <- dim(x)[2]
  # result <- ifelse(y[]==0, 0, x[]/y[])
  result <- x[]/y[]
  one.mat <- matrix(rep(1, n1*n2), nrow = n1, ncol = n2)
  result.rst <- raster(one.mat, xmn = 0, xmx = 1, ymn = 0, ymx = 1, 
                       crs = NA)
  result.rst[] <- result
  return(result.rst)
}

# my division square function on raster
myrstdivide2fun <- function(x, y, na.rm){
  n1 <- dim(x)[1]
  n2 <- dim(x)[2]
  # result <- ifelse(y[]==0, 0, x[]/(y[]^2))
  result <- x[]/(y[]^2)
  one.mat <- matrix(rep(1, n1*n2), nrow = n1, ncol = n2)
  result.rst <- raster(one.mat, xmn = 0, xmx = 1, ymn = 0, ymx = 1, 
                       crs = NA)
  result.rst[] <- result
  return(result.rst)
}

# my numerical solution with homogenization function with Rcpp matrix (without advection and death)
mynumericalPDEhomofun.Rcpp_mat <- function(n, epsilon, initial.vec, stream.vec,
                                           diffu.vec, growth.vec, 
                                           length_of_grid, time, step, delta_t, 
                                           min.value = 0, max.value = 1){
  # m <- n/epsilon
  # 
  # index.n <- c(1:(n*n))
  # s1.n <- rep(c(1:n), each = n)
  # s2.n <- rep(c(1:n), n)
  # 
  # index.m <- c(1:(m*m))
  # s1.m <- rep(c(1:m), each = m)
  # s2.m <- rep(c(1:m), m)
  
  # original raster and rescaling raster
  # onematrix.NbyN <- matrix(rep(1, n*n), nrow = n, ncol = n)
  # origin.rst <- raster(onematrix.NbyN, xmn = 0, xmx = 1, ymn = 0, ymx = 1, 
  #                      crs = NA)
  # rescaling.rst <- aggregate(origin.rst, fact = epsilon, na.rm = TRUE, 
  #                            fun = myharmonicmeanfun)
  # 
  # neighbormatrix <- neighborhood(rescaling.rst)
  
  diffu.rst <- origin.rst
  diffu.rst[] <- diffu.vec
  
  diffu.bar <- aggregate(diffu.rst, fact = epsilon, na.rm = TRUE, 
                         fun = myharmonicmeanfun)
  
  # growth rate and homogenization
  growth.rst <- origin.rst
  growth.rst[] <- growth.vec
  
  growth.bar <- diffu.bar * aggregate(myrstdividefun(growth.rst, diffu.rst), 
                                      fact = epsilon, 
                                      na.rm = TRUE, 
                                      fun = mymeanfun)
  
  growth2.bar <- diffu.bar * aggregate(myrstdivide2fun(growth.rst, diffu.rst), 
                                       fact = epsilon, 
                                       na.rm = TRUE, 
                                       fun = mymeanfun)
  
  # initial value for PDE and its rescaling
  mu0 <- initial.vec
  mu.rst <- origin.rst
  mu.rst[] <- mu0
  
  c0 <- rescaling.rst
  c0 <- aggregate(mu.rst * diffu.rst, fact = epsilon, na.rm = TRUE, fun = mymeanfun)
  # c0[] <- raster::extract(mu.rst * diffu.rst, SpatialPoints(rescaling.rst))
  
  # will be the same if set cmat <- c0[], since matrix multiplication will convert it to matrix
  cmat <- as.matrix(c0[])
  
  # rate coef for the function
  diffu.rc <- diffu.bar[]
  growth.rc <- growth.bar[]
  growth2.rc <- growth2.bar[]
  
  H1.mat <- propagator_1st_Rcpp(neighbor = neighbormatrix, 
                                diffu_rate = diffu.rc, 
                                growth_rate = growth.rc, 
                                growth_rate2 = growth2.rc,
                                dx = epsilon*length_of_grid, 
                                dy = epsilon*length_of_grid, 
                                dt = delta_t)
  
  # H1.mat <- propagator_1st_Rcpp(neighbor = neighbormatrix, 
  #                               diffu_rate = diffu.rc, 
  #                               growth_rate = growth.rc, 
  #                               growth_rate2 = growth2.rc,
  #                               dx = epsilon/n, dy = epsilon/n, dt = delta_t)
  # H1.sp <- Matrix(H1.mat, sparse = TRUE)
  
  # H2.mat <- propagator.2nd(neighbor = neighbormatrix, diffu.rate = diffu.rc,
  #                          growth.rate = growth.rc, growth.rate2 = growth2.rc,
  #                          dx = epsilon/n, dy = epsilon/n, dt = delta_t)
  
  H2.vec <- propagator_2nd_vec_Rcpp(neighbor = neighbormatrix, 
                                    diffu_rate = diffu.rc, 
                                    growth_rate = growth.rc, 
                                    growth_rate2 = growth2.rc,
                                    dx = epsilon*length_of_grid, 
                                    dy = epsilon*length_of_grid, 
                                    dt = delta_t)
  
  # H2.vec <- propagator_2nd_vec_Rcpp(neighbor = neighbormatrix, 
  #                                   diffu_rate = diffu.rc, 
  #                                   growth_rate = growth.rc, 
  #                                   growth_rate2 = growth2.rc,
  #                                   dx = epsilon/n, dy = epsilon/n, dt = delta_t)
  
  # set up result matrix
  cmat <- cmat * (1-stream.vec)
  c.rst <- rescaling.rst
  c.rst[] <- as.vector(cmat)
  mu.rst <- disaggregate(c.rst, epsilon)/1
  mu.rst[] <- mu.rst[]/diffu.rst[]
  result <- matrix(mu.rst[], ncol = 1)
  result.c <- matrix(c.rst[], ncol = 1)
  
  # count step, to identify the 1st step
  step_count <- 1
  
  for (t in 2:(time*step+1)) {
    
    # cmat <- H1.mat %*% cmat + H2.mat %*% cmat^2
    cmat <- H1.mat %*% cmat + H2.vec * cmat^2
    # cmat <- H1.sp %*% cmat + H2.vec * cmat^2
    cmat <- cmat * (1-stream.vec)
    
    if(step_count==step){
      # obtain mu
      c.rst <- rescaling.rst
      c.rst[] <- as.vector(cmat)
      mu.rst <- disaggregate(c.rst, epsilon)/1
      # mu.rst[] <- ifelse(diffu.rst[]==0, 0, mu.rst[]/diffu.rst[])
      # mu.rst[] <- mu.rst[]/diffu.rst[]*(1-stream.vec)
      mu.rst[] <- mu.rst[]/diffu.rst[]
      result <- cbind(result, mu.rst[])
      result.c <- cbind(result.c, c.rst[])
      step_count <- 1
    } else {
      step_count <- step_count+1
    }
  }
  
  # if(sum(result[,-c(1,2,3,4)] == min.value) > 0){
  #   warning("There is value equal to min value.")
  # }
  # 
  # if(sum(result[,-c(1,2,3,4)] < min.value) > 0){
  #   warning("There is value less than min value.")
  # }
  # 
  # if(sum(result[,-c(1,2,3,4)] == max.value) > 0){
  #   warning("There is value equal to max value.")
  # }
  # 
  # if(sum(result[,-c(1,2,3,4)] > max.value) > 0){
  #   warning("There is value greater than max value.")
  # }
  
  # return(list(result.mu = result, result.c = result.c, 
  #             diffu.bar = diffu.bar[], 
  #             growth.bar = growth.bar[],
  #             growth2.bar =growth2.bar[],
  #             H1 = H1.mat,
  #             H2 = H2.mat))
  
  return(list(result.mu = result, result.c = result.c, 
              diffu.bar = diffu.bar[], 
              growth.bar = growth.bar[],
              growth2.bar =growth2.bar[],
              H1 = H1.mat,
              H2 = H2.vec))
  
}


# my numerical solution of PDE (with diffusion and logistic growth) 
# $\frac{\partial}{\partial t}\mu=\frac{\partial^2}{\partial \bs^2}(D\mu)+(G\mu)(1-\mu)$
mynumericalPDEfun <- function(n, initial.vec, stream.vec,
                              diffu.vec, growth.vec, 
                              length_of_grid, time, step, delta_t, 
                              min.value = 0, max.value = 1){
  
  # diffusion rate
  diffu.rc <- diffu.vec
  
  # growth rate
  growth.rc <- growth.vec
  
  H1.mat <- propagator.FP.1st(neighbor = neighbormatrix, 
                              diffu.rate = diffu.rc, 
                              growth.rate = growth.rc,
                              dx = length_of_grid, dy = length_of_grid, dt = delta_t)
  
  # H1.mat <- propagator.FP.1st(neighbor = neighbormatrix, 
  #                             diffu.rate = diffu.rc, 
  #                             growth.rate = growth.rc,
  #                             dx = 1/n, dy = 1/n, dt = delta_t)
  
  H2.mat <- propagator.FP.2nd(neighbor = neighbormatrix, 
                              diffu.rate = diffu.rc, 
                              growth.rate = growth.rc,
                              dx = length_of_grid, dy = length_of_grid, dt = delta_t)
  
  # H2.mat <- propagator.FP.2nd(neighbor = neighbormatrix, 
  #                             diffu.rate = diffu.rc, 
  #                             growth.rate = growth.rc,
  #                             dx = 1/n, dy = 1/n, dt = delta_t)
  
  # initial value for PDE
  mu0 <- initial.vec
  mumat <- as.matrix(mu0)
  
  # set up result 
  result <- matrix(mu0, ncol = 1)
  
  # count step, to identify the 1st step
  step_count <- 1
  
  for (t in 2:(time*step+1)) {
    
    mumat <- H1.mat %*% mumat + H2.mat %*% mumat^2
    mumat <- mumat * (1-stream.vec)
    
    if(step_count==step){
      # obtain mu
      mu <- as.vector(mumat)
      # mu <- mu * (1-stream.vec)
      result <- cbind(result, mu)
      step_count <- 1
    } else {
      step_count <- step_count+1
    }
  }
  
  # if(sum(result[,-c(1,2,3,4)] == min.value) > 0){
  #   warning("There is value equal to min value.")
  # }
  # 
  # if(sum(result[,-c(1,2,3,4)] < min.value) > 0){
  #   warning("There is value less than min value.")
  # }
  # 
  # if(sum(result[,-c(1,2,3,4)] == max.value) > 0){
  #   warning("There is value equal to max value.")
  # }
  # 
  # if(sum(result[,-c(1,2,3,4)] > max.value) > 0){
  #   warning("There is value greater than max value.")
  # }
  
  return(list(result.mu = result,
              diffu = diffu.rc, 
              growth = growth.rc,
              H1 = H1.mat,
              H2 = H2.mat))
}

# my numerical solution of PDE (with diffusion and exponential growth) 
# $\frac{\partial}{\partial t}\mu=\frac{\partial^2}{\partial \bs^2}(D\mu)+G\mu$
mynumericalPDEfun_expgr <- function(n, initial.vec, stream.vec,
                                    diffu.vec, growth.vec, 
                                    length_of_grid, time, step, delta_t, 
                                    min.value = 0, max.value = 1){
  
  # diffusion rate
  diffu.rc <- diffu.vec
  
  # growth rate
  growth.rc <- growth.vec
  
  H1.mat <- propagator.FP.1st(neighbor = neighbormatrix, 
                              diffu.rate = diffu.rc, 
                              growth.rate = growth.rc,
                              dx = length_of_grid, dy = length_of_grid, dt = delta_t)
  
  # H1.mat <- propagator.FP.1st(neighbor = neighbormatrix, 
  #                             diffu.rate = diffu.rc, 
  #                             growth.rate = growth.rc,
  #                             dx = 1/n, dy = 1/n, dt = delta_t)
  
  # initial value for PDE
  mu0 <- initial.vec
  mumat <- as.matrix(mu0)
  
  # set up result 
  result <- matrix(mu0, ncol = 1)
  
  # count step, to identify the 1st step
  step_count <- 1
  
  for (t in 2:(time*step+1)) {
    
    mumat <- H1.mat %*% mumat
    mumat <- mumat * (1-stream.vec)
    
    if(step_count==step){
      # obtain mu
      mu <- as.vector(mumat)
      # mu <- mu * (1-stream.vec)
      result <- cbind(result, mu)
      step_count <- 1
    } else {
      step_count <- step_count+1
    }
  }
  
  # if(sum(result[,-c(1,2,3,4)] == min.value) > 0){
  #   warning("There is value equal to min value.")
  # }
  # 
  # if(sum(result[,-c(1,2,3,4)] < min.value) > 0){
  #   warning("There is value less than min value.")
  # }
  # 
  # if(sum(result[,-c(1,2,3,4)] == max.value) > 0){
  #   warning("There is value equal to max value.")
  # }
  # 
  # if(sum(result[,-c(1,2,3,4)] > max.value) > 0){
  #   warning("There is value greater than max value.")
  # }
  
  return(list(result.mu = result,
              diffu = diffu.rc, 
              growth = growth.rc,
              H1 = H1.mat))
}

# my numerical solution with homogenization function with Rcpp matrix (with exponential growth)
mynumericalPDEhomofun_expgr.Rcpp_mat <- function(n, epsilon, initial.vec, stream.vec,
                                                 diffu.vec, growth.vec, 
                                                 length_of_grid, time, step, delta_t, 
                                                 min.value = 0, max.value = 1){
  # m <- n/epsilon
  # 
  # index.n <- c(1:(n*n))
  # s1.n <- rep(c(1:n), each = n)
  # s2.n <- rep(c(1:n), n)
  # 
  # index.m <- c(1:(m*m))
  # s1.m <- rep(c(1:m), each = m)
  # s2.m <- rep(c(1:m), m)
  
  # original raster and rescaling raster
  # onematrix.NbyN <- matrix(rep(1, n*n), nrow = n, ncol = n)
  # origin.rst <- raster(onematrix.NbyN, xmn = 0, xmx = 1, ymn = 0, ymx = 1, 
  #                      crs = NA)
  # rescaling.rst <- aggregate(origin.rst, fact = epsilon, na.rm = TRUE, 
  #                            fun = myharmonicmeanfun)
  # 
  # neighbormatrix <- neighborhood(rescaling.rst)
  
  diffu.rst <- origin.rst
  diffu.rst[] <- diffu.vec
  
  diffu.bar <- aggregate(diffu.rst, fact = epsilon, na.rm = TRUE, 
                         fun = myharmonicmeanfun)
  
  # growth rate and homogenization
  growth.rst <- origin.rst
  growth.rst[] <- growth.vec
  
  growth.bar <- diffu.bar * aggregate(myrstdividefun(growth.rst, diffu.rst), 
                                      fact = epsilon, 
                                      na.rm = TRUE, 
                                      fun = mymeanfun)
  
  growth2.bar <- diffu.bar * aggregate(myrstdivide2fun(growth.rst, diffu.rst), 
                                       fact = epsilon, 
                                       na.rm = TRUE, 
                                       fun = mymeanfun)
  
  # initial value for PDE and its rescaling
  mu0 <- initial.vec
  mu.rst <- origin.rst
  mu.rst[] <- mu0
  
  c0 <- rescaling.rst
  c0 <- aggregate(mu.rst * diffu.rst, fact = epsilon, na.rm = TRUE, fun = mymeanfun)
  # c0[] <- raster::extract(mu.rst * diffu.rst, SpatialPoints(rescaling.rst))
  
  # will be the same if set cmat <- c0[], since matrix multiplication will convert it to matrix
  cmat <- as.matrix(c0[])
  
  # rate coef for the function
  diffu.rc <- diffu.bar[]
  growth.rc <- growth.bar[]
  growth2.rc <- growth2.bar[]
  
  H1.mat <- propagator_1st_Rcpp(neighbor = neighbormatrix, 
                                diffu_rate = diffu.rc, 
                                growth_rate = growth.rc, 
                                growth_rate2 = growth2.rc,
                                dx = epsilon*length_of_grid, 
                                dy = epsilon*length_of_grid, 
                                dt = delta_t)
  
  # H1.mat <- propagator_1st_Rcpp(neighbor = neighbormatrix, 
  #                               diffu_rate = diffu.rc, 
  #                               growth_rate = growth.rc, 
  #                               growth_rate2 = growth2.rc,
  #                               dx = epsilon/n, dy = epsilon/n, dt = delta_t)
  
  # H1.sp <- Matrix(H1.mat, sparse = TRUE)
  
  # H2.mat <- propagator.2nd(neighbor = neighbormatrix, diffu.rate = diffu.rc,
  #                          growth.rate = growth.rc, growth.rate2 = growth2.rc,
  #                          dx = epsilon/n, dy = epsilon/n, dt = delta_t)
  
  # set up result matrix
  cmat <- cmat * (1-stream.vec)
  c.rst <- rescaling.rst
  c.rst[] <- as.vector(cmat)
  mu.rst <- disaggregate(c.rst, epsilon)/1
  mu.rst[] <- mu.rst[]/diffu.rst[]
  result <- matrix(mu.rst[], ncol = 1)
  result.c <- matrix(c.rst[], ncol = 1)
  
  # count step, to identify the 1st step
  step_count <- 1
  
  for (t in 2:(time*step+1)) {
    
    # cmat <- H1.mat %*% cmat + H2.mat %*% cmat^2
    cmat <- H1.mat %*% cmat
    # cmat <- H1.sp %*% cmat + H2.vec * cmat^2
    cmat <- cmat * (1-stream.vec)
    
    if(step_count==step){
      # obtain mu
      c.rst <- rescaling.rst
      c.rst[] <- as.vector(cmat)
      mu.rst <- disaggregate(c.rst, epsilon)/1
      # mu.rst[] <- ifelse(diffu.rst[]==0, 0, mu.rst[]/diffu.rst[])
      # mu.rst[] <- mu.rst[]/diffu.rst[]*(1-stream.vec)
      mu.rst[] <- mu.rst[]/diffu.rst[]
      result <- cbind(result, mu.rst[])
      result.c <- cbind(result.c, c.rst[])
      step_count <- 1
    } else {
      step_count <- step_count+1
    }
  }
  
  # if(sum(result[,-c(1,2,3,4)] == min.value) > 0){
  #   warning("There is value equal to min value.")
  # }
  # 
  # if(sum(result[,-c(1,2,3,4)] < min.value) > 0){
  #   warning("There is value less than min value.")
  # }
  # 
  # if(sum(result[,-c(1,2,3,4)] == max.value) > 0){
  #   warning("There is value equal to max value.")
  # }
  # 
  # if(sum(result[,-c(1,2,3,4)] > max.value) > 0){
  #   warning("There is value greater than max value.")
  # }
  
  # return(list(result.mu = result, result.c = result.c, 
  #             diffu.bar = diffu.bar[], 
  #             growth.bar = growth.bar[],
  #             growth2.bar =growth2.bar[],
  #             H1 = H1.mat,
  #             H2 = H2.mat))
  
  return(list(result.mu = result, result.c = result.c, 
              diffu.bar = diffu.bar[], 
              growth.bar = growth.bar[],
              H1 = H1.mat))
  
}