## Sub-ordinate functions required to perform Minimum Divergence Estimation (MDE) of parameters of a stochastic system
## Written by: Alexandre Breuzé, David Price, and Olivier Restif*
## Contact: or226@cam.ac.uk
## Date: 21st March, 2017

# Required R packages
library(expm)
library(foreach)
library(doParallel)
library(Rcpp)
require(powell)

# Substitute individual parameter values by name
# - x and rep are named vectors.
replace.par <- function(x,rep)
{
  if(length(rep)==0) return(x)
  pos <- sapply(names(rep),function(n) which(names(x)==n))
  replace(x,pos,rep)
}

# Labels for the 9 moments (expectation and variances)
mom.names <- c('E.NB','E.NL','E.NS','V.NB','V.NL','V.NS','V.NB.NL','V.NB.NS','V.NL.NS')



# Calculate the KL divergence between two, n-dimensional, multivariate normal distributions 
# with mean vectors mu0 and mu1, and cov matrices cov0 and cov1
KL.div <- function(mu0,cov0,mu1,cov1){
  k <- length(mu0)
  inv.C1 <- solve(cov1)
  as.numeric(sum(diag(inv.C1 %*% cov0)) + (mu1-mu0) %*% inv.C1 %*% (mu1-mu0) - k + log(det(cov1)/det(cov0)))/2
}

# Calculate the KL divergence from a predicted distribution to an observed sample, both characterised by their vectors of 9 moments
# - pred: predicted moments
# - obs: observed moments
# - sub: vector specifying which organs to include, e.g. c(1,2) for B and L. Default 1:3
KL.div.M2M <- function(pred,obs,sub=1:3){
  # Predicted moments
  mu0 <- pred[sub]
  cov0 <- matrix(pred[c(4,7,8,7,5,9,8,9,6)],3)[sub,sub]
  # Observed moments
  mu1 <- obs[sub]
  cov1 <- matrix(obs[c(4,7,8,7,5,9,8,9,6)],3)[sub,sub]
  # KL divergence
  KL.div(mu0,cov0,mu1,cov1)
}


# -------------------------- General Version: N Organs, Matrix Exponential Solution to Integral -------------------------------
# A generalised version for calculating moments of a Radial, N organ system.
# Author: Alexandrè Breuze
# Date: July 2016
# INPUTS:
# N: Number of organs
# t: observation time
# par: Parameters of system, in order c("m12","m13","m14",...,"m21","m31","m41",...,"k2","k3",...,"r2","r3",...) : NOTE: Code (currently) does not allow r/k in organ 1.
# M0: Initial state of the system, in order (E.O1, E.O2, ... E.ON, V.O1, V.O2, ..., V.ON, V.O1.O2, V.O1.O3, ..., V.O(N-1).ON), where Oi is organ i.


WITS.moment.sol.N.radial <- function(N,t,par,M0,met="Higham08.b"){with(as.list(par),{
  # First moment: M.1(t) = (E[n1],E[n2],...,E[nN]) is solution of M.1'(t) = A * M.1
  # We take n1=nB
  A <- matrix (0 , N, N)
  A[1,2:N] <- par[N:(2*N-2)]
  A[2:N,1] <- par[1:(N-1)]
  for (i in 1:N-1){
    A[i+1,i+1] <- -par[i+N-1]-par[i+2*N-2]+par[i+3*N-3]
    # A[i+1,i+1]=paste("-",par[i+N-1],"-",par[i+2*N-2],"+",par[i+3*N-3], sep="")
  }
  A[1,1] <- -sum(par[1:(N-1)])
  # A[1,1]= paste("-",paste(par[1:(N-1)], collapse = "-"), sep="")
  
  # Solution
  M.1 <- expAtv(A,M0[1:N],t)$eAtv
  
  
  ## Second moment: M.2(t) = (E(n1?),E(n2?),..,E(nN?),E(n1,n2),..,E(nN-1,nN)) is solution of M.2' = B*M.1 + C*M.2
  B <- matrix(0, N*(N+1)/2, N)
  B[1,2:N] <- par[N:(2*N-2)]
  B[2:N,1] <- par[1:(N-1)]
  B[(N+1):(2*N-1),1] <- -par[1:(N-1)]
  for (i in 1:(N-1)){
    B[i+1,i+1] <- par[i+N-1]+par[i+2*N-2]+par[i+3*N-3]
    B[i+N,i+1] <- -par[i+N-1]
  }
  B[1,1] <- sum(par[1:N-1])
  
  C <- matrix(0, N*(N+1)/2, N*(N+1)/2)
  
  #E(n1 2)
  C[1,(N+1):(2*N-1)] <- 2*par[N:(2*N-2)] 
  
  #E(ni 2)
  for (i in 1:N-1){  
    C[i+1,i+1] <- 2*(-par[i+N-1]-par[i+2*N-2]+par[i+3*N-3]) 
  }
  
  #  Else error replacement length zero
  for (i in 2:N){
    C[i,N+i-1] <- 2*par[i-1]
  }
  
  #E(n1ni)
  C[(N+1):(2*N-1),1] <- par[1:(N-1)]       
  for (i in 1:(N-1)){
    C[i+N,i+1] <- par[i+N-1]
    C[i+N,i+N] <- -par[i+N-1]-par[i+2*N-2]+par[i+3*N-3]-sum(par[1:N-1])
  }
  
  #E(n1n2)
  ################
  C[(N+1),(2*N):(3*N-3)] <- par[(N+1):(2*N-2)]
  
  #E(n1ni)
  m <- N
  for(i in 2:(N-1)){
    n <- 0
    for(j in 2:i){
      C[(N+i),(2*N+i-2+n)] <- par[N+j-2]
      n <- n+N-j-1
    }
    if (i<N-1){
      m <- m+N-i
      C[(N+i),(N+m):(N+m+N-i-2)] <- par[(N+i):(2*N-2)]
    }
  }
  
  #E(ninj) 
  i <- 2
  m <- 0
  l <- 0
  while(i<N)
  {
    for (j in 1:(N-i)){
      C[(2*N-1+j+m),(N+i-1)] <- par[j+l+1]
      C[(2*N-1+j+m),(2*N-1+j+m)] <-  -par[i+N-2]-par[i+2*N-3]+par[i+3*N-4]-par[i+j+N-2]-par[i+j+2*N-3]+par[i+j+3*N-4]
    }
    
    for (j in 2:(N-i+1)){
      C[2*N-2+j+m,N+j+l] <- par[l+1]
    }
    m <- m+N-i
    i <- i+1
    l <- l+1
  }
  
  C[1,1] <- -2*sum(par[1:(N-1)])
  
  
  #### Here, we calculate the moments via the computation of a matrix exponential (cf article 'Computing integrals using the matrix exponential')
  
  D <- matrix(0,N*(N+1)/2,N)
  K <- matrix(0,N*(N+3)/2,N*(N+3)/2)
  
  ## K contains the A,B,C matrices
  K[1:(N*(N+1)/2),1:(N*(N+1)/2)] <- C[1:(N*(N+1)/2),1:(N*(N+1)/2)]
  K[1:(N*(N+1)/2),(N*(N+1)/2+1):(N*(N+3)/2)] <- B[1:(N*(N+1)/2),1:N]
  K[(N*(N+1)/2+1):(N*(N+3)/2),(N*(N+1)/2+1):(N*(N+3)/2)] <- A[1:N,1:N]
  
  ## We take the exponential of K 
  M <- expm(t*K,m=met)
  
  ## D represents the matrix solution of the second moment equation
  D[1:(N*(N+1)/2),1:N] <- M[1:(N*(N+1)/2),(N*(N+1)/2+1):(N*(N+3)/2)]
  
  ## Solution
  M.2 <- D %*% M0[1:N]
  M.2 <- as.numeric(expm(t*C,m=met) %*% M0[(N+1):(N*(N+3)/2)] + M.2)
  
  ## Solutions
  return(c(M.1,M.2))
})}


