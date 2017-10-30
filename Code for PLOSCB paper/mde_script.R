## Example script to perform Minimum Divergence Estimation (MDE) of parameters of a stochastic system
##
## The example given below reproduces the results of the naive treatment group, in "An Efficient Moments-Based Inference Method for Within-Host Bacterial Infection Dynamics", by Price et al. (2017)
## Sub-ordinate functions are in file: mde_funcs.R
##
## Written by: Alexandre Breuzé, David Price, and Olivier Restif*
## Contact: or226@cam.ac.uk
## Date: 21st March, 2017


source("mde_funcs.R")

# Moments at 0.5 and 6 hours corresponding to the Naïve group, adjusted for observation noise using simulation-based approach as detailed in the manuscript
obs.mom.1 <- c(6.2977113, 12.6324056, 4.3383687, 13.8912109, 27.0637181, 4.4629803, -0.2500392, 0.7818082, 1.9076971)
obs.mom.2 <- c(0.6360223, 4.3245888, 5.6004023, 1.8920138, 15.4520982, 21.7673993, 0.1313821, 0.9310497, 5.1778312)


# The parameters listed in "par.values" will be updated, the remainder will stay as defined in par.to.est. 
# For simplicity, those being estimated are set to 1. Set e_L=e_S=0 as these are fixed at 0 throughout.
# The starting points for each of the remaining parameters are randomly generated in [0,1]

par.estimates <- vector("numeric", length=7)
par.values <- c(cL=runif(1),cS=runif(1),kL=runif(1),kS=runif(1),rL=runif(1),rS=runif(1))
par.to.est <- c(cL=1,cS=1,eL=0,eS=0,kL=1,kS=1,rL=1,rS=1)

set.seed(1)


optim.powell.est <- NULL
# Estimate inoculum
init <- c(xi=0.25,par.values)

# Assume no inoculum lost
# init <- c(par.values)


optim.powell.est <- powell(init, function(par.ext){
  if(min(par.ext)<0) return(1E100)
  names(par.ext) <- c("x.i",names(par.values))
  par <- par.ext[-1]
  # Update parameters being estimated
  par.i <- replace.par(par.to.est,par)
  # Estimate proportion of inoculum lost. Otherwise, comment out this line
  inoc <- par.ext[1]*200   ## * 200 to change scale
  # Assume that inoculum is poisson distributed -- so mean/var in first compartment are equal
  init.mom.i <- c(inoc,0,0, inoc,0,0, 0,0,0)
  # Evalute moments at 0.5 and 6 hours for 3 compartment model with parameters and initial conditions
  mom.i.1 <- WITS.moment.sol.N.radial(3,0.5,par.i,init.mom.i)
  mom.i.2 <- WITS.moment.sol.N.radial(3,6,par.i,init.mom.i)
  # Evaluate divergence between observed moments, and those predicted for par.i and mom.i
  # Note that the blood (1st compartment) is not used at t=6hrs
  div.1 <- KL.div.M2M(mom.i.1,obs.mom.1,1:3) 
  div.2 <- KL.div.M2M(mom.i.2,obs.mom.2,2:3)
  # Objective: minimise the sum of divergences
  div.1+div.2
},control=list(rhoend=1E-5, maxit=1e+07))


par.estimates <- optim.powell.est$par
obs.div <-  optim.powell.est$value
par.estimates[1] <- par.estimates[1] *200
# Estimated parameters and observed divergence
c(par.estimates, obs.div)

