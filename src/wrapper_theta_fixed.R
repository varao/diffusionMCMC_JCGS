library('ggplot2')
library('tidyverse')
library('reshape2')
library('coda')


source('./sder.R')

hyp_sin <- 2   # 1 in the sine drift, 2 is the hyperbolic drift

if(hyp_sin == 1) { 
    source('./ea1_func.R')
    ea_func <- ea1_func
    drift <- drift1
 } else { 
    source('./ea2_func.R')
    ea_func <- ea2_func
    drift <- drift2
 }

# Generate synthetic data
set.seed(1)
x0 <- 0; tEnd <- 20; nObs <- 9
sdv <- .5
sample_size = 15000
Prt = 50
step_size=.01

# Theta is fixed to 1
theta=1

dat <- data_gen(x0 = 0, t0 = 0, sd = sdv, tEnd = tEnd, nObs = nObs,drift=function(x) drift(x,theta))
tDat <- dat$tDat
xDat <- dat$xDat


# Run HMC sampler with theta fixed to 1
time_spent1 <- system.time(rslt_hmc <- sde_hmc(reps = sample_size, x0=x0, t0=0,
                                         tEnd = tEnd, xDat = xDat, tDat = tDat,
                                         func = ea_func,
                                         sdv = sdv,
                                         theta = 1, M=100, eps=.2, L=5,
                                         ))
print("time spent is")
print(time_spent1)
timeCost_hmc <- time_spent1["elapsed"]
hmcObsMat <- rslt_hmc$xObsMat
es_hmc <- effectiveSize(t(hmcObsMat))
ess_hmc <- es_hmc / timeCost_hmc

# Run Euler-Maruyama pPMCMC   
time_spent2 <- system.time(rslt_eul <- 
                          pmcmc_Eul(reps = sample_size, P=Prt, t0=0, stepSize=step_size, tEnd=tEnd, xDat=xDat, 
                          tDat=tDat, sdv = sdv, theta=1, func=drift))
rslt_eul <- rslt_eul$xObsMat
print("time spent is")
print(time_spent2)
timeCost_eul <- time_spent2["elapsed"]
es_eul <- effectiveSize(t(rslt_eul))
ess_eul <- es_eul / timeCost_eul

# Summarize results:
if(hyp_sin == 1) { 
   sin_mn <- rbind(xDat, rowMeans(rslt_hmc$xObsMat), rowMeans(rslt_eul))
   sin_mn
} else {
   hyp_mn <- rbind(xDat, rowMeans(rslt_hmc$xObsMat), rowMeans(rslt_eul))
   hyp_mn
}
