source('./sder.R')


if(hyp_sin == 1) { 
  source('./ea1_func.R')
  ea_func <- ea1_func 
  drft <- drift1
} else { 
  source('./ea2_func.R')
  ea_func <- ea2_func
  drft <- drift2
}

x0 <- 0; tEnd <- 20; nObs <- 9
sdv <- .5
sample_size = 5000
Prt = 50
step_size=.01
theta=1


set.seed(1)

dat <- data_gen(x0 = 0, t0 = 0, sd = sdv, tEnd = tEnd, nObs = nObs,drift=function(x) drift(x,theta))
tDat <- dat$tDat
xDat <- dat$xDat

# Run HMC sampler (update theta)
time_spend1 <- system.time(rslt_hmc <- sde_hmc(reps = sample_size, x0=x0, t0=0,
                                         tEnd = tEnd, xDat = xDat, tDat = tDat,
                                         func = ea_func,
                                         sdv = sdv,
                                         theta = NULL, M=100, eps=.2, L=5,
                                         ))

# Run Euler pMCMC sampler (update theta)
time_spend2 <- system.time(rslt_eul <- 
                          pmcmc_Eul(reps = sample_size, P=Prt, t0=0, stepSize=0.01, tEnd=tEnd, 
                                    xDat=xDat, tDat=tDat, sdv = sdv, theta=NULL, func=drft))


############
# Summarize results
timeCost_hmc <- time_spend1["elapsed"]
timeCost_eul <- time_spend2["elapsed"]
      
hmc_theta <- rslt_hmc$theta
es_hmc    <- effectiveSize(hmc_theta)
ess_hmc   <- es_hmc / timeCost_hmc
      
es_eul    <- effectiveSize(rslt_eul$theta)
ess_eul   <- es_eul / timeCost_eul
    
Idx       <- seq(1, sample_size, by = 50)
ksRslt    <- ks.test(rslt_hmc$theta[Idx], rslt_eul$theta[Idx])
pValue    <- ksRslt$p.value
ksDist    <- ksRslt$statistic
    
dfm_cmp   <- tibble(p=Prt, T=tEnd, N=nObs, tHMC=timeCost_hmc, esHMC=es_hmc, essHMC=ess_hmc, 
                                           tEul=timeCost_eul, esEul=es_eul, essEul=ess_eul, 
                                           pVal=pValue, ks=ksDist)
dfm_cmp
