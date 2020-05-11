#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#  
#   This is the fixed version of before I remove all Qi's functions

hmc = function(U, grad_U, epsilon = 0.1, L = 10, M = 1, current_q, path_status, func)
{
  if (length(current_q) == 0) {
    stop("Please provide current q for HMC.")
  }

  q = current_q
  p = rnorm(length(q), 0, sd = sqrt(M))  # independent standard normal variates
  current_p = p
  
  # Make a half step for momentum at the beginning
  
  p = p - epsilon * grad_U(q, path_status, func) / 2
  
  # Alternate full steps for position and momentum
  
  for (i in 1:L)
  {
    # Make a full step for the position
    
    q = q + epsilon * p / M
    
    # Make a full step for the momentum, except at end of trajectory
    
    if (i != L) p = p - epsilon * grad_U(q, path_status, func)
  }
  
  # Make a half step for momentum at the end.
  
  p = p - epsilon * grad_U(q, path_status, func) / 2
  
  if (any(is.na(p)) | any(is.na(q))) {
    #print("NA produced")
    return(current_q)
  }
  
  if (any(abs(p) > 1E50) | any(abs(q) > 1E50)) {
    #print("The component of p or q is too large")
    return(current_q)
  }
  
  
  # Negate momentum at end of trajectory to make the proposal symmetric
  
  p = -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  
  current_U = U(current_q, path_status, func)
  current_K = sum(current_p^2) / (2 * M)
  proposed_U = U(q, path_status, func)
  proposed_K = sum(p^2) / (2 * M)
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position

  acc <- (current_U-proposed_U+current_K-proposed_K)

  # Occasionally, because of very nearby Poisson events, the leapfrog dynamics for a finite epsilon will produce weird results.
  # reject those too
  if(is.nan(acc) || log(runif(1)) > acc)
  {
    return (current_q)  # reject
  }
  else
  {
    return (q)  # accept
  }
}

poisson_process <- function(t0, tEnd, lambda) {
  # Homogeneous Possion Process
  
  duration <- tEnd - t0
  
  nSample <- rpois(1, lambda * duration)
  
  return(sort(duration * runif(nSample)))
  
}

brownian_motion <- function(tVec, t0, tEnd, x0) {
  # Standardard Brownian Motion with mean 0 and variance 1
  
  t <- c(t0, tVec, tEnd)
  
  # Check whether t is ascending
  if (is.unsorted(t)) t <- sort(t)
  
  m <- length(t)
  increment <- sqrt(diff(t)) * rnorm(m - 1)
  W <- cumsum(c(x0, increment))
  return(W)
}

brownian_bridge_mod <- function(tVec, t0,  tEnd, x0, xEnd, step_size = .01) {
  # Standardard Brownian Bridge
  # t is discreate time point
  # start and end parameter correspond to the first and last time point in t
  
  tLen <- length(tVec)
  if (tLen == 0) {
    return(NULL)
  }
  
  time  <- c(t0, tVec, tEnd)
  
  flag <- 0
  
  if (max(abs(diff(time))) > step_size) {
    # time grid is too coarse, need to be implemented in a finer grid
    
    flag <- 1
    tIdx <- rep(0, length(tVec))
    
    if (t0 == tVec[1]) {
      timeSeg <- c(t0, tVec[1])
    } else {
      timeSeg <- seq(t0, tVec[1], by = step_size)
    }
    
    tIdx[1] <- length(timeSeg)
    
    if (tLen >= 2) {
      for(i in 2:tLen) {
        if (tVec[i-1] == tVec[i]) {
          stop("time overlaps")
        }
        timeSeg <- c(timeSeg, seq(tVec[i - 1], tVec[i], by = step_size)[-1])
        tIdx[i] <- length(timeSeg)
      }
    }
    
    if (tEnd == tVec[tLen]) {
      timeSeg <- c(timeSeg, tEnd)
    } else {
      timeSeg <- c(timeSeg, seq(tVec[tLen], tEnd, by = step_size)[-1])
    }
  } else {
    timeSeg = time
  }
  
  
  m <- length(timeSeg)
  W <- brownian_motion(tVec = timeSeg[-c(1, m)], t0 = t0, tEnd = tEnd, x0 = 0)
  B <- W  + x0 + (timeSeg - timeSeg[1]) / (timeSeg[m] - timeSeg[1]) * (xEnd - x0 - W[m])
  
  if (flag) {
    return(B[tIdx])
  } else {
    return(B[-c(1, m)])
  }
}

brownian_bridge <- function(tVec, t0,  tEnd, x0, xEnd) {
  # Standardard Brownian Bridge
  # t is discreate time point
  # start and end parameter correspond to the first and last time point in t
  
  if (length(tVec) == 0) {
    return(NULL)
  }
  
  t <- c(t0, tVec, tEnd)
  
  m <- length(t)
  W <- brownian_motion(tVec = tVec, t0 = t0, tEnd = tEnd, x0 = 0)
  B <- W  + x0 + (t - t[1]) / (t[m] - t[1]) * (xEnd - x0 - W[m])
  
  return(B[-c(1, m)])
}


brownian_bridge_impute <-
  function(tVec, t0, tVecNew, tEnd, x0, xVec, xEnd) {
    
    if (length(tVecNew) == 0) {
      return(xVecNew = NULL)
    }
    
    if (length(tVec) == 0) {
      return(xVecNew = brownian_bridge(tVec = tVecNew, t0 = t0, tEnd = tEnd, x0 = x0, xEnd = xEnd))
    }
    
    xVecNew <- NULL
    
    tLen <- length(tVec)
    
    for (i in 1:tLen) {
      if (i == 1) {
        xVecNew <- brownian_bridge(tVec = tVecNew[tVecNew < tVec[i]], t0 = t0,
                                   tEnd = tVec[i], x0 = x0, xEnd = xVec[i])
      } else {
        xVecNew <- c(xVecNew,
                     brownian_bridge(tVec = tVecNew[tVecNew >= tVec[i - 1] & tVecNew < tVec[i]],
                                     t0 = tVec[i - 1],
                                     tEnd = tVec[i], x0 = xVec[i - 1],
                                     xEnd = xVec[i]))
      }
    }
    
    xVecNew <- c(xVecNew,
                 brownian_bridge(tVec = tVecNew[tVecNew >= tVec[tLen]],
                                 t0 = tVec[tLen], tEnd = tEnd, x0 = xVec[i], xEnd = xEnd))
    
    if (length(xVecNew) != length(tVecNew)) {
      stop("The length of xVecNew doesn't agree with the length of tVecNew")
    }
    
    return(xVecNew)
  }

get_idx <- function(tVec, tDat) {
  # if (length(tVec) == 0 | length(tDat) == 0) {
  #   print("tVec or tDat is NULL")
  # } 
  #tVecLen <- length(tVec)
  tComb <- c(tVec, tDat)
  tCombSort <- order(tComb)
    
  return(tCombSort)
}


all_lglk <- function(t0, tVec, tEnd, x0, xVec, xEnd,
                     xDat, tDat, xObs, sdv, idx, func) {
  # Calculate log-likelihood the whole path
  # regardless of tVec, xVec is NULL or not, these terms always has its contribution
  tLen <- length(tVec)
  lglk <- func$H1_lglk(xEnd = xEnd, tEnd = tEnd, x0 = x0, t0 = t0)

  Mpoisson = func$Mpoisson
  sdv = func$sdv

  # contribution of the poisson process term
  lglk <- lglk + tLen * log(Mpoisson) - Mpoisson*(tEnd-t0)  #+ func$poisson_lglk(tVec = tVec, Mpoisson = Mpoisson, t0 = t0, tEnd = tEnd)
  
  if (length(xDat) == 0) {
    # No observations
    if (tLen > 0) {
      # contribution of brownian motion
      lglk <- lglk + func$BB_lglk(tVec = tVec, xVec = xVec, t0 = t0, tEnd = tEnd, x0 = x0, xEnd = xEnd)
    
      # contribution of the path likelihood term
      lglk <- lglk + func$path_lglk(xVec = xVec)
    }
    return(lglk = lglk)
  } else {
    # with observation
    if (tLen == 0) {
      lglk <- lglk + func$BB_lglk(tVec = tDat, xVec = xObs, t0 = t0, tEnd = tEnd, x0 = x0, xEnd = xEnd)
      lglk <- lglk + func$obs_lglk(xObs = xObs, xDat = xDat, .sd = sdv)
      return(lglk)
    } else {
      lglk <- lglk + func$BB_lglk(tVec = c(tVec, tDat)[idx], xVec = c(xVec, xObs)[idx],
                                  t0 = t0, tEnd = tEnd, x0 = x0, xEnd = xEnd)
      lglk <- lglk + func$path_lglk(xVec = xVec)
      lglk <- lglk + func$obs_lglk(xObs = xObs, xDat = xDat, .sd = sdv)
      return(lglk)
    }
  }
}

all_grad <- function(t0, tVec, tEnd, x0, xVec, xEnd, xDat, tDat, xObs,
                     sdv, idx, func) {
  
  tLen <- length(tVec)
  
  if (length(xDat) == 0) {
    grad <- rep(0, tLen + 1)
    
    # The contribution just from h1
    grad[tLen + 1] <- func$H1_grad(xEnd = xEnd, tEnd = tEnd, x0 = x0, t0 = t0)
    
    if (tLen > 0) {
      
      # contribution of brownian motion
      grad <- grad + func$BB_grad(tVec = tVec, xVec = xVec, t0 = t0, tEnd = tEnd, x0 = x0, xEnd = xEnd)
      
      # contribution of the path likelihood term
      grad[1:tLen] <- grad[1:tLen] + func$path_grad(xVec = xVec)
    }
    return(grad = grad)

  } else {
    datLen <- length(xDat)
    
    if (tLen == 0) {
      grad <- rep(0, datLen + 1)
      grad[datLen + 1] <- func$H1_grad(xEnd = xEnd, tEnd = tEnd, x0 = x0, t0 = t0)
      grad <- grad +  func$BB_grad(tVec = tDat, xVec = xObs, t0 = t0, tEnd = tEnd, x0 = x0, xEnd = xEnd)
      grad[1:datLen] <- grad[1:datLen] + func$obs_grad(xObs = xObs, xDat = xDat)
    } else {

      combLen <- tLen + datLen
      grad <- rep(0, combLen  + 1)
      grad[c(idx, combLen+1)] <- func$BB_grad(tVec = c(tVec, tDat)[idx], xVec = c(xVec, xObs)[idx],
                                  t0 = t0, tEnd = tEnd, x0 = x0, xEnd = xEnd)
      
      grad[combLen + 1] <- grad[combLen+1] + func$H1_grad(xEnd = xEnd, tEnd = tEnd, x0 = x0, t0 = t0)
      
      grad[1:tLen] <- grad[1:tLen] + func$path_grad(xVec = xVec)
      grad[(tLen+1):combLen] <- grad[(tLen+1):combLen] + func$obs_grad(xObs = xObs, xDat = xDat)
    }
    return(grad = grad)
  }
}

U <- function(current_q, path_status, func) {
  qLen <- length(current_q)
  if (length(path_status$xDat) == 0) {
    return(-all_lglk(t0 = path_status$t0,
                     tVec = path_status$tVec,
                     tEnd = path_status$tEnd,
                     x0 = path_status$x0,
                     xVec = current_q[-qLen],
                     xDat = NULL,
                     xEnd = current_q[qLen],
                     func = func))
  } else {
    tLen <- length(path_status$tVec)
    if (tLen == 0) {
      return(-all_lglk(t0 = path_status$t0,
                       tVec = path_status$tVec,
                       tEnd = path_status$tEnd,
                       x0 = path_status$x0,
                       xVec = NULL,
                       xEnd = current_q[qLen],
                       xDat = path_status$xDat,
                       tDat = path_status$tDat,
                       xObs = current_q[-qLen],
                       sdv = path_status$sdv,
                       idx = NULL,
                       func = func))
    } else {
      return(-all_lglk(t0 = path_status$t0,
                       tVec = path_status$tVec,
                       tEnd = path_status$tEnd,
                       x0 = path_status$x0,
                       xVec = current_q[1:tLen],
                       xEnd = current_q[qLen],
                       xDat = path_status$xDat,
                       tDat = path_status$tDat,
                       xObs = current_q[-qLen][-(1:tLen)],
                       sdv = path_status$sdv,
                       idx = path_status$idx,
                       func = func))
    }
  }
}

grad_U <- function(current_q, path_status, func = func) {
  qLen <- length(current_q)
  if (length(path_status$xDat) == 0) {
    return(-all_grad(t0 = path_status$t0,
                     tVec = path_status$tVec,
                     tEnd = path_status$tEnd,
                     x0 = path_status$x0,
                     xVec = current_q[-qLen],
                     xEnd = current_q[qLen],
                     xDat = NULL,
                     func = func))
  } else {
    tLen <- length(path_status$tVec)
    if (tLen == 0) {
      return(-all_grad(t0 = path_status$t0,
                       tVec = path_status$tVec,
                       tEnd = path_status$tEnd,
                       x0 = path_status$x0,
                       xVec = NULL,
                       xEnd = current_q[qLen],
                       xDat = path_status$xDat,
                       tDat = path_status$tDat,
                       xObs = current_q[-qLen],
                       sdv = path_status$sdv,
                       idx = NULL,
                       func = func))
    } else {
      return(-all_grad(t0 = path_status$t0,
                       tVec = path_status$tVec,
                       tEnd = path_status$tEnd,
                       x0 = path_status$x0,
                       xVec = current_q[1:tLen],
                       xEnd = current_q[qLen],
                       xDat = path_status$xDat,
                       tDat = path_status$tDat,
                       xObs = current_q[-qLen][-(1:tLen)],
                       sdv = path_status$sdv,
                       idx = path_status$idx,
                       func = func))
    }
  }
}

# This function updates both the path and the Poisson grid
path_hmc <- function(t0, tVec, tEnd, x0, xVec, xEnd,
                     xDat = NULL, tDat = NULL, xObs = NULL,
                     epsilon, L, M, func) {

  Mpoisson = func$Mpoisson
  sdv = func$sdv

  path_status <- list(t0 = t0,
                      tVec = tVec,
                      tEnd = tEnd,
                      x0 = x0,
                      xDat = xDat,
                      tDat = tDat,
                      xObs = xObs,
                      Mpoisson = Mpoisson,
                      sdv = sdv,
                      idx = NULL)
  
  tLen <- length(tVec)
  dLen <- length(tDat)
  
  # without observations and Poisson grid
  if (dLen == 0 && tLen == 0) {
      xEnd <- func$bias_sample(tEnd = tEnd, t0 = t0, x0 = x0)
  } else {
    
      path_status$idx <- get_idx(tVec = tVec, tDat = tDat)
      pathNew <- hmc(U, grad_U, epsilon = epsilon, L = L, M = M,
                     current_q = c(xVec, xObs, xEnd), path_status = path_status, func = func)
      
      newLen <- length(pathNew)
      
      xEnd <- pathNew[newLen]

      xVec <- xObs <- NULL
      if(tLen > 0) xVec <- pathNew[1:tLen]
      if(dLen > 0) xObs <- pathNew[(tLen+1):(newLen-1)]
  }
    
    
  tVecNew <- poisson_process(t0 = t0, tEnd = tEnd, lambda = Mpoisson)
    
  if (length(tVecNew) == 0) {
    xVec <- NULL
    tVec <- NULL
  } else {
    
    tComb <- c(tVec, tDat)
    idx <- order(tComb)
    xVecNew <- brownian_bridge_impute(tVec = tComb[idx], t0 = t0,  tVecNew = tVecNew, tEnd = tEnd,
                                      x0 = x0, xVec = pathNew[-newLen][idx], xEnd = xEnd)
    
    poissonIdx <- runif(length(xVecNew)) > func$phi(xVecNew) / Mpoisson
    
    xVec <- xVecNew[poissonIdx]
    tVec <- tVecNew[poissonIdx]
  }
  return(list(tVec = tVec, xVec = xVec, xObs = xObs, xEnd = xEnd))
}


theta_lglk <- function(theta,
                       x0, t0,
                       xEnd, tEnd,
                       xVec, tVec,
                       xObs, tObs,
                       sdv,
                       func_gen) {
  
  
  
  func <- func_gen(const = 1, theta = theta, sdv = sdv)
  Mpoisson <- func$Mpoisson
  L <- func$L
  A_func <- func$A_func
  
  lglk <- A_func(xEnd) - A_func(x0)   

  #H1 <- function(x) exp(func$H1_lglk(xEnd = x, tEnd = tEnd, x0 = x0, t0 = t0))

  # contribution of the poisson process term
  lglk <- lglk - (L/2 + Mpoisson) * (tEnd - t0)
  
  # contribution of the path likelihood term  (i.e. the Poisson thinning)
  lglk <- lglk + func$path_lglk(xVec = xVec)
  
  # contribution of the path instantiated at observations
  del_t = diff(tObs)
  del_x = diff(xObs)
  lglk <- lglk + sum(dnorm(del_x, 0, del_t, log=T))

  return(lglk = lglk)
  
}

new_theta <- function(theta,
                      x0, t0,
                      xEnd, tEnd,
                      xVec, tVec,
                      xObs, tObs,
                      sdv,
                      func) {
  
  theta_new <- func(1,theta,sdv)$MHprop()   # Propose from the prior
                       
  #  includes prior
  diff_lglk <- theta_lglk(theta_new, x0 = x0, t0 = t0, xEnd = xEnd, tEnd = tEnd, xVec = xVec, tVec = tVec, xObs=xObs, tObs=tObs, sdv = sdv, func = func) -
    theta_lglk(theta, x0 = x0, t0 = t0, xEnd = xEnd, tEnd = tEnd, xVec = xVec, tVec = tVec, xObs=xObs, tObs=tObs, sdv = sdv, func = func) 
  
  if (log(runif(1)) < diff_lglk) {
    return(theta_new)
  } else {
    return(theta)
  }
}

## Vinayak: Use this (cleaned up sde_hmc_old)
sde_hmc <- function(reps = 1000,
                    x0 = 0, t0 = 0, tEnd = 10, xDat = NULL, tDat = NULL,
                    func = ea2_func,
                    sdv = 0.05,
                    epsilon = 0.2, L = 5, M = 100,
                    theta = NULL,
                    updt_data=F, init_scl = 1) {
  
  # initialize theta if needed
  if(is.null(theta)) {
    updt_theta = T 
    theta = rexp(1, rate = 1)
    thetaVec <- rep(NA, reps)
    thetaVec[1] <- theta
  } else {
    updt_theta = F
    thetaVec <- theta
  }
  
  # initialize xEnd and xObs
  if (length(xDat) == 0) {
    # without observation
    xObsMat = NULL 
    xObs    = NULL
    xEnd <- 0
    updt_data <- F
  } else {
    # with observations
    xObsMat = matrix(0,nrow=length(xDat), ncol = reps)
    xObs <- xDat + init_scl * sdv*rnorm(length(xDat))
    xEnd <- xDat[length(xDat)]
  }
  
  phi <- func(theta = theta, sdv = sdv)$phi
  Mpoisson <- (theta^2 + theta) / 2
  tVec <- NULL
  xVec <- NULL
  
  # variables to store variable length sufficient statistics from every MCMC iteration
  curPois <- 1
  maxLen <- 30 * reps
  increment <- reps
  time <- rep(NA, maxLen)
  omega <- rep(NA, maxLen)
  chainIdx <- rep(NA, maxLen)
  pointType <- rep(NA, maxLen)
  x0Vec <- rep(NA, reps)
  xTVec <- rep(NA, reps)

  # Begin MCMC
  for(i in 1:reps) {
    if((i %% 200) == 0) print(i)
    # This function updates both the path and the Poisson grid
    newPathRslt <- path_hmc(t0 = t0, tVec = tVec, tEnd = tEnd,
                            x0 = x0, xVec = xVec, xEnd = xEnd,
                            xDat = xDat, tDat = tDat, xObs = xObs,
                            epsilon = epsilon, L = L, M = M,
                            func = func(theta = theta, sdv = sdv))
    
    xVec <- newPathRslt$xVec
    tVec <- newPathRslt$tVec
    xEnd <- newPathRslt$xEnd
    xObs <- newPathRslt$xObs
    pathLen <- length(xVec) + length(xObs) + 2
    
    if(updt_theta) {  
      theta <- new_theta(theta = theta,
                       x0 = x0, t0 = t0,
                       sdv = sdv,
                       xEnd = xEnd, tEnd = tEnd,
                       xVec = xVec, tVec = tVec,
                       xObs = xObs, tObs = tDat,
                       func = func)
      thetaVec[i] <- theta
    }
    
    # update xDat (update data for debugging)
    if(updt_data) xDat <- xObs + rnorm(length(xObs), sd = sdv)
    
    
    # Reallocate more space to store MCMC results
    if (curPois > maxLen) {
      time <- c(time, rep(NA, increment))
      omega <- c(omega, rep(NA, increment))
      chainIdx <- c(chainIdx, rep(NA, increment))
      pointType <- c(pointType, rep(NA, increment))
      maxLen <- maxLen + increment
    }
    #########

    
    # Store results from this MCMC iteration
    if (length(xDat) > 0) xObsMat[,i] = xObs
    
    x0Vec[i] <- x0
    xTVec[i] <- xEnd
   
    idx <- get_idx(tVec = tVec, tDat = tDat)

    time[curPois:(curPois + pathLen - 1)] <- c(t0, c(tVec, tDat)[idx], tEnd)
    omega[curPois:(curPois + pathLen - 1)] <- c(x0, c(xVec, xObs)[idx], xEnd)
    chainIdx[curPois:(curPois + pathLen - 1)] <- rep(i,pathLen)
    pointType[curPois:(curPois + pathLen - 1)] <- c(0, as.numeric(idx > length(tVec)), 0)
    
    curPois <- curPois + pathLen
    #########
  }
  rslt <- list(dfm = data.frame(time = time[1:(curPois - 1)],
                                omega = omega[1:(curPois - 1)],
                                pointType = pointType[1:(curPois - 1)],
                                chainIdx = chainIdx[1:(curPois - 1)]),
               x0 = x0Vec, xT = xTVec,
               theta = thetaVec, xObsMat = xObsMat
               )
  
  return(rslt)
}



EA1 <- function(reps = 1000, evalPoint = NULL, x0 = 0, tEnd = 10) {
  
  t0 <- 0
  Mpoisson <- 9 / 8
  mu0 <- 0
  sigma02 <- .1
  
  samplePath <- numeric()
  count <- 0
  nTotal <- 0
  
  evalTrue <- (length(evalPoint) != 0)
  
  
  
  curPois <- 1
  maxLen <- 30 * reps
  increment <- reps
  time <- rep(NA, maxLen)
  omega <- rep(NA, maxLen)
  chainIdx <- rep(NA, maxLen)
  
  
  while (count < reps){
    
    # x0 <- rnorm(1, mean = 0, sd = sqrt(0.1))
    
    nTotal <- nTotal + 1
    
    tVec <- poisson_process(t0, tEnd, Mpoisson)
    
    #x0 <- rnorm(1, mean = mu0, sd = sqrt(sigma02) )
    
    xEnd <- ea1_func(1)$bias_sample(tEnd, t0 = t0, x0 = x0)
    
    xVec <- brownian_bridge(tVec = tVec, t0 = t0, tEnd = tEnd, x0 = x0, xEnd = xEnd)
    
    if (curPois > maxLen) {
      time <- c(time, rep(NA, increment))
      omega <- c(omega, rep(NA, increment))
      chainIdx <- c(chainIdx, rep(NA, increment))
      maxLen <- maxLen + increment
    }
    
    if (length(xVec) == 0) {
      count <- count + 1
      
      if (evalTrue) {
        evalVal <-  brownian_bridge(tVec = evalPoint,
                                    t0 = t0,
                                    tEnd = tEnd, x0 =x0, xEnd = xEnd)
        
        
        pathLen <- 3
        time[curPois:(curPois + pathLen - 1)] <- c(t0, evalPoint, tEnd)
        omega[curPois:(curPois + pathLen - 1)] <- c(x0, evalVal, xEnd)
        chainIdx[curPois:(curPois + pathLen - 1)] <- rep(count, pathLen)
        curPois <- curPois + pathLen
        
        
      } else {
        pathLen <- 2
        time[curPois:(curPois + pathLen - 1)] <- c(t0, tEnd)
        omega[curPois:(curPois + pathLen - 1)] <- c(x0, xEnd)
        chainIdx[curPois:(curPois + pathLen - 1)] <- rep(count, pathLen)
        curPois <- curPois + pathLen
        # samplePath <- rbind(samplePath,
        #                     data.frame(time = c(t0, tEnd),
        #                                  omega = c(x0, xEnd),
        #                                  chainIdx = rep(count, 2)))
      }
      
      next
    }
    
    weight <- ea1_func(1)$phi(xVec) / Mpoisson
    
    if (!any(runif(length(weight)) < weight)) {
      count <- count + 1
      tCur <- c(t0 ,tVec, tEnd)
      xCur <- c(x0, xVec, xEnd)
      
      if (evalTrue) {
        tNew <- sort(c(tCur, evalPoint))
        evalIdx <- which(tNew == evalPoint)
        evalVal <- brownian_bridge(tVec = evalPoint,
                                   t0 = tNew[evalIdx - 1],
                                   tEnd = tNew[evalIdx + 1],
                                   x0 = xCur[evalIdx - 1],
                                   xEnd = xCur[evalIdx])
        xCur <- c(xCur[1:(evalIdx - 1)], evalVal, xCur[evalIdx:length(xCur)])
        tCur <- tNew
      }
      
      
      pathLen <- length(tCur)
      time[curPois:(curPois + pathLen - 1)] <- tCur
      omega[curPois:(curPois + pathLen - 1)] <- xCur
      chainIdx[curPois:(curPois + pathLen - 1)] <- rep(count, pathLen)
      curPois <- curPois + pathLen
    }
  }
  samplePath <- data.frame(time = time[1:(curPois - 1)], omega = omega[1:(curPois - 1)], chainIdx = chainIdx[1:(curPois - 1)])
  print(paste0("The acceptance rate is", count / nTotal))
  return(samplePath)
}

EA2 <- function(reps = 1000, evalPoint = NULL, x0 = 0, tEnd = 10, theta = 1) {
  
  t0 <- 0
  Mpoisson <-  optimize(ea2_func(const = 1, theta = theta)$phi, interval = c(0, tEnd), maximum = TRUE)$object
  
  mu0 <- 0
  sigma02 <- .1
  
  samplePath <- numeric()
  count <- 0
  nTotal <- 0
  
  evalTrue <- (length(evalPoint) != 0)
  
  
  
  curPois <- 1
  maxLen <- 30 * reps
  increment <- reps
  time <- rep(NA, maxLen)
  omega <- rep(NA, maxLen)
  chainIdx <- rep(NA, maxLen)
  
  
  while (count < reps){
    
    # x0 <- rnorm(1, mean = 0, sd = sqrt(0.1))
    
    nTotal <- nTotal + 1
    
    tVec <- poisson_process(t0, tEnd, Mpoisson)
    
    #x0 <- rnorm(1, mean = mu0, sd = sqrt(sigma02) )
    
    xEnd <- ea2_func(const = 1, theta = theta)$bias_sample(tEnd = tEnd, t0 = t0, x0 = x0)
    
    xVec <- brownian_bridge(tVec = tVec, t0 = t0, tEnd = tEnd, x0 = x0, xEnd = xEnd)
    
    if (curPois > maxLen) {
      time <- c(time, rep(NA, increment))
      omega <- c(omega, rep(NA, increment))
      chainIdx <- c(chainIdx, rep(NA, increment))
      maxLen <- maxLen + increment
    }
    
    if (length(xVec) == 0) {
      count <- count + 1
      
      if (evalTrue) {
        evalVal <-  brownian_bridge(tVec = evalPoint,
                                    t0 = t0,
                                    tEnd = tEnd, x0 =x0, xEnd = xEnd)
        
        
        pathLen <- 3
        time[curPois:(curPois + pathLen - 1)] <- c(t0, evalPoint, tEnd)
        omega[curPois:(curPois + pathLen - 1)] <- c(x0, evalVal, xEnd)
        chainIdx[curPois:(curPois + pathLen - 1)] <- rep(count, pathLen)
        curPois <- curPois + pathLen
        
        
      } else {
        pathLen <- 2
        time[curPois:(curPois + pathLen - 1)] <- c(t0, tEnd)
        omega[curPois:(curPois + pathLen - 1)] <- c(x0, xEnd)
        chainIdx[curPois:(curPois + pathLen - 1)] <- rep(count, pathLen)
        curPois <- curPois + pathLen
      }
      
      next
    }
    
    weight <- ea2_func(const = 1, theta = theta)$phi(xVec) / Mpoisson
    
    if (!any(runif(length(weight)) < weight)) {
      count <- count + 1
      tCur <- c(t0 ,tVec, tEnd)
      xCur <- c(x0, xVec, xEnd)
      
      if (evalTrue) {
        tNew <- sort(c(tCur, evalPoint))
        evalIdx <- which(tNew == evalPoint)
        evalVal <- brownian_bridge(tVec = evalPoint,
                                   t0 = tNew[evalIdx - 1],
                                   tEnd = tNew[evalIdx + 1],
                                   x0 = xCur[evalIdx - 1],
                                   xEnd = xCur[evalIdx])
        xCur <- c(xCur[1:(evalIdx - 1)], evalVal, xCur[evalIdx:length(xCur)])
        tCur <- tNew
      }
      
      
      pathLen <- length(tCur)
      time[curPois:(curPois + pathLen - 1)] <- tCur
      omega[curPois:(curPois + pathLen - 1)] <- xCur
      chainIdx[curPois:(curPois + pathLen - 1)] <- rep(count, pathLen)
      curPois <- curPois + pathLen
    }
  }
  samplePath <- data.frame(time = time[1:(curPois - 1)], omega = omega[1:(curPois - 1)], chainIdx = chainIdx[1:(curPois - 1)])
  print(paste0("The acceptance rate is", count / nTotal))
  return(samplePath)
}

ea_path <- function(t0 = 0,
                    tEnd = 10,
                    x0 = 0,
                    xEnd = NULL,
                    Mpoisson = 9/8,
                    maxIter = 50000,
                    func = ea1_func(1)) {
  
  
  if (is.null(xEnd)) {
    # xEnd is not provided
    xEnd <- func$bias_sample(tEnd, t0 = t0, x0 = x0)
  }
  
  for(i in 1:maxIter) {
    
    tVec <- poisson_process(t0, tEnd, Mpoisson)
    
    if (length(tVec) == 0) {
      pathLen <- 2
      return(list(time = c(t0, tEnd), omega = c(x0, xEnd), iterTime = i))
    }
    
    tVec <- tVec + t0
    
    xVec <- brownian_bridge(tVec = tVec, t0 = t0, tEnd = tEnd, x0 = x0, xEnd = xEnd)
    
    weight <- func$phi(xVec) / Mpoisson
    
    if (!any(runif(length(weight)) < weight)) {
      return(list(time = c(t0 ,tVec, tEnd), omega =  c(x0, xVec, xEnd), iterTime = i))
    }
  }
  
  print('max iteration reached')
  return(list(time = c(t0, tEnd), omega = c(x0, xEnd), iterTime = maxIter))
}

ea_seg <- function(reps = 1000, tEnd = 10, x0 = 0, segLen = 5, Mpoisson = 9 / 8, func = ea1_func(1)) {
  
  t0 <- 0
  
  samplePath <- numeric()
  
  
  
  curPois <- 1
  maxLen <- 3 * ceiling(Mpoisson * (tEnd - t0)) * reps
  increment <- reps
  time <- rep(NA, maxLen)
  omega <- rep(NA, maxLen)
  chainIdx <- rep(NA, maxLen)
  iterTime <- 0
  
  if (tEnd >= 2 * segLen) {
    nSeg <- (tEnd - .00001) %/% segLen
    tSeq <- (1:nSeg) * segLen
    tSeqLen <- length(tSeq)
    
    xEnd <- func$bias_sample(tEnd = tEnd, t0 = t0, x0 = x0)
    
    xSeq <- brownian_bridge(tVec = tSeq, t0 = t0, tEnd = tEnd, x0 = x0, xEnd = xEnd)
  }
  
  for (i in 1:reps) {
    
    if (curPois > maxLen) {
      time <- c(time, rep(NA, increment))
      omega <- c(omega, rep(NA, increment))
      chainIdx <- c(chainIdx, rep(NA, increment))
      maxLen <- maxLen + increment
    }
    
    
    if (tEnd < 2 * segLen) {
      
      rslt <- ea_path(t0 = t0,
                      tEnd = tEnd,
                      x0 = x0,
                      xEnd = NULL,
                      Mpoisson = Mpoisson,
                      maxIter = 50000,
                      func = func)
      pathLen <- length(rslt$time)
      time[curPois:(curPois + pathLen - 1)] <- rslt$time
      omega[curPois:(curPois + pathLen - 1)] <- rslt$omega
      iterTime <- iterTime + rslt$iterTime
      chainIdx[curPois:(curPois + pathLen - 1)] <- rep(i, pathLen)
      curPois <- curPois + pathLen
    } else {
      for(j in 1:tSeqLen) {
        if (j == 1) {
          rslt <- ea_path(t0 = t0,
                          tEnd = tSeq[j],
                          x0 = x0,
                          xEnd = xSeq[j],
                          Mpoisson = Mpoisson,
                          maxIter = 50000,
                          func = func)
          timeTmp <- rslt$time
          omegaTmp <- rslt$omega
          iterTimeTmp <- rslt$iterTime
        } else {
          rslt <- ea_path(t0 = tSeq[j - 1],
                          tEnd = tSeq[j],
                          x0 = xSeq[j - 1],
                          xEnd = xSeq[j],
                          Mpoisson = Mpoisson,
                          maxIter = 50000,
                          func = func)
          curPathLen <- length(timeTmp)
          xSeq[j - 1] <- brownian_bridge(tVec = tSeq[j - 1],
                                         t0 = timeTmp[curPathLen - 1],
                                         tEnd = rslt$time[2],
                                         x0 = omegaTmp[curPathLen - 1],
                                         xEnd = rslt$omega[2])
          timeTmp <- c(timeTmp, rslt$time[-1])
          omegaTmp <- c(omegaTmp, rslt$omega[-1])
          iterTimeTmp <- iterTimeTmp + rslt$iterTime
        }
      }
      
      xEnd <- func$bias_sample(tEnd = tEnd, t0 = t0, x0 = x0)
      
      rslt <- ea_path(t0 = tSeq[tSeqLen],
                      tEnd = tEnd,
                      x0 = xSeq[tSeqLen],
                      xEnd = xEnd,
                      Mpoisson = Mpoisson,
                      maxIter = 50000,
                      func = func)
      
      curPathLen <- length(timeTmp)
      xSeq[tSeqLen] <- brownian_bridge(tVec = tSeq[tSeqLen],
                                       t0 = timeTmp[curPathLen - 1],
                                       tEnd = rslt$time[2],
                                       x0 = omegaTmp[curPathLen - 1],
                                       xEnd = rslt$omega[2])
      timeTmp <- c(timeTmp, rslt$time[-1])
      omegaTmp <- c(omegaTmp, rslt$omega[-1])
      iterTimeTmp <- iterTimeTmp + rslt$iterTime
      
      pathLen <- length(timeTmp)
      time[curPois:(curPois + pathLen - 1)] <- timeTmp
      omega[curPois:(curPois + pathLen - 1)] <- omegaTmp
      iterTime <- iterTime + iterTimeTmp
      chainIdx[curPois:(curPois + pathLen - 1)] <- rep(i, pathLen)
      curPois <- curPois + pathLen
      
    }
  }
  return(list(time = time[1:(curPois - 1)],
              omega = omega[1:(curPois - 1)],
              chainIdx = chainIdx[1:(curPois - 1)],
              iterTime = iterTime))
}

ea_seg2 <- function(reps = 1000, tEnd = 10, segLen = 5, x0 = 0, Mpoisson = 9 / 8, func = ea1_func(1)) {
  
  t0 <- 0
  
  samplePath <- numeric()
  
  curPois <- 1
  maxLen <- 3 * ceiling(Mpoisson * (tEnd - t0)) * reps
  increment <- reps
  time <- rep(NA, maxLen)
  omega <- rep(NA, maxLen)
  chainIdx <- rep(NA, maxLen)
  iterTime <- 0
  
  if (tEnd >= 2 * segLen) {
    nSeg <- (tEnd - .00001) %/% segLen
    tSeq <- (1:nSeg) * segLen
    tSeqLen <- length(tSeq)
    
    xEnd <- func$bias_sample(tEnd = tEnd, t0 = t0, x0 = x0)
    
    xSeq <- brownian_bridge(tVec = tSeq, t0 = t0, tEnd = tEnd, x0 = x0, xEnd = xEnd)
  }
  
  for (i in 1:reps) {
    
    if (curPois > maxLen) {
      time <- c(time, rep(NA, increment))
      omega <- c(omega, rep(NA, increment))
      chainIdx <- c(chainIdx, rep(NA, increment))
      maxLen <- maxLen + increment
    }
    
    
    if (tEnd < 2 * segLen) {
      
      rslt <- ea_path(t0 = t0,
                      tEnd = tEnd,
                      x0 = x0,
                      xEnd = NULL,
                      Mpoisson = Mpoisson,
                      maxIter = 50000,
                      func = func)
      pathLen <- length(rslt$time)
      time[curPois:(curPois + pathLen - 1)] <- rslt$time
      omega[curPois:(curPois + pathLen - 1)] <- rslt$omega
      iterTime <- iterTime + rslt$iterTime
      chainIdx[curPois:(curPois + pathLen - 1)] <- rep(i, pathLen)
      curPois <- curPois + pathLen
    } else {
      for(j in 1:tSeqLen) {
        if (j == 1) {
          xSeq[j] <- func$bias_sample(tEnd = tSeq[j], t0 = t0, x0 = x0)
          rslt <- ea_path(t0 = t0,
                          tEnd = tSeq[j],
                          x0 = x0,
                          xEnd = xSeq[j],
                          Mpoisson = Mpoisson,
                          maxIter = 50000,
                          func = func)
          timeTmp <- rslt$time
          omegaTmp <- rslt$omega
          iterTimeTmp <- rslt$iterTime
        } else {
          xSeq[j] <- func$bias_sample(tEnd = tSeq[j], t0 = tSeq[j - 1], x0 = xSeq[j - 1])
          rslt <- ea_path(t0 = tSeq[j - 1],
                          tEnd = tSeq[j],
                          x0 = xSeq[j - 1],
                          xEnd = xSeq[j],
                          Mpoisson = Mpoisson,
                          maxIter = 50000,
                          func = func)
          curPathLen <- length(timeTmp)
          timeTmp <- c(timeTmp, rslt$time[-1])
          omegaTmp <- c(omegaTmp, rslt$omega[-1])
          iterTimeTmp <- iterTimeTmp + rslt$iterTime
        }
      }
      
      xEnd <- func$bias_sample(tEnd = tEnd, t0 = tSeq[tSeqLen], x0 = xSeq[tSeqLen])
      
      rslt <- ea_path(t0 = tSeq[tSeqLen],
                      tEnd = tEnd,
                      x0 = xSeq[tSeqLen],
                      xEnd = xEnd,
                      Mpoisson = Mpoisson,
                      maxIter = 50000,
                      func = func)
      
      curPathLen <- length(timeTmp)
      timeTmp <- c(timeTmp, rslt$time[-1])
      omegaTmp <- c(omegaTmp, rslt$omega[-1])
      iterTimeTmp <- iterTimeTmp + rslt$iterTime
      
      pathLen <- length(timeTmp)
      time[curPois:(curPois + pathLen - 1)] <- timeTmp
      omega[curPois:(curPois + pathLen - 1)] <- omegaTmp
      iterTime <- iterTime + iterTimeTmp
      chainIdx[curPois:(curPois + pathLen - 1)] <- rep(i, pathLen)
      curPois <- curPois + pathLen
      
    }
  }
  return(list(time = time[1:(curPois - 1)],
              omega = omega[1:(curPois - 1)],
              chainIdx = chainIdx[1:(curPois - 1)],
              iterTime = iterTime))
}


timegrid <- function(t0 = 0, tEnd = 10, stepSize = 0.01) {
  
  return(data.frame(time = seq(t0, tEnd, by = stepSize)))
}

timegrid_impute <- function(rslt, dfm) {
  rslt %>% group_by(chainIdx) %>%
    do(data.frame(time = dfm$time,
                  omega = brownian_bridge_impute(tVec = .$time[-c(1, length(.$time))], t0 = .$time[1],
                                                 tVecNew = dfm$time, tEnd = .$time[length(.$time)],
                                                 x0 = .$omega[1], xVec = .$omega[-c(1, length(.$time))],
                                                 xEnd = .$omega[length(.$time)]))) %>% select(time, omega, chainIdx)
  
}


get_val <- function(rslt, tEnd = 10, val = 5) {
  dfm <- timegrid(t0 = 0, tEnd = tEnd, stepSize = 0.01)
  rslt_impute <- timegrid_impute(rslt = rslt, dfm = dfm)
  return(rslt_impute[rslt_impute$time == val, ])
}

ea_plot <- function(rslt) {
  dfm <- data.frame(time = rslt$time, omega = rslt$omega, chainIdx = rslt$chainIdx)
  ggplot(dfm) + geom_path(aes(x = time, y = omega, group = chainIdx))
}

sde_euler1 <- function(t0 = 0, tEnd = 10, x0 = 0, stepSize = 0.01) {
  # Simulate periodical drift vis Euler–Maruyama method
  
  pathLen <- round((tEnd - t0) / stepSize) + 1
  omega <- rep(0, pathLen)
  omega[1] <- x0
  for(i in 2:pathLen) {
    omega[i] <- omega[i - 1] + sin(omega[i - 1]) * stepSize + rnorm(1, mean = 0, sd = sqrt(stepSize))
  }
  
  return(omega = omega)
}

euler1 <- function(reps = 1000, t0 = 0, tEnd = 10, x0 = 0, stepSize = 0.01) {
  
  timeGrid <- seq(t0, tEnd, by = stepSize)
  pathLen <- length(timeGrid)
  time <- rep(timeGrid, times = reps)
  omega <- rep(NA, times = reps * pathLen)
  chainIdx <- rep(1:reps, each = pathLen)
  
  for(i in 1:reps) {
    omega[(((i - 1) * pathLen + 1) : (i * pathLen))] <-
      sde_euler1(t0 = t0, tEnd = tEnd, x0 = x0, stepSize = stepSize)
  }
  
  return(data.frame(time = time, omega = omega, chainIdx = chainIdx))
}

sde_euler2 <- function(t0 = 0, tEnd = 10, x0 = 0, stepSize = 0.01, theta = 1) {
  # Simulate hyperbolic process vis Euler–Maruyama method
  
  pathLen <- round((tEnd - t0) / stepSize) + 1
  omega <- rep(0, pathLen)
  omega[1] <- x0
  for(i in 2:pathLen) {
    omega[i] <- omega[i - 1] - theta * omega[i - 1] / sqrt(1 + omega[i - 1]^2) * stepSize +
      rnorm(1, mean = 0, sd = sqrt(stepSize))
  }
  
  return(omega = omega)
}

euler2 <- function(reps = 1000, t0 = 0, tEnd = 10, x0 = 0, stepSize = 0.01, theta = 1) {
  
  timeGrid <- seq(t0, tEnd, by = stepSize)
  pathLen <- length(timeGrid)
  time <- rep(timeGrid, times = reps)
  omega <- rep(NA, times = reps * pathLen)
  chainIdx <- rep(1:reps, each = pathLen)
  tend <- rep(NA, reps)

  
  for(i in 1:reps) {
    tmp <- sde_euler2(t0 = t0, tEnd = tEnd, x0 = x0, stepSize = stepSize, theta = theta)
    omega[(((i - 1) * pathLen + 1) : (i * pathLen))]  <- tmp
    tend[i] = tmp[length(tmp)]
  }
  
  return(list(tend = tend, df = data.frame(time = time, omega = omega, chainIdx = chainIdx)))
}

get_path <- function(pathRslt, choosenIdx, Amat, tDatIdx) {
  idx <- choosenIdx
  acesChain <- NULL
  tDatLen <- length(tDatIdx)
  for (i in nrow(Amat):1){
    acesChain <- c(acesChain, Amat[i, idx])
    idx <- Amat[i, idx]
  }
  
  acesChain <- rev(acesChain)
  path <- pathRslt[1:tDatIdx[1], acesChain[1]]
  
  for (i in 2:(tDatLen - 1)) {
    path <- c(path, pathRslt[(tDatIdx[i - 1] + 1):tDatIdx[i], acesChain[i]])
  }
  
  path <-  c(path, pathRslt[(tDatIdx[tDatLen - 1] + 1):tDatIdx[tDatLen], choosenIdx])
  return(path)
}

grid_euler <- function(time = NULL, start = 0, func) {
  # Simulate SDE vis Euler–Maruyama method
  
  pathLen <- length(time)
  omega <- rep(0, pathLen)
  omega[1] <- start
  timeDiff <- diff(time)
  for(i in 2:pathLen) {
    omega[i] <- omega[i - 1] + func(omega[i - 1]) +
      rnorm(1, 0, sd = sqrt(timeDiff[i - 1]))
  }
  
  return(omega)
}

particle_euler <- function(startVec, time, func) {
  return(sapply(startVec, function(x) grid_euler(time, start = x, func)))
}

grid_gen <- function(tDat = NULL, t0 = 0, tEnd = 10, stepSize = 0.01) {
  
  tLen <- length(tDat)
  timeTmp <- seq(t0, tDat[1], by = stepSize)
  dfm <- data.frame(time = timeTmp,
                    class = c(rep("sim", times = length(timeTmp) - 1),  "obs"))
  
  for(i in 2:tLen) {
    timeTmp <- seq(tDat[i - 1], tDat[i], by = stepSize)[-1]
    dfm <- rbind(dfm, data.frame(time = timeTmp,
                                 class = c(rep("sim", times = length(timeTmp) - 1),  "obs")))
    
  }
  
  timeTmp <- seq(tDat[tLen], tEnd, by = stepSize)[-1]
  
  dfm <- rbind(dfm, data.frame(time = timeTmp,
                               class = rep("sim", times = length(timeTmp))))
  return(dfm)
}

log_weight_cal <- function(xObs = 0, xDat = 0, sd = 1) {
  
  if (length(xDat) == 1) {
    return(dnorm(xObs, mean = xDat, sd = sd, log = TRUE))
  } else {
    log_sum <- apply(xObs, 2, function(x) {
      sum(dnorm(x, mean = xDat, sd = sd, log = TRUE))
    })
    return(log_sum)
  }
}

data_gen <- function(x0 = 0, t0 = 0, tEnd = 10, sd = 0.05, nObs = 9, step_size = 0.1, drift) {
  tDat <- unique(round(seq(t0, tEnd, length.out = nObs + 2), digits = 1))
  tLen <- length(tDat)
  tDat[tLen-1] = tEnd-step_size
  xDat <- rep(x0, tLen)
  tDatDiff <- diff(tDat)

  grd <- step_size
  xp <- x0
  for(i in 2:tLen) {
    tgrid <- seq(tDat[i - 1], tDat[i], by=grd)
    for(j in 1:length(tgrid)) {
      xp <- xp + drift(xp) * grd + rnorm(1, sd = sqrt(grd))
    }  
    xDat[i] <- rnorm(1, mean = xp, sd = sd)
  }
  return(data.frame(tDat = tDat[-c(1, tLen)], xDat = xDat[-c(1, tLen)]))
}


data_gen_sin <- function(x0 = 0, t0 = 0, tEnd = 10, sd = 0.05, nObs = 9, step_size = 0.1) {
  tDat <- unique(round(seq(t0, tEnd, length.out = nObs + 2), digits = 1))
  tLen <- length(tDat)
  tDat[tLen-1] = tEnd-step_size
  xDat <- rep(x0, tLen)
  tDatDiff <- diff(tDat)

  grd <- step_size
  xp <- x0
  for(i in 2:tLen) {
    tgrid <- seq(tDat[i - 1], tDat[i], by=grd)
    for(j in 1:length(tgrid)) {
      xp <- xp + sin(xp) * grd + rnorm(1, sd = sqrt(grd))
    }  
    xDat[i] <- rnorm(1, mean = xp, sd = sd)
  }
  return(data.frame(tDat = tDat[-c(1, tLen)], xDat = xDat[-c(1, tLen)]))
}

data_gen_hyper <- function(x0 = 0, t0 = 0, tEnd = 10, sd = 0.05, nObs = 9, theta = 1, step_size = 0.1) {
  tDat <- unique(round(seq(t0, tEnd, length.out = nObs + 2), digits = 1))
  tLen <- length(tDat)
  tDat[tLen-1] = tEnd-step_size
  xDat <- rep(x0, tLen)
  tDatDiff <- diff(tDat)
  for(i in 2:tLen) {
    xDat[i] <- xDat[i - 1] - theta * xDat[i - 1] / sqrt(1 + xDat[i - 1]^2) * tDatDiff[i - 1] +
      rnorm(1, sd = sqrt(tDatDiff[i - 1]))
    rnorm(1, sd = sd)
  }
  return(data.frame(tDat = tDat[-c(1, tLen)], xDat = xDat[-c(1, tLen)]))
}

loc_chain <- function(chainIdx = 0, chainIdxVec = 0) {
  idx <- which(chainIdxVec == chainIdx)
  if (length(idx) == 0) stop("the chain number is not in the chain list")
  if (length(idx) > 1) stop("multiple chains match the chain number")
  
  return(idx)
}

log_trick_sum <- function(logWeightDfm = NULL) {
  if (length(logWeightDfm) == 0) stop("Input logWeightDfm")
  sum(apply(logWeightDfm, 1, function(x) {
    x_max <- max(x)
    return(x_max + log(mean(exp(x - x_max))))
  }))
}


sim_Eul <- function(time = NULL, start = 0, func) {
  # Simulate SDE vis Euler–Maruyama method
  
  pathLen <- length(time)
  dm      <- length(start)
  omega <- matrix(0, dm, pathLen)
  omega[,1] <- start
  timeDiff <- diff(time)
  for(i in 2:pathLen) {
    omega[,i] <- omega[,i - 1] + func(omega[,i - 1])*timeDiff[i - 1] +
      rnorm(dm, 0, sd = sqrt(timeDiff[i - 1]))
  }
  
  return(omega)
}


pf_Eul <- function(P = 100, xDat = NULL, tDat = NULL, x0 = 0, t0 = 0, tEnd = 10, sdv, stepSize, func) {
  # Particle filtering algorithm for Euler-Maruyama scheme
  # P: the number of particles
  # xDat: the observed realization
  # tDat: the time of observed realization
  # Assume tDat is strictly contained within (t0, tEnd)
  
  xDat <- c(NA, xDat, NA)
  tVec <- c(t0, tDat, tEnd) - t0
  tLen = length(tVec)
  logW <- rep(0, P) 
  logZ <- rep(0, tLen) 
  maxWt <- rep(0, tLen) 
  omega <- matrix(x0, nrow=P, ncol = tLen)
  omegaNew   <- rep(0,P)
  Amat  <- matrix(0,  nrow=P, ncol = tLen)
  
  path <- rep(x0, tLen)  # includes endpoints

  for (i in 2:tLen) {
    if((i %% 50) == 0) print(i)

    prp  <- sim_Eul(time=seq(tVec[i-1], tVec[i], by = stepSize), omega[,i-1], func=func)
    omegaNew  <- prp[,ncol(prp)]

    if(i == tLen || is.na(xDat[i])) {   # No obs at tEnd, and is.na is for imputation
      logW <- rep(0,P)
    } else {
      logW <- dnorm(omegaNew, xDat[i], sd = sdv, log=T)
    }

    maxWt[i] <-  max(logW)
    logW <- logW - maxWt[i]

    Amat[,i] <- sample(1:P, size = P, replace = TRUE, prob = exp(logW))
    omega[,i] <- omegaNew[Amat[,i]] 

    logZ[i] = log(sum(exp(logW))) + maxWt[i] 
  }
    
  AIdx <- Amat[1, tLen]  
  path[-1]   <- omega[AIdx, tLen]
  
  for (i in (tLen-1):1){
    path[i] <- omega[AIdx,i]
    AIdx    <- Amat[AIdx,i] 
  }
  
  return(list(path = data.frame(time = tVec, omega = path, xDat = xDat), 
              logProb = sum(logZ), est = data.frame(x = omega[,tLen-1], w = exp(logW)) ))
}

pmcmc_Eul <- function(reps = 1000, P = 100, x0 = 0, t0 = 0, tEnd = 10, xDat = NULL, tDat = NULL, 
                       stepSize = 0.01, sdv = 1, theta=NULL, func, updt_data=FALSE) {
  
  if (length(tDat) == 0) {
    stop("Observation must be provided")
  }
  
  tDatLen <- length(tDat)
  xObsMat = matrix(0,nrow=length(xDat), ncol = reps)

  if(is.null(theta)) {
    updt_theta = T 
    theta_n = func(NULL, NULL, dr='prop')
  } else {
    updt_theta = F
    theta_n = theta
  }
  theta_vec   = rep(theta_n, reps)
  
  init <- pf_Eul(P, xDat, tDat, x0, t0, tEnd, sdv, stepSize, function(x) func(x,theta_n))
  Z_o  <- init$logProb 
  path_o <- init$path$omega 

  paths <- matrix(0, reps, length(path_o))
  paths[1,] <- path_o
  
  for (i in 2:reps) {
    if((i %% 50) == 0) print(i)

    if(updt_theta) {
      theta_n <- func(NULL, NULL, dr='prop')   # Proposing from the prior
    } else {
      theta_n = theta_vec[i-1]
    }

      prop <- pf_Eul(P, xDat, tDat, x0, t0, tEnd, sdv, stepSize,  function(x) func(x,theta_n))
      Z_n <- prop$logProb    
      acc <- Z_n - Z_o

    
    if (log(runif(1)) < acc) {
      xObsMat[,i] = prop$path$omega[2:(tDatLen+1)]
      paths[i,]      <- prop$path$omega
      theta_vec[i]   <- theta_n
      Z_o            <- Z_n
    } else {
      xObsMat[,i] = xObsMat[,i-1]
      paths[i,] <- paths[i-1,]
      theta_vec[i] <- theta_vec[i-1]
    }
    # update xDat (update data for debugging)
    if(updt_data) {
      xDat <- xObsMat[,i] + rnorm(length(xDat), sd = sdv)
      stop('You need to update Z_o as well, havent implemented this!')
    }
    ## Need to recalculate
    
  }
  rownames(xObsMat) <- tDat
  
  return(list(paths = paths, theta = theta_vec, xObsMat=xObsMat))
}
  
pf_Fearn <- function(P = 100, xDat = NULL, tDat = NULL, x0 = 0, t0 = 0, tEnd = 10, theta, sdv, func) {
  # Particle filtering algorithm for Fearnhead scheme
  # P: the number of particles
  # xDat: the observed realization
  # tDat: the time of observed realization
  # Assume tDat is strictly contained within (t0, tEnd)

  func <- func(theta=theta, sdv=sdv)
  
  xDat <- c(NA, xDat, NA)
  tVec <- c(t0, tDat, tEnd) - t0
  tLen = length(tVec)
  logW <- rep(0, P) 
  logZ <- rep(0, tLen) 
  maxWt <- rep(0, tLen) 
  omega <- matrix(x0, nrow=P, ncol = tLen)
  Amat  <- matrix(0,  nrow=P, ncol = tLen)
  
  path <- rep(x0, tLen)  # includes endpoints

  for (i in 2:tLen) {
    tDelta <- tVec[i] - tVec[i-1]
    omegaNew   <- rnorm(P,omega[,i-1], sqrt(tDelta))
    
    if(is.na(xDat[i])) {          # is.na is for imputation
      logW <- rep(0,P)
    } else if(i == tLen) {   # No obs at tEnd
      logW <- func$H1_lglk(omegaNew, tEnd, x0, t0) 
    } else {
      logW <- dnorm(omegaNew, xDat[i], sd = sdv, log=T)
    }
    #logW <- logW + func$A_func(omegaNew) - func$A_func(omega[,i-1]) 
    for(p in 1:P) {
      tVecNew <- poisson_process(t0 = tVec[i-1], tEnd = tVec[i], lambda = func$Mpoisson)
      if (length(tVecNew) > 0) {
      
        xVecNew <- brownian_bridge(tVec = tVecNew, t0 = tVec[i-1], tEnd = tVec[i], 
                                   x0 = omega[p,i-1], xEnd = omegaNew[p])
         logW[p] = logW[p] + sum(log(func$Mpoisson-func$phi(xVecNew)))  - length(tVecNew)*log(func$Mpoisson)
      }
    }

    maxWt[i] <-  max(logW)
    logW <- logW - maxWt[i]
    
    Amat[,i] <- sample(1:P, size = P, replace = TRUE, prob = exp(logW))
    omega[,i] <- omegaNew[Amat[,i]] 

    logZ[i] = log(sum(exp(logW))) + maxWt[i] 
  }
    
  AIdx <- Amat[1, tLen]  
  path[-1]   <- omega[AIdx, tLen]
  
  for (i in (tLen-1):1){
    path[i] <- omega[AIdx,i]
    AIdx    <- Amat[AIdx,i] 
  }
  if(any(is.na(xDat))) {        
    omeg_ret <- omega
  } else omeg_ret = NULL
  
  return(list(path = data.frame(time = tVec, omega = path, xDat = xDat), 
              logProb = sum(logZ), est = data.frame(x = omega[,tLen-1], w = exp(logW)), omeg_ret = omeg_ret ))
}
  
pf_Fearn_old <- function(P = 100, xDat = NULL, tDat = NULL, x0 = 0, t0 = 0, tEnd = 10, theta, sdv, func) {
  # Particle filtering algorithm for Fearnhead scheme
  # P: the number of particles
  # xDat: the observed realization
  # tDat: the time of observed realization
  # Assume tDat is strictly contained within (t0, tEnd)

  func <- func(theta=theta, sdv=sdv)
  
  xDat <- c(NA, xDat, NA)
  tVec <- c(t0, tDat, tEnd) - t0
  tLen = length(tVec)
  logW <- rep(0, P) 
  logZ <- rep(0, tLen) 
  omegaNew <- rep(0, P) 
  maxWt <- rep(0, tLen) 
  omega <- matrix(x0, nrow=P, ncol = tLen)
  Amat  <- matrix(0,  nrow=P, ncol = tLen)
  
  path <- rep(x0, tLen)  # includes endpoints
    

  for (i in tLen:1) {

    if(i==tLen) {   # assume also no obs at tEnd
      for(p in 1:P) omegaNew[p] = func$bias_sample(tEnd, t0, x0)
      # logW <- func$H1_lglk(omegaNew, tEnd, x0, t0) # + func$A_func(omegaNew) - func$A_func(x0) 
    } else {

      for(p in 1:P) {
      
        tVecNew <- poisson_process(t0 = tVec[i], tEnd = tVec[i+1], lambda = func$Mpoisson)
        xVecNew <- brownian_bridge(tVec = c(tVec[i],tVecNew), t0 = t0, tEnd = tVec[i+1], 
                                   x0 = x0, xEnd = omega[p,i+1])
        omegaNew[p] <- xVecNew[1]
        xVecNew     <- xVecNew[-1]
        if (length(tVecNew) > 0) {
          logW[p] = logW[p] + sum(log(func$Mpoisson-func$phi(xVecNew)))  - length(tVecNew)*log(func$Mpoisson)
        }
      }
      if(!is.na(xDat[i]) && i > 1) {         # is.na is for imputation 
        logW <- dnorm(omegaNew, xDat[i], sd = sdv, log=T)
      }
    }

    maxWt[i] <-  max(logW)
    logW <- logW - maxWt[i]
    
    Amat[,i] <- sample(1:P, size = P, replace = TRUE, prob = exp(logW))
    omega[,i] <- omegaNew[Amat[,i]] 

    logZ[i] = log(sum(exp(logW))) + maxWt[i] 
  }
    
  path[1]   <- omega[1, 1]
  AIdx <- Amat[1, 1]  
  
  for (i in 2:tLen){
    path[i] <- omega[AIdx,i]
    AIdx    <- Amat[AIdx,i] 
  }
  
  return(list(path = data.frame(time = tVec, omega = path, xDat = xDat), 
              logProb = sum(logZ), est = data.frame(x = omega[,tLen-1], w = exp(logW)) ))
}

pmcmcFearn <- function(reps = 2000, P = 100, x0 = 0, t0 = 0, tEnd = 10, xDat = NULL, tDat = NULL, theta=1, sdv = 1, func) {

  init <- pf_Fearn(P, xDat, tDat, x0, t0, tEnd, theta, sdv, func)
  Z_o  <- init$logProb
  path_o <- init$path$omega 

  paths <- matrix(0, reps, length(path_o))
  paths[1,] <- path_o

  for(i in 2:reps) {
    if((i %% 100) == 0) print(i)
    prop <- pf_Fearn(P, xDat, tDat, x0, t0, tEnd, theta, sdv, func)
    Z_n <- prop$logProb
    if(log(runif(1)) <= Z_n - Z_o) {
      paths[i,] <- prop$path$omega 
    } else {
      paths[i,] <- paths[i-1,]
    }
  }
  paths = as.data.frame(paths)
  colnames(paths) <- init$path$time
  return(paths)
}

  
FearnPredict <- function(P = 100, t0, x0, xDat, tPred = NULL, theta=1, sdv=.1,  func) {
  # Particle filtering algorithm for Fearnhead scheme
  # P: the number of particles
  
  func <- func(theta=theta, sdv=sdv)
  tVec  <- tPred - t0
  tLen  <- length(tVec)
  logW  <- matrix(0, P, tLen) 
  xPred <- matrix(0, P, tLen) 
  err   <- rep(0,tLen)
  pred  <- rep(0,tLen)

  print(func$Mpoisson)

  for (i in 1:tLen) {
    xPred[,i]   <- rnorm(P,x0, sqrt(tVec[i]))
    
    #logW <- func$H1_lglk(omegaNew, tEnd, x0, t0) 
    #logW <- logW + func$A_func(omegaNew) - func$A_func(omega[,i-1]) 
    for(p in 1:P) {
      tVecNew <- poisson_process(t0 = 0, tEnd = tVec[i], lambda = func$Mpoisson)
      if (length(tVecNew) > 0) {
        xVecNew <- brownian_bridge(tVec = tVecNew, t0 = 0, tEnd = tVec[i], 
                                   x0 = x0[p], xEnd = xPred[p,i])
         logW[p,i] = logW[p,i] + sum(log(func$Mpoisson-func$phi(xVecNew)))  - length(tVecNew)*log(func$Mpoisson)
      }
    }
    wt = exp(logW[,i])
    wt = wt/sum(wt)
    pred[i] = sum(wt*xPred[,i])
    err[i] = abs(xDat[i]- pred[i])
  }
  return(tibble(pred = pred, err = err))
}

  
FearnPredict_old <- function(P = 100, x0, xDat, tPred = NULL, t0 = 0, func) {
  # Particle filtering algorithm for Fearnhead scheme
  # P: the number of particles
  # tPred: the time of observed realization
  # Assume tPred is strictly contained within (t0, tEnd)
  
  tVec <- c(t0, tPred) - t0
  tLen = length(tVec)
  logW <- rep(0, P) 
  logZ <- rep(0, tLen) 
  maxWt <- rep(0, tLen) 
  omega <- matrix(x0, nrow=P, ncol = tLen)
  #omega[,1] <- x0
  Amat  <- matrix(0,  nrow=P, ncol = tLen)
  
  path <- rep(0, tLen)  # includes endpoints

  for (i in 2:tLen) {
    tDelta <- tVec[i] - tVec[i-1]
    omegaNew   <- rnorm(P, omega[,i-1], sqrt(tDelta))
    
    logW <- rep(0,P)
    logW <- logW + func$A_func(omegaNew) - func$A_func(omega[,i-1]) 
    for(p in 1:P) {
      tVecNew <- poisson_process(t0 = tVec[i-1], tEnd = tVec[i], lambda = func$Mpoisson)
      if (length(tVecNew) > 0) {
      
        xVecNew <- brownian_bridge(tVec = tVecNew, t0 = tVec[i-1], tEnd = tVec[i], 
                                   x0 = omega[p,i-1], xEnd = omegaNew[p])
         logW[p] = logW[p] + sum(func$Mpoisson-func$phi(xVecNew))
      }
    }

    maxWt[i] <-  max(logW)
    logW <- logW - maxWt[i]
    
  Amat[,i] <- sample(1:P, size = P, replace = TRUE, prob = exp(logW))
  omega[,i] <- omegaNew[Amat[,i]] 

  logZ[i] = log(sum(exp(logW))) + maxWt[i] 
  }
  
AIdx <- Amat[1, tLen]  
path[-1]   <- omega[AIdx, tLen]

pred_err <- 0
for (i in (tLen-1):1){
  path[i]  <- omega[AIdx,i]
  AIdx     <- Amat[AIdx,i] 
  pred_err <- pred_err + (path[i]-xDat[i-1])^2
}

return(list(path = data.frame(time = tVec, omega = path), 
            logProb = sum(logZ), pred_err = pred_err))
}
