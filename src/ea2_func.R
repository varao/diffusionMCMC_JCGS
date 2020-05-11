ea2_func <- function(const = 1, theta = 1, sdv) {

  Mpoisson_g <- (const^2 * theta^2 + const * theta) / 2
  L_g <- -theta
  logPrior <- dexp(x = theta, rate = 1, log = TRUE)
  

  A_func <- function(x, .theta = theta) {
    return(theta-theta*sqrt(1+x^2))
  }

  H1_lglk <- function(xEnd, tEnd, x0, t0, .const = const, .theta = theta) {
    # Calculate the log-likelihood h1(x_end)
    # Follow the example of EA2: dX_t = -theta X_t / sqrt(1 + X_t^2)dt + dB_t
    
    lglk = -.const * .theta * sqrt(1 + xEnd^2) - (xEnd - x0)^2 / (2 * (tEnd - t0)) +
      .const * .theta *  sqrt(1 + x0^2)
    return(lglk = lglk)
  }
  
  H1_grad <- function(xEnd, tEnd, x0, t0, .const = const, .theta = theta) {
    # Calculate the gradient of h1(x_end)
    # Follow the example of EA2: dX_t = -theta X_t / sqrt(1 + X_t^2)dt + dB_t
    
    grad = -.const * .theta * xEnd / sqrt(1 + xEnd^2)  - (xEnd - x0) / (tEnd - t0)
    return(grad = grad)
  }
  
  BB_lglk <- function(tVec, xVec, t0, tEnd, x0, xEnd) {
    # Calculate the log-likelihood of the Brownian bridge
    # starts at (t0, x0) and ends at (tEnd, xEnd)
    # Suppose x0 and xEnd is fixed
    # x_i ~ Normail(mu_i, sigma_i)
    # mu_i = (t_i - t_0) / (t_{i + 1} - t_0) x_{i + 1} + (t_{i + 1} - t_i)  / (t_{i + 1} - t_0} x_0
    # sigma_i^2 = (t_{i + 1} - t_i)(t_i - t_0) / (t_{i + 1} - t0)
    
    tLen <- length(tVec)
    
    if (tLen == 0) {
      # The Path just has begining and end point
      return(lglk = 0)
    }
    
    if (tLen == 1){
      # The Path just has 1 point in the middle
      mu <- (tVec - t0) / (tEnd - t0) * xEnd + (tEnd - tVec) / (tEnd - t0) * x0
      sigma2 <- (tEnd - tVec) * (tVec - t0) / (tEnd - t0)
      lglk <- dnorm(xVec, mean = mu, sd = sqrt(sigma2), log = TRUE)
      return(lglk = lglk)
    }
    
    # calculate t_{i + 1} - t_i in vector, i <= tLen
    tDiff <- c(tVec[-1], tEnd) - tVec
    
    # calculate mu
    mu <- (tVec - t0) / (c(tVec[-1], tEnd) - t0) * c(xVec[-1], xEnd) +
      tDiff / (c(tVec[-1], tEnd) - t0) * x0
    
    # calculate sigma
    sigma2 <- tDiff * (tVec - t0) / (c(tVec[-1], tEnd) - t0)
    
    # calculate log-likelihood of the path
    lglk <- sum(dnorm(xVec, mean = mu, sd = sqrt(sigma2), log = TRUE))
    
    return(lglk = lglk)
  }
  
  BB_grad <- function(tVec, xVec, t0, tEnd, x0, xEnd) {
    # Calculate the log-likelihood and gradient of the Brownian bridge
    # starts at (t0, x0) and ends at (tEnd, xEnd)
    # Suppose x0 and xEnd is fixed
    # x_i ~ Normail(mu_i, sigma_i)
    # mu_i = (t_i - t_0) / (t_{i + 1} - t_0) x_{i + 1} + (t_{i + 1} - t_i)  / (t_{i + 1} - t_0} x_0
    # sigma_i^2 = (t_{i + 1} - t_i)(t_i - t_0) / (t_{i + 1} - t0)
    
    
    tLen <- length(tVec)
    
    if (tLen == 0) {
      # The Path just has begining and end point
      return(list(lglk = 0, grad = 0))
    }
    
    if (tLen == 1){
      # The Path just has 1 point in the middle
      mu <- (tVec - t0) / (tEnd - t0) * xEnd + (tEnd - tVec) / (tEnd - t0) * x0
      sigma2 <- (tEnd - tVec) * (tVec - t0) / (tEnd - t0)
      partialDerivXi <- -(xVec - mu) / sigma2
      grad <- c(partialDerivXi, 0)
      grad <- grad + c(0, -partialDerivXi* (tVec - t0) / (tEnd - t0))
      return(grad = grad)
    }
    
    # calculate t_{i + 1} - t_i in vector, i <= tLen
    tDiff <- c(tVec[-1], tEnd) - tVec
    
    # calculate mu
    mu <- (tVec - t0) / (c(tVec[-1], tEnd) - t0) * c(xVec[-1], xEnd) +
      tDiff / (c(tVec[-1], tEnd) - t0) * x0
    
    # calculate sigma
    sigma2 <- tDiff * (tVec - t0) / (c(tVec[-1], tEnd) - t0)
    
    # The gradient of the log-likelihood is the summation of two parts
    # The first part is the partial derivative of Normal(mu_i, sigma_i) on x_i
    # It is -(x_i - mu_i) / sigma^2 for i = 1, 2, ...., tLen
    # It is 0 for x0 and xEnd
    
    partialDerivXi <- -(xVec - mu) / sigma2
    grad <- c(partialDerivXi, 0)
    
    # The second part is the partial derivative of Normal(mu_i, sigma_i) on x_{i + 1}
    # It is (x_i - mu_i) / sigma^2 * (ti - t0) / (t_{i + 1}  - t0), for i = 2, 3, ..., tLen, end
    # It is 0 for x0 and x1
    
    grad <- grad + c(0, -partialDerivXi* (tVec - t0) / (c(tVec[-1], tEnd) - t0))
    
    return(grad = grad)
  }
  
  phi <- function(x, .const = const, .theta = theta, .L = L_g) {
    if (length(x) == 0) print("NULL comes from Phi")
    return(0.5 * (.const^2 * .theta^2 * x^2 / (1 + x^2) -
                    .const * .theta / (1 + x^2)^1.5 - .const * .L))
    
  }
  
  path_lglk <- function(xVec,
                        .Mpoisson = Mpoisson_g,
                        .const = const,
                        .phi = phi) {
    # Calculate the log-likelihood of Multiply(1 - phi(x)/M)
    # M is the upper bound of phi(x) on (t0, tEnd)
    
    if(length(xVec) == 0)  {
      return(lglk = 0)
    } else {
      thn <- .Mpoisson - .phi(xVec)
      thn[thn < 0] <- 0    # occasional numerical issues for large xVec produce annoying warnings.
      lglk = sum(log(thn))   
      return(lglk = lglk)
    }
  }
  
  path_grad <- function(xVec,
                        .Mpoisson = Mpoisson_g,
                        .const = const,
                        .theta = theta, .phi=phi) {
    # Calculate the gradient of log(M - phi(x))
    # M is the upper bound of phi(x) on (t0, tEnd)
    
    if(length(xVec) == 0)  {
      return(grad = 0)
    }
    
    grad = -0.5 * (2 * .const^2 * .theta^2 * xVec / (1 + xVec^2) ^ 2 + 3 * .const * .theta * xVec / (1 + xVec^2)^2.5) /
      (.Mpoisson - .phi(xVec))
    
    return(grad = grad)
  }
  
  
  
  bias_sample <- function(tEnd, t0, x0, .const = const, .theta = theta) {
    # Sample from exp(-cos(u) + cos(x) - (u - x)^2 / 2T)
    # Use reject sampling to sample from this distribution
    # propose from a normal density: exp(- (u - x) ^2 / 2T)
    # i.e. q ~ N(x, T)
    # if (length(x0)  == 0) x0 <- rnorm(1, mean = 0, sd = sqrt(.1))
    
    while (TRUE) {
      u <- rnorm(1, mean = x0, sd = sqrt(tEnd - t0))
      #if (log(runif(1)) < .const * .theta -.const * .theta * (sqrt(1 + u^2))) {
      if (log(runif(1)) <  -.const * .theta * (sqrt(1 + u^2))) {
        return(u)
      }
    }
  }
  
  obs_lglk <- function(xObs, xDat, .sd = sdv) {
    # Assume normal distribution with known variance
    # x_i ~ Normal(mu_i, sigma^2)
    
    if (length(xObs) == 0) {
      stop("obsKernel doesn't accept NULL input")
    }
    
    lglk <- sum(-.5 * log(2 * pi * .sd^2) - (xDat - xObs)^2 / (2 * .sd^2))
    
    return(lglk = lglk)
  }
  
  obs_grad <- function(xObs, xDat, .sd = sdv) {
    # Assume normal distribution with known variance
    # x_i ~ Normal(mu_i, sigma^2)
    
    if (length(xObs) == 0) {
      stop("obsKernel doesn't accept NULL input")
    }
    
    grad <- -(xObs - xDat) / .sd^2
    
    return(grad = grad)
    
  }
  
  
  poisson_lglk <- function(tVec, Mpoisson = Mpoisson_g,
                           t0, tEnd) {
    # Calculate the log-likelihood of Poisson process
    # lglk = M^n * exp(-MT)
    
    n <- length(tVec)
    lglk <- n * log(Mpoisson) - Mpoisson * (tEnd - t0)
    return(lglk = lglk)
  }
  
  MHprop <- function() return(rexp(1))
  
  return(list(H1_lglk = H1_lglk, H1_grad = H1_grad,
              BB_lglk = BB_lglk, BB_grad = BB_grad,
              path_lglk = path_lglk, path_grad = path_grad,
              phi = phi,
              bias_sample = bias_sample,
              obs_grad = obs_grad, obs_lglk = obs_lglk,
              #poisson_lglk = poisson_lglk,
              const = const, Mpoisson = Mpoisson_g, sdv = sdv, L = L_g, A_func = A_func,
              theta = theta, logPrior = logPrior, MHprop = MHprop))
}

drift2 <- function(x, theta, dr = T)  {
  if(dr==T) {                 # Ugly hack
    return(-theta*x/sqrt(1+x*x))
  } else if(dr == 'prop') {
    return(rexp(1,1))
  } else {
    return(dexp(x = theta, rate = 1, log = TRUE))
  }
}
