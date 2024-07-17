SIRICC_VDA <- function(x0, yT, cT, window_length, mu, B0, noise, N_scale,
                       state_pop) {
  B         <- B0[1:2, 1:2] # only want covariance for beta and gamma
  Binv      <- solve(B) # invert covariance matrix
  
  test <- optim(fn = costFx, par = x0, yT = yT, 
                cT = cT, 
                s = window_length,
                mu = mu, Binv = Binv, noise = noise, N_scale = N_scale,
                state_pop = state_pop,
                method = "Nelder-Mead")
  
  test$estimate        <- test$par 
  names(test$estimate) <- c("beta", "gamma", "N", "C0")
  
  preds <- deSolve::ode(y     = c(C = cT[1]) ,
                        times = 0:(length(yT)+0),
                        func  = SIRICC_ode,
                        parms = test$estimate)
  new_G <- SIRICC_ode(t     = 0, 
                      C     = preds[-1,2],
                      parms = test$estimate)[[1]]
  
  temp <- data.frame(fitG  = new_G, 
                     fitC  = preds[-1,2],
                     times = 1:(length(yT)+0))
  p <- ggplot(temp)+
    geom_line(data.frame(trueG = yT, trueC = cT), 
              mapping =aes( x = trueC, y = trueG, color = "truth"))+
    geom_line(mapping = aes(x = fitC, y = fitG, color = "fit"))
  
  return(list(params = test$estimate, fig = p, initC = temp$fitC[window_length]))
}


costFx <- function(x0, yT, cT, s, mu, Binv, noise, N_scale, state_pop) {
  
  B <- Binv
  R <- noise
  
  beta   <- x0[1]
  gamma  <- x0[2]
  N      <- x0[3]
  Cstart <- x0[4]
  mu0    <- mu
  
  # Jx = 1/2*(x0 - mu0)'* B^(-1)*(x0 - mu0)
  # see line 3 -- B was pre-inverted
  Jx = 1/2*(x0[1:2] - mu0[1:2]) %*% B %*% (x0[1:2] - mu0[1:2])
  
  parms <- c(beta  = x0[1], 
             gamma = x0[2],
             N     = x0[3],
             C0    = x0[4])
  
  times  <- 0:length(yT)
  tC_mat <- deSolve::ode(y = cT[1], times = times, func = SIRICC_ode, parms = parms)
  new_G  <- SIRICC_ode(0, tC_mat[-1,2], parms)[[1]]
  
  diff  <- sum((1/diag(R) * (new_G - yT)^2))
  Jx  <- Jx + 1/2*diff;

  # penalize if effective population is smaller than final cumulative number
  if (N < cT[length(cT)]) { 
    Jx  <- 2*Jx
  }

  # penalize if parameters aren't realistic
  if (beta <= 0 | gamma <= 0 | beta/gamma > 20){
    Jx  <- 2*Jx
  }

  # penalize if effective population is smaller than starting population
  if (N > state_pop){
    Jx  <- 2*Jx
  }

  # penalize if initial condition is larger than effective population
  if (Cstart > N){
    Jx  <- 2*Jx
  }

  return(Jx)
}




# This function calculates the ICC curve for the SIR model, and includes
# the initial condition as C0 = C(0) = N (1 - kappa).
# The value of dy is set to 0 in regions where x >= N or if dy < 0.
# Created 5/24/2019 by J. Lega
# Updated 7/14/2024 by H. Biegel
#
# e.g., parms = c(beta = 0.1, gamma = 0.01, N = 1000, C0 = 100)
# x = c(C = 50)
SIRICC_ode <- function(t, C, parms, ...){
  with(as.list(c(parms, C)), {

    beta  <- parms["beta"]
    gamma <- parms["gamma"]
    N     <- parms["N"]
    C0    <- parms["C0"]
    dx    <- (beta*C + N*gamma*log(abs(N-C)/(N-C0))) * (1-C/N)
    dx[which(N-C <= 0)] <- 0
    
    return(list(dx))
  })
}

quiet_SIRICC_VDA <- purrr::quietly(SIRICC_VDA)


