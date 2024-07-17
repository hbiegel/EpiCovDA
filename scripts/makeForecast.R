makeForecast <- function(state_abbr        = "AZ", 
                         forecast_date     = "06-03-2020", 
                         forecast_window   = 30, 
                         historical_window = 30, 
                         r.seed            = 101,
                         path              = "./") {
  
  # load smoothed data
  state_dats <- loadAndSmoothData(state_abbr        = state_abbr,
                                  forecast_date     = forecast_date,
                                  forecast_window   = forecast_window,
                                  historical_window = historical_window,
                                  path              = path)
  
  # identify smoothed data to use for forecasting
  G <- state_dats$smoothed |>
    filter(status == "available") |>
    pull(G)
  
  # create pseudo observation
  # approximately normal with mean = variance given by smoothed incidence value
  # normal approximation justified by poisson distribution for large N
  noise <- diag(G)
  set.seed(r.seed)
  
  nG <- pmax(G + MASS::mvrnorm(n = 1, mu = rep(0,length(G)), Sigma = noise), 0)
  
  # Define the noisy cumulative values as the sum of the noisy indicence
  # recenter around truth
  nC <- cumsum(nG) 
  
  
  # load state population 
  # will use to set initial search for effective population
  state_pop <- read.csv(paste0(path,
                               "/param_data/state_pops/state_pop",
                               state_abbr, ".csv")) |> 
    as.numeric()
  
  # load prior information
  p_params <- loadEpiGroParams(path = path)
  N_scale  <- mean(p_params$Ns)
  # covariance matrix
  B0       <- cov(p_params)
  
  # average parameters
  mu_beta  <- mean(p_params$betas)
  mu_gamma <- mean(p_params$gammas)
  mu_N     <- mean(p_params$Ns)
  mu       <- c(mu_beta, mu_gamma, mu_N)
  
  # Set initial conditions
  P0 <- c(mu_beta, mu_gamma, state_pop/3, nG[1])
  
  # Assimilate real data using K observations
  fitParams <- quiet_SIRICC_VDA(x0 = P0, 
                                yT = nG, 
                                cT = nC,
                                window_length = historical_window, 
                                mu = mu, 
                                B0 = B0, 
                                noise = noise, 
                                N_scale = N_scale,
                                state_pop = state_pop)$result
  
  preds <- predictSIRICC(initC = fitParams$initC,
                         times = seq(0, forecast_window, 1),
                         parms = fitParams$params) |>
    mutate(date = as.Date(forecast_date, format = "%m-%d-%Y") + Time) |>
    mutate(date = format(date, "%m-%d-%Y")) |>
    mutate(status = "prediction")
  
  fit <- predictSIRICC(initC = nC[1],
                         times = seq(0, historical_window, 1),
                         parms = fitParams$params) |>
    mutate(date = as.Date(forecast_date, format = "%m-%d-%Y") + Time - historical_window) |>
    mutate(date = format(date, "%m-%d-%Y")) |>
    mutate(status = "fit") |>
    mutate(date = as.Date(date, format = "%m-%d-%Y"))
  
  
  out_data <- left_join(
    state_dats$raw |>
      mutate(date = as.Date(date)) |>
      mutate(date = format(date, "%m-%d-%Y")),
    preds |> select(date, predictedIncidence = G, predictedC = C),
    by = "date"
  ) |>
    mutate(smoothedIncidence = smooth_data(positiveIncrease))|> 
    mutate(date = as.Date(date, format = "%m-%d-%Y"))
  
  
  plt <- ggplot(out_data |>
                  filter(date <= as.Date(forecast_date, format = "%m-%d-%Y") + forecast_window + 30)) +
    geom_line(mapping = aes(x = date, y = positiveIncrease, color = "Raw Data")) +
    geom_line(mapping = aes(x = date, y = smoothedIncidence, color = "Time-Averaged"), linewidth = 1,
              col = "black")+
    geom_point(out_data |>
                 filter(status %in% c("available")), mapping = aes(x = date, y = positiveIncrease, 
                                                                   color = "Used for Fitting"))+
    geom_line(data = fit, mapping = aes(x = date, y = G, color = "Model Fit"), linewidth = 1) +
    geom_line(mapping = aes(x = date, y = predictedIncidence, color = "Predicted"), linewidth = 1,
              linetype = "dashed") +
    scale_color_manual(values = c("Predicted" = "deeppink2",
                                  "Time-Averaged Data" = "black",
                                  "Raw Data" = "black",
                                  "Model Fit" = "deeppink2",
                                  "Used for Fitting" = "deepskyblue3"))+
    theme_bw() + 
    labs(x = "Date", y = "Incident Cases", color = element_blank()) +
    theme(legend.position = "top")
  
  return(list(plt = plt, dat = out_data, 
              forecast_date = as.Date(forecast_date, format = "%m-%d-%Y")))
}


predictSIRICC <- function(initC = 100, 
                          times = seq(0, 30, 1), 
                          parms = c(beta = 0.2, gamma = 0.1, N = 1000, C0 = 50)) {
  preds <- deSolve::ode(y     = c(C = initC) ,
                        times = times,
                        func  = SIRICC_ode,
                        parms = parms)
  new_G <- SIRICC_ode(t     = 0, 
                      C     = preds[-1,2],
                      parms = parms)[[1]]
  
  return(data.frame(Time = times[-1],
                    G = new_G,
                    C = preds[-1, 2]))
}


