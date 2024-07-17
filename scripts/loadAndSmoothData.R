loadAndSmoothData <- function(state_abbr, 
                              forecast_date, 
                              forecast_window = 30, 
                              historical_window = 30, 
                              path = "./",
                              smoothing_window = 7) {
  forecast_date <- as.Date(forecast_date, format = "%m-%d-%Y")
  state_data <- read.csv(paste0(path, "data/state_", state_abbr, ".csv")) |>
    mutate(status = case_when(
      date < forecast_date - historical_window + 1~ "historical",
      date <= forecast_date ~ "available",
      date <= forecast_date + forecast_window  ~ "forecast window",
      date > forecast_date + forecast_window  ~ "future beyone forecast window"
    ))
  
  # case incidence data
  G     <- state_data$positiveIncrease[which(state_data$status %in% 
                                               c("historical", "available"))]
  smG   <- zoo::rollmean(G, k = smoothing_window, na.pad = FALSE, 
                         align = "center")
  smSmG <- pmax(smooth_data(G),1)

  # death incidence data
  dG <- state_data$deathIncrease[which(state_data$status %in% 
                                         c("historical", "available"))]
  smSmdG <- pmax(smooth_data(dG),1)
  
  # make cumulative
  C  <- cumsum(smSmG) + state_data$positive[1] - smSmG[1]
  DC <- cumsum(smSmdG) + state_data$death[1] - smSmdG[1]
  
  smoothed_data <- data.frame(date = state_data$date[which(state_data$status %in%
                                                             c("historical", "available"))],
                              G = smSmG, # smoothed case incidence
                              C = C, # cumulative cases
                              DG = smSmdG, # smoothed death incidence
                              DC = DC,
                              status = state_data$status[which(state_data$status %in%
                                                               c("historical", "available"))]) # cumulative deaths 
  
 
  return(list(smoothed = smoothed_data,
              raw      = state_data))
}

smooth_data <- function(x, smoothing_window = 7){

  half_window <- floor((smoothing_window-1)/2)
  smoothX <- zoo::rollmean(x, k = smoothing_window, 
                           na.pad = TRUE, align = "center")
  smoothSmoothX <- zoo::rollmean(smoothX, k = smoothing_window, 
                                 na.pad = TRUE, align = "center")
  
  for (j in 1:(smoothing_window - 1)) {
    minIndFront <- max(j - half_window, 1)
    maxIndFront <- min(j + half_window, length(x))
    
    minIndBack <- max(length(x) - j - half_window, 1)
    maxIndBack <- min(length(x) - j + half_window, length(x))
    smoothSmoothX[j] <- mean(x[minIndFront:maxIndFront])
    smoothSmoothX[length(x)-j+1] <- mean(x[minIndBack:maxIndBack])
  }
  return(smoothSmoothX)
}

