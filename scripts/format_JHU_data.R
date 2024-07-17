format_JHU_data <- function(all_data, state_fip){
  state_dat <- all_data |> 
    filter(FIPS == state_fip) |>
    mutate(state_abbr = usdata::state2abbr(Province_State)) |>
    arrange(Date) |>
    select(date         = Date, 
           state        = Province_State, 
           state_abbr   = state_abbr,
           positive     = Confirmed, 
           death        = Deaths, 
           hospitalized = People_Hospitalized) |>
    mutate(positiveIncrease = c(0, positive[-1] - positive[-length(positive)]),
           deathIncrease    = c(0, death[-1] - death[-length(death)]),
           hospitalIncrease = c(0, hospitalized[-1] - 
                                  hospitalized[-length(hospitalized)])) |>
    filter(date != min(date))
  
  # check if missing any dates
  minDate = min(state_dat$date)
  maxDate = max(state_dat$date)
  
  missingDates <- setdiff(minDate:maxDate, state_dat$date)
  if ( length(missingDates) > 0 ) {
    print(sprintf("Missing dates for %s.",unique(state_dat$state_abbr)))
  }
  
  return(state_dat)
}