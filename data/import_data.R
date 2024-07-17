# update JHU data

library(tidyverse)

source("scripts/format_JHU_data.R")

# path to data on github
JHU_path <- paste0("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/",
                   "master/csse_covid_19_data/csse_covid_19_daily_reports_us/")

# specify dates to fit epi curve
dates_of_interest <- (as.Date("04-13-2020", format = "%m-%d-%Y") + (0:365)) |>
  format("%m-%d-%Y")

# import data from github
all_data <- map_df(dates_of_interest, .f = function(dt){
  read.csv(paste0(JHU_path,dt,".csv"))
})  |>
  mutate(Date = as.Date(Date))
  

# specify locations of interest
# all 50 states plus DC
all_fips <- c(1, 2, 4:6, 8:13, 15:42, 44:51, 53:56)

# split data into state files with daily increase (e.g., incidence) included
temp <- map(all_fips, .f = function(fip){
  state_dat  <- format_JHU_data(all_data, fip) |>
    mutate(date = as.character(date))
  save_name  <- sprintf("data/state_%s.csv", unique(state_dat$state_abbr))
  write.csv(state_dat,save_name,row.names=FALSE)
})



