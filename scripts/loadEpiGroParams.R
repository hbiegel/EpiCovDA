loadEpiGroParams <- function(path = "./") {
  
  data <- read.csv(paste0(path,"param_data/state_parameters_case_counts_4-27.csv"))
  
  # states with at least 1000 cases before April 1st
  # MI,WA,NY,CA,NJ,PA,CT,FL,IL,GA,LA,OH,NC,MA,TN,CO,TX,WI,KY,MD,IN,SC,AZ
  early_state_rows = c(23,49,35,5,32,10,15,11,
                       20,44,6,39,7,36,28,
                       45,50,18,21,16,42,4)
  
  
  params_out <- data[early_state_rows, ] |>
    select(betas = beta_opt, 
           gammas = gamma_opt,
           Ns = N_opt)
  
  return(params_out)
}
