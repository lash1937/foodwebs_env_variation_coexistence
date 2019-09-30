## Code to run the intertidal simulation and partition coexistence

#### Source simulation functions ####

source("intertidal_model/intertidal_model_functions.R")
source("intertidal_model/intertidal_coexistence_partition.R")


#### Load packages ####

library(tidyverse)


#### Generate larval supply scenarios ####

larval_scenarios_for_full_run <- tribble (
  ~scenario, ~species, ~supply_level,
  1, "B", "high",
  1, "C", "high",
  1, "L", "high",
  1, "P", "high",
  2, "B", "low",
  2, "C", "low",
  2, "L", "low",
  2, "P", "low",
  3, "B", "high",
  3, "C", "low",
  3, "L", "low",
  3, "P", "low",
  4, "B", "low",
  4, "C", "high",
  4, "L", "low",
  4, "P", "low",
  5, "B", "low",
  5, "C", "low",
  5, "L", "high",
  5, "P", "low",
  6, "B", "low",
  6, "C", "low",
  6, "L", "low",
  6, "P", "high"
)

# Based on values in Forde & Doak (2004) Table 1
larval_supply_full <- tribble (
  ~species, ~mean_recruit, ~variance_recruit, ~mean_supply_level, ~variance_supply_level,
  "B", 90000, sqrt(4.6*10^9), "high", "low",
  "B", 90000, sqrt(3.24*10^10), "high", "med",
  "B", 90000, sqrt(9.9*10^10), "high", "high",
  
  "B", 50000, sqrt(1.41*10^10), "med", "low",
  "B", 50000, sqrt(1*10^10), "med", "med",
  "B", 50000, sqrt(3*10^10), "med", "high",
  
  "B", 6000, sqrt(2.025*10^7), "low", "low",
  "B", 6000, sqrt(1.44*10^8), "low", "med",
  "B", 6000, sqrt(4.41*10^8), "low", "high",
  
  "C", 70000, sqrt(2.75*10^9), "high", "low",
  "C", 70000, sqrt(1.96*10^10), "high", "med",
  "C", 70000, sqrt(6*10^10), "high", "high",
  
  "C", 30000, sqrt(5.1*10^8), "med", "low",
  "C", 30000, sqrt(3.6*10^9), "med", "med",
  "C", 30000, sqrt(1.1*10^10), "med", "high",
  
  "C", 6000, sqrt(2.025*10^7), "low", "low",
  "C", 6000, sqrt(1.44*10^8), "low", "med",
  "C", 6000, sqrt(4.41*10^8), "low", "high",
  
  "L", 3000, sqrt(3.8*10^6), "high", "low", 
  "L", 3000, sqrt(2.8*10^7), "high", "med", 
  "L", 3000, sqrt(8.1*10^7), "high", "high", 
  
  "L", 2400, sqrt(2.4*10^6), "med", "low",
  "L", 2400, sqrt(1.7*10^7), "med", "med",
  "L", 2400, sqrt(5.2*10^7), "med", "high",
  
  "L", 200, sqrt(16900), "low", "low",
  "L", 200, sqrt(1.2*10^5), "low", "med",
  "L", 200, sqrt(3.6*10^5), "low", "high",
  
  "P", 6873, sqrt(3.02*10^7), "high", "low",
  "P", 6873, sqrt(1.2*10^8), "high", "med",
  "P", 6873, sqrt(1.4*10^8), "high", "high",
  
  "P", 3800, sqrt(9.2*10^6), "med", "low",
  "P", 3800, sqrt(3.7*10^7), "med", "med",
  "P", 3800, sqrt(4.2*10^7), "med", "high",
  
  "P", 727, sqrt(3.4*10^5), "low", "low",
  "P", 727, sqrt(1.3*10^6), "low", "med",
  "P", 727, sqrt(1.5*10^6), "low", "high",
)

larval_scenarios_for_full_run %>%
  rename(mean_supply_level = "supply_level") %>%
  mutate(variance_supply_level = "low") %>%
  left_join(larval_supply_full) %>%
  select(-mean_supply_level, -variance_supply_level) %>%
  gather(key = variable, value = value, mean_recruit, variance_recruit) %>%
  unite(temp, species, variable) %>%
  spread(temp, value) -> larval_scenarios_input # table always need it to be named this


### Run simulation and partition ####

# Likely need to run in batches:
csv_names <- c("001_100", "101_200", "201_300",
               "301_400", "401_500")

# Set number of simulations
number_simulations <- 1
  # Was set at 100 for the manuscript, but set lower here for example
  # number_simulations <- 100

for (i in 1:5) {
  print(paste0("i = ", i))
  sim_output_list <- vector(mode = "list", length = 6)
  names(sim_output_list) <- c("high", "low", "balanus_high", "chthamalus_high",
                              "limpets_high", "pisaster_high")
  
  for (k in 1:length(sim_output_list)) {
    print(paste0("k = ", k))
    sim_output_list[[k]] <- do.larval.supply.simulation(k = k, n_sim = number_simulations)
  }
  
  sim_output_list %>%
    bind_rows(.id = "larval_scenario") %>%
    write_csv(path = paste0("final_larval_maintext_results_", csv_names[i], ".csv"))
}

#### Read in & summarize output ####

# modify as needed in your local environment:
output_files <- list.files()
sim_output_list_all <- lapply(X = output_files, FUN = read_csv)
names(sim_output_list_all) <- paste0("rep_", csv_names)

# get data frame of all results
sim_output_list_all %>%
  bind_rows(.id = "rep") %>%
  # Clarify unique ids
  mutate(rep_num = as.numeric(str_sub(rep, start = -7, end = -5)) - 1) %>%
  mutate(simulation_loop = simulation_loop + rep_num) %>%
  select(-rep_num, -rep) -> sim_output_df

# summarize output
sim_output_df %>%
  group_by(larval_scenario, coexistence_partition, species) %>%
  summarise(mean_cs = mean(coexistence_strength),
            sd_cs = sd(coexistence_strength),
            n_cs = n()) %>%
  mutate(se_cs = sd_cs/sqrt(n_cs))

