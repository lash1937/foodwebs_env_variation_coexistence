## Functions to partition coexistence of species in the intertidal model
# Multiple ways to consider coexistence are considered below, including:
# - consider the influence of variation in predator recruitment only
# - consider the influence of variation in adult predator abundance only
# - consider both the influence of variation in predator recruitment and adult abundance

# As well as code to run the coexistence partitioning across multiple levels
# of larval supply to barnacles, limpets, and sea stars.


#### Function to run the first step of coexistence partitioning ####

do.intertidal.rbar <- function(years_set = 100) {
  
  # 1. run to eqilibrium
  
  fd_results_1 <- do.intertidal.simulation(years_set = years_set)
  
  # 2. low density growth rate calculation
  
  # balanus
  # run without balanus adults or recruits -> equilibrium
  fd_results_ldr_b_absent <- do.intertidal.simulation(years_set = years_set, B_1 = 0)
  # set adult balanus population size very small (~1/1000 of usual size) and run using 
  # densities from last time step
  fd_results_ldr_b_invade <- 
    do.intertidal.simulation(years_set = years_set, 
                             B_1 = 1,
                             C_1 = fd_results_ldr_b_absent$chthamalus_dalli[nrow(fd_results_ldr_b_absent)],
                             L_1 = fd_results_ldr_b_absent$limpets[nrow(fd_results_ldr_b_absent)],
                             W_1 = fd_results_ldr_b_absent$whelks[nrow(fd_results_ldr_b_absent)],
                             P_1 = fd_results_ldr_b_absent$pisaster_ochraceus[nrow(fd_results_ldr_b_absent)],
                             total_1 = fd_results_ldr_b_absent$free_space[nrow(fd_results_ldr_b_absent)]
    )
  # calculate low density growth rates
  gr_b_invade <- do.growth.rates(results = fd_results_ldr_b_invade, col_nums = c(2:6))
  
  
  # chthamalus
  fd_results_ldr_c_absent <- do.intertidal.simulation(years_set = years_set, C_1 = 0)
  fd_results_ldr_c_invade <- 
    do.intertidal.simulation(years_set = years_set, 
                             B_1 = fd_results_ldr_c_absent$balanus_glandula[nrow(fd_results_ldr_c_absent)],
                             C_1 = 1,
                             L_1 = fd_results_ldr_c_absent$limpets[nrow(fd_results_ldr_c_absent)],
                             W_1 = fd_results_ldr_c_absent$whelks[nrow(fd_results_ldr_c_absent)],
                             P_1 = fd_results_ldr_c_absent$pisaster_ochraceus[nrow(fd_results_ldr_c_absent)],
                             total_1 = fd_results_ldr_c_absent$free_space[nrow(fd_results_ldr_c_absent)]
    )
  # low density growthrates
  gr_c_invade <- do.growth.rates(results = fd_results_ldr_c_invade, col_nums = c(2:6))
  
  
  # limpets
  fd_results_ldr_l_absent <- do.intertidal.simulation(years_set = years_set, L_1 = 0)
  fd_results_ldr_l_invade <- 
    do.intertidal.simulation(years_set = years_set, 
                             B_1 = fd_results_ldr_l_absent$balanus_glandula[nrow(fd_results_ldr_l_absent)],
                             C_1 = fd_results_ldr_l_absent$chthamalus_dalli[nrow(fd_results_ldr_l_absent)],
                             L_1 = 1,
                             W_1 = fd_results_ldr_l_absent$whelks[nrow(fd_results_ldr_l_absent)],
                             P_1 = fd_results_ldr_l_absent$pisaster_ochraceus[nrow(fd_results_ldr_l_absent)],
                             total_1 = fd_results_ldr_l_absent$free_space[nrow(fd_results_ldr_l_absent)]
    )
  # low density growthrates
  gr_l_invade <- do.growth.rates(results = fd_results_ldr_l_invade, col_nums = c(2:6))
  
  
  # whelks
  fd_results_ldr_w_absent <- do.intertidal.simulation(years_set = years_set, W_1 = 0)
  fd_results_ldr_w_invade <- 
    do.intertidal.simulation(years_set = years_set, 
                             B_1 = fd_results_ldr_w_absent$balanus_glandula[nrow(fd_results_ldr_w_absent)],
                             C_1 = fd_results_ldr_w_absent$chthamalus_dalli[nrow(fd_results_ldr_w_absent)],
                             L_1 = fd_results_ldr_w_absent$limpets[nrow(fd_results_ldr_w_absent)],
                             W_1 = 1,
                             P_1 = fd_results_ldr_w_absent$pisaster_ochraceus[nrow(fd_results_ldr_w_absent)],
                             total_1 = fd_results_ldr_w_absent$free_space[nrow(fd_results_ldr_w_absent)]
    )
  # low density growthrates
  gr_w_invade <- do.growth.rates(results = fd_results_ldr_w_invade, col_nums = c(2:6))
  
  
  # sea stars (pisaster)
  fd_results_ldr_p_absent <- do.intertidal.simulation(years_set = years_set, P_1 = 0)
  fd_results_ldr_p_invade <- 
    do.intertidal.simulation(years_set = years_set, 
                             B_1 = fd_results_ldr_p_absent$balanus_glandula[nrow(fd_results_ldr_p_absent)],
                             C_1 = fd_results_ldr_p_absent$chthamalus_dalli[nrow(fd_results_ldr_p_absent)],
                             L_1 = fd_results_ldr_p_absent$limpets[nrow(fd_results_ldr_p_absent)],
                             W_1 = fd_results_ldr_p_absent$whelks[nrow(fd_results_ldr_p_absent)],
                             P_1 = 1,
                             total_1 = fd_results_ldr_p_absent$free_space[nrow(fd_results_ldr_p_absent)]
    )
  # low density growthrates
  gr_p_invade <- do.growth.rates(results = fd_results_ldr_p_invade, col_nums = c(2:6))
  
  # get coexistence strengths (r_bar)
  # r_bar = mean low density growth rate of i - (1/2) (mean )
  r_bar_b <- mean(gr_b_invade[, "balanus_glandula"], na.rm=TRUE) - (1/2)*(mean(gr_b_invade[, "chthamalus_dalli"], na.rm=TRUE) + 
                                                                            mean(gr_b_invade[, "limpets"], na.rm=TRUE))
  r_bar_c <- mean(gr_c_invade[, "chthamalus_dalli"], na.rm=TRUE) - (1/2)*(mean(gr_c_invade[, "balanus_glandula"], na.rm=TRUE) +
                                                                            mean(gr_c_invade[, "limpets"], na.rm=TRUE))
  r_bar_l <- mean(gr_l_invade[, "limpets"], na.rm=TRUE) - (1/2)*(mean(gr_l_invade[, "chthamalus_dalli"], na.rm=TRUE) + 
                                                                   mean(gr_l_invade[, "balanus_glandula"], na.rm=TRUE))
  
  return(list(r_bar_result = tibble(r_bar = c(r_bar_b, r_bar_c, r_bar_l),
                                    species = c("balanus_glandula", "chthamalus_dalli", "limpets")),
              var_average = tibble(var_b_average = mean(fd_results_1$larvae.B),
                                   var_c_average = mean(fd_results_1$larvae.C),
                                   var_l_average = mean(fd_results_1$larvae.L),
                                   var_p_average = mean(fd_results_1$larvae.P)),
              # also return results in case want to look at each
              # results from overall run
              results_1 = fd_results_1,
              # results from absence/invasion of each species
              results_2_b_invade = bind_rows(fd_results_ldr_b_absent, fd_results_ldr_b_invade),
              results_2_c_invade = bind_rows(fd_results_ldr_c_absent, fd_results_ldr_c_invade),
              results_2_l_invade = bind_rows(fd_results_ldr_l_absent, fd_results_ldr_l_invade),
              results_2_p_invade = bind_rows(fd_results_ldr_p_absent, fd_results_ldr_p_invade)
              
  )
  )
}



#### Function to run coexistence partitioning just on the variation in recruitment rates ####

do.intertidal.all.average <- function(var_P_input, var_B_input, var_C_input, var_L_input) {
  
  # low density growth rate calculation
  
  # balanus
  # run without balanus adults or recruits -> equilibrium
  fd_results_ldr_b_absent <- do.intertidal.simulation(years_set = 100, B_1 = 0,
                                                      var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input)
  # set adult balanus population size very small (~1/1000 of usual size) and run using 
  # densities from last time step
  fd_results_ldr_b_invade <- 
    do.intertidal.simulation(years_set = 100, 
                             B_1 = 1,
                             C_1 = fd_results_ldr_b_absent$chthamalus_dalli[nrow(fd_results_ldr_b_absent)],
                             L_1 = fd_results_ldr_b_absent$limpets[nrow(fd_results_ldr_b_absent)],
                             W_1 = fd_results_ldr_b_absent$whelks[nrow(fd_results_ldr_b_absent)],
                             P_1 = fd_results_ldr_b_absent$pisaster_ochraceus[nrow(fd_results_ldr_b_absent)],
                             total_1 = fd_results_ldr_b_absent$free_space[nrow(fd_results_ldr_b_absent)],
                             var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input
    )
  # calculate low density growth rates
  gr_b_invade <- do.growth.rates(results = fd_results_ldr_b_invade, col_nums = c(2:6))
  
  
  # chthamalus
  fd_results_ldr_c_absent <- do.intertidal.simulation(years_set = 100, C_1 = 0,
                                                      var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input)
  fd_results_ldr_c_invade <- 
    do.intertidal.simulation(years_set = 100, 
                             B_1 = fd_results_ldr_c_absent$balanus_glandula[nrow(fd_results_ldr_c_absent)],
                             C_1 = 1,
                             L_1 = fd_results_ldr_c_absent$limpets[nrow(fd_results_ldr_c_absent)],
                             W_1 = fd_results_ldr_c_absent$whelks[nrow(fd_results_ldr_c_absent)],
                             P_1 = fd_results_ldr_c_absent$pisaster_ochraceus[nrow(fd_results_ldr_c_absent)],
                             total_1 = fd_results_ldr_c_absent$free_space[nrow(fd_results_ldr_c_absent)],
                             var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input
    )
  # low density growthrates
  gr_c_invade <- do.growth.rates(results = fd_results_ldr_c_invade, col_nums = c(2:6))
  
  
  # limpets
  fd_results_ldr_l_absent <- do.intertidal.simulation(years_set = 100, L_1 = 0,
                                                      var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input)
  fd_results_ldr_l_invade <- 
    do.intertidal.simulation(years_set = 100, 
                             B_1 = fd_results_ldr_l_absent$balanus_glandula[nrow(fd_results_ldr_l_absent)],
                             C_1 = fd_results_ldr_l_absent$chthamalus_dalli[nrow(fd_results_ldr_l_absent)],
                             L_1 = 1,
                             W_1 = fd_results_ldr_l_absent$whelks[nrow(fd_results_ldr_l_absent)],
                             P_1 = fd_results_ldr_l_absent$pisaster_ochraceus[nrow(fd_results_ldr_l_absent)],
                             total_1 = fd_results_ldr_l_absent$free_space[nrow(fd_results_ldr_l_absent)],
                             var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input
    )
  # low density growthrates
  gr_l_invade <- do.growth.rates(results = fd_results_ldr_l_invade, col_nums = c(2:6))
  
  
  # whelks
  fd_results_ldr_w_absent <- do.intertidal.simulation(years_set = 100, W_1 = 0,
                                                      var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input)
  fd_results_ldr_w_invade <- 
    do.intertidal.simulation(years_set = 100, 
                             B_1 = fd_results_ldr_w_absent$balanus_glandula[nrow(fd_results_ldr_w_absent)],
                             C_1 = fd_results_ldr_w_absent$chthamalus_dalli[nrow(fd_results_ldr_w_absent)],
                             L_1 = fd_results_ldr_w_absent$limpets[nrow(fd_results_ldr_w_absent)],
                             W_1 = 1,
                             P_1 = fd_results_ldr_w_absent$pisaster_ochraceus[nrow(fd_results_ldr_w_absent)],
                             total_1 = fd_results_ldr_w_absent$free_space[nrow(fd_results_ldr_w_absent)],
                             var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input
    )
  # low density growthrates
  gr_w_invade <- do.growth.rates(results = fd_results_ldr_w_invade, col_nums = c(2:6))
  
  
  # sea stars (pisaster)
  fd_results_ldr_p_absent <- do.intertidal.simulation(years_set = 100, P_1 = 0,
                                                      var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input)
  fd_results_ldr_p_invade <- 
    do.intertidal.simulation(years_set = 100, 
                             B_1 = fd_results_ldr_p_absent$balanus_glandula[nrow(fd_results_ldr_p_absent)],
                             C_1 = fd_results_ldr_p_absent$chthamalus_dalli[nrow(fd_results_ldr_p_absent)],
                             L_1 = fd_results_ldr_p_absent$limpets[nrow(fd_results_ldr_p_absent)],
                             W_1 = fd_results_ldr_p_absent$whelks[nrow(fd_results_ldr_p_absent)],
                             P_1 = 1,
                             total_1 = fd_results_ldr_p_absent$free_space[nrow(fd_results_ldr_p_absent)],
                             var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input
    )
  # low density growthrates
  gr_p_invade <- do.growth.rates(results = fd_results_ldr_p_invade, col_nums = c(2:6))
  
  # get coexistence strengths (r_bar)
  # r_bar = mean low density growth rate of i - (1/2) (mean )
  r_bar_b <- mean(gr_b_invade[, "balanus_glandula"], na.rm=TRUE) - (1/2)*(mean(gr_b_invade[, "chthamalus_dalli"], na.rm=TRUE) + 
                                                                            mean(gr_b_invade[, "limpets"], na.rm=TRUE))
  r_bar_c <- mean(gr_c_invade[, "chthamalus_dalli"], na.rm=TRUE) - (1/2)*(mean(gr_c_invade[, "balanus_glandula"], na.rm=TRUE) +
                                                                            mean(gr_c_invade[, "limpets"], na.rm=TRUE))
  r_bar_l <- mean(gr_l_invade[, "limpets"], na.rm=TRUE) - (1/2)*(mean(gr_l_invade[, "chthamalus_dalli"], na.rm=TRUE) + 
                                                                   mean(gr_l_invade[, "balanus_glandula"], na.rm=TRUE))
  
  return(list(r_bar_result = tibble(r_bar = c(r_bar_b, r_bar_c, r_bar_l),
                                    species = c("balanus_glandula", "chthamalus_dalli", "limpets")),
              var_average = tibble(var_b_average = mean(fd_results_1$larvae.B),
                                   var_c_average = mean(fd_results_1$larvae.C),
                                   var_l_average = mean(fd_results_1$larvae.L),
                                   var_p_average = mean(fd_results_1$larvae.P)),
              # also return results in case want to look at each
              # results from overall run
              results_1 = fd_results_1,
              # results from absence/invasion of each species
              results_2_b_invade = bind_rows(fd_results_ldr_b_absent, fd_results_ldr_b_invade),
              results_2_c_invade = bind_rows(fd_results_ldr_c_absent, fd_results_ldr_c_invade),
              results_2_l_invade = bind_rows(fd_results_ldr_l_absent, fd_results_ldr_l_invade),
              results_2_p_invade = bind_rows(fd_results_ldr_p_absent, fd_results_ldr_p_invade)
              
  )
  )
}





#### Function to run coexistence partitioning on +/- variation in predator abundance overall ####

do.intertidal.predator.removal <- function(var_P_input, var_B_input, var_C_input, var_L_input,
                                           P_avg_input, W_avg_input, years_set = 50) {
  
  # low density growth rate calculation
  
  # balanus
  # run without balanus adults or recruits -> equilibrium
  fd_results_ldr_b_absent <- do.intertidal.simulation(years_set = years_set, B_1 = 0,
                                                      var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                                                      P_avg = P_avg_input, W_avg = W_avg_input)
  # set adult balanus population size very small (~1/1000 of usual size) and run using 
  # densities from last time step
  fd_results_ldr_b_invade <- 
    do.intertidal.simulation(years_set = years_set, 
                             B_1 = 1,
                             C_1 = fd_results_ldr_b_absent$chthamalus_dalli[nrow(fd_results_ldr_b_absent)],
                             L_1 = fd_results_ldr_b_absent$limpets[nrow(fd_results_ldr_b_absent)],
                             W_1 = fd_results_ldr_b_absent$whelks[nrow(fd_results_ldr_b_absent)],
                             P_1 = fd_results_ldr_b_absent$pisaster_ochraceus[nrow(fd_results_ldr_b_absent)],
                             total_1 = fd_results_ldr_b_absent$free_space[nrow(fd_results_ldr_b_absent)],
                             var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                             P_avg = P_avg_input, W_avg = W_avg_input
    )
  # calculate low density growth rates
  gr_b_invade <- do.growth.rates(results = fd_results_ldr_b_invade, col_nums = c(2:6))
  
  
  # chthamalus
  fd_results_ldr_c_absent <- do.intertidal.simulation(years_set = years_set, C_1 = 0,
                                                      var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                                                      P_avg = P_avg_input, W_avg = W_avg_input)
  fd_results_ldr_c_invade <- 
    do.intertidal.simulation(years_set = years_set, 
                             B_1 = fd_results_ldr_c_absent$balanus_glandula[nrow(fd_results_ldr_c_absent)],
                             C_1 = 1,
                             L_1 = fd_results_ldr_c_absent$limpets[nrow(fd_results_ldr_c_absent)],
                             W_1 = fd_results_ldr_c_absent$whelks[nrow(fd_results_ldr_c_absent)],
                             P_1 = fd_results_ldr_c_absent$pisaster_ochraceus[nrow(fd_results_ldr_c_absent)],
                             total_1 = fd_results_ldr_c_absent$free_space[nrow(fd_results_ldr_c_absent)],
                             var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                             P_avg = P_avg_input, W_avg = W_avg_input
    )
  # low density growthrates
  gr_c_invade <- do.growth.rates(results = fd_results_ldr_c_invade, col_nums = c(2:6))
  
  
  # limpets
  fd_results_ldr_l_absent <- do.intertidal.simulation(years_set = years_set, L_1 = 0,
                                                      var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                                                      P_avg = P_avg_input, W_avg = W_avg_input)
  fd_results_ldr_l_invade <- 
    do.intertidal.simulation(years_set = years_set, 
                             B_1 = fd_results_ldr_l_absent$balanus_glandula[nrow(fd_results_ldr_l_absent)],
                             C_1 = fd_results_ldr_l_absent$chthamalus_dalli[nrow(fd_results_ldr_l_absent)],
                             L_1 = 1,
                             W_1 = fd_results_ldr_l_absent$whelks[nrow(fd_results_ldr_l_absent)],
                             P_1 = fd_results_ldr_l_absent$pisaster_ochraceus[nrow(fd_results_ldr_l_absent)],
                             total_1 = fd_results_ldr_l_absent$free_space[nrow(fd_results_ldr_l_absent)],
                             var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                             P_avg = P_avg_input, W_avg = W_avg_input
    )
  # low density growthrates
  gr_l_invade <- do.growth.rates(results = fd_results_ldr_l_invade, col_nums = c(2:6))
  
  
  # whelks
  fd_results_ldr_w_absent <- do.intertidal.simulation(years_set = years_set, W_1 = 0,
                                                      var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                                                      P_avg = P_avg_input, W_avg = W_avg_input)
  fd_results_ldr_w_invade <- 
    do.intertidal.simulation(years_set = years_set, 
                             B_1 = fd_results_ldr_w_absent$balanus_glandula[nrow(fd_results_ldr_w_absent)],
                             C_1 = fd_results_ldr_w_absent$chthamalus_dalli[nrow(fd_results_ldr_w_absent)],
                             L_1 = fd_results_ldr_w_absent$limpets[nrow(fd_results_ldr_w_absent)],
                             W_1 = 1,
                             P_1 = fd_results_ldr_w_absent$pisaster_ochraceus[nrow(fd_results_ldr_w_absent)],
                             total_1 = fd_results_ldr_w_absent$free_space[nrow(fd_results_ldr_w_absent)],
                             var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                             P_avg = P_avg_input, W_avg = W_avg_input
    )
  # low density growthrates
  gr_w_invade <- do.growth.rates(results = fd_results_ldr_w_invade, col_nums = c(2:6))
  
  
  # sea stars (pisaster)
  fd_results_ldr_p_absent <- do.intertidal.simulation(years_set = years_set, P_1 = 0,
                                                      var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                                                      P_avg = P_avg_input, W_avg = W_avg_input)
  fd_results_ldr_p_invade <- 
    do.intertidal.simulation(years_set = years_set, 
                             B_1 = fd_results_ldr_p_absent$balanus_glandula[nrow(fd_results_ldr_p_absent)],
                             C_1 = fd_results_ldr_p_absent$chthamalus_dalli[nrow(fd_results_ldr_p_absent)],
                             L_1 = fd_results_ldr_p_absent$limpets[nrow(fd_results_ldr_p_absent)],
                             W_1 = fd_results_ldr_p_absent$whelks[nrow(fd_results_ldr_p_absent)],
                             P_1 = 1,
                             total_1 = fd_results_ldr_p_absent$free_space[nrow(fd_results_ldr_p_absent)],
                             var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                             P_avg = P_avg_input, W_avg = W_avg_input
    )
  # low density growthrates
  gr_p_invade <- do.growth.rates(results = fd_results_ldr_p_invade, col_nums = c(2:6))
  
  # get coexistence strengths (r_bar)
  # r_bar = mean low density growth rate of i - (1/2) (mean )
  r_bar_b <- mean(gr_b_invade[, "balanus_glandula"], na.rm=TRUE) - (1/2)*(mean(gr_b_invade[, "chthamalus_dalli"], na.rm=TRUE) + 
                                                                            mean(gr_b_invade[, "limpets"], na.rm=TRUE))
  r_bar_c <- mean(gr_c_invade[, "chthamalus_dalli"], na.rm=TRUE) - (1/2)*(mean(gr_c_invade[, "balanus_glandula"], na.rm=TRUE) +
                                                                            mean(gr_c_invade[, "limpets"], na.rm=TRUE))
  r_bar_l <- mean(gr_l_invade[, "limpets"], na.rm=TRUE) - (1/2)*(mean(gr_l_invade[, "chthamalus_dalli"], na.rm=TRUE) + 
                                                                   mean(gr_l_invade[, "balanus_glandula"], na.rm=TRUE))
  
  return(list(r_bar_result = tibble(r_bar = c(r_bar_b, r_bar_c, r_bar_l),
                                    species = c("balanus_glandula", "chthamalus_dalli", "limpets")),
              var_average = tibble(var_b_average = mean(fd_results_1$larvae.B),
                                   var_c_average = mean(fd_results_1$larvae.C),
                                   var_l_average = mean(fd_results_1$larvae.L),
                                   var_p_average = mean(fd_results_1$larvae.P)),
              # also return results in case want to look at each
              # results from overall run
              results_1 = fd_results_1,
              # results from absence/invasion of each species
              results_2_b_invade = bind_rows(fd_results_ldr_b_absent, fd_results_ldr_b_invade),
              results_2_c_invade = bind_rows(fd_results_ldr_c_absent, fd_results_ldr_c_invade),
              results_2_l_invade = bind_rows(fd_results_ldr_l_absent, fd_results_ldr_l_invade),
              results_2_p_invade = bind_rows(fd_results_ldr_p_absent, fd_results_ldr_p_invade),
              results_2_w_invade = bind_rows(fd_results_ldr_w_absent, fd_results_ldr_w_invade)
              
  )
  )
}



#### Functions to run coexistence partitioning on +/- variation in predator abundance overall and allow variation in larval supply ####

do.intertidal.rbar.with.larval.var <- function(B.mean_input, B.stdev_input, 
                                               C.mean_input, C.stdev_input,
                                               L.mean_input, L.stdev_input, 
                                               P.mean_input, P.stdev_input) {
  
  # 1. run to eqilibrium
  
  fd_results_1 <- do.intertidal.simulation(years_set = 100,
                                           B.mean = B.mean_input, B.stdev = B.stdev_input,
                                           C.mean = C.mean_input, C.stdev = C.stdev_input,
                                           L.mean = L.mean_input, L.stdev = L.stdev_input,
                                           P.mean = P.mean_input, P.stdev = P.stdev_input)
  
  # 2. low density growth rate calculation
  
  # balanus
  # run without balanus adults or recruits -> equilibrium
  fd_results_ldr_b_absent <- do.intertidal.simulation(years_set = 50, B_1 = 0,
                                                      B.mean = B.mean_input, B.stdev = B.stdev_input, 
                                                      C.mean = C.mean_input, C.stdev = C.stdev_input,
                                                      L.mean = L.mean_input, L.stdev = L.stdev_input, 
                                                      P.mean = P.mean_input, P.stdev = P.stdev_input)
  # set adult balanus population size very small (~1/1000 of usual size) and run using 
  # densities from last time step
  fd_results_ldr_b_invade <- 
    do.intertidal.simulation(years_set = 50, 
                             B_1 = 1,
                             C_1 = mean(fd_results_ldr_b_absent$chthamalus_dalli[(nrow(fd_results_ldr_b_absent)-5):nrow(fd_results_ldr_b_absent)]),
                             L_1 = mean(fd_results_ldr_b_absent$limpets[(nrow(fd_results_ldr_b_absent)-5):nrow(fd_results_ldr_b_absent)]),
                             W_1 = mean(fd_results_ldr_b_absent$whelks[(nrow(fd_results_ldr_b_absent)-5):nrow(fd_results_ldr_b_absent)]),
                             P_1 = mean(fd_results_ldr_b_absent$pisaster_ochraceus[(nrow(fd_results_ldr_b_absent)-5):nrow(fd_results_ldr_b_absent)]),
                             total_1 = mean(fd_results_ldr_b_absent$free_space[(nrow(fd_results_ldr_b_absent)-5):nrow(fd_results_ldr_b_absent)]),
                             B.mean = B.mean_input, B.stdev = B.stdev_input, 
                             C.mean = C.mean_input, C.stdev = C.stdev_input,
                             L.mean = L.mean_input, L.stdev = L.stdev_input, 
                             P.mean = P.mean_input, P.stdev = P.stdev_input
    )
  # calculate low density growth rates
  gr_b_invade <- do.growth.rates(results = fd_results_ldr_b_invade, col_nums = c(2:6))
  
  
  # chthamalus
  fd_results_ldr_c_absent <- do.intertidal.simulation(years_set = 50, C_1 = 0,
                                                      B.mean = B.mean_input, B.stdev = B.stdev_input, 
                                                      C.mean = C.mean_input, C.stdev = C.stdev_input,
                                                      L.mean = L.mean_input, L.stdev = L.stdev_input, 
                                                      P.mean = P.mean_input, P.stdev = P.stdev_input)
  fd_results_ldr_c_invade <- 
    do.intertidal.simulation(years_set = 50, 
                             B_1 = mean(fd_results_ldr_c_absent$balanus_glandula[(nrow(fd_results_ldr_c_absent)-5):nrow(fd_results_ldr_c_absent)]),
                             C_1 = 1,
                             L_1 = mean(fd_results_ldr_c_absent$limpets[(nrow(fd_results_ldr_c_absent)-5):nrow(fd_results_ldr_c_absent)]),
                             W_1 = mean(fd_results_ldr_c_absent$whelks[(nrow(fd_results_ldr_c_absent)-5):nrow(fd_results_ldr_c_absent)]),
                             P_1 = mean(fd_results_ldr_c_absent$pisaster_ochraceus[(nrow(fd_results_ldr_c_absent)-5):nrow(fd_results_ldr_c_absent)]),
                             total_1 = mean(fd_results_ldr_c_absent$free_space[(nrow(fd_results_ldr_c_absent)-5):nrow(fd_results_ldr_c_absent)]),
                             B.mean = B.mean_input, B.stdev = B.stdev_input, 
                             C.mean = C.mean_input, C.stdev = C.stdev_input,
                             L.mean = L.mean_input, L.stdev = L.stdev_input, 
                             P.mean = P.mean_input, P.stdev = P.stdev_input
    )
  # low density growthrates
  gr_c_invade <- do.growth.rates(results = fd_results_ldr_c_invade, col_nums = c(2:6))
  
  
  # limpets
  fd_results_ldr_l_absent <- do.intertidal.simulation(years_set = 50, L_1 = 0,
                                                      B.mean = B.mean_input, B.stdev = B.stdev_input, 
                                                      C.mean = C.mean_input, C.stdev = C.stdev_input,
                                                      L.mean = L.mean_input, L.stdev = L.stdev_input, 
                                                      P.mean = P.mean_input, P.stdev = P.stdev_input)
  fd_results_ldr_l_invade <- 
    do.intertidal.simulation(years_set = 50, 
                             B_1 = mean(fd_results_ldr_l_absent$balanus_glandula[(nrow(fd_results_ldr_l_absent)-5):nrow(fd_results_ldr_l_absent)]),
                             C_1 = mean(fd_results_ldr_l_absent$chthamalus_dalli[(nrow(fd_results_ldr_l_absent)-5):nrow(fd_results_ldr_l_absent)]),
                             L_1 = 1,
                             W_1 = mean(fd_results_ldr_l_absent$whelks[(nrow(fd_results_ldr_l_absent)-5):nrow(fd_results_ldr_l_absent)]),
                             P_1 = mean(fd_results_ldr_l_absent$pisaster_ochraceus[(nrow(fd_results_ldr_l_absent)-5):nrow(fd_results_ldr_l_absent)]),
                             total_1 = mean(fd_results_ldr_l_absent$free_space[(nrow(fd_results_ldr_l_absent)-5):nrow(fd_results_ldr_l_absent)]),
                             B.mean = B.mean_input, B.stdev = B.stdev_input, 
                             C.mean = C.mean_input, C.stdev = C.stdev_input,
                             L.mean = L.mean_input, L.stdev = L.stdev_input, 
                             P.mean = P.mean_input, P.stdev = P.stdev_input
    )
  # low density growthrates
  gr_l_invade <- do.growth.rates(results = fd_results_ldr_l_invade, col_nums = c(2:6))
  
  
  # whelks
  fd_results_ldr_w_absent <- do.intertidal.simulation(years_set = 50, W_1 = 0,
                                                      B.mean = B.mean_input, B.stdev = B.stdev_input, 
                                                      C.mean = C.mean_input, C.stdev = C.stdev_input,
                                                      L.mean = L.mean_input, L.stdev = L.stdev_input, 
                                                      P.mean = P.mean_input, P.stdev = P.stdev_input)
  fd_results_ldr_w_invade <- 
    do.intertidal.simulation(years_set = 50, 
                             B_1 = mean(fd_results_ldr_w_absent$balanus_glandula[(nrow(fd_results_ldr_w_absent)-5):nrow(fd_results_ldr_w_absent)]),
                             C_1 = mean(fd_results_ldr_w_absent$chthamalus_dalli[(nrow(fd_results_ldr_w_absent)-5):nrow(fd_results_ldr_w_absent)]),
                             L_1 = mean(fd_results_ldr_w_absent$limpets[(nrow(fd_results_ldr_w_absent)-5):nrow(fd_results_ldr_w_absent)]),
                             W_1 = 1,
                             P_1 = mean(fd_results_ldr_w_absent$pisaster_ochraceus[(nrow(fd_results_ldr_w_absent)-5):nrow(fd_results_ldr_w_absent)]),
                             total_1 = mean(fd_results_ldr_w_absent$free_space[(nrow(fd_results_ldr_w_absent)-5):nrow(fd_results_ldr_w_absent)]),
                             B.mean = B.mean_input, B.stdev = B.stdev_input, 
                             C.mean = C.mean_input, C.stdev = C.stdev_input,
                             L.mean = L.mean_input, L.stdev = L.stdev_input, 
                             P.mean = P.mean_input, P.stdev = P.stdev_input
    )
  # low density growthrates
  gr_w_invade <- do.growth.rates(results = fd_results_ldr_w_invade, col_nums = c(2:6))
  
  
  # sea stars (pisaster)
  fd_results_ldr_p_absent <- do.intertidal.simulation(years_set = 50, P_1 = 0,
                                                      B.mean = B.mean_input, B.stdev = B.stdev_input, 
                                                      C.mean = C.mean_input, C.stdev = C.stdev_input,
                                                      L.mean = L.mean_input, L.stdev = L.stdev_input, 
                                                      P.mean = P.mean_input, P.stdev = P.stdev_input)
  fd_results_ldr_p_invade <- 
    do.intertidal.simulation(years_set = 50, 
                             B_1 = mean(fd_results_ldr_p_absent$balanus_glandula[(nrow(fd_results_ldr_p_absent)-5):nrow(fd_results_ldr_p_absent)]),
                             C_1 = mean(fd_results_ldr_p_absent$chthamalus_dalli[(nrow(fd_results_ldr_p_absent)-5):nrow(fd_results_ldr_p_absent)]),
                             L_1 = mean(fd_results_ldr_p_absent$limpets[(nrow(fd_results_ldr_p_absent)-5):nrow(fd_results_ldr_p_absent)]),
                             W_1 = mean(fd_results_ldr_p_absent$whelks[(nrow(fd_results_ldr_p_absent)-5):nrow(fd_results_ldr_p_absent)]),
                             P_1 = 1,
                             total_1 = mean(fd_results_ldr_p_absent$free_space[(nrow(fd_results_ldr_p_absent)-5):nrow(fd_results_ldr_p_absent)]),
                             B.mean = B.mean_input, B.stdev = B.stdev_input, 
                             C.mean = C.mean_input, C.stdev = C.stdev_input,
                             L.mean = L.mean_input, L.stdev = L.stdev_input, 
                             P.mean = P.mean_input, P.stdev = P.stdev_input
    )
  # low density growthrates
  gr_p_invade <- do.growth.rates(results = fd_results_ldr_p_invade, col_nums = c(2:6))
  
  # get coexistence strengths (r_bar)
  # r_bar = mean low density growth rate of i - (1/2) (mean )
  r_bar_b <- mean(gr_b_invade[, "balanus_glandula"], na.rm=TRUE) - (1/2)*(mean(gr_b_invade[, "chthamalus_dalli"], na.rm=TRUE) + 
                                                                            mean(gr_b_invade[, "limpets"], na.rm=TRUE))
  r_bar_c <- mean(gr_c_invade[, "chthamalus_dalli"], na.rm=TRUE) - (1/2)*(mean(gr_c_invade[, "balanus_glandula"], na.rm=TRUE) +
                                                                            mean(gr_c_invade[, "limpets"], na.rm=TRUE))
  r_bar_l <- mean(gr_l_invade[, "limpets"], na.rm=TRUE) - (1/2)*(mean(gr_l_invade[, "chthamalus_dalli"], na.rm=TRUE) + 
                                                                   mean(gr_l_invade[, "balanus_glandula"], na.rm=TRUE))
  
  return(list(r_bar_result = tibble(r_bar = c(r_bar_b, r_bar_c, r_bar_l),
                                    species = c("balanus_glandula", "chthamalus_dalli", "limpets")),
              var_average = tibble(var_b_average = mean(fd_results_1$larvae.B),
                                   var_c_average = mean(fd_results_1$larvae.C),
                                   var_l_average = mean(fd_results_1$larvae.L),
                                   var_p_average = mean(fd_results_1$larvae.P)),
              # also return results in case want to look at each
              # results from overall run
              results_1 = fd_results_1,
              ## results from absence/invasion of each species
              results_2_b_invade = bind_rows(fd_results_ldr_b_absent, fd_results_ldr_b_invade),
              results_2_c_invade = bind_rows(fd_results_ldr_c_absent, fd_results_ldr_c_invade),
              results_2_l_invade = bind_rows(fd_results_ldr_l_absent, fd_results_ldr_l_invade),
              results_2_p_invade = bind_rows(fd_results_ldr_p_absent, fd_results_ldr_p_invade),
              results_2_w_invade = bind_rows(fd_results_ldr_w_absent, fd_results_ldr_w_invade)
              
  )
  )
}



do.intertidal.predator.removal.with.larval.var <- function(var_P_input, var_B_input, var_C_input, 
                                                           var_L_input, P_avg_input, W_avg_input,
                                                           B.mean_input, B.stdev_input, C.mean_input, C.stdev_input,
                                                           L.mean_input, L.stdev_input, P.mean_input, P.stdev_input) {
  
  # 1. run to eqilibrium
  
  fd_results_1 <- do.intertidal.simulation(years_set = 100,
                                           var_P = var_P_input, var_B = var_B_input, 
                                           var_C = var_C_input, var_L = var_L_input,
                                           P_avg = P_avg_input, W_avg = W_avg_input,
                                           B.mean = B.mean_input, B.stdev = B.stdev_input,
                                           C.mean = C.mean_input, C.stdev = C.stdev_input,
                                           L.mean = L.mean_input, L.stdev = L.stdev_input,
                                           P.mean = P.mean_input, P.stdev = P.stdev_input)
  
  # 2. low density growth rate calculation
  
  # balanus
  # run without balanus adults or recruits -> equilibrium
  fd_results_ldr_b_absent <- do.intertidal.simulation(years_set = 50, B_1 = 0,
                                                      var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                                                      P_avg = P_avg_input, W_avg = W_avg_input,
                                                      B.mean = B.mean_input, B.stdev = B.stdev_input, 
                                                      C.mean = C.mean_input, C.stdev = C.stdev_input,
                                                      L.mean = L.mean_input, L.stdev = L.stdev_input, 
                                                      P.mean = P.mean_input, P.stdev = P.stdev_input)
  # set adult balanus population size very small (~1/1000 of usual size) and run using 
  # densities from last time step
  fd_results_ldr_b_invade <- 
    do.intertidal.simulation(years_set = 50, 
                             B_1 = 1,
                             C_1 = mean(fd_results_ldr_b_absent$chthamalus_dalli[(nrow(fd_results_ldr_b_absent)-5):nrow(fd_results_ldr_b_absent)]),
                             L_1 = mean(fd_results_ldr_b_absent$limpets[(nrow(fd_results_ldr_b_absent)-5):nrow(fd_results_ldr_b_absent)]),
                             W_1 = mean(fd_results_ldr_b_absent$whelks[(nrow(fd_results_ldr_b_absent)-5):nrow(fd_results_ldr_b_absent)]),
                             P_1 = mean(fd_results_ldr_b_absent$pisaster_ochraceus[(nrow(fd_results_ldr_b_absent)-5):nrow(fd_results_ldr_b_absent)]),
                             total_1 = mean(fd_results_ldr_b_absent$free_space[(nrow(fd_results_ldr_b_absent)-5):nrow(fd_results_ldr_b_absent)]),
                             var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                             P_avg = P_avg_input, W_avg = W_avg_input,
                             B.mean = B.mean_input, B.stdev = B.stdev_input, 
                             C.mean = C.mean_input, C.stdev = C.stdev_input,
                             L.mean = L.mean_input, L.stdev = L.stdev_input, 
                             P.mean = P.mean_input, P.stdev = P.stdev_input
    )
  # calculate low density growth rates
  gr_b_invade <- do.growth.rates(results = fd_results_ldr_b_invade, col_nums = c(2:6))
  
  
  # chthamalus
  fd_results_ldr_c_absent <- do.intertidal.simulation(years_set = 50, C_1 = 0,
                                                      var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                                                      P_avg = P_avg_input, W_avg = W_avg_input,
                                                      B.mean = B.mean_input, B.stdev = B.stdev_input, 
                                                      C.mean = C.mean_input, C.stdev = C.stdev_input,
                                                      L.mean = L.mean_input, L.stdev = L.stdev_input, 
                                                      P.mean = P.mean_input, P.stdev = P.stdev_input)
  fd_results_ldr_c_invade <- 
    do.intertidal.simulation(years_set = 50, 
                             B_1 = mean(fd_results_ldr_c_absent$balanus_glandula[(nrow(fd_results_ldr_c_absent)-5):nrow(fd_results_ldr_c_absent)]),
                             C_1 = 1,
                             L_1 = mean(fd_results_ldr_c_absent$limpets[(nrow(fd_results_ldr_c_absent)-5):nrow(fd_results_ldr_c_absent)]),
                             W_1 = mean(fd_results_ldr_c_absent$whelks[(nrow(fd_results_ldr_c_absent)-5):nrow(fd_results_ldr_c_absent)]),
                             P_1 = mean(fd_results_ldr_c_absent$pisaster_ochraceus[(nrow(fd_results_ldr_c_absent)-5):nrow(fd_results_ldr_c_absent)]),
                             total_1 = mean(fd_results_ldr_c_absent$free_space[(nrow(fd_results_ldr_c_absent)-5):nrow(fd_results_ldr_c_absent)]),
                             var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                             P_avg = P_avg_input, W_avg = W_avg_input,
                             B.mean = B.mean_input, B.stdev = B.stdev_input, 
                             C.mean = C.mean_input, C.stdev = C.stdev_input,
                             L.mean = L.mean_input, L.stdev = L.stdev_input, 
                             P.mean = P.mean_input, P.stdev = P.stdev_input
    )
  # low density growthrates
  gr_c_invade <- do.growth.rates(results = fd_results_ldr_c_invade, col_nums = c(2:6))
  
  
  # limpets
  fd_results_ldr_l_absent <- do.intertidal.simulation(years_set = 50, L_1 = 0,
                                                      var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                                                      P_avg = P_avg_input, W_avg = W_avg_input,
                                                      B.mean = B.mean_input, B.stdev = B.stdev_input, 
                                                      C.mean = C.mean_input, C.stdev = C.stdev_input,
                                                      L.mean = L.mean_input, L.stdev = L.stdev_input, 
                                                      P.mean = P.mean_input, P.stdev = P.stdev_input)
  fd_results_ldr_l_invade <- 
    do.intertidal.simulation(years_set = 50, 
                             B_1 = mean(fd_results_ldr_l_absent$balanus_glandula[(nrow(fd_results_ldr_l_absent)-5):nrow(fd_results_ldr_l_absent)]),
                             C_1 = mean(fd_results_ldr_l_absent$chthamalus_dalli[(nrow(fd_results_ldr_l_absent)-5):nrow(fd_results_ldr_l_absent)]),
                             L_1 = 1,
                             W_1 = mean(fd_results_ldr_l_absent$whelks[(nrow(fd_results_ldr_l_absent)-5):nrow(fd_results_ldr_l_absent)]),
                             P_1 = mean(fd_results_ldr_l_absent$pisaster_ochraceus[(nrow(fd_results_ldr_l_absent)-5):nrow(fd_results_ldr_l_absent)]),
                             total_1 = mean(fd_results_ldr_l_absent$free_space[(nrow(fd_results_ldr_l_absent)-5):nrow(fd_results_ldr_l_absent)]),
                             var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                             P_avg = P_avg_input, W_avg = W_avg_input,
                             B.mean = B.mean_input, B.stdev = B.stdev_input, 
                             C.mean = C.mean_input, C.stdev = C.stdev_input,
                             L.mean = L.mean_input, L.stdev = L.stdev_input, 
                             P.mean = P.mean_input, P.stdev = P.stdev_input
    )
  # low density growthrates
  gr_l_invade <- do.growth.rates(results = fd_results_ldr_l_invade, col_nums = c(2:6))
  
  
  # whelks
  fd_results_ldr_w_absent <- do.intertidal.simulation(years_set = 50, W_1 = 0,
                                                      var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                                                      P_avg = P_avg_input, W_avg = W_avg_input,
                                                      B.mean = B.mean_input, B.stdev = B.stdev_input, 
                                                      C.mean = C.mean_input, C.stdev = C.stdev_input,
                                                      L.mean = L.mean_input, L.stdev = L.stdev_input, 
                                                      P.mean = P.mean_input, P.stdev = P.stdev_input)
  fd_results_ldr_w_invade <- 
    do.intertidal.simulation(years_set = 50, 
                             B_1 = mean(fd_results_ldr_w_absent$balanus_glandula[(nrow(fd_results_ldr_w_absent)-5):nrow(fd_results_ldr_w_absent)]),
                             C_1 = mean(fd_results_ldr_w_absent$chthamalus_dalli[(nrow(fd_results_ldr_w_absent)-5):nrow(fd_results_ldr_w_absent)]),
                             L_1 = mean(fd_results_ldr_w_absent$limpets[(nrow(fd_results_ldr_w_absent)-5):nrow(fd_results_ldr_w_absent)]),
                             W_1 = 1,
                             P_1 = mean(fd_results_ldr_w_absent$pisaster_ochraceus[(nrow(fd_results_ldr_w_absent)-5):nrow(fd_results_ldr_w_absent)]),
                             total_1 = mean(fd_results_ldr_w_absent$free_space[(nrow(fd_results_ldr_w_absent)-5):nrow(fd_results_ldr_w_absent)]),
                             var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                             P_avg = P_avg_input, W_avg = W_avg_input,
                             B.mean = B.mean_input, B.stdev = B.stdev_input, 
                             C.mean = C.mean_input, C.stdev = C.stdev_input,
                             L.mean = L.mean_input, L.stdev = L.stdev_input, 
                             P.mean = P.mean_input, P.stdev = P.stdev_input
    )
  # low density growthrates
  gr_w_invade <- do.growth.rates(results = fd_results_ldr_w_invade, col_nums = c(2:6))
  
  
  # sea stars (pisaster)
  fd_results_ldr_p_absent <- do.intertidal.simulation(years_set = 50, P_1 = 0,
                                                      var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                                                      P_avg = P_avg_input, W_avg = W_avg_input,
                                                      B.mean = B.mean_input, B.stdev = B.stdev_input, 
                                                      C.mean = C.mean_input, C.stdev = C.stdev_input,
                                                      L.mean = L.mean_input, L.stdev = L.stdev_input, 
                                                      P.mean = P.mean_input, P.stdev = P.stdev_input)
  fd_results_ldr_p_invade <- 
    do.intertidal.simulation(years_set = 50, 
                             B_1 = mean(fd_results_ldr_p_absent$balanus_glandula[(nrow(fd_results_ldr_p_absent)-5):nrow(fd_results_ldr_p_absent)]),
                             C_1 = mean(fd_results_ldr_p_absent$chthamalus_dalli[(nrow(fd_results_ldr_p_absent)-5):nrow(fd_results_ldr_p_absent)]),
                             L_1 = mean(fd_results_ldr_p_absent$limpets[(nrow(fd_results_ldr_p_absent)-5):nrow(fd_results_ldr_p_absent)]),
                             W_1 = mean(fd_results_ldr_p_absent$whelks[(nrow(fd_results_ldr_p_absent)-5):nrow(fd_results_ldr_p_absent)]),
                             P_1 = 1,
                             total_1 = mean(fd_results_ldr_p_absent$free_space[(nrow(fd_results_ldr_p_absent)-5):nrow(fd_results_ldr_p_absent)]),
                             var_P = var_P_input, var_B = var_B_input, var_C = var_C_input, var_L = var_L_input,
                             P_avg = P_avg_input, W_avg = W_avg_input,
                             B.mean = B.mean_input, B.stdev = B.stdev_input, 
                             C.mean = C.mean_input, C.stdev = C.stdev_input,
                             L.mean = L.mean_input, L.stdev = L.stdev_input, 
                             P.mean = P.mean_input, P.stdev = P.stdev_input
    )
  # low density growthrates
  gr_p_invade <- do.growth.rates(results = fd_results_ldr_p_invade, col_nums = c(2:6))
  
  # get coexistence strengths (r_bar)
  # r_bar = mean low density growth rate of i - (1/2) (mean )
  r_bar_b <- mean(gr_b_invade[, "balanus_glandula"], na.rm=TRUE) - (1/2)*(mean(gr_b_invade[, "chthamalus_dalli"], na.rm=TRUE) + 
                                                                            mean(gr_b_invade[, "limpets"], na.rm=TRUE))
  r_bar_c <- mean(gr_c_invade[, "chthamalus_dalli"], na.rm=TRUE) - (1/2)*(mean(gr_c_invade[, "balanus_glandula"], na.rm=TRUE) +
                                                                            mean(gr_c_invade[, "limpets"], na.rm=TRUE))
  r_bar_l <- mean(gr_l_invade[, "limpets"], na.rm=TRUE) - (1/2)*(mean(gr_l_invade[, "chthamalus_dalli"], na.rm=TRUE) + 
                                                                   mean(gr_l_invade[, "balanus_glandula"], na.rm=TRUE))
  
  return(list(r_bar_result = tibble(r_bar = c(r_bar_b, r_bar_c, r_bar_l),
                                    species = c("balanus_glandula", "chthamalus_dalli", "limpets")),
              var_average = tibble(var_b_average = mean(fd_results_1$larvae.B),
                                   var_c_average = mean(fd_results_1$larvae.C),
                                   var_l_average = mean(fd_results_1$larvae.L),
                                   var_p_average = mean(fd_results_1$larvae.P)),
              # also return results in case want to look at each
              # results from overall run
              results_1 = fd_results_1,
              # results from absence/invasion of each species
              results_2_b_invade = bind_rows(fd_results_ldr_b_absent, fd_results_ldr_b_invade),
              results_2_c_invade = bind_rows(fd_results_ldr_c_absent, fd_results_ldr_c_invade),
              results_2_l_invade = bind_rows(fd_results_ldr_l_absent, fd_results_ldr_l_invade),
              results_2_p_invade = bind_rows(fd_results_ldr_p_absent, fd_results_ldr_p_invade),
              results_2_w_invade = bind_rows(fd_results_ldr_w_absent, fd_results_ldr_w_invade)
              
  )
  )
}




#### Function to run simulations across larval supply rates ####


do.larval.supply.simulation <- function(k, # indexing number for which larval scenario to use
                                        n_sim = 50) {
  
  simulation_loop_output_tmp <- vector(mode = "list", length = n_sim)
  
  for (i in 1:length(simulation_loop_output_tmp)) {
    
    print(i)
    
    B.mean_input <- larval_scenarios_input$B_mean_recruit[k]
    B.stdev_input <- larval_scenarios_input$B_variance_recruit[k]
    C.mean_input <- larval_scenarios_input$C_mean_recruit[k]
    C.stdev_input <- larval_scenarios_input$C_variance_recruit[k]
    L.mean_input <- larval_scenarios_input$L_mean_recruit[k]
    L.stdev_input <- larval_scenarios_input$L_variance_recruit[k]
    P.mean_input <- larval_scenarios_input$P_mean_recruit[k]
    P.stdev_input <- larval_scenarios_input$P_variance_recruit[k]
    
    # run model to equilibrium, get low density growth rates for invader/resident combinations
    fd_tmp_1 <- do.intertidal.rbar.with.larval.var(B.mean_input = B.mean_input, B.stdev_input = B.stdev_input, 
                                                   C.mean_input = C.mean_input, C.stdev_input = C.stdev_input,
                                                   L.mean_input = L.mean_input, L.stdev_input = L.stdev_input, 
                                                   P.mean_input = P.mean_input, P.stdev_input = P.stdev_input)
    
    # get long term averages:
    fd_tmp_1_var_C <- mean(fd_tmp_1$results_1$larvae.C)
    fd_tmp_1_var_B <- mean(fd_tmp_1$results_1$larvae.B)
    fd_tmp_1_var_L <- mean(fd_tmp_1$results_1$larvae.L)
    fd_tmp_1_var_P <- mean(fd_tmp_1$results_1$larvae.P)
    fd_tmp_1_P_avg <- mean(fd_tmp_1$results_1$pisaster_ochraceus)
    fd_tmp_1_W_avg <- mean(fd_tmp_1$results_1$whelks)
    
    # run model to equilibrium and get long term and low density growth rates
    fd_tmp_3a <- do.intertidal.predator.removal.with.larval.var(var_B_input = fd_tmp_1_var_B, 
                                                                var_C_input = fd_tmp_1_var_C, 
                                                                var_L_input = fd_tmp_1_var_L,
                                                                var_P_input = fd_tmp_1_var_P,
                                                                P_avg_input = fd_tmp_1_P_avg,
                                                                W_avg_input = fd_tmp_1_W_avg,
                                                                B.mean = B.mean_input, B.stdev = B.stdev_input, 
                                                                C.mean = C.mean_input, C.stdev = C.stdev_input,
                                                                L.mean = L.mean_input, L.stdev = L.stdev_input, 
                                                                P.mean = P.mean_input, P.stdev = P.stdev_input)
    
    # set variation in competitor recruitment to average (constant)
    fd_tmp_3b <- do.intertidal.predator.removal.with.larval.var(var_B_input = fd_tmp_1_var_B, 
                                                                var_C_input = fd_tmp_1_var_C, 
                                                                var_L_input = fd_tmp_1_var_L,
                                                                var_P_input = NULL,
                                                                P_avg_input = NULL,
                                                                W_avg_input = NULL,
                                                                B.mean = B.mean_input, B.stdev = B.stdev_input, 
                                                                C.mean = C.mean_input, C.stdev = C.stdev_input,
                                                                L.mean = L.mean_input, L.stdev = L.stdev_input, 
                                                                P.mean = P.mean_input, P.stdev = P.stdev_input)
    
    # set variation in predator ABUNDANCE to average (constant)
    fd_tmp_3c <- do.intertidal.predator.removal.with.larval.var(var_B_input = NULL, 
                                                                var_C_input = NULL, 
                                                                var_L_input = NULL,
                                                                var_P_input = fd_tmp_1_var_P,
                                                                P_avg_input = fd_tmp_1_P_avg,
                                                                W_avg_input = fd_tmp_1_W_avg,
                                                                B.mean = B.mean_input, B.stdev = B.stdev_input, 
                                                                C.mean = C.mean_input, C.stdev = C.stdev_input,
                                                                L.mean = L.mean_input, L.stdev = L.stdev_input, 
                                                                P.mean = P.mean_input, P.stdev = P.stdev_input)
    
    list(r_bar = fd_tmp_1$r_bar_result,
         delta_0 = fd_tmp_3a$r_bar_result,
         delta_p = fd_tmp_3b$r_bar_result,
         delta_c = fd_tmp_3c$r_bar_result) %>%
      bind_rows(.id = "id") -> fd_tmp_3_df
    
    fd_tmp_3_df  %>%
      spread(key = "id", value = r_bar) %>%
      mutate(delta_cp = r_bar - (delta_0 + delta_c + delta_p)) -> simulation_loop_output_tmp[[i]]
  }
  
  # check if any are null
  for(j in 1:length(simulation_loop_output_tmp)) {
    if (is.null(simulation_loop_output_tmp[[i]])) {
      simulation_loop_output_tmp[[i]] <- tibble(species = c("balanus_glandula", "chthamalus_dalli", "limpets"),
                                                delta_0 = rep(NA, 3),
                                                delta_c = rep(NA, 3),
                                                delta_p = rep(NA, 3),
                                                r_bar = rep(NA, 3),
                                                delta_cp = rep(NA, 3))
    }
  }
  
  simulation_loop_output_tmp %>%
    map(gather, delta_0:delta_cp, key = "coexistence_partition", value = "coexistence_strength") %>%
    map(mutate, coexistence_partition = factor(coexistence_partition, levels = c("r_bar", 
                                                                                 "delta_0", "delta_c", 
                                                                                 "delta_p", "delta_cp"))) %>%
    bind_rows(.id = "simulation_loop") %>%
    return()
  
}










