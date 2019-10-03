# Implementation of the classic diamond model.
# Equations and parameterizations from Vasseur and Fox (2007) Ecology Letters
# with temporal variation in mortality rates
# Implemented with tree species mortality rates through time (e.g. from Vasseur and Fox appendix)

library(deSolve)

# ----------------------------------------------------------------------------------------------------
# Function to run model with both species for overall dynamics

VassFox_Cvar <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dP = - (M_P * P) + ( ( (J_P * P) * ( (O_P_C1 * C1) + ( (1 - O_P_C1) * C2) ) ) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) ) 
    dC1 = - (M_C1 * C1) + ( (O_C1_R * J_C1 * C1 * R) / (R + R_0_1) ) - ( (O_P_C1 * J_P * P * C1) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
    dC2 = - (M_C2 * C2) + ( (O_C2_R * J_C2 * C2 * R) / (R + R_0_2) ) - ( ((1-O_P_C1) * J_P * P * C2) / ( (O_P_C1 * C1) + ((1 - O_P_C1) * C2) + C_0) )
    dR = r * R * (1 - (R / K)) - ( (O_C1_R * J_C1 * C1 * R) / (R + R_0_1) ) - ( (O_C2_R * J_C2 * C2 * R) / (R + R_0_2) )
    
    return(list(c(dP, dC1, dC2, dR)))
    
  })
}

# ----------------------------------------------------------------------------------------------------
# parameters
# resource intrinsic rate of growth
r = 1.0
# resource carrying capacity
K = 1.0
# consumer 1 ingestion rate
J_C1 = 0.8036
# consumer 2 ingestion rate
J_C2 = 0.7
# predator ingestion rate
J_P = 0.4
# predator mortality rate
M_P = 0.08
# half saturation constant
R_0_1 = 0.16129
R_0_2 = 0.9
C_0 = 0.5
# preference coefficient
O_P_C1 = 0.92
O_C1_R = 1.0
O_C2_R = 0.98

# strength of env. on mortality rate
sigma=0.55

# cross-correlation of C1 and C2
rho=0

# number of timesteps to run the model 
time  <- 5000 

# ----------------------------------------------------------------------------------------------------
# looping over multiple runs

runs <- 5 # example with 5 runs for computational efficiency. All figures were create from 500 runs

# set up matrices to hold results
C1_final_mechanisms <- matrix(data=NA, nrow=runs, ncol=5)
C2_final_mechanisms <- matrix(data=NA, nrow=runs, ncol=5)
colnames(C1_final_mechanisms ) <- c("C1_r_bar", "C1_delta_0", "C1_delta_P", "C1_delta_E", "C1_delta_EP")
colnames(C2_final_mechanisms ) <- c("C2_r_bar", "C2_delta_0", "C2_delta_P", "C2_delta_E", "C2_delta_EP")

for (run_loop in 1:runs) {
  
  # ----------------------------------------------------------------------------------------------------
  # starting conditions
  State <- c(P = 1, C1 = 1, C2 = 1, R = 1) # starting parameters
  
  # Calculate temporal variation in mortality timeseries
  z <- matrix(data = rnorm(n = 2*time, mean = 0, sd = 1), nrow = 2)
  cholesky <- matrix(data = c(sigma^2, rho * (sigma ^2), rho * (sigma ^2), sigma ^2), 
                     nrow = 2)
  g <- cholesky %*% z
  M_C1_temp <- 0.4 * exp(g[1,])
  M_C2_temp <- 0.2 * exp(g[2,])
  
  results <- matrix(data=NA, nrow=time, ncol=4)
  results[1,] <- State
  
  for (t in 2:time) {
    # medial consumer 1 mortality rate
    M_C1 <- M_C1_temp[t]
    # medial consumer 2 mortality rate
    M_C2 <- M_C2_temp[t]
    
    pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    
    # Udate state variables to output from last timestep
    State <- c(P = results[t-1,1], C1 = results[t-1,2], C2 = results[t-1,3], R = results[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    
    # Update results matrix
    results[t,] <- c(VF_out[2,2], VF_out[2,3], VF_out[2,4], VF_out[2,5])
    
  }
  
  # ----------------------------------------------------------------------------------------------------
  # Low density growth rate calculations
  
  # C1 as the invader, C2 as the resident species
  State <- c(P = 1, C1 = 0, C2 = 1, R = 1) # starting parameters
  
  # invasion time
  invade_start_time <- time/2
  
  # Run to equilibrium
  C1_invade_C2_resident <- matrix(data=NA, nrow=time, ncol=4)
  C1_invade_C2_resident[1,] <- State
  
  for (t in 2:time) {
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    
    pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    
    # Udate state variables to output from last timestep
    State <- c(P = C1_invade_C2_resident[t-1,1], C1 = C1_invade_C2_resident[t-1,2], 
               C2 = C1_invade_C2_resident[t-1,3], R = C1_invade_C2_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    
    # Update results matrix
    C1_invade_C2_resident[t,] <- c(VF_out[2,2], VF_out[2,3], VF_out[2,4], VF_out[2,5])
    
  }
  
  # now invade C1
  C1_ldgr <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  C2_resident <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  counter <- 1
  
  for (t in invade_start_time:time) {
    
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    
    pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    
    # Udate state variables to output from last timestep
    State <- c(P = C1_invade_C2_resident[t-1,1], C1 = 0.001, 
               C2 = C1_invade_C2_resident[t-1,3], R = C1_invade_C2_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    C1_ldgr[counter] <- log (VF_out[2,3]/.001)
    C2_resident[counter] <- log (VF_out[2,4]/VF_out[1,4])
    
    counter <- counter + 1
  }
  
  # C2 as the invader, C1 as the resident
  State <- c(P = 1, C1 = 1, C2 = 0, R = 1) # starting parameters
  
  # invasion time
  invade_start_time <- time/2
  
  # Run to equilibrium
  C2_invade_C1_resident <- matrix(data=NA, nrow=time, ncol=4)
  C2_invade_C1_resident[1,] <- State
  
  for (t in 2:time) {
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    
    pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    
    # Udate state variables to output from last timestep
    State <- c(P = C2_invade_C1_resident[t-1,1], C1 = C2_invade_C1_resident[t-1,2], 
               C2 = C2_invade_C1_resident[t-1,3], R = C2_invade_C1_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    
    # Update results matrix
    C2_invade_C1_resident[t,] <- c(VF_out[2,2], VF_out[2,3], VF_out[2,4], VF_out[2,5])
    
  }
  
  # now invade C2
  C2_ldgr <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  C1_resident <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  counter <- 1
  
  for (t in invade_start_time:time) {
    
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    
    pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    
    # Udate state variables to output from last timestep
    State <- c(P = C2_invade_C1_resident[t-1,1], C1 = C2_invade_C1_resident[t-1,2], 
               C2 = .001, R = C2_invade_C1_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    C2_ldgr[counter] <- log (VF_out[2,4]/.001)
    C1_resident[counter] <- log (VF_out[2,3]/VF_out[1,3])
    
    counter <- counter + 1
  }
  
  # calculate r_bar
  C1_r_bar <- mean(C1_ldgr)-mean(C2_resident)
  C2_r_bar <- mean(C2_ldgr)-mean(C1_resident)
  
  # ----------------------------------------------------------------------------------------------------
  # Partitioning coexistence mechanisms
  # set up non-fluctuating conditions
  avg_M_C1 <- mean(M_C1_temp[invade_start_time:time])
  avg_M_C2 <- mean(M_C2_temp[invade_start_time:time])
  
  avg_predator_C1_invade_C2_resident <- mean(C1_invade_C2_resident[invade_start_time:time,1])
  avg_predator_C2_invade_C1_resident <- mean(C2_invade_C1_resident[invade_start_time:time,1])
  
  # ----------------------------------------------------------------------------------------------------
  # calculate delta_0 
  
  # C1 as the invader
  State <- c(P = avg_predator_C1_invade_C2_resident, C1 = 0, C2 = 1, R = 1) # starting parameters
  
  # Run to equilibrium
  epsilon_0_C1_invade_C2_resident <- matrix(data=NA, nrow=time, ncol=4)
  epsilon_0_C1_invade_C2_resident[1,] <- State
  
  M_C1 <- avg_M_C1 
  M_C2 <- avg_M_C2 
  pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
  
  for (t in 2:time) {
    
    # Udate state variables to output from last timestep
    State <- c(P = avg_predator_C1_invade_C2_resident, C1 = epsilon_0_C1_invade_C2_resident[t-1,2], 
               C2 = epsilon_0_C1_invade_C2_resident[t-1,3], R = epsilon_0_C1_invade_C2_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    
    # Update results matrix
    epsilon_0_C1_invade_C2_resident[t,] <- c(avg_predator_C1_invade_C2_resident, VF_out[2,3], VF_out[2,4], VF_out[2,5])
    
  }
  
  
  # now invade C1
  C1_epsilon_0 <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  C2_resident_epsilon_0 <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  counter <- 1
  
  M_C1 <- avg_M_C1 
  M_C2 <- avg_M_C2 
  pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
  
  for (t in invade_start_time:time) {
    
    pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    
    # Udate state variables to output from last timestep
    State <- c(P = epsilon_0_C1_invade_C2_resident[t-1,1], C1 = 0.001, 
               C2 = epsilon_0_C1_invade_C2_resident[t-1,3], R = epsilon_0_C1_invade_C2_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    C1_epsilon_0[counter] <- log (VF_out[2,3]/.001)
    C2_resident_epsilon_0[counter] <- log (VF_out[2,4]/VF_out[1,4])
    
    counter <- counter + 1
  }
  
  # C2 as the invader
  State <- c(P = avg_predator_C2_invade_C1_resident, C1 = 1, C2 = 0, R = 1) # starting parameters
  
  # Run to equilibrium
  epsilon_0_C2_invade_C1_resident <- matrix(data=NA, nrow=time, ncol=4)
  epsilon_0_C2_invade_C1_resident[1,] <- State
  
  M_C1 <- avg_M_C1 
  M_C2 <- avg_M_C2 
  pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
  
  for (t in 2:time) {
    
    # Udate state variables to output from last timestep
    State <- c(P = avg_predator_C2_invade_C1_resident, C1 = epsilon_0_C2_invade_C1_resident[t-1,2], 
               C2 = epsilon_0_C2_invade_C1_resident[t-1,3], R = epsilon_0_C2_invade_C1_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    
    # Update results matrix
    epsilon_0_C2_invade_C1_resident[t,] <- c(avg_predator_C2_invade_C1_resident, VF_out[2,3], VF_out[2,4], VF_out[2,5])
    
  }
  
  
  # now invade C2
  C2_epsilon_0 <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  C1_resident_epsilon_0 <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  counter <- 1
  
  M_C1 <- avg_M_C1 
  M_C2 <- avg_M_C2 
  pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
  
  for (t in invade_start_time:time) {
    
    # Udate state variables to output from last timestep
    State <- c(P = epsilon_0_C2_invade_C1_resident[t-1,1], C1 = epsilon_0_C2_invade_C1_resident[t-1,2], 
               C2 = 0.001, R = epsilon_0_C2_invade_C1_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    C2_epsilon_0[counter] <- log (VF_out[2,4]/.001)
    C1_resident_epsilon_0[counter] <- log (VF_out[2,3]/VF_out[1,3])
    
    counter <- counter + 1
  }
  
  C1_delta_0 <- mean(C1_epsilon_0)-mean(C2_resident_epsilon_0)
  C2_delta_0 <- mean(C2_epsilon_0)-mean(C1_resident_epsilon_0)
  
  # ----------------------------------------------------------------------------------------------------
  # calculate delta_E, with predator abundances set at their average 
  
  # C1 as the invader
  State <- c(P = avg_predator_C1_invade_C2_resident, C1 = 0, C2 = 1, R = 1) # starting parameters
  
  # Run to equilibrium
  epsilon_E_C1_invade_C2_resident <- matrix(data=NA, nrow=time, ncol=4)
  epsilon_E_C1_invade_C2_resident[1,] <- State
  
  for (t in 2:time) {
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    
    pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    
    # Udate state variables to output from last timestep
    State <- c(P = avg_predator_C1_invade_C2_resident, C1 = epsilon_E_C1_invade_C2_resident[t-1,2], 
               C2 = epsilon_E_C1_invade_C2_resident[t-1,3], R = epsilon_E_C1_invade_C2_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    
    # Update results matrix
    epsilon_E_C1_invade_C2_resident[t,] <- c(avg_predator_C1_invade_C2_resident, VF_out[2,3], VF_out[2,4], VF_out[2,5])
    
  }
  
  
  # now invade C1
  C1_epsilon_E <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  C2_resident_epsilon_E <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  counter <- 1
  
  for (t in invade_start_time:time) {
    
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    
    pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    
    # Udate state variables to output from last timestep
    State <- c(P = epsilon_E_C1_invade_C2_resident[t-1,1], C1 = 0.001, 
               C2 = epsilon_E_C1_invade_C2_resident[t-1,3], R = epsilon_E_C1_invade_C2_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    C1_epsilon_E[counter] <- log (VF_out[2,3]/.001)
    C2_resident_epsilon_E[counter] <- log (VF_out[2,4]/VF_out[1,4])
    
    counter <- counter + 1
  }
  
  # C2 as the invader
  State <- c(P = avg_predator_C2_invade_C1_resident, C1 = 1, C2 = 0, R = 1) # starting parameters
  
  # Run to equilibrium
  epsilon_E_C2_invade_C1_resident <- matrix(data=NA, nrow=time, ncol=4)
  epsilon_E_C2_invade_C1_resident[1,] <- State
  
  for (t in 2:time) {
    
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    
    pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    
    # Udate state variables to output from last timestep
    State <- c(P = avg_predator_C2_invade_C1_resident, C1 = epsilon_E_C2_invade_C1_resident[t-1,2], 
               C2 = epsilon_E_C2_invade_C1_resident[t-1,3], R = epsilon_E_C2_invade_C1_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    
    # Update results matrix
    epsilon_E_C2_invade_C1_resident[t,] <- c(avg_predator_C2_invade_C1_resident, VF_out[2,3], VF_out[2,4], VF_out[2,5])
    
  }
  
  
  # now invade C2
  C2_epsilon_E <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  C1_resident_epsilon_E <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  counter <- 1
  
  for (t in invade_start_time:time) {
    
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    
    pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
    
    
    # Udate state variables to output from last timestep
    State <- c(P = epsilon_E_C2_invade_C1_resident[t-1,1], C1 = epsilon_E_C2_invade_C1_resident[t-1,2], 
               C2 = 0.001, R = epsilon_E_C2_invade_C1_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    C2_epsilon_E[counter] <- log (VF_out[2,4]/.001)
    C1_resident_epsilon_E[counter] <- log (VF_out[2,3]/VF_out[1,3])
    
    counter <- counter + 1
  }
  
  C1_delta_E <- mean(C1_epsilon_E)-mean(C2_resident_epsilon_E) - C1_delta_0
  C2_delta_E <- mean(C2_epsilon_E)-mean(C1_resident_epsilon_E) - C2_delta_0
  
  # ----------------------------------------------------------------------------------------------------
  # calculate delta_P 
  # predator population size varies, but consumer mortality remains constant
  
  # C1 as the invader
  State <- c(P = 1, C1 = 0, C2 = 1, R = 1) # starting parameters
  
  # Run to equilibrium
  epsilon_P_C1_invade_C2_resident <- matrix(data=NA, nrow=time, ncol=4)
  epsilon_P_C1_invade_C2_resident[1,] <- State
  
  M_C1 <- avg_M_C1 
  M_C2 <- avg_M_C2 
  pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
  
  for (t in 2:time) {
    
    # Udate state variables to output from last timestep
    State <- c(P = epsilon_P_C1_invade_C2_resident[t-1,1], C1 = epsilon_P_C1_invade_C2_resident[t-1,2], 
               C2 = epsilon_P_C1_invade_C2_resident[t-1,3], R = epsilon_P_C1_invade_C2_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    
    # Update results matrix
    epsilon_P_C1_invade_C2_resident[t,] <- c(VF_out[2,2], VF_out[2,3], VF_out[2,4], VF_out[2,5])
    
  }
  
  
  # now invade C1
  C1_epsilon_P <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  C2_resident_epsilon_P <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  counter <- 1
  
  M_C1 <- avg_M_C1 
  M_C2 <- avg_M_C2 
  pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
  
  for (t in invade_start_time:time) {
    
    # Udate state variables to output from last timestep
    State <- c(P = epsilon_P_C1_invade_C2_resident[t-1,1], C1 = 0.001, 
               C2 = epsilon_P_C1_invade_C2_resident[t-1,3], R = epsilon_P_C1_invade_C2_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    C1_epsilon_P[counter] <- log (VF_out[2,3]/.001)
    C2_resident_epsilon_P[counter] <- log (VF_out[2,4]/VF_out[1,4])
    
    counter <- counter + 1
  }
  
  # C2 as the invader
  State <- c(P = 1, C1 = 1, C2 = 0, R = 1) # starting parameters
  
  # Run to equilibrium
  epsilon_P_C2_invade_C1_resident <- matrix(data=NA, nrow=time, ncol=4)
  epsilon_P_C2_invade_C1_resident[1,] <- State
  
  M_C1 <- avg_M_C1 
  M_C2 <- avg_M_C2 
  pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
  
  for (t in 2:time) {
    
    # Udate state variables to output from last timestep
    State <- c(P = epsilon_P_C2_invade_C1_resident[t-1,1], C1 = epsilon_P_C2_invade_C1_resident[t-1,2], 
               C2 = epsilon_P_C2_invade_C1_resident[t-1,3], R = epsilon_P_C2_invade_C1_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    
    # Update results matrix
    epsilon_P_C2_invade_C1_resident[t,] <- c(VF_out[2,2], VF_out[2,3], VF_out[2,4], VF_out[2,5])
    
  }
  
  
  # now invade C2
  C2_epsilon_P <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  C1_resident_epsilon_P <- matrix(data=NA, nrow=(time-invade_start_time), ncol=1)
  counter <- 1
  
  M_C1 <- avg_M_C1 
  M_C2 <- avg_M_C2 
  pars <- c(r, K, J_C1, J_C2, J_P, M_C1, M_C2, M_P, R_0_1, R_0_2, C_0, O_P_C1, O_C1_R, O_C2_R, sigma, rho)
  
  for (t in invade_start_time:time) {
    
    # Udate state variables to output from last timestep
    State <- c(P = epsilon_P_C2_invade_C1_resident[t-1,1], C1 = epsilon_P_C2_invade_C1_resident[t-1,2], 
               C2 = 0.001, R = epsilon_P_C2_invade_C1_resident[t-1,4])
    
    # run ODE solver
    VF_out <- as.data.frame(ode(func = VassFox_Cvar, y = State, parms = pars, times = seq(0,1)), 
                            events = list(func = eventfun))
    C2_epsilon_P[counter] <- log (VF_out[2,4]/.001)
    C1_resident_epsilon_P[counter] <- log (VF_out[2,3]/VF_out[1,3])
    
    counter <- counter + 1
  }
  
  C1_delta_P <- mean(C1_epsilon_P)-mean(C2_resident_epsilon_P) - C1_delta_0
  C2_delta_P <- mean(C2_epsilon_P)-mean(C1_resident_epsilon_P) - C2_delta_0
  
  # ----------------------------------------------------------------------------------------------------
  # calculate delta_EP 
  C1_delta_EP <- C1_r_bar - (C1_delta_0 + C1_delta_P + C1_delta_E)
  C2_delta_EP <- C2_r_bar - (C2_delta_0 + C2_delta_P + C2_delta_E)
  
  # ----------------------------------------------------------------------------------------------------
  # Update final results vector
  C1_results <- c(C1_r_bar, C1_delta_0, C1_delta_P, C1_delta_E, C1_delta_EP)
  C2_results <- c(C2_r_bar, C2_delta_0, C2_delta_P, C2_delta_E, C2_delta_EP)
  
  C1_final_mechanisms[run_loop,] <- C1_results
  C2_final_mechanisms[run_loop,] <- C2_results
  
  # print how far we've come
   print(run_loop)

}

# ----------------------------------------------------------------------------------------------------
# to plot run plotting_mechanisms.R

