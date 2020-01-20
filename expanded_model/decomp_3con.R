#mechanistic decomposition of 3 con 1 pred model

rm(list = ls())
library(deSolve)
library(here)
library(foreach)
library(doSNOW)


# ----------------------------------------------------------------------------------------------------
# Model of 3 cons

VassFox_Cvar <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    all_Cons <- (O_P_C1 * C1) + (O_P_C2 * C2) + (1 - (O_P_C1 + O_P_C2)) * C3
    
    dP = -(M_P * P) + (((J_P * P) * all_Cons) / (all_Cons + C_0))
    
    dC1 = -(M_C1 * C1) + ((O_C1_R * J_C1 * C1 * R) / (R + R_0_1)) - ((O_P_C1 * J_P * P * C1) /
                                                                       (all_Cons + C_0))
    
    dC2 = -(M_C2 * C2) + ((O_C2_R * J_C2 * C2 * R) / (R + R_0_2)) - ((O_P_C2 * J_P * P * C2) / (all_Cons +
                                                                                                  C_0))
    
    dC3 = -(M_C3 * C3) + ((O_C3_R * J_C3 * C3 * R) / (R + R_0_3)) - (((1 - (O_P_C1 +
                                                                              O_P_C2)) * J_P * P * C3) / (all_Cons + C_0))
    
    dR = r * R * (1 - (R / K)) - ((O_C1_R * J_C1 * C1 * R) / (R + R_0_1)) - ((O_C2_R * J_C2 * C2 * R) / (R + R_0_2)) - ((O_C3_R * J_C3 * C3 * R) / (R + R_0_3))
    return(list(c(dP, dC1, dC2, dC3, dR)))
    
  })
}

#estimated parameters
#parameter set 1
est<-c(0.262571  ,  0.862126  ,  0.227814  ,  0.881397 ,   0.065556 ,   0.533164) 
  

# parameters
# resource intrinsic rate of growth
r <- 1.0
# resource carrying capacity
K <- 1.0
# consumer 1 ingestion rate
J_C1 <- 0.8036
# consumer 2 ingestion rate
J_C2 <- 0.7
# consumer 3 ingestion rate (estimated)
J_C3 <- est[2]
# predator ingestion rate
J_P <- 0.4
# predator mortality rate
M_P <- 0.08
# half saturation constant
R_0_1 <- 0.16129
R_0_2 <- 0.9
# half saturation constant (estimated)
R_0_3 <- est[3]
C_0 <- 0.5
# preference coefficient (estimated)

O_P_C1 <- est[4]#set 3
O_P_C2 <- est[5] #set3


O_C1_R <- 1.0
O_C2_R <- 0.98
#estimated

O_C3_R <- est[6]
time <- 5000

# Run to equilibrium
runtoequal <- function(State_int,
                       M_C1_temp,
                       M_C2_temp,
                       M_C3_temp) {
  mat <- matrix(data = NA,
                nrow = time,
                ncol = 5)
  mat[1,] <- State_int
  
  for (t in 2:time) {
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    M_C3 <- M_C3_temp[t]
    
    
    pars <-
      c(
        r = r,
        K = K,
        J_C1 = J_C1,
        J_C2 = J_C2,
        J_C3 = J_C3,
        J_P = J_P,
        M_C1 = M_C1,
        M_C2 = M_C2,
        M_C3 = M_C3,
        M_P = M_P,
        R_0_1 = R_0_1,
        R_0_2 = R_0_2,
        R_0_3 = R_0_3,
        C_0 = C_0,
        O_P_C1 = O_P_C1,
        O_P_C2 = O_P_C2,
        O_C1_R = O_C1_R,
        O_C2_R = O_C2_R,
        O_C3_R = O_C3_R
      )
    
    # Udate state variables to output from last timestep
    State <- c(
      P = mat[t - 1, 1],
      C1 = mat[t - 1, 2],
      C2 = mat[t - 1, 3],
      C3 = mat[t - 1, 4],
      R = mat[t - 1, 5]
    )
    
    # run ODE solver
    VF <-
      as.data.frame(
        ode(
          func = VassFox_Cvar,
          y = State,
          parms = pars,
          times = seq(0, 1)
        ),
        events = list(func = eventfun)
      )
    
    # Update results matrix
    mat[t,] <- c(VF[2, 2], VF[2, 3], VF[2, 4], VF[2, 5], VF[2, 6])
    
    
  }
  return(mat)
  
}


#sets pred to mean value
runtoequal_set <-
  function(State_int,
           M_C1_temp,
           M_C2_temp,
           M_C3_temp,
           P) {
    mat <- matrix(data = NA,
                  nrow = time,
                  ncol = 5)
    mat[1,] <- State_int
    
    for (t in 2:time) {
      M_C1 <- M_C1_temp[t]
      M_C2 <- M_C2_temp[t]
      M_C3 <- M_C3_temp[t]
      
      pars <-
        c(
          r = r,
          K = K,
          J_C1 = J_C1,
          J_C2 = J_C2,
          J_C3 = J_C3,
          J_P = J_P,
          M_C1 = M_C1,
          M_C2 = M_C2,
          M_C3 = M_C3,
          M_P = M_P,
          R_0_1 = R_0_1,
          R_0_2 = R_0_2,
          R_0_3 = R_0_3,
          C_0 = C_0,
          O_P_C1 = O_P_C1,
          O_P_C2 = O_P_C2,
          O_C1_R = O_C1_R,
          O_C2_R = O_C2_R,
          O_C3_R = O_C3_R
        )
      
      # Udate state variables to output from last timestep
      State <- c(
        P = P,
        C1 = mat[t - 1, 2],
        C2 = mat[t - 1, 3],
        C3 = mat[t - 1, 4],
        R = mat[t - 1, 5]
      )
      
      
      # run ODE solver
      VF <-
        as.data.frame(
          ode(
            func = VassFox_Cvar,
            y = State,
            parms = pars,
            times = seq(0, 1)
          ),
          events = list(func = eventfun)
        )
      
      # Update results matrix
      mat[t,] <- c(P, VF[2, 3], VF[2, 4], VF[2, 5], VF[2, 6])
      
    }
    return(mat)
    
  }

#calculates what happens when C1 invades
invade_C1 <- function(mat, M_C1_temp, M_C2_temp, M_C3_temp) {
  invade_start_time <- time / 2
  
  
  C1_ldgr <- matrix(
    data = NA,
    nrow = (time - invade_start_time),
    ncol = 1
  )
  C2_resident <-
    matrix(
      data = NA,
      nrow = (time - invade_start_time),
      ncol = 1
    )
  
  C3_resident <-  matrix(
    data = NA,
    nrow = (time - invade_start_time),
    ncol = 1
  )
  
  counter <- 1
  
  for (t in invade_start_time:time) {
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    M_C3 <- M_C3_temp[t]
    
    pars <-
      c(
        r = r,
        K = K,
        J_C1 = J_C1,
        J_C2 = J_C2,
        J_C3 = J_C3,
        J_P = J_P,
        M_C1 = M_C1,
        M_C2 = M_C2,
        M_C3 = M_C3,
        M_P = M_P,
        R_0_1 = R_0_1,
        R_0_2 = R_0_2,
        R_0_3 = R_0_3,
        C_0 = C_0,
        O_P_C1 = O_P_C1,
        O_P_C2 = O_P_C2,
        O_C1_R = O_C1_R,
        O_C2_R = O_C2_R,
        O_C3_R = O_C3_R
      )
    
    # Udate state variables to output from last timestep
    State <- c(
      P = mat[t - 1, 1],
      C1 = 0.001,
      C2 = mat[t - 1, 3],
      C3 = mat[t - 1, 4],
      R = mat[t - 1, 5]
    )
    
    # run ODE solver
    VF_out <-
      as.data.frame(
        ode(
          func = VassFox_Cvar,
          y = State,
          parms = pars,
          times = seq(0, 1)
        ),
        events = list(func = eventfun)
      )
    C1_ldgr[counter] <- log (VF_out[2, 3] / .001)
    C2_resident[counter] <- log (VF_out[2, 4] / VF_out[1, 4])
    C3_resident[counter] <- log (VF_out[2, 5] / VF_out[1, 5])
    counter <- counter + 1
  }
  
  
  
  result <-
    data.frame(C1_ldgr = C1_ldgr,
               C2_resident = C2_resident,
               C3_resident = C3_resident)
  return(result)
  
}

#calculate what happens when C2 invades
invade_C2 <- function(mat, M_C1_temp, M_C2_temp, M_C3_temp) {
  invade_start_time <- time / 2
  
  
  C2_ldgr <- matrix(
    data = NA,
    nrow = (time - invade_start_time),
    ncol = 1
  )
  C1_resident <-
    matrix(
      data = NA,
      nrow = (time - invade_start_time),
      ncol = 1
    )
  C3_resident <-
    matrix(
      data = NA,
      nrow = (time - invade_start_time),
      ncol = 1
    )
  
  
  counter <- 1
  
  for (t in invade_start_time:time) {
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    M_C3 <- M_C3_temp[t]
    
    
    pars <-
      c(
        r = r,
        K = K,
        J_C1 = J_C1,
        J_C2 = J_C2,
        J_C3 = J_C3,
        J_P = J_P,
        M_C1 = M_C1,
        M_C2 = M_C2,
        M_C3 = M_C3,
        M_P = M_P,
        R_0_1 = R_0_1,
        R_0_2 = R_0_2,
        R_0_3 = R_0_3,
        C_0 = C_0,
        O_P_C1 = O_P_C1,
        O_P_C2 = O_P_C2,
        O_C1_R = O_C1_R,
        O_C2_R = O_C2_R,
        O_C3_R = O_C3_R
      )
    
    # run ODE solver
    State <- c(
      P = mat[t - 1, 1],
      C1 = mat[t - 1, 2],
      C2 = .001,
      C3 = mat[t - 1, 4],
      R = mat[t - 1, 5]
    )
    
    # run ODE solver
    VF_out <-
      as.data.frame(
        ode(
          func = VassFox_Cvar,
          y = State,
          parms = pars,
          times = seq(0, 1)
        ),
        events = list(func = eventfun)
      )
    C2_ldgr[counter] <- log (VF_out[2, 4] / .001)
    C1_resident[counter] <- log (VF_out[2, 3] / VF_out[1, 3])
    C3_resident[counter] <- log (VF_out[2, 5] / VF_out[1, 5])
    counter <- counter + 1
  }
  
  result <-
    data.frame(C2_ldgr = C2_ldgr,
               C1_resident = C1_resident,
               C3_resident = C3_resident)
  return(result)
  
}

#calculate what happens when C3 invades
invade_C3 <- function(mat, M_C1_temp, M_C2_temp, M_C3_temp) {
  invade_start_time <- time / 2
  
  C3_ldgr <- matrix(
    data = NA,
    nrow = (time - invade_start_time),
    ncol = 1
  )
  C1_resident <-
    matrix(
      data = NA,
      nrow = (time - invade_start_time),
      ncol = 1
    )
  C2_resident <-
    matrix(
      data = NA,
      nrow = (time - invade_start_time),
      ncol = 1
    )
  
  
  counter <- 1
  
  for (t in invade_start_time:time) {
    M_C1 <- M_C1_temp[t]
    M_C2 <- M_C2_temp[t]
    M_C3 <- M_C3_temp[t]
    
    
    pars <-
      c(
        r = r,
        K = K,
        J_C1 = J_C1,
        J_C2 = J_C2,
        J_C3 = J_C3,
        J_P = J_P,
        M_C1 = M_C1,
        M_C2 = M_C2,
        M_C3 = M_C3,
        M_P = M_P,
        R_0_1 = R_0_1,
        R_0_2 = R_0_2,
        R_0_3 = R_0_3,
        C_0 = C_0,
        O_P_C1 = O_P_C1,
        O_P_C2 = O_P_C2,
        O_C1_R = O_C1_R,
        O_C2_R = O_C2_R,
        O_C3_R = O_C3_R
      )
    
    # run ODE solver
    State <- c(
      P = mat[t - 1, 1],
      C1 = mat[t - 1, 2],
      C2 = mat[t - 1, 3],
      C3 = .001,
      R = mat[t - 1, 5]
    )
    
    # run ODE solver
    VF_out <-
      as.data.frame(
        ode(
          func = VassFox_Cvar,
          y = State,
          parms = pars,
          times = seq(0, 1)
        ),
        events = list(func = eventfun)
      )
    C3_ldgr[counter] <- log (VF_out[2, 5] / .001)
    C1_resident[counter] <- log (VF_out[2, 3] / VF_out[1, 3])
    C2_resident[counter] <- log (VF_out[2, 4] / VF_out[1, 4])
    counter <- counter + 1
  }
  
  result <-
    data.frame(C3_ldgr = C3_ldgr,
               C1_resident = C1_resident,
               C2_resident = C2_resident)
  return(result)
  
}


runall <- function(sigma, rho, time) {
  #evaluate all combinations of invaders, wrapped up in a function so that it can be called in parallel
  # ----------------------------------------------------------------------------------------------------
  # starting conditions
  State <- c(
    P = 1,
    C1 = 1,
    C2 = 1,
    C3 = 1,
    R = 1
  ) # starting parameters
  
  # Calculate temporal variation in mortality timeseries
  z <- matrix(data = rnorm(
    n = 2 * time,
    mean = 0,
    sd = 1
  ), nrow = 2)
  
  cholesky <-
    matrix(data = c(sigma ^ 2, rho * (sigma ^ 2), rho * (sigma ^ 2), sigma ^
                      2),
           nrow = 2)
  
  g <- cholesky %*% z
  M_C1_temp <- 0.4 * exp(g[1,])
  M_C2_temp <- 0.2 * exp(g[2,])
  
  #reroll the random matrix
  z <- matrix(data = rnorm(
    n = 2 * time,
    mean = 0,
    sd = 1
  ), nrow = 2)
  
  g <- cholesky %*% z
  
  #growth rate estimated
  
  M_C3_temp <- est[1] * exp(g[1,])
  
  invade_start_time <- time / 2
  
  State <- c(
    P = 1,
    C1 = 0,
    C2 = 1,
    C3 = 1,
    R = 1
  ) # starting parameters
  
  C1_invade_C2_resident_C3_resident <-
    runtoequal(State, M_C1_temp, M_C2_temp, M_C3_temp)
  # now invade C1
  invadeC1 <-
    invade_C1(C1_invade_C2_resident_C3_resident,
              M_C1_temp,
              M_C2_temp,
              M_C3_temp)
  C1_ldgr <- invadeC1$C1_ldgr
  C2_resident_C1_invade <- invadeC1$C2_resident
  C3_resident_C1_invade <- invadeC1$C3_resident
  
  # C2 as the invader
  State <- c(
    P = 1,
    C1 = 1,
    C2 = 0,
    C3 = 1,
    R = 1
  ) # starting parameters
  
  C2_invade_C1_resident_C3_resident <-
    runtoequal(State, M_C1_temp, M_C2_temp, M_C3_temp)
  # now invade C2
  invadeC2 <-
    invade_C2(C2_invade_C1_resident_C3_resident,
              M_C1_temp,
              M_C2_temp,
              M_C3_temp)
  C2_ldgr <- invadeC2$C2_ldgr
  
  C1_resident_C2_invade <- invadeC2$C1_resident
  C3_resident_C2_invade <- invadeC2$C3_resident
  
  # C3 as the invader
  State <- c(
    P = 1,
    C1 = 1,
    C2 = 1,
    C3 = 0,
    R = 1
  ) # starting parameters
  
  C3_invade_C1_resident_C2_resident <-
    runtoequal(State, M_C1_temp, M_C2_temp, M_C3_temp)
  # now invade C3
  invadeC3 <-
    invade_C3(C3_invade_C1_resident_C2_resident,
              M_C1_temp,
              M_C2_temp,
              M_C3_temp)
  C3_ldgr <- invadeC3$C3_ldgr
  
  C1_resident_C3_invade <- invadeC3$C1_resident
  C2_resident_C3_invade <- invadeC3$C1_resident
  
  
  
  
  # calculate r_bar notice division by 2,
  # ----------------------------------------------------------------------------------------------------
  
  
  C1_r_bar <-
    mean(C1_ldgr) - (mean(C2_resident_C1_invade) + mean(C3_resident_C1_invade)) /
    2
  C2_r_bar <-
    mean(C2_ldgr) - (mean(C1_resident_C2_invade) + mean(C3_resident_C2_invade)) /
    2
  C3_r_bar <-
    mean(C3_ldgr) - (mean(C1_resident_C3_invade) + mean(C2_resident_C3_invade)) /
    2
  # ----------------------------------------------------------------------------------------------------
  # Partitioning coexistence mechanisms
  # set up non-fluctuating conditions
  avg_M_C1 <- mean(M_C1_temp[invade_start_time:time])
  avg_M_C2 <- mean(M_C2_temp[invade_start_time:time])
  avg_M_C3 <- mean(M_C3_temp[invade_start_time:time])
  
  
  avg_predator_C1_invade_C2_resident_C3_resident <-
    mean(C1_invade_C2_resident_C3_resident[invade_start_time:time, 1])
  
  avg_predator_C2_invade_C1_resident_C3_resident <-
    mean(C2_invade_C1_resident_C3_resident[invade_start_time:time, 1])
  
  avg_predator_C3_invade_C1_resident_C2_resident <-
    mean(C3_invade_C1_resident_C2_resident[invade_start_time:time, 1])
  
  # ----------------------------------------------------------------------------------------------------
  # calculate delta_0
  
  # C1 as the invader
  State <-
    c(
      P = avg_predator_C1_invade_C2_resident_C3_resident,
      C1 = 0,
      C2 = 1,
      C3 = 1,
      R = 1
    ) # starting parameters
  
  
  epsilon_0_C1_invade_C2_resident_C3_resident <-
    runtoequal_set(
      State,
      rep(avg_M_C1, length(M_C1_temp)),
      rep(avg_M_C2, length(M_C2_temp)),
      rep(avg_M_C3, length(M_C3_temp)),
      avg_predator_C1_invade_C2_resident_C3_resident
    )
  
  # now invade C1
  invadeC1 <-
    invade_C1(
      epsilon_0_C1_invade_C2_resident_C3_resident,
      rep(avg_M_C1, length(M_C1_temp)),
      rep(avg_M_C2, length(M_C2_temp)),
      rep(avg_M_C2, length(M_C2_temp))
    )
  
  
  C1_epsilon_0 <- invadeC1$C1_ldgr
  C2_resident_epsilon_0_invade_C1 <- invadeC1$C2_resident
  C3_resident_epsilon_0_invade_C1 <- invadeC1$C3_resident
  
  # C2 as the invader
  State <-
    c(
      P = avg_predator_C2_invade_C1_resident_C3_resident,
      C1 = 1,
      C2 = 0,
      C3 = 1,
      R = 1
    ) # starting parameters
  
  epsilon_0_C2_invade_C1_resident_C3_resident <-
    runtoequal_set(
      State,
      rep(avg_M_C1, length(M_C1_temp)),
      rep(avg_M_C2, length(M_C2_temp)),
      rep(avg_M_C3, length(M_C3_temp)),
      avg_predator_C2_invade_C1_resident_C3_resident
    )
  
  # now invade C2
  invadeC2 <-
    invade_C2(
      epsilon_0_C2_invade_C1_resident_C3_resident,
      rep(avg_M_C1, length(M_C1_temp)),
      rep(avg_M_C2, length(M_C2_temp)),
      rep(avg_M_C3, length(M_C3_temp))
    )
  
  C2_epsilon_0 <- invadeC2$C2_ldgr
  C1_resident_epsilon_0_invade_C2 <- invadeC2$C1_resident
  C3_resident_epsilon_0_invade_C2 <- invadeC2$C3_resident
  
  
  
  # C3 as the invader
  State <-
    c(
      P = avg_predator_C3_invade_C1_resident_C2_resident,
      C1 = 1,
      C2 = 1,
      C3 = 0,
      R = 1
    ) # starting parameters
  
  epsilon_0_C3_invade_C1_resident_C2_resident <-
    runtoequal_set(
      State,
      rep(avg_M_C1, length(M_C1_temp)),
      rep(avg_M_C2, length(M_C2_temp)),
      rep(avg_M_C3, length(M_C3_temp)),
      avg_predator_C3_invade_C1_resident_C2_resident
    )
  
  # now invade C3
  invadeC3 <-
    invade_C3(
      epsilon_0_C3_invade_C1_resident_C2_resident,
      rep(avg_M_C1, length(M_C1_temp)),
      rep(avg_M_C2, length(M_C2_temp)),
      rep(avg_M_C3, length(M_C3_temp))
    )
  
  C3_epsilon_0 <- invadeC3$C3_ldgr
  C1_resident_epsilon_0_invade_C3 <- invadeC3$C1_resident
  C2_resident_epsilon_0_invade_C3 <- invadeC3$C2_resident
  
  
  
  # note the division by 2 because we are comparing to 2 other cons now rather then just 1
  # ----------------------------------------------------------------------------------------------------
  
  C1_delta_0 <-
    mean(C1_epsilon_0) - (mean(C2_resident_epsilon_0_invade_C1) + mean(C3_resident_epsilon_0_invade_C1)) /
    2
  C2_delta_0 <-
    mean(C2_epsilon_0) - (mean(C1_resident_epsilon_0_invade_C2) + mean(C3_resident_epsilon_0_invade_C2)) /
    2
  C3_delta_0 <-
    mean(C3_epsilon_0) - (mean(C1_resident_epsilon_0_invade_C3) + mean(C2_resident_epsilon_0_invade_C3)) /
    2
  
  # ----------------------------------------------------------------------------------------------------
  # calculate delta_E
  
  # C1 as the invader
  State <-
    c(
      P = avg_predator_C1_invade_C2_resident_C3_resident,
      C1 = 0,
      C2 = 1,
      C3 = 1,
      R = 1
    ) # starting parameters
  
  
  epsilon_E_C1_invade_C2_resident_C3_resident <-
    runtoequal_set(
      State,
      M_C1_temp,
      M_C2_temp,
      M_C3_temp,
      avg_predator_C1_invade_C2_resident_C3_resident
    )
  
  # now invade C1
  invadeC1 <-
    invade_C1(epsilon_E_C1_invade_C2_resident_C3_resident,
              M_C1_temp,
              M_C2_temp,
              M_C3_temp)
  
  
  C1_epsilon_E <- invadeC1$C1_ldgr
  
  C2_resident_epsilon_E_invade_C1 <- invadeC1$C2_resident
  C3_resident_epsilon_E_invade_C1 <- invadeC1$C3_resident
  
  
  # C2 as the invader
  State <-
    c(
      P = avg_predator_C2_invade_C1_resident_C3_resident,
      C1 = 1,
      C2 = 0,
      C3 = 1,
      R = 1
    )
  
  epsilon_E_C2_invade_C1_resident_C3_resident <-
    runtoequal_set(
      State,
      M_C1_temp,
      M_C2_temp,
      M_C3_temp,
      avg_predator_C2_invade_C1_resident_C3_resident
    )
  # now invade C2
  invadeC2 <-
    invade_C2(epsilon_E_C2_invade_C1_resident_C3_resident,
              M_C1_temp,
              M_C2_temp,
              M_C3_temp)
  C2_epsilon_E <- invadeC2$C2_ldgr
  
  C1_resident_epsilon_E_invade_C2 <- invadeC2$C1_resident
  C3_resident_epsilon_E_invade_C2 <- invadeC2$C3_resident
  
  
  # C3 as the invader
  State <-
    c(
      P = avg_predator_C3_invade_C1_resident_C2_resident,
      C1 = 1,
      C2 = 1,
      C3 = 0,
      R = 1
    )
  
  epsilon_E_C3_invade_C1_resident_C2_resident <-
    runtoequal_set(
      State,
      M_C1_temp,
      M_C2_temp,
      M_C3_temp,
      avg_predator_C3_invade_C1_resident_C2_resident
    )
  # now invade C3
  invadeC3 <-
    invade_C3(epsilon_E_C3_invade_C1_resident_C2_resident,
              M_C1_temp,
              M_C2_temp,
              M_C3_temp)
  C3_epsilon_E <- invadeC3$C3_ldgr
  
  C1_resident_epsilon_E_invade_C3 <- invadeC3$C1_resident
  C2_resident_epsilon_E_invade_C3 <- invadeC3$C1_resident
  
  
  
  
  C1_delta_E <-
    mean(C1_epsilon_E) - (mean(C2_resident_epsilon_E_invade_C1) + mean(C3_resident_epsilon_E_invade_C1)) /
    2 - C1_delta_0
  C2_delta_E <-
    mean(C2_epsilon_E) - (mean(C1_resident_epsilon_E_invade_C2) + mean(C3_resident_epsilon_E_invade_C2)) /
    2 - C2_delta_0
  C3_delta_E <-
    mean(C3_epsilon_E) - (mean(C1_resident_epsilon_E_invade_C3) + mean(C2_resident_epsilon_E_invade_C3)) /
    2 - C3_delta_0
  
  
  
  
  # ----------------------------------------------------------------------------------------------------
  # calculate delta_P
  # predator population size varies, but mortality remains constant
  
  # C1 as the invader
  State <- c(
    P = 1,
    C1 = 0,
    C2 = 1,
    C3 = 1,
    R = 1
  )
  
  epsilon_P_C1_invade_C2_resident_C3_resident <-
    runtoequal(State,
               rep(avg_M_C1, length(M_C1_temp)),
               rep(avg_M_C2, length(M_C2_temp)),
               rep(avg_M_C3, length(M_C3_temp)))
  
  
  # now invade C1
  invadeC1 <-
    invade_C1(
      epsilon_P_C1_invade_C2_resident_C3_resident,
      rep(avg_M_C1, length(M_C1_temp)),
      rep(avg_M_C2, length(M_C2_temp)),
      rep(avg_M_C3, length(M_C3_temp))
    )
  
  C1_epsilon_P <- invadeC1$C1_ldgr
  
  C2_resident_epsilon_P_invade_C1 <- invadeC1$C2_resident
  C3_resident_epsilon_P_invade_C1 <- invadeC1$C3_resident
  
  
  
  # C2 as the invader
  State <- c(
    P = 1,
    C1 = 1,
    C2 = 0,
    C3 = 1,
    R = 1
  ) # starting parameters
  
  epsilon_P_C2_invade_C1_resident_C3_resident <-
    runtoequal(State,
               rep(avg_M_C1, length(M_C1_temp)),
               rep(avg_M_C2, length(M_C2_temp)),
               rep(avg_M_C3, length(M_C3_temp)))
  
  # now invade C2
  invadeC2 <-
    invade_C2(
      epsilon_P_C2_invade_C1_resident_C3_resident ,
      rep(avg_M_C1, length(M_C1_temp)),
      rep(avg_M_C2, length(M_C2_temp)),
      rep(avg_M_C3, length(M_C3_temp))
    )
  
  
  C2_epsilon_P <- invadeC2$C2_ldgr
  
  C1_resident_epsilon_P_invade_C2 <- invadeC2$C1_resident
  C3_resident_epsilon_P_invade_C2 <- invadeC2$C3_resident
  
  
  
  # C3 as the invader
  State <- c(
    P = 1,
    C1 = 1,
    C2 = 1,
    C3 = 0,
    R = 1
  ) # starting parameters
  
  epsilon_P_C3_invade_C1_resident_C2_resident <-
    runtoequal(State,
               rep(avg_M_C1, length(M_C1_temp)),
               rep(avg_M_C2, length(M_C2_temp)),
               rep(avg_M_C3, length(M_C3_temp)))
  
  # now invade C2
  invadeC3 <-
    invade_C3(
      epsilon_P_C3_invade_C1_resident_C2_resident ,
      rep(avg_M_C1, length(M_C1_temp)),
      rep(avg_M_C2, length(M_C2_temp)),
      rep(avg_M_C3, length(M_C3_temp))
    )
  
  
  C3_epsilon_P <- invadeC3$C3_ldgr
  
  C1_resident_epsilon_P_invade_C3 <- invadeC3$C1_resident
  C2_resident_epsilon_P_invade_C3 <- invadeC3$C2_resident
  
  
  # same issue with dealing with 3 terms
  # ----------------------------------------------------------------------------------------------------
  
  C1_delta_P <-
    mean(C1_epsilon_P) - (mean(C2_resident_epsilon_P_invade_C1) + mean(C3_resident_epsilon_P_invade_C1)) /
    2 - C1_delta_0
  C2_delta_P <-
    mean(C2_epsilon_P) - (mean(C1_resident_epsilon_P_invade_C2) + mean(C3_resident_epsilon_E_invade_C2)) /
    2 - C2_delta_0
  
  C3_delta_P <-
    mean(C3_epsilon_P) - (mean(C1_resident_epsilon_P_invade_C3) + mean(C1_resident_epsilon_E_invade_C3)) /
    2 - C3_delta_0
  
  
  # ----------------------------------------------------------------------------------------------------
  # calculate delta_EP
  C1_delta_EP <- C1_r_bar - (C1_delta_0 + C1_delta_P + C1_delta_E)
  C2_delta_EP <- C2_r_bar - (C2_delta_0 + C2_delta_P + C2_delta_E)
  C3_delta_EP <- C3_r_bar - (C3_delta_0 + C3_delta_P + C3_delta_E)
  
  # ----------------------------------------------------------------------------------------------------
  # Update final results vector
  C1_results <-
    c(C1_r_bar, C1_delta_0, C1_delta_P, C1_delta_E, C1_delta_EP)
  C2_results <-
    c(C2_r_bar, C2_delta_0, C2_delta_P, C2_delta_E, C2_delta_EP)
  C3_results <-
    c(C3_r_bar, C3_delta_0, C3_delta_P, C3_delta_E, C3_delta_EP)
  
  
  done <- NULL
  #record results
  done$C1_final_mechanisms <- C1_results
  done$C2_final_mechanisms <- C2_results
  done$C3_final_mechanisms <- C3_results
  
  return(done)
  
}

# strength of env. on mortality rate
sigma <- 0.55

# cross-correlation of C1 and C2
rho <- 0

#set up for parallization, change the number based on how many cores your machine has. Most personal laptops have 4 or 8. If you don't know how many you have its probably 4. If you have a fancy cluster crank this up to like 100.
cluster = makeCluster(100, type = "SOCK")
registerDoSNOW(cluster)
# looping over multiple runs
runs <- 500

#call the function that will do all of the calculations and store the results
results <-
  foreach (
    run_loop = 1:runs,
    .packages = c("deSolve") ,
    .combine = 'rbind'
  ) %dopar% {
    d = runall(sigma, rho, time)
  }


#format things nicely and save
C1_final_mechanisms <- t(as.data.frame(results[, 1]))
colnames(C1_final_mechanisms) <-
  c("C1_r_bar",
    "C1_delta_0",
    "C1_delta_P",
    "C1_delta_E",
    "C1_delta_EP")

C2_final_mechanisms <-  t(as.data.frame(results[, 2]))
colnames(C2_final_mechanisms) <-
  c("C2_r_bar",
    "C2_delta_0",
    "C2_delta_P",
    "C2_delta_E",
    "C2_delta_EP")


C3_final_mechanisms <-  t(as.data.frame(results[, 3]))
colnames(C3_final_mechanisms) <-
  c("C3_r_bar",
    "C3_delta_0",
    "C3_delta_P",
    "C3_delta_E",
    "C3_delta_EP")


write.csv(
  C1_final_mechanisms,
  file = here("C1_final_mechanisms_3con_param1.csv"),
  row.names = FALSE
)
write.csv(
  C2_final_mechanisms,
  file = here("C2_final_mechanisms_3con_param1.csv"),
  row.names = FALSE
)
write.csv(
  C3_final_mechanisms,
  file = here("C3_final_mechanisms_3con_param1.csv"),
  row.names = FALSE
)

#release the cores we used, this is important! Always remember to run this step if you have created a cluster, or else your cores won't be released until you restart the R session.
stopCluster(cluster)
