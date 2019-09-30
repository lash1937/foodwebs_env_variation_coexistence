

# Estimate parameters for model with 3Cs
# This program was run twice to get the parameter values given in the supplemental table.

rm(list = ls())
library(deSolve)
library(here)
library(foreach)
library(doSNOW)
library(ggplot2)
library(cowplot)
library(reshape)
library(tidyverse)
library(R.utils)
library(DEoptim)
library(lubridate)
library(compiler)

optimize_function <- function(pars) {
  #exapanded model with 3C
  VassFox_2Cvar <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      all_Cons <- (O_P_C1 * C1) + (O_P_C2 * C2) + (1 - (O_P_C1 + O_P_C2)) * C3
      
      dP = -(M_P * P) + (((J_P * P) * all_Cons) / (all_Cons + C_0))
      
      dC1 = -(M_C1 * C1) + ((O_C1_R * J_C1 * C1 * R) / (R + R_0_1)) - ((O_P_C1 * J_P * P * C1) /
                                                                         (all_Cons + C_0))
      
      dC2 = -(M_C2 * C2) + ((O_C2_R * J_C2 * C2 * R) / (R + R_0_2)) - ((O_P_C2 * J_P * P * C2) / (all_Cons +
                                                                                                    C_0))
      
      dC3 = -(M_C3 * C3) + ((O_C3_R * J_C3 * C3 * R) / (R + R_0_3)) - (((1 - (
        O_P_C1 + O_P_C2
      )) * J_P * P * C3) / (all_Cons + C_0))
      
      dR = r * R * (1 - (R / K)) - ((O_C1_R * J_C1 * C1 * R) / (R + R_0_1)) - ((O_C2_R * J_C2 * C2 * R) / (R + R_0_2)) - ((O_C3_R * J_C3 * C3 * R) / (R + R_0_3))
      return(list(c(dP, dC1, dC2, dC3, dR)))
      
    })
  }
  
  # strength of env. on mortality rate
  # if you want to simulate a system with env. flux use sigma <-0.55.
  #sigma <- 0.55
  sigma <- 0
  # cross-correlation of C1 and C2
  ro <- 0
  
  z <- matrix(data = rnorm(
    n = 2 * 5000,
    mean = 0,
    sd = 1
  ), nrow = 2)
  cholesky <-
    matrix(data = c(sigma ^ 2, ro * (sigma ^ 2), ro * (sigma ^ 2), sigma ^
                      2),
           nrow = 2)
  g <- cholesky %*% z
  M_C1_temp <- 0.4 * exp(g[1, ])
  M_C2_temp <- 0.2 * exp(g[2, ])
  
  #redraw fluxs for C3, rather then make a 3x matrix, its easier to do this.
  z <- matrix(data = rnorm(
    n = 2 * 5000,
    mean = 0,
    sd = 1
  ), nrow = 2)
  g <- cholesky %*% z
  
  M_C3_temp <- pars[1] * exp(g[1, ])
  
  
  # ----------------------------------------------------------------------------------------------------
  # parameters
  # resource intrinsic rate of growth
  r <- 1.0
  # resource carrying capacity
  K <- 1.0
  # consumer 1 ingestion rate
  J_C1 <- 0.8036
  # consumer 2 ingestion rate
  J_C2 <- 0.7
  # consumer 3 ingestion rate
  J_C3 <- pars[2]
  # predator ingestion rate
  J_P <- 0.4
  # predator mortality rate
  M_P <- 0.08
  
  # half saturation constant
  R_0_1 <- 0.16129
  R_0_2 <- 0.9
  R_0_3 <- pars[3]
  
  
  C_0 <- 0.5
  # preference coefficient (need to estimate these because the preference for C3 is just 1-(C1+C2) )
  O_P_C1 <- pars[4]
  O_P_C2 <- pars[5]
  
  O_C1_R <- 1.0
  O_C2_R <- 0.98
  O_C3_R <- pars[6]

  time <- 5000
  
  
  # make sure that preferences are fractions and that C3 doesnt have a zero preference
  if (O_P_C1 + O_P_C2 >= 1) {
    return(Inf)
  }
  
  
  State <-
    c(
      P = 1,
      C1 = 1,
      C2 = 1,
      C3 = 1,
      R = 1
    ) # starting parameters
  
  mat <- matrix(data = NA,
                nrow = time,
                ncol = 5)
  mat[1, ] <- State
  
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
    
    # Update state variables to output from last timestep
    State <- c(
      P = mat[t - 1, 1],
      C1 = mat[t - 1, 2],
      C2 = mat[t - 1, 3],
      C3 = mat[t - 1, 4],
      R = mat[t - 1, 5]
    )
    
    
    an.error.occured <- FALSE
    
    
    #if the eval produces errorrs or takes a long time don't bother just throw an error
    tryCatch({
      VF <- withTimeout({
        as.data.frame(ode(
          func = VassFox_2Cvar,
          y = State,
          parms = pars,
          times = seq(0, 1)
        ))
      },
      timeout = 2,
      events = list(func = eventfun))
    },
    error = function(e) {
      an.error.occured <<- TRUE
    },
    TimeoutException = function(e) {
      an.error.occured <<- TRUE
    })
    
    
    #if the call to the ode solver produced and error this is a bad parameter set
    if (an.error.occured) {
      return(Inf)
    }
    
    # print(VF)
    
    # Update results matrix
    mat[t, ] <- c(VF[2, 2], VF[2, 3], VF[2, 4], VF[2, 5], VF[2, 6])
    
    
  }
  
  
  dat <-
    data.frame(
      pred1 = mat[, 1],
      con1 = mat[, 2],
      con2 = mat[, 3],
      con3 = mat[, 4],
      res = mat[, 5]
    )
  suppressMessages({
    d <- melt(dat)
  })
  
  d$time <- c(rep(seq(1, 5000), 5))
  
  
  #if something went wrong and the parameter set does work
  if (any(is.na(d$value))) {
    return(Inf)
  }
  
  d <- d[d$time > 1000,]
  
  
  #anytime a population crashes below this threshhold it get penalized in the scoring function
  Crash_penality <- sum(d$value < 0.001)
  
  
  return(as.numeric(Crash_penality))
  
}

#attempt to see if compiling makes things any faster
optimize_function_cmp <- cmpfun(optimize_function)
#compiling the function does seem to make it run a little faster, either way it doesnt hurt.

#lower and upper bounds on parameters that we are estimating.
lower <- c(rep(0, 6))
high <- c(rep(1, 6))


#function will print top scoring individual in each generation and its parameters. It will run for 200 iterations maximum, but it normally gets to zero in 30 or 40 iterations.
opt_out <- DEoptim(
  optimize_function_cmp,
  lower = lower,
  upper = high,
  DEoptim.control(
    NP = 100,
    itermax = 200,
    strategy = 3,
    c = 0.05,
    parallelType = 1,
    packages = c("deSolve", "reshape", "tidyverse", "R.utils","lubridate")
  )
)

#you may need to restart you R session to kill zombie threads if you abort the optimization.
