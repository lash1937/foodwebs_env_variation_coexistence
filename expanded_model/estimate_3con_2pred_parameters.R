#mplementation of the Vasseur and Fox Model with temporal variation in mortality rates
# Compare positive, no, and negative autocorrelation
# For tree species mortality rates through time (e.g. from Vasseur and Fox appendix)
# estimate parameters for 3 model

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
library(compiler)




optimize_function<-function(pars){

  pars_from_3_system<-c(0.262571  ,  0.862126  ,  0.227814  ,  0.881397 ,   0.065556 ,   0.533164) 
      
      VassFox_2P_3C <- function(Time, State, Pars) {
        with(as.list(c(State, Pars)), {
          
          #variables to make math easier to read
          
          all_c_p1 <-
            ((O_P1_C1 * C1) +  (O_P1_C2 * C2) + (1 - (O_P1_C1 + O_P1_C2)) * C3)
          
          
          all_c_p2 <-
            ((O_P2_C1 * C1) +  (O_P2_C2 * C2) + (1 - (O_P2_C1 + O_P2_C2)) * C3)
          
          both_preds_eat_C1 <-
            ((O_P1_C1 * J_P1 * P1 * C1) / (all_c_p1 + C_0)) + ((O_P2_C1 * J_P2 * P2 * C1) / (all_c_p2 + C_0_P2))
          
          both_preds_eat_C2 <-
            ((O_P1_C2 * J_P1 * P1 * C2) / (all_c_p1  + C_0)) + ((O_P2_C2 * J_P2 * P2 * C2) / (all_c_p2 + C_0_P2))
          
          both_preds_eat_C3 <-
            (((1 - (
              O_P1_C1 + O_P1_C2
            )) * J_P1 * P1 * C3) / (all_c_p1 + C_0)) + (((1 - (
              O_P2_C1 + O_P2_C2
            )) * J_P2 * P2 * C3) / (all_c_p2 + C_0_P2))
          
          
          
          dP1 = -(M_P1 * P1) + (((J_P1 * P1) * (
            (O_P1_C1 * C1) + (O_P1_C2 * C2) + (1 - (O_P1_C1 + O_P1_C2)) * C3
          )) / (all_c_p1  + C_0))
          
          
          dP2 = -(M_P2 * P2) + (((J_P2 * P2) * (
            (O_P2_C1 * C1) + (O_P2_C2 * C2) + (1 - (O_P2_C1 + O_P2_C2)) * C3
          )) / (all_c_p2 + C_0_P2))
          
          
          dC1 = -(M_C1 * C1) + ((O_C1_R * J_C1 * C1 * R) / (R + R_0_1)) - both_preds_eat_C1
          dC2 = -(M_C2 * C2) + ((O_C2_R * J_C2 * C2 * R) / (R + R_0_2)) - both_preds_eat_C2
          dC3 = -(M_C3 * C3) + ((O_C3_R * J_C3 * C3 * R) / (R + R_0_3)) - both_preds_eat_C3
          
          dR = r * R * (1 - (R / K)) - ((O_C1_R * J_C1 * C1 * R) / (R + R_0_1)) - ((O_C2_R * J_C2 * C2 * R) / (R + R_0_2)) -  ((O_C3_R * J_C3 * C3 * R) / (R + R_0_3))
          
          return(list(c(dP1, dP2, dC1, dC2, dC3, dR)))
          
        })
        
      }
  # strength of env. on mortality rate
  #sigma <- 0.55
  sigma<-0
  # cross-correlation of C1 and C2
  ro <- 0
  
  z <- matrix(data = rnorm(n = 2*5000, mean = 0, sd = 1), nrow = 2)
  cholesky <- matrix(data = c(sigma^2, ro * (sigma ^2), ro * (sigma ^2), sigma ^2), 
                     nrow = 2)
  g <- cholesky %*% z
  M_C1_temp <- 0.4 * exp(g[1,])
  M_C2_temp <- 0.2 * exp(g[2,])
  
  #redraw fluxs for C3, rather then make a 3x matrix, its easier to do this.
  z <- matrix(data = rnorm(n = 2*5000, mean = 0, sd = 1), nrow = 2)
  g <- cholesky %*% z
  
  M_C3_temp <- pars_from_3_system[1] * exp(g[1,])
  

  
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
  #
  J_C3<- pars_from_3_system[2]
  # predator1 ingestion rate
  J_P1 <- 0.4
  # predator1 mortality rate
  M_P1 <- 0.08
  # predator2 ingestion rate
  J_P2 <- pars[1]
  # predator2 mortality rate
  M_P2 <- pars[2]
  
  # half saturation constant
  R_0_1 <- 0.16129
  R_0_2 <- 0.9
  R_0_3 <- pars_from_3_system[3]
  
  
  C_0 <- 0.5
  # preference coefficient of pred 1
  O_P1_C1 <- pars_from_3_system[4]
  O_P1_C2 <- pars_from_3_system[5]
  
  
  O_P2_C1 <- pars[3]
  O_P2_C2 <- pars[4]
  
  C_0_P2<-pars[5]
  O_C1_R <- 1.0
  O_C2_R <- 0.98
  O_C3_R <- pars_from_3_system[6]
  time <- 5000
  
  if(O_P1_C1 + O_P1_C2 >.95 ){
    return(Inf)
  }
  
  if(O_P2_C1 + O_P2_C2 >.95 ){
    return(Inf)
  }
  
  State <-
    c(P1 =1,
      P2 =1,
      C1 = 1,
      C2 = 1,
      C3 =1,
      R = 1) # starting parameters
  
  mat <- matrix(data = NA,
                nrow = time,
                ncol = 6)
  mat[1,] <- State
  
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
        J_P1 = J_P1,
        J_P2 = J_P2,
        M_C1 = M_C1,
        M_C2 = M_C2,
        M_C3 = M_C3,
        M_P1 = M_P1,
        M_P2 = M_P2,
        R_0_1 = R_0_1,
        R_0_2 = R_0_2,
        R_0_3 = R_0_3,
        C_0 = C_0,
	      C_0_P2 = C_0_P2,
        O_P1_C1 = O_P1_C1,
        O_P1_C2 = O_P1_C2,
        O_P2_C1 = O_P2_C1,
        O_P2_C2 = O_P2_C2,
        O_C1_R = O_C1_R,
        O_C2_R = O_C2_R,
        O_C3_R = O_C3_R
      )
    
    # Udate state variables to output from last timestep
    State <- c(
      P1 = mat[t - 1, 1],
      P2 = mat[t - 1, 2],
      C1 = mat[t - 1, 3],
      C2 = mat[t - 1, 4],
      C3 = mat[t - 1, 5],
      R = mat[t - 1, 6]
    )
    
    
    an.error.occured <- FALSE
   
    
    #uf the eval produces errorrs or takes a long time don't bother just throw an error
    tryCatch({
      VF <- withTimeout( {
        as.data.frame(
          ode(
            func = VassFox_2P_3C,
            y = State,
            parms = pars,
            times = seq(0, 1)
          ) ) },
        timeout = 2,
        events = list(func = eventfun) )
    }, 
    error = function(e) {
      an.error.occured <<- TRUE
    },
    TimeoutException = function(e) { an.error.occured <<- TRUE
    }
    
    )
    
    
    #if the call to the ode solver produced and error this is a bad parameter set
    if (an.error.occured) {
      return(Inf)
    }
    
    #print(VF)
    
    # Update results matrix
    mat[t,] <- c(VF[2, 2], VF[2, 3], VF[2, 4], VF[2, 5], VF[2,6],VF[2,7])
    
    #print(mat[t,])
  }
  
  
  dat <-
    data.frame(
      pred1 = mat[, 1],
      pred2 = mat[, 2],
      con1 = mat[, 3],
      con2 = mat[, 4],
      con3 = mat[, 5],
      res = mat[, 6]
    )
  suppressMessages({
    d <- melt(dat)
  })
  
  d$time <- c(rep(seq(1, 5000), 6))
  
  if(any(is.na(d$value))){
    return(Inf)
  }
  
  d <- d[d$time > 1000, ]
  
 
  #anytime a population crashes below this threshhold it get penalized in the scoring function
  Crash_penality <- sum(d$value < 0.001)
  
  
  return(as.numeric(Crash_penality))

  
}





#attempt to see if compling makes things any faster
#system.time(optimize_function(rep(0.2,4)))
optimize_function_cmp<-cmpfun(optimize_function)
#system.time(optimize_function_cmp(rep(0.2,4)))
#print(x)
lower <- c(0,0,0.05,0.05,0.01)
high <- c(rep(1, 5))



opt_out<-DEoptim(
  optimize_function_cmp,
  lower = lower,
  upper = high,
  DEoptim.control(
    NP=100,
    itermax=500,
    strategy=3,
    c=0.05,
    parallelType = 1,
    packages = c("deSolve", "reshape", "tidyverse", "R.utils")
 )
)
