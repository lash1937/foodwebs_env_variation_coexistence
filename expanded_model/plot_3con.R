#plots 3con model, variation can be turned off and on with the sigma parameter.

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


# this function is the same as the optimize function except it returns the data rather then scoring the data
optimize_function <- function(pars) {
  #This is a 3C model, I just didn't change the labels from the 2C model
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
  sigma <- 0.55
  #sigma <- 0
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
  # preference coefficient
  O_P_C1 <- pars[4]
  O_P_C2 <- pars[5]
  
  O_C1_R <- 1.0
  O_C2_R <- 0.98
  O_C3_R <- pars[6]
  time <- 5000
  
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
  return(dat)
  
  
}


#these parameters were estimated using estimate_3con_parameters.R, just comment out whichever one you dont want to draw

#parameter set 1
est <-
  c(0.097798 ,   0.670232  ,  0.297716  ,  0.332333,    0.001743  ,  0.926310)

#parameter set 2
#est<-c(0.142901 ,   0.917076  ,  0.328642,    0.262048 ,   0.012992,    0.919742)


#call the function to get dyanmics
p <- optimize_function(est)

#store everything for printing
dat <-
  data.frame(
    pred1 = p[, 1],
    con1 = p[, 2],
    con2 = p[, 3],
    con3 = p[, 4],
    res = p[, 5]
  )
suppressMessages({
  d <- melt(dat)
})

#add time to data
d$time <- c(rep(seq(1, 5000), 5))

#just take the last 4000 points
d <- d[d$time > 999,]

#make a plot for each species
g <-
  ggplot(d[d$variable == "pred1", ], aes(x = time, y = value)) + geom_line(aes(y = value), size = 1) + ylab("predator\ndensity") +
  scale_y_continuous(expand = c(0, .001), limits = c(0, max(d$value))) +
  scale_x_continuous(expand = c(0, .001)) + theme(plot.margin = unit(c(0, 1.2, 0, 0), "cm"))
g1 <-
  ggplot(d[d$variable == "con1", ], aes(x = time, y = value)) + geom_line(aes(y = value), size = 1) + ylab("competitor 1\ndensity") +
  scale_y_continuous(expand = c(0, .001), limits = c(0, max(d$value))) +
  scale_x_continuous(expand = c(0, .001)) + theme(plot.margin = unit(c(0, 1.2, 0, 0), "cm"))
g2 <-
  ggplot(d[d$variable == "con2", ], aes(x = time, y = value)) + geom_line(aes(y = value), size = 1) + ylab("competitor 2\ndensity") +
  scale_y_continuous(expand = c(0, .001), limits = c(0, max(d$value))) +
  scale_x_continuous(expand = c(0, .001)) + theme(plot.margin = unit(c(0, 1.2, 0, 0), "cm"))
g3 <-
  ggplot(d[d$variable == "con3", ], aes(x = time, y = value)) + geom_line(aes(y = value), size = 1) + ylab("competitor 3\ndensity") +
  scale_y_continuous(expand = c(0, .001), limits = c(0, max(d$value))) +
  scale_x_continuous(expand = c(0, .001)) + theme(plot.margin = unit(c(0, 1.2, 0, 0), "cm"))
g4 <-
  ggplot(d[d$variable == "res", ], aes(x = time, y = value)) + geom_line(aes(y = value), size = 1) + ylab("resource\ndensity") +
  scale_y_continuous(expand = c(0, .001), limits = c(0, max(d$value))) +
  scale_x_continuous(expand = c(0, .001)) + theme(plot.margin = unit(c(0, 1.2, 0, 0), "cm"))


#arrange each plot
pp <- plot_grid(
  g,
  g1,
  g2,
  g3,
  g4,
  labels = c("A", "B", "C", "D", "E"),
  align = 'vh',
  hjust = -1,
  nrow = 5
)


#save the plot
save_plot(
  "3con_dynamics.pdf",
  plot = pp,
  base_height = 12,
  units = "in"
)