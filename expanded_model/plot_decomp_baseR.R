#plots results from mechanistic decomposition using base R functions rather then ggplot so that it looks like the other figures. It leaves a box on the left to draw in the system with some graphics editing package
#you can take out the empty box plots by hashing them out and changing the par argument to make room for 3 plots rather than 4.

library(here)

rm(list = ls())



#load in files
C1_final_mechanisms <-
  read.csv(file = "C1_final_mechanisms_3con_param1.csv")
C2_final_mechanisms <-
  read.csv(file = "C2_final_mechanisms_3con_param1.csv")
C3_final_mechanisms <-
  read.csv(file = "C3_final_mechanisms_3con_param1.csv")


#set min and max of plots
max_y <- .21
min_y <- -.12


C1_r_bar <- mean(C1_final_mechanisms[, 1])

C1_delta_0 <- mean(C1_final_mechanisms[, 2])

C1_delta_P <- mean(C1_final_mechanisms[, 3])

C1_delta_E <- mean(C1_final_mechanisms[, 4])

C1_delta_EP <- mean(C1_final_mechanisms[, 5])



C2_r_bar <- mean(C2_final_mechanisms[, 1])

C2_delta_0 <- mean(C2_final_mechanisms[, 2])

C2_delta_P <- mean(C2_final_mechanisms[, 3])

C2_delta_E <- mean(C2_final_mechanisms[, 4])

C2_delta_EP <- mean(C2_final_mechanisms[, 5])


C3_r_bar <- mean(C3_final_mechanisms[, 1])

C3_delta_0 <- mean(C3_final_mechanisms[, 2])

C3_delta_P <- mean(C3_final_mechanisms[, 3])

C3_delta_E <- mean(C3_final_mechanisms[, 4])

C3_delta_EP <- mean(C3_final_mechanisms[, 5])


C1_results_sigma1 <-
  c(C1_r_bar, C1_delta_0, C1_delta_P, C1_delta_E, C1_delta_EP)

C2_results_sigma1 <-
  c(C2_r_bar, C2_delta_0, C2_delta_P, C2_delta_E, C2_delta_EP)

C3_results_sigma1 <-
  c(C3_r_bar, C3_delta_0, C3_delta_P, C3_delta_E, C3_delta_EP)

# and calculate standard deviation

sdC1_r_bar <- sd(C1_final_mechanisms[, 1])

sdC1_delta_0 <- sd(C1_final_mechanisms[, 2])

sdC1_delta_P <- sd(C1_final_mechanisms[, 3])

sdC1_delta_E <- sd(C1_final_mechanisms[, 4])

sdC1_delta_EP <- sd(C1_final_mechanisms[, 5])



sdC2_r_bar <- sd(C2_final_mechanisms[, 1])

sdC2_delta_0 <- sd(C2_final_mechanisms[, 2])

sdC2_delta_P <- sd(C2_final_mechanisms[, 3])

sdC2_delta_E <- sd(C2_final_mechanisms[, 4])

sdC2_delta_EP <- sd(C2_final_mechanisms[, 5])


sdC3_r_bar <- sd(C3_final_mechanisms[, 1])

sdC3_delta_0 <- sd(C3_final_mechanisms[, 2])

sdC3_delta_P <- sd(C3_final_mechanisms[, 3])

sdC3_delta_E <- sd(C3_final_mechanisms[, 4])

sdC3_delta_EP <- sd(C3_final_mechanisms[, 5])


sdC1_results_sigma1 <-
  c(sdC1_r_bar,
    sdC1_delta_0,
    sdC1_delta_P,
    sdC1_delta_E,
    sdC1_delta_EP)

sdC2_results_sigma1 <-
  c(sdC2_r_bar,
    sdC2_delta_0,
    sdC2_delta_P,
    sdC2_delta_E,
    sdC2_delta_EP)

sdC3_results_sigma1 <-
  c(sdC3_r_bar,
    sdC3_delta_0,
    sdC3_delta_P,
    sdC3_delta_E,
    sdC3_delta_EP)


# ----------------------------------------------------------------------------------------------------

# Plot results

quartz(width = 6, height = 6)

pdf(file = "sigma_variation_results_3con_2pred.pdf",
    width = 6,
    height = 6)



par(
  mfrow = c(2, 4),
  oma = c(4, 2, 1.5, 1),
  mar = c(1.5, .75, 3, 5)
)


barplot(
  0,
  type = 'n',
  axes = FALSE,
  ann = FALSE,
  ylim = c(min_y, max_y),
  width = 0
)
box(which = "plot", lty = "solid")
text(x = -.8, y = .2, "A)")
mtext(expression("System"),
      side = 3,
      outer = FALSE,
      adj = 0.5)

par(mar = c(1.5, .15, 3, 2))


barplot(
  C1_results_sigma1,
  ylim = c(min_y, max_y),
  xlab = "",
  ylab = "Growth Rate When Rare",
  
  col = c("#000000", "#0072B2", "#56B4E9", "#CC79A7", "#E69F00"),
  yaxt = "n",
  besides = FALSE,
  axisnames = FALSE
)

abline(h = 0)

axis(
  side = 2,
  at = seq(-0.1, 0.4, by = .1),
  labels = TRUE,
  tick = TRUE
)

axis(
  side = 1,
  at = c(.7, 1.9, 3.1, 4.3, 5.5),
  lab = c(
    "a" = expression(bar("r")[i] - bar("r")[r]) ,
    
    "b" = expression(Delta[i] ^
                       0),
    
    "c" = expression(Delta[i] ^
                       P),
    
    "d" = expression(Delta[i] ^
                       E),
    
    "e" = expression(Delta[i] ^
    {
      E * P
    })
  )
)





box(which = "plot", lty = "solid")

title(
  xlab = "",
  ylab = "Growth Rate When Rare",
  outer = TRUE,
  line = -8.5,
  cex.lab = 1.5
)
mtext(
  expression("Competitor 1"),
  side = 3,
  outer = FALSE,
  adj = 0.5
)

text(x = .4, y = .2, "B)")

arrows(
  x0 = c(.7, 1.9, 3.1, 4.3, 5.5),
  y0 = C1_results_sigma1 - sdC1_results_sigma1,
  
  x1 = c(.7, 1.9, 3.1, 4.3, 5.5),
  y1 = C1_results_sigma1 + sdC1_results_sigma1,
  length = .05,
  
  angle = 90,
  col = c("black"),
  code = 3
)



barplot(
  C2_results_sigma1,
  ylim = c(min_y, max_y),
  xlab = "",
  ylab = c(""),
  
  col = c("#000000", "#0072B2", "#56B4E9", "#CC79A7", "#E69F00"),
  yaxt = "n"
)

abline(h = 0)

axis(
  side = 2,
  at = seq(-0.1, 0.4, by = .1),
  labels = FALSE,
  tick = TRUE
)
axis(
  side = 1,
  at = c(.7, 1.9, 3.1, 4.3, 5.5),
  lab = c(
    "a" = expression(bar("r")[i] - bar("r")[r]) ,
    
    "b" = expression(Delta[i] ^
                       0),
    
    "c" = expression(Delta[i] ^
                       P),
    
    "d" = expression(Delta[i] ^
                       E),
    
    "e" = expression(Delta[i] ^
    {
      E * P
    })
  )
)



axis(
  side = 2,
  at = seq(-0.15, 0.4, by = .1),
  lab = c("", "", "")
)

box(which = "plot", lty = "solid")

arrows(
  x0 = c(.7, 1.9, 3.1, 4.3, 5.5),
  y0 = C2_results_sigma1 - sdC2_results_sigma1,
  
  x1 = c(.7, 1.9, 3.1, 4.3, 5.5),
  y1 = C2_results_sigma1 + sdC2_results_sigma1,
  length = .05,
  
  angle = 90,
  col = c("black"),
  code = 3
)

text(x = .4, y = .2, "C)")


mtext(
  expression("Competitor 2"),
  side = 3,
  outer = FALSE,
  adj = 0.5
)


#--- comp 3

barplot(
  C3_results_sigma1,
  ylim = c(min_y, max_y),
  xlab = "",
  ylab = c(""),
  
  col = c("#000000", "#0072B2", "#56B4E9", "#CC79A7", "#E69F00"),
  yaxt = "n"
)

abline(h = 0)
axis(
  side = 2,
  at = seq(-0.1, 0.4, by = .1),
  labels = FALSE,
  tick = TRUE
)

axis(
  side = 1,
  at = c(.7, 1.9, 3.1, 4.3, 5.5),
  lab = c(
    "a" = expression(bar("r")[i] - bar("r")[r]) ,
    
    "b" = expression(Delta[i] ^
                       0),
    
    "c" = expression(Delta[i] ^
                       P),
    
    "d" = expression(Delta[i] ^
                       E),
    
    "e" = expression(Delta[i] ^
    {
      E * P
    })
  )
)



axis(
  side = 2,
  at = seq(-0.1, 0.4, by = .1),
  lab = c("", "", "")
)

box(which = "plot", lty = "solid")

arrows(
  x0 = c(.7, 1.9, 3.1, 4.3, 5.5),
  y0 = C3_results_sigma1 - sdC3_results_sigma1,
  
  x1 = c(.7, 1.9, 3.1, 4.3, 5.5),
  y1 = C3_results_sigma1 + sdC3_results_sigma1,
  length = .05,
  
  angle = 90,
  col = c("black"),
  code = 3
)

text(x = .4, y = .2, "D)")



mtext(
  expression("Competitor 3"),
  side = 3,
  outer = FALSE,
  adj = 0.5
)

#----other param set


#change this if you want to plot the both of the 2 con and 1 pred model.



C1_final_mechanisms <-
  read.csv(file = "C1_final_mechanisms_3con_2pred.csv")
C2_final_mechanisms <-
  read.csv(file = "C2_final_mechanisms_3con_2pred.csv")
C3_final_mechanisms <-
  read.csv(file = "C3_final_mechanisms_3con_2pred.csv")




C1_r_bar <- mean(C1_final_mechanisms[, 1])

C1_delta_0 <- mean(C1_final_mechanisms[, 2])

C1_delta_P <- mean(C1_final_mechanisms[, 3])

C1_delta_E <- mean(C1_final_mechanisms[, 4])

C1_delta_EP <- mean(C1_final_mechanisms[, 5])



C2_r_bar <- mean(C2_final_mechanisms[, 1])

C2_delta_0 <- mean(C2_final_mechanisms[, 2])

C2_delta_P <- mean(C2_final_mechanisms[, 3])

C2_delta_E <- mean(C2_final_mechanisms[, 4])

C2_delta_EP <- mean(C2_final_mechanisms[, 5])


C3_r_bar <- mean(C3_final_mechanisms[, 1])

C3_delta_0 <- mean(C3_final_mechanisms[, 2])

C3_delta_P <- mean(C3_final_mechanisms[, 3])

C3_delta_E <- mean(C3_final_mechanisms[, 4])

C3_delta_EP <- mean(C3_final_mechanisms[, 5])


C1_results_sigma1 <-
  c(C1_r_bar, C1_delta_0, C1_delta_P, C1_delta_E, C1_delta_EP)

C2_results_sigma1 <-
  c(C2_r_bar, C2_delta_0, C2_delta_P, C2_delta_E, C2_delta_EP)

C3_results_sigma1 <-
  c(C3_r_bar, C3_delta_0, C3_delta_P, C3_delta_E, C3_delta_EP)

# and calculate standard deviation

sdC1_r_bar <- sd(C1_final_mechanisms[, 1])

sdC1_delta_0 <- sd(C1_final_mechanisms[, 2])

sdC1_delta_P <- sd(C1_final_mechanisms[, 3])

sdC1_delta_E <- sd(C1_final_mechanisms[, 4])

sdC1_delta_EP <- sd(C1_final_mechanisms[, 5])



sdC2_r_bar <- sd(C2_final_mechanisms[, 1])

sdC2_delta_0 <- sd(C2_final_mechanisms[, 2])

sdC2_delta_P <- sd(C2_final_mechanisms[, 3])

sdC2_delta_E <- sd(C2_final_mechanisms[, 4])

sdC2_delta_EP <- sd(C2_final_mechanisms[, 5])


sdC3_r_bar <- sd(C3_final_mechanisms[, 1])

sdC3_delta_0 <- sd(C3_final_mechanisms[, 2])

sdC3_delta_P <- sd(C3_final_mechanisms[, 3])

sdC3_delta_E <- sd(C3_final_mechanisms[, 4])

sdC3_delta_EP <- sd(C3_final_mechanisms[, 5])


sdC1_results_sigma1 <-
  c(sdC1_r_bar,
    sdC1_delta_0,
    sdC1_delta_P,
    sdC1_delta_E,
    sdC1_delta_EP)

sdC2_results_sigma1 <-
  c(sdC2_r_bar,
    sdC2_delta_0,
    sdC2_delta_P,
    sdC2_delta_E,
    sdC2_delta_EP)

sdC3_results_sigma1 <-
  c(sdC3_r_bar,
    sdC3_delta_0,
    sdC3_delta_P,
    sdC3_delta_E,
    sdC3_delta_EP)


# ----------------------------------------------------------------------------------------------------

# Plot results
par(mar = c(1.5, .75, 3, 5))


barplot(
  0,
  type = 'n',
  axes = FALSE,
  ann = FALSE,
  ylim = c(min_y, max_y),
  width = 0
)
box(which = "plot", lty = "solid")
text(x = -.8, y = .2, "E)")
mtext(expression("System"),
      side = 3,
      outer = FALSE,
      adj = 0.5)


par(mar = c(1.5, .15, 3, 2))

x <-
  barplot(
    C1_results_sigma1,
    ylim = c(min_y, max_y),
    xlab = "",
    ylab = c("Growth Rate When Rare"),
    
    col = c("#000000", "#0072B2", "#56B4E9", "#CC79A7", "#E69F00"),
    yaxt = "n"
  )

abline(h = 0)

axis(
  side = 2,
  at = seq(-0.1, 0.4, by = .1),
  labels = TRUE,
  tick = TRUE
)

axis(
  side = 1,
  at = c(.7, 1.9, 3.1, 4.3, 5.5),
  lab = c(
    "a" = expression(bar("r")[i] - bar("r")[r]) ,
    
    "b" = expression(Delta[i] ^
                       0),
    
    "c" = expression(Delta[i] ^
                       P),
    
    "d" = expression(Delta[i] ^
                       E),
    
    "e" = expression(Delta[i] ^
    {
      E * P
    })
  )
)





box(which = "plot", lty = "solid")

mtext(
  expression("Competitor 1"),
  side = 3,
  outer = FALSE,
  adj = 0.5
)


text(x = .4, y = .2, "F)")

arrows(
  x0 = c(.7, 1.9, 3.1, 4.3, 5.5),
  y0 = C1_results_sigma1 - sdC1_results_sigma1,
  
  x1 = c(.7, 1.9, 3.1, 4.3, 5.5),
  y1 = C1_results_sigma1 + sdC1_results_sigma1,
  length = .05,
  
  angle = 90,
  col = c("black"),
  code = 3
)



barplot(
  C2_results_sigma1,
  ylim = c(min_y, max_y),
  xlab = "",
  ylab = c(""),
  
  col = c("#000000", "#0072B2", "#56B4E9", "#CC79A7", "#E69F00"),
  yaxt = "n"
)

abline(h = 0)

axis(
  side = 2,
  at = seq(-0.1, 0.4, by = .1),
  labels = FALSE,
  tick = TRUE
)

axis(
  side = 1,
  at = c(.7, 1.9, 3.1, 4.3, 5.5),
  lab = c(
    "a" = expression(bar("r")[i] - bar("r")[r]) ,
    
    "b" = expression(Delta[i] ^
                       0),
    
    "c" = expression(Delta[i] ^
                       P),
    
    "d" = expression(Delta[i] ^
                       E),
    
    "e" = expression(Delta[i] ^
    {
      E * P
    })
  )
)



axis(
  side = 2,
  at = seq(-0.1, 0.4, by = .1),
  lab = c("", "", "")
)

box(which = "plot", lty = "solid")

arrows(
  x0 = c(.7, 1.9, 3.1, 4.3, 5.5),
  y0 = C2_results_sigma1 - sdC2_results_sigma1,
  
  x1 = c(.7, 1.9, 3.1, 4.3, 5.5),
  y1 = C2_results_sigma1 + sdC2_results_sigma1,
  length = .05,
  
  angle = 90,
  col = c("black"),
  code = 3
)

text(x = .4, y = .2, "G)")

mtext(
  "Mechanistic Partitioning",
  side = 1,
  outer = TRUE,
  adj = .65,
  line = 1.25
)



mtext(
  expression("Competitor 2"),
  side = 3,
  outer = FALSE,
  adj = 0.5
)


#--- comp 3

barplot(
  C3_results_sigma1,
  ylim = c(min_y, max_y),
  xlab = "",
  ylab = c(""),
  
  col = c("#000000", "#0072B2", "#56B4E9", "#CC79A7", "#E69F00"),
  yaxt = "n"
)

abline(h = 0)
axis(
  side = 2,
  at = seq(-0.1, 0.4, by = .1),
  labels = FALSE,
  tick = TRUE
)

axis(
  side = 1,
  at = c(.7, 1.9, 3.1, 4.3, 5.5),
  lab = c(
    "a" = expression(bar("r")[i] - bar("r")[r]) ,
    
    "b" = expression(Delta[i] ^
                       0),
    
    "c" = expression(Delta[i] ^
                       P),
    
    "d" = expression(Delta[i] ^
                       E),
    
    "e" = expression(Delta[i] ^
    {
      E * P
    })
  )
)



axis(
  side = 2,
  at = seq(-0.1, 0.4, by = .1),
  lab = c("", "", "")
)

box(which = "plot", lty = "solid")

arrows(
  x0 = c(.7, 1.9, 3.1, 4.3, 5.5),
  y0 = C3_results_sigma1 - sdC3_results_sigma1,
  
  x1 = c(.7, 1.9, 3.1, 4.3, 5.5),
  y1 = C3_results_sigma1 + sdC3_results_sigma1,
  length = .05,
  
  angle = 90,
  col = c("black"),
  code = 3
)

text(x = .4, y = .2, "H)")


mtext(
  expression("Competitor 3"),
  side = 3,
  outer = FALSE,
  adj = 0.5
)


dev.off()