# Dominic_liveCoding.R
# Day 1: Modelling for fundamental insight into disease transmission
# SACEMA-UNITID Epidemiological Modelling Workshop
# University of Narobi, Nairobi, Kenya, 8-11 April 2019
#
# Emanuel Dominic, 2019
#
# Code to reproduce Figure 1 of:
# Dushoff et al. (2004) Dynamical resonance can account for seasonality of influenza epidemics.
# PNAS 101(48): 16915–16916

rm(list=ls())

# Libraries ---------------------------------------------------------------
library(deSolve)

# Model functions ---------------------------------------------------------

sirs <- function(t,y,parms){
  # with function gives access to param values within the local environment
  with(as.list(c(y,parms)),{
    
    betat <- beta.t(beta0, beta1, t)
    
    dSdt <- (N-I-S)/LL - betat*I*S/N
    dIdt <- betat*I*S/N - I/DD #don't need to specify dRdt because population is fixed
    
    output <- c(dSdt, dIdt)
    return(list(output))
  })
}

beta.t <- function(beta0, beta1, t){
  beta0 * (1+ beta1*cos(2*pi*t))
}

endemicEq <- function(parms){
  with(as.list(parms),{
    S.star <- N/(beta0 * DD)
    I.star <- (N-S.star)/(1+LL/DD)
    return(c(S = S.star, I = I.star))
  })
}

R0.t <- function(beta0, beta1, t, D){
  beta.t(beta0, beta1, t) * D
}

intrinsicT <- function(parms){
  with(as.list(parms),{
    t <- 2*pi*sqrt(DD*LL/(beta0*DD - 1))
    return(t)
  })
}
# Variables and parameters ------------------------------------------------

fig1A <- c(N = 500000,
           LL = 4,
           DD = 0.02,
           beta0 = 500,
           beta1 = 0.04)

fig1B <- c(N = 500000,
           LL = 8,
           DD = 0.025,
           beta0 = 400,
           beta1 = 0.04)

init1A <- endemicEq(fig1A)
init1A #does not aad up to 500000 because R is omitted

init1B <- endemicEq(fig1B)
init1B #does not add up to 500000 because R is omitted

time.out <- seq(0,20, 0.05)


# Figure 1a: numerical integration and visualisation ----------------------

?lsoda

sirs.out.fig1A <- lsoda(y = init1A,
                        times = time.out,
                        func = sirs,
                        parms = fig1A)
head(sirs.out.fig1A)
sirs.out.fig1A.df <- data.frame(sirs.out.fig1A)

par(mfcol = c(2,1), mar=c(4,4,1,0)) #par combine multiple plots in 1 graph, mar to set margins

plot(x = sirs.out.fig1A.df$time, y = sirs.out.fig1A.df$I,
     xlab = "Time (years)",
     ylab = "Number of infected",
     xlim = c(10,20),
     ylim = c(0,4000),
     col = "red",
     lwd = 2,
     type = "l",
     yaxt = "none") #turn of y axis

axis(side = 2,seq(0,4000,4000/8), las=1)#las = 1 means ticks are horizontal


with(as.list(fig1A), {
  S.star.t <- N / R0.t(beta0, beta1, time.out, DD)
  I.star.t <- (N - S.star.t)/(1+LL/DD)
  lines(x= sirs.out.fig1A.df$time, y = I.star.t, lwd = 2)
})



# Figure 1b: numerical integration and visualisation ----------------------

sirs.out.fig1B <- lsoda(y = init1B,
                        times = time.out,
                        func = sirs,
                        parms = fig1B)
head(sirs.out.fig1B)
sirs.out.fig1B.df <- data.frame(sirs.out.fig1B)

plot(x = sirs.out.fig1B.df$time, y = sirs.out.fig1B.df$I,
     xlab = "Time (years)",
     ylab = "Number of infected",
     xlim = c(10,20),
     ylim = c(0,4000),
     col = "red",
     lwd = 2,
     type = "l",
     yaxt = "none")

axis(side = 2,seq(0,4000,4000/8), las=1)#las = 1 means ticks are horizontal

with(as.list(fig1B), {
  S.star.t <- N / R0.t(beta0, beta1, time.out, DD)
  I.star.t <- (N - S.star.t)/(1+LL/DD)
  lines(x= time.out, y = I.star.t, lwd = 2)
})


intrinsicT(fig1A)
intrinsicT(fig1B)
