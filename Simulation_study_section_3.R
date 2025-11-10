# Cleaning objects from the workplace
rm(list=ls())

# Packages (may first require installations: install.packages())
library(tidyverse) # for data manipulation
library(tidyr) # for data manipulation
library(scales) # to manipulate the internal scaling infrastructure used by ggplot2
library(RColorBrewer) # for colors and palettes
library(gridExtra) # to show multiple plots
library(cowplot) # to show multiple plots
library(deSolve) # to solve a system of ordinary differential equations
library(boot) # for the inverse logit function 'inv.logit'

setwd("C:/Users/wb9928/OneDrive - Princeton University/Desktop/OMEGA/Data_and_codes")

# Execute functions from external R script
source("Functions.R")

Sys.setlocale("LC_TIME", "English")

nA <- 100
nB <- 100

q0 <- 10/100

I0w <- 1E-7
I0m <- q0/(1-q0)*I0w

init <- c("SA"=nA-I0w-I0m, "IAw"=I0w, "IAm"=I0m, "RA"=0,
          "SB"=nB-I0w, "IBw"=I0w, "IBm"=0, "RB"=0)

tM <- 10
tf <- tM + 40
times <- seq(0,tf,1)

bAw <- bBw <- bAm <- bBm <- 0.15
gAw <- gBw <- 0.2
gAm <- gBm <- gAw/2

s <- -(gAm-gAw) # true selection

mA <- mB <- 0

parms <- c("bAw"=bAw, "bBw"=bBw, "bAm"=bAm, "bBm"=bBm,
           "gAw"=gAw, "gBw"=gBw, "gAm"=gAm, "gBm"=gBm,
           "mA"=mA, "mB"=mB)

simul_init <- simulate(init, times=times[times<=tM], parms=parms, approx.WM=FALSE)

init2 <- simul_init[nrow(simul_init),2:9] %>% unlist

m <- c(0,0.01,0.1,0.5)
lnQ0 <- 0:3
sigma <- 0.25

simul_data <- data.frame()
estim_rel_err <- data.frame()

for(i in m){
  
  parms[["mA"]] <- parms[["mB"]] <- i
  
  for(j in lnQ0){
    
    init2[["IBm"]] <- inv.logit(logit(init2[["IAm"]]/(init2[["IAw"]]+init2[["IAm"]]))-j)*init2[["IBw"]]
    
    simul <- simulate(init2, times=times[times>=tM], parms=parms, approx.WM=FALSE)
    
    simul$logit_q_obs <- simul$logit_q + rnorm(nrow(simul), 0, sigma)
    simul$logit_qA_obs <- simul$logit_qA + rnorm(nrow(simul), 0, sigma)
    simul$logit_qB_obs <- simul$logit_qB + rnorm(nrow(simul), 0, sigma)
    
    simul_data <- simul_data %>%
      rbind(data.frame("time"=simul$time-tM,
                       "logit_q"=c(simul$logit_q, simul$logit_qA, simul$logit_qB),
                       "logit_q_obs"=c(simul$logit_q_obs, simul$logit_qA_obs, simul$logit_qB_obs),
                       "pop"=rep(c("Metapopulation", "Pop. A", "Pop. B"), each=nrow(simul)),
                       "m"=i, "lnQ0"=j))
    
    reg <- lm(logit_q_obs~time, data=simul)
    regA <- lm(logit_qA_obs~time, data=simul)
    regB <- lm(logit_qB_obs~time, data=simul)
    
    CI95 <- confint(reg, level=0.95)[2,]
    CI95A <- confint(regA, level=0.95)[2,]
    CI95B <- confint(regB, level=0.95)[2,]
    
    estim_rel_err <- estim_rel_err %>%
      rbind(cbind(
        data.frame("rel_err"=c(reg$coefficients[["time"]], regA$coefficients[["time"]], regB$coefficients[["time"]]),
                   "lower"=c(CI95[1], CI95A[1], CI95B[1]), "upper"=c(CI95[2], CI95A[2], CI95B[2]))/s-1,
        "pop"=c("Metapopulation", "Pop. A", "Pop. B"), "m"=i, "lnQ0"=j))
  }
}

pd <- position_dodge(width=0.4)
colors <- c('#FDB4C5', '#BF9035', '#818332', '#104261')

(Fig1 <- plot_grid(
  
  simul_data %>%
    mutate(lnQ0=factor(lnQ0, levels=lnQ0, labels=paste("ln(\U1D4AC(0)) =", lnQ0))) %>%
    ggplot(aes(x=time, col=as.factor(m))) +
    facet_grid(pop~lnQ0) +
    labs(x="Time", y="Variant logit-frequency") +
    geom_line(aes(y=logit_q), linewidth=0.3) +
    geom_point(aes(y=logit_q_obs), cex=0.5) +
    scale_color_manual(values=colors) +
    theme_bw() + theme(legend.position='none'),

  ggplot(estim_rel_err, aes(x=lnQ0, col=as.factor(m))) +
    geom_hline(yintercept=0, lty="dashed", linewidth=0.3) +
    facet_grid(pop~.) +
    labs(x="Initial log-differentiation, ln(\U1D4AC(0))",
         y="Relative error in selection estimation",
         col=expression(mu), fill=expression(mu)) +
    geom_segment(aes(y=lower, yend=upper), linewidth=0.5, position=pd) +
    geom_point(aes(y=rel_err, fill=as.factor(m)), shape=21, col="black", size=1.5, position=pd) +
    scale_y_continuous(labels=percent_format(), breaks=seq(-0.6,0.6,0.2)) +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    theme_bw(),
  ncol=2, rel_widths=c(1,0.9)))

################################################################################

# Pop. B is 10 x smaller than pop. A

nA <- 100
nB <- 10

init <- c("SA"=nA-I0w-I0m, "IAw"=I0w, "IAm"=I0m, "RA"=0,
          "SB"=nB-I0w, "IBw"=I0w, "IBm"=0, "RB"=0)

simul_init <- simulate(init, times=times[times<=tM], parms=parms, approx.WM=FALSE)

init2 <- simul_init[nrow(simul_init),2:9] %>% unlist

simul_data <- data.frame()
estim_rel_err <- data.frame()

for(i in m){
  
  parms[["mA"]] <- parms[["mB"]] <- i
  
  for(j in lnQ0){
    
    init2[["IBm"]] <- inv.logit(logit(init2[["IAm"]]/(init2[["IAw"]]+init2[["IAm"]]))-j)*init2[["IBw"]]
    
    simul <- simulate(init2, times=times[times>=tM], parms=parms, approx.WM=FALSE)
    
    simul$logit_q_obs <- simul$logit_q + rnorm(nrow(simul), 0, sigma)
    simul$logit_qA_obs <- simul$logit_qA + rnorm(nrow(simul), 0, sigma)
    simul$logit_qB_obs <- simul$logit_qB + rnorm(nrow(simul), 0, sigma)
    
    simul_data <- simul_data %>%
      rbind(data.frame("time"=simul$time-tM,
                       "logit_q"=c(simul$logit_q, simul$logit_qA, simul$logit_qB),
                       "logit_q_obs"=c(simul$logit_q_obs, simul$logit_qA_obs, simul$logit_qB_obs),
                       "pop"=rep(c("Metapopulation", "Pop. A", "Pop. B"), each=nrow(simul)),
                       "m"=i, "lnQ0"=j))
    
    reg <- lm(logit_q_obs~time, data=simul)
    regA <- lm(logit_qA_obs~time, data=simul)
    regB <- lm(logit_qB_obs~time, data=simul)
    
    CI95 <- confint(reg, level=0.95)[2,]
    CI95A <- confint(regA, level=0.95)[2,]
    CI95B <- confint(regB, level=0.95)[2,]
    
    estim_rel_err <- estim_rel_err %>%
      rbind(cbind(
        data.frame("rel_err"=c(reg$coefficients[["time"]], regA$coefficients[["time"]], regB$coefficients[["time"]]),
                   "lower"=c(CI95[1], CI95A[1], CI95B[1]), "upper"=c(CI95[2], CI95A[2], CI95B[2]))/s-1,
        "pop"=c("Metapopulation", "Pop. A", "Pop. B"), "m"=i, "lnQ0"=j))
  }
}

pd <- position_dodge(width=0.3)

(Fig2 <- plot_grid(
  
  simul_data %>%
    mutate(lnQ0=factor(lnQ0, levels=lnQ0, labels=paste("ln(\U1D4AC(0)) =", lnQ0))) %>%
    ggplot(aes(x=time, col=as.factor(m))) +
    facet_grid(pop~lnQ0) +
    labs(x="Time", y="Variant logit-frequency") +
    geom_line(aes(y=logit_q), linewidth=0.3) +
    geom_point(aes(y=logit_q_obs), cex=0.5) +
    scale_color_manual(values=colors) +
    theme_bw() + theme(legend.position='none'),
  
  ggplot(estim_rel_err, aes(x=lnQ0, col=as.factor(m))) +
    geom_hline(yintercept=0, lty="dashed", linewidth=0.3) +
    facet_grid(pop~.) +
    labs(x="Initial log-differentiation, ln(\U1D4AC(0))",
         y="Relative error in selection estimation",
         col=expression(mu), fill=expression(mu)) +
    geom_segment(aes(y=lower, yend=upper), linewidth=0.5, position=pd) +
    geom_point(aes(y=rel_err, fill=as.factor(m)), shape=21, col="black", size=1.5, position=pd) +
    scale_y_continuous(labels=percent_format(), breaks=seq(-0.6,0.6,0.2)) +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    theme_bw(),
  ncol=2, rel_widths=c(1,0.9)))

################################################################################

# Final plot

plot_grid(
  ggdraw() + draw_label(expression('A)'~N^B/N^A~'='~1), size=15, hjust=0, x=0.05),
  Fig1,
  ggdraw() + draw_label(expression('B)'~N^B/N^A~'='~0.1), size=15, hjust=0, x=0.05),
  Fig2,
  ncol=1, rel_heights=rep(c(0.1,0.9))
) # Landscape 9 x 8

ggsave("Relative_errors_selection.png", width=9.5, height=9, dpi=1000)
