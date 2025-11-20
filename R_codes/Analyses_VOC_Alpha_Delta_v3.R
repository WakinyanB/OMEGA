rm(list=ls()) # Cleaning objects from the workplace

# Packages
library(tidyverse)
library(readODS)
library(plyr)
library(lubridate)
library(boot)
library(cowplot)
library(scales)
library(lme4)
library(ggpubr)
library(RColorBrewer)
library(grid)
library(deSolve)
library(ggtext)

setwd("...")

source("Functions.R")

regions <- c("East Midlands", "East of England", "London", "North East", "North West",
             "South East", "South West", "West Midlands", "Yorkshire and Humber")

colors <- c("#e7c459","#5fbb55","#da6967","#4c97cc",
            "#bc67da","#daa167","#da67ab","#67c3da","#8a96ac")

mytheme <- theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=11),
        axis.title.y=element_text(size=16), axis.text.y=element_text(size=13),
        plot.title=element_text(hjust=0.5, size=20, face='italic'),
        legend.background=element_rect(color='black', linewidth=0.1),
        legend.title=element_text(size=11), legend.text=element_text(size=11))

mytheme2 <- theme_bw() +
  theme(axis.title.x=element_text(size=15), axis.text.x=element_text(size=13),
        axis.title.y=element_text(size=15), axis.text.y=element_text(size=11),
        legend.position='none')

# Import and process data

## Alpha variant

Alpha <- read_ods("Data/VOC_202012_01_Technical_Briefing_5_Data_England.ods",
                  sheet=6, col_names=TRUE)
Alpha$week <- dmy(Alpha$week)

Alpha_All <- Alpha %>% ddply(~week, function(X){
  return(apply(X[,colnames(X) %in%
                   c("n_Confirmed S-gene","n_Confirmed SGTF","n_Total")],2,sum))
})

Alpha$q <- Alpha$`percent_Confirmed SGTF`/100
Alpha_All$q <- Alpha_All$`n_Confirmed SGTF`/Alpha_All$n_Total

## Delta variant

Delta <- readRDS("Data/Delta_data.rds")
Delta <- Delta[-which(Delta$phecname==""),]
Delta$date <- ymd(Delta$date)
Delta$week <- floor_date(Delta$date, unit="week", week_start=1)

Delta2 <- Delta %>% ddply(~date+phecname, function(X){
  return(data.frame("N"=nrow(X),
                    "NDelta"=sum(X$Delta)))
  })

Delta2_All <- Delta %>% ddply(~date, function(X){
  return(data.frame("N"=nrow(X), "NDelta"=sum(X$Delta)))})

Delta2$q <- Delta2$NDelta/Delta2$N
Delta2_All$q <- Delta2_All$NDelta/Delta2_All$N

# Compute selection coefficient

s_Alpha <- Alpha %>% subset(q>=0.1 & q<=0.99) %>% ddply(~Region, function(X){
  te <- min(X$week)
  reg <- lm(logit(q)~week, X)
  slope <- reg$coefficients["week"]
  CI95 <- confint(reg, level=0.95)[2,]
  return(data.frame("te"=te, "slope"=slope, "CI_2.5"=CI95[1], "CI_97.5"=CI95[2]))
})

fit_Alpha_All <- Alpha_All %>% subset(q>=0.1 & q<=0.99) %>% lm(logit(q)~week, .)
CI95 <- confint(fit_Alpha_All, level=0.95)[2,]
s_Alpha_All <- data.frame("te"=min(subset(Alpha_All, q>=0.1 & q<=0.99)$week),
                          "slope"=fit_Alpha_All$coefficients["week"],
                          "CI_2.5"=CI95[1], "CI_97.5"=CI95[2])

s_Delta <- Delta2 %>% subset(q>=0.1 & q<=0.99) %>% ddply(~phecname, function(X){
  te <- min(X$date)
  reg <- lm(logit(q)~date, X)
  slope <- reg$coefficients["date"]
  CI95 <- confint(reg, level=0.95)[2,]
  return(data.frame("te"=te, "slope"=slope, "CI_2.5"=CI95[1], "CI_97.5"=CI95[2]))
})

fit_Delta_All <- Delta2_All %>% subset(q>=0.1 & q<=0.99) %>% lm(logit(q)~date, .)
CI95 <- confint(fit_Delta_All, level=0.95)[2,]
s_Delta_All <- data.frame("te"=min(subset(Delta2_All, q>=0.1 & q<=0.99)$date),
                          "slope"=fit_Delta_All$coefficients["date"],
                          "CI_2.5"=CI95[1], "CI_97.5"=CI95[2])
set.seed(1996)
jitter <- runif(9,-2,2)

# Fig. 1. Growth of the Alpha and Delta variants in England

plot_grid(
  plot_grid(
    
    Alpha %>% subset(q>=0.1 & q<=0.99) %>%
      ggplot(aes(x=week, y=logit(q), col=Region)) +
      geom_point(aes(size=n_Total)) +
      labs(title='Alpha', y="Variant logit-frequency", col="Region", size="Sample size") +
      scale_color_manual(values=colors) +
      scale_x_date(date_labels = "%Y %b.") +
      scale_y_continuous(limits=c(-2.5,5.5), breaks=seq(-2,6,2)) +
      scale_size_continuous(range=c(1,4)) +
      geom_smooth(method='lm', se=FALSE, linewidth=0.6) +
      mytheme +
      theme(legend.position=c(0.25,0.8)) +
      guides(col='none'),
    
    Delta2 %>% subset(q>=0.1 & q<=0.99) %>%
      ggplot(aes(x=date, y=logit(q), col=phecname)) +
      geom_point(aes(size=N)) +
      labs(title='Delta', y="Variant logit-frequency", size="Sample size") +
      scale_x_date(date_labels = "%Y %b.") +
      scale_y_continuous(limits=c(-2.5,5.5), breaks=seq(-2,6,2)) +
      scale_color_manual(values=colors) +
      scale_size_continuous(range=c(0.5,2), breaks=c(500,1000,2000)) +
      geom_smooth(method='lm', se=FALSE, linewidth=0.6) +
      mytheme + theme(axis.title.y=element_blank(), legend.position=c(0.25,0.8)) +
      guides(col='none'),
    
    ggplot(s_Alpha, aes(x=te+jitter, y=slope, col=Region)) +
      geom_rect(aes(xmin=ymd("2020-10-20"), xmax=ymd("2020-12-25"),
                    ymin=s_Alpha_All$CI_2.5, ymax=s_Alpha_All$CI_97.5), col=NA, fill="grey90") +
      geom_hline(yintercept=s_Alpha_All$slope) +
      geom_point(size=3.5, pch=18) +
      geom_segment(aes(y=`CI_2.5`, yend=`CI_97.5`), linewidth=0.75) +
      labs(x="\nEmergence date", y=expression(Estimated~slope~(day^{-1}))) +
      scale_x_date(limits=ymd(c("2020-10-20","2020-12-25")),
                   breaks=ymd(c("2020-11-15","2020-12-15")),
                   expand=c(0,0),
                   date_labels = "%b. %d") +
      scale_y_continuous(limits=c(0.05,0.15), breaks=seq(0.05,0.15,0.02)) +
      scale_color_manual(values=colors) +
      mytheme2,
    
    ggplot(s_Delta, aes(x=te+jitter, y=slope, col=phecname)) +
      geom_rect(aes(xmin=ymd("2021-04-01"), xmax=ymd("2021-05-20"),
                    ymin=s_Delta_All$CI_2.5, ymax=s_Delta_All$CI_97.5), col=NA, fill="grey90") +
      geom_hline(yintercept=s_Delta_All$slope) +
      geom_point(size=3.5, pch=18) +
      geom_segment(aes(y=`CI_2.5`, yend=`CI_97.5`), linewidth=0.75) +
      labs(x="\nEmergence date", y=expression(Estimated~slope~(day^{-1}))) +
      scale_x_date(limits=ymd(c("2021-04-01","2021-05-20")),
                   breaks=ymd(c("2021-04-15","2021-05-15")),
                   expand=c(0,0),
                   date_labels = "%b. %d") +
      scale_y_continuous(limits=c(0.05,0.15), breaks=seq(0.05,0.15,0.02)) +
      scale_color_manual(values=colors) +
      mytheme2 + theme(axis.title.y=element_blank()),
    
    ncol=2, nrow=2, align='v', labels=LETTERS[1:4], label_size=20),
  
  get_legend(ggplot(Alpha, aes(x=week, y=q, col=Region)) +
               geom_point(cex=3) + geom_line(cex=0.7) +
               scale_color_manual(values=colors) + theme_minimal() +
               theme(legend.title=element_text(size=12),
                     legend.text=element_text(size=10))),
  ncol=2, rel_widths=c(0.8,0.2)) # 10.5 x 8

# ggsave("Fig1.png", width=10.5, height=8, dpi=690)

# Compute spatial differentiation

regions_alpha <- c("South East", "London", "East of England",
                   "South West", "North East", "West Midlands",
                   "East Midlands", "North West", "Yorkshire and Humber")

Diff_alpha <- Alpha %>% subset(q>=0.1 & q<=0.99) %>%
  ddply(~week, function(X){
    if(nrow(X)>1){
      tab <- data.frame()
      X <- X[match(regions_alpha, X$Region),] %>% drop_na
      for(i in 1:(nrow(X)-1)){
        for(j in (i+1):nrow(X)){
          tab <- tab %>% rbind(data.frame("lnQ"=logit(X$q[i])-logit(X$q[j]),
                                          "regions"=paste(X$Region[i], X$Region[j], sep="-")))
        }
      }
      return(tab)
    }
  }) %>% ddply(~regions, function(X){
  
  pair <- unique(X$regions)
  
  if(pair == "North East-West Midlands"){  # changing the direction of one pair
    X$lnQ <- -X$lnQ
    X$regions <- strsplit(pair, "-")[[1]] %>% rev %>% paste(collapse="-")
    return(X)
  }else{
    return(X)
  }
})

colnames(Diff_alpha)[1] <- "date"

Diff_alpha$time <- NA

for(i in 1:nrow(Diff_alpha)){
  Diff_alpha$time[i] <- Diff_alpha$date[i]-
    min(subset(Diff_alpha, regions==Diff_alpha$regions[i])$date) %>% as.numeric
}

Diff_alpha_summary <- Diff_alpha %>% ddply(~time, function(X){
  if(nrow(X)>=10){
    return(c("mean"=mean(X$lnQ), quantile(X$lnQ, probs=c(0.1,0.5,0.9))))
  }
})

regions_delta <- c("North West", "London", "South West",
                   "South East", "East of England", "East Midlands",
                   "West Midlands", "North East", "Yorkshire and Humber")

Diff_delta <- Delta2 %>% subset(q>=0.1 & q<0.99) %>%
  ddply(~date, function(X){
    if(nrow(X)>1){
      tab <- data.frame()
      X <- X[match(regions_delta, X$phecname),] %>% drop_na
      for(i in 1:(nrow(X)-1)){
        for(j in (i+1):nrow(X)){
          tab <- tab %>% rbind(data.frame("lnQ"=logit(X$q[i])-logit(X$q[j]),
                                          "regions"=paste(X$phecname[i], X$phecname[j], sep="-")))
        }
      }
      return(tab)
    }
  }) %>% ddply(~regions, function(X){  # changing the direction of some pairs
    
    pair <- unique(X$regions)
    
    if(pair %in% c("North West-London", "North West-South West", "London-South West",
                   "South East-East of England", "West Midlands-North East")){
      X$lnQ <- -X$lnQ
      X$regions <- strsplit(pair, "-")[[1]] %>% rev %>% paste(collapse="-")
      return(X)
      }else{
        return(X)
        }
    })

Diff_delta$time <- NA

for(i in 1:nrow(Diff_delta)){
  Diff_delta$time[i] <- Diff_delta$date[i]-
    min(subset(Diff_delta, regions==Diff_delta$regions[i])$date) %>% as.numeric
}

Diff_delta_summary <- Diff_delta %>% ddply(~time, function(X){
  if(nrow(X)>=10){
    return(c("mean"=mean(X$lnQ), quantile(X$lnQ, probs=c(0.1,0.5,0.9))))
  }
})

Delta2.week <- Delta %>% ddply(~week+phecname, function(X){ # Aggregated per week
  return(data.frame("N"=nrow(X),
                    "NDelta"=sum(X$Delta)))
})

Delta2.week$q <- Delta2.week$NDelta/Delta2.week$N

Diff_delta.week <- Delta2.week %>% subset(q>=0.1 & q<=0.99) %>%
  ddply(~week, function(X){
    if(nrow(X)>1){
      tab <- data.frame()
      X <- X[match(regions_delta, X$phecname),] %>% drop_na
      for(i in 1:(nrow(X)-1)){
        for(j in (i+1):nrow(X)){
          tab <- tab %>% rbind(data.frame("lnQ"=logit(X$q[i])-logit(X$q[j]),
                                          "regions"=paste(X$phecname[i], X$phecname[j], sep="-")))
        }
      }
      return(tab)
    }
  }) %>% ddply(~regions, function(X){ # changing the direction of some pairs
    
    pair <- unique(X$regions)
    
    if(pair %in% c("North West-London", "North West-South West", "London-South West",
                   "South East-East of England", "West Midlands-North East")){
      X$lnQ <- -X$lnQ
      X$regions <- strsplit(pair, "-")[[1]] %>% rev %>% paste(collapse="-")
      return(X)
    }else{
      return(X)
    }
  })

colnames(Diff_delta.week)[1] <- "date"

abbrev <- c("East Midlands"="EM", "East of England"="EE", "London"="LDN",
            "North East"="NE", "North West"="NW", "South East"="SE",
            "South West"="SW", "West Midlands"="WM", "Yorkshire and Humber"="YH")

Diff_delta.week$time <- NA

for(i in 1:nrow(Diff_delta.week)){
  Diff_delta.week$time[i] <- Diff_delta.week$date[i]-
    min(subset(Diff_delta.week, regions==Diff_delta.week$regions[i])$date) %>% as.numeric
}

Diff_delta_summary.week <- Diff_delta.week %>% ddply(~time, function(X){
  if(nrow(X)>=10){
    return(c("mean"=mean(X$lnQ), quantile(X$lnQ, probs=c(0.1,0.5,0.9))))
  }
})

levels_alpha <- paste(rep(regions_alpha, each=9), regions_alpha, sep="-")

Diff_alpha$regions <- factor(Diff_alpha$regions, levels=levels_alpha,
                             labels=levels_alpha %>% stringr::str_replace_all(abbrev))

levels_delta <- paste(rep(regions_delta, each=9), regions_delta, sep="-")

i1 <- which(levels_delta=="South West-North West")
i2 <- which(levels_delta=="South West-London")
levels_delta[c(i1,i2)] <- levels_delta[c(i2,i1)]

Diff_delta$regions <- factor(Diff_delta$regions, levels=levels_delta,
                             labels=levels_delta %>% stringr::str_replace_all(abbrev))

Diff_delta.week$regions <- factor(Diff_delta.week$regions, levels=levels_delta,
                                  labels=levels_delta %>% stringr::str_replace_all(abbrev))

i1 <- which(levels_delta=="South West-North West")
i2 <- which(levels_delta=="South West-London")
levels_delta[c(i1,i2)] <- levels_delta[c(i2,i1)]

Diff_alpha$regions_colored <- factor(html_style(Diff_alpha$regions),
                                     levels=html_style(levels(Diff_alpha$regions)))

Diff_delta$regions_colored <- factor(html_style(Diff_delta$regions),
                                     levels=html_style(levels(Diff_delta$regions)))

Diff_delta.week$regions_colored <- factor(html_style(Diff_delta.week$regions),
                                          levels=html_style(levels(Diff_delta.week$regions)))

lnQ_limits <- c(max(abs(c(Diff_delta$lnQ, Diff_delta.week$lnQ))))*c(-1,1)


legend1 <- get_legend(
  ggplot(data.frame(x=1:2, y=1:2), aes(x=x, y=y, col='trajectory')) +
    geom_line(linewidth=0.5, alpha=0.5) +
    scale_color_manual(values='black') +
    theme_minimal() + theme(legend.title=element_blank(), legend.text=element_text(size=8))
  )

legend2 <- get_legend(
  ggplot(data.frame(x=1:2, y=1:2, ymax=2:3, id='median (80% interval)'), aes(x=x)) +
    geom_line(aes(y=y, col=id), cex=0.5) +
    geom_ribbon(aes(ymin=y, ymax=ymax, fill=id), col=NA, alpha=0.2) +
    scale_fill_manual(values='blue') + scale_color_manual(values='blue') +
    theme(legend.title=element_blank(), legend.text=element_text(size=8))
)

# Fig. 5. The decay of spatial differentiation over time for the Alpha and Delta variants in England

plot_grid(
  
  ggdraw()+draw_label('Alpha', size=18, hjust=0.5, fontface='italic'),
  ggdraw()+draw_label('Delta', size=18, hjust=0.5, fontface='italic'),
  
  ggplot(Diff_alpha, aes(x=date, y=regions_colored, fill=lnQ)) +
    geom_tile(col='black', size=0.2) +
    scale_x_date(expand=c(0,0), date_labels= "%Y %b.", breaks=c(ymd("2020-12-1", "2021-1-1"))) +
    scale_y_discrete(position = "right") +
    scale_fill_distiller(palette = "RdBu", direction = -1, limits=lnQ_limits) +
    labs(fill="ln(\U1D4AC)") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=10), axis.title.x=element_blank(),
          axis.text.y.right=element_markdown(size=5), axis.title.y=element_blank(),
          legend.title=element_text(size=15), panel.grid = element_blank(),
          legend.position=c(0.1,0.7)),
  
  ggplot(Diff_delta.week, aes(x=date, y=regions_colored, fill=lnQ)) +
    geom_tile(col='black', size=0.2) +
    scale_x_date(expand=c(0,0), date_labels="%Y %b.", breaks=c(ymd("2021-4-1", "2021-5-1", "2021-6-1"))) +
    scale_fill_distiller(palette="RdBu", direction=-1, limits=lnQ_limits) +
    scale_y_discrete(position="right") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=10), axis.title.x=element_blank(),
          axis.text.y.right=element_markdown(size=5), axis.title.y=element_blank(),
          panel.grid=element_blank(), legend.position='none'),
  
  ggplot() +
    geom_ribbon(data=Diff_alpha_summary, aes(x=time, ymin=`10%`, ymax=`90%`, col='trajectory'),
                col=NA, fill='blue', alpha=0.2) +
    geom_line(data=Diff_alpha, aes(x=time, y=lnQ, group=regions), linewidth=0.2, alpha=0.4) +
    geom_line(data=Diff_alpha_summary, aes(x=time, y=`50%`), cex=0.5, col='blue') +
    geom_hline(yintercept = 0, col='red', linewidth=0.4) +
    labs(x='Elapsed time (days)', y="ln(\U1D4AC), spatial log-differentiation\n") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), limits=c(-3.5,3.5), breaks=seq(-3,3,1)) +
    theme_bw() +
    theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12)) +
    annotation_custom(grob=legend1, x=3, xmax=17, ymin=-1.75, ymax=-1.75) +
    annotation_custom(grob=legend2, x=3, xmax=30, ymin=-2.5, ymax=-2.5),
    
  ggplot() +
    geom_ribbon(data=Diff_delta_summary.week, aes(x=time, ymin=`10%`, ymax=`90%`),
                col=NA, fill='blue', alpha=0.2) +
    geom_line(data=Diff_delta.week, aes(x=time, y=lnQ, group=regions), linewidth=0.2, alpha=0.4) +
    geom_line(data=Diff_delta_summary, aes(x=time, y=`50%`), cex=0.5, col='blue') +
    geom_hline(yintercept = 0, col='red', linewidth=0.4) +
    labs(x='Elapsed time (days)', y="ln(\U1D4AC), spatial log-differentiation\n") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), limits=c(-3.5,3.5), breaks=seq(-3,3,1)) +
    theme_bw() +
    theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12)),
  
  ncol=2, align='v', labels=c('', '', LETTERS[1:4]),
  rel_heights=c(0.075,0.45,0.55)) # 9 x 7

# ggsave("Fig5.png", width=9, height=7, dpi=700)

# Fig. S1. Delays between emergence of the SARS-CoV-2 Alpha and Delta variants across regions of England

Alpha_first_date <- Alpha %>% subset(q>=0.1) %>%
  ddply(~Region, function(X){return(c('date'=min(X$week)))})

Delta_first_date <- Delta2 %>% subset(q>=0.1) %>%
  ddply(~phecname, function(X){return(c('date'=min(X$date)))})

date_diff <- rbind(
  data.frame('Variant'='Alpha',
             'date_diff'=Alpha_first_date$date %>% combn(2, function(X){return(diff(X))})),
  data.frame('Variant'='Delta',
             'date_diff'=Delta_first_date$date %>% combn(2, function(X){return(diff(X))})))

ggplot(date_diff, aes(x=abs(date_diff), fill=Variant)) +
  labs(x="\nInterval between emergence dates (days)", y="Count\n") +
  geom_density(aes(y=after_stat(count)*5), bounds=c(0, Inf), color=NA, alpha=0.3) +
  geom_bar(position = 'identity', col='black', alpha=0.75) +
  scale_x_continuous(breaks=seq(0,60,10), limits=c(-0.5,62)) +
  scale_y_continuous(breaks=seq(0,10,2)) +
  scale_fill_manual(values=c('#d95f02', "#1b9e77")) +
  theme_minimal() +
  theme(legend.position=c(0.87,0.8)) # 6 x 4

# Fig. S8. Approximations of the early dynamics of the spatial differentiation

## A-B (Theory)

bAw <- bBw <- bAm <- bBm <- 0.15
gAw <- gBw <- 0.1
gAm <- gBm <- 0.05

mA <- mB <- 1E-2

parms <- c("bAw"=bAw, "bBw"=bBw, "bAm"=bAm, "bBm"=bBm,
           "gAw"=gAw, "gBw"=gBw, "gAm"=gAm, "gBm"=gAm,
           "mA"=mA, "mB"=mB)

nA <- nB <- 100

I0w <- 1E-3
qA0 <- 1/100

times <- seq(0,80,0.1)

lnQ0 <- 0:3

simul1 <- lnQ0 %>% lapply(function(x){
  
  qB0 <- inv.logit(logit(qA0)-x)
  
  I0Am <- qA0/(1-qA0)*I0w
  I0Bm <- qB0/(1-qB0)*I0w
  
  init <- c("SA"=nA-I0w-I0Am, "IAw"=I0w, "IAm"=I0Am, "RA"=0,
            "SB"=nB-I0w-I0Bm, "IBw"=I0w, "IBm"=I0Bm, "RB"=0)
  
  return(simulate(init, times, parms, approx.WM=FALSE) %>%
           cbind("lnQ0"=x, "mu"=paste("Migration prob. =", parms[["mA"]])))
}) %>% bind_rows

parms[["mA"]] <- parms[["mA"]] <- 0.1

simul2 <- lnQ0 %>% lapply(function(x){
  
  qB0 <- inv.logit(logit(qA0)-x)
  
  I0Am <- qA0/(1-qA0)*I0w
  I0Bm <- qB0/(1-qB0)*I0w
  
  init <- c("SA"=nA-I0w-I0Am, "IAw"=I0w, "IAm"=I0Am, "RA"=0,
            "SB"=nB-I0w-I0Bm, "IBw"=I0w, "IBm"=I0Bm, "RB"=0)
  
  return(simulate(init, times, parms, approx.WM=FALSE) %>%
           cbind("lnQ0"=x, "mu"=paste("Migration prob. =", parms[["mA"]])))
}) %>% bind_rows

reg1 <- simul1 %>% ddply(~lnQ0, function(X){
  return(data.frame("slope"=(X$lnQ[X$time==3]-X$lnQ[X$time==0])/3, "mu"=X$mu[1]))
})

reg2 <- simul2 %>% ddply(~lnQ0, function(X){
  return(data.frame("slope"=(X$lnQ[X$time==3]-X$lnQ[X$time==0])/3, "mu"=X$mu[1]))
})

reg <- rbind(reg1, reg2)
simul <- rbind(simul1, simul2)

levels_mu <- reg$mu %>% unique %>% rev

reg$mu <- factor(reg$mu, levels=levels_mu)
simul$mu <- factor(simul$mu, levels=levels_mu)

# C (Data)

Diff_alpha_reg <- Diff_alpha %>% ddply(~regions, function(X){
  reg <- lm(lnQ~time, data=subset(X, time<=35))
  return(c(reg$coefficients, confint(reg)) %>%
           set_names(c("intercept","slope","intercept_CI_lower","slope_CI_lower",
                       "intercept_CI_upper","slope_CI_upper")))
})

Diff_alpha_reg$Variant <- "Alpha"

Diff_delta_reg <- Diff_delta %>% ddply(~regions, function(X){
  reg <- lm(lnQ~time, data=subset(X, time<=35))
  return(c(reg$coefficients, confint(reg)) %>%
           set_names(c("intercept","slope","intercept_CI_lower","slope_CI_lower",
                       "intercept_CI_upper","slope_CI_upper")))
})

Diff_delta_reg$Variant <- "Delta"

reg_reg_alpha <- lm(intercept~slope, Diff_alpha_reg)
reg_reg_delta <- lm(intercept~slope, Diff_delta_reg)

## Final plot

purples <- c("pink","#BF87B3","#8255a1","#1F007F")

plot_grid(
  ggdraw()+draw_label("Theory", size=16, fontface="italic"),
  ggdraw()+draw_label("Data", size=16, fontface="italic"),
  plot_grid(
    plot_grid(
      NULL,
       simul %>% ggplot(aes(x=time, y=lnQ, col=as.factor(lnQ0))) +
        facet_wrap(~mu) +
        geom_line(cex=1.5, alpha=0.4) +
        geom_abline(data=reg, aes(slope=slope, intercept=lnQ0, col=as.factor(lnQ0))) +
        labs(x="Time", y="Spatial log-differentiation, ln(\U1D4AC)\n",
             col="Initial spatial\nlog-differentiation") +
        scale_x_continuous(expand=c(0,0)) +
        scale_color_manual(values=purples) +
          theme_bw() +
        theme(axis.title.x=element_text(size=9), axis.text.x=element_text(size=8),
              axis.title.y=element_text(size=9), axis.text.y=element_text(size=8),
              legend.title=element_text(size=8), legend.text=element_text(size=8),
              legend.position='none'),
      NULL, ncol=3, rel_widths=c(0.03,0.9,0.07)),
  
    plot_grid(
      NULL,
      ggplot(reg, aes(x=slope, y=lnQ0, group=mu, pch=mu)) +
        geom_hline(yintercept=0, linewidth=0.5, col='gray60') +
        geom_vline(xintercept=0, linewidth=0.5, col='gray60') +
        geom_line(linewidth=0.2, lty='dashed') +
        geom_point(aes(col=as.factor(lnQ0)), cex=3) +
        geom_segment(x=-0.075, xend=-0.175, y=1, yend=1,
                     arrow=arrow(length=unit(0.25,"cm"), type="closed"),
                     linewidth=0.25) +
        annotate(geom='text', label='Migration', x=-0.12, y=0.8, size=3) +
        labs(x="\nSlope (early decay rate)", y="Intercept (initial value)\n", pch="Migration prob.") +
        scale_color_manual(values=purples) +
        scale_shape_manual(values=c(15,19), labels=c(0.1,0.01)) +
        theme_minimal() +
        theme(axis.title.x=element_text(size=9), axis.text.x=element_text(size=8),
              axis.title.y=element_text(size=9), axis.text.y=element_text(size=8),
              legend.title=element_text(size=8), legend.text=element_text(size=8),
              legend.position='top', legend.direction='horizontal',
              legend.box.just='left', legend.justification=c(0, 1)) +
        guides(col='none'),
      NULL, ncol=3, rel_widths=c(0.15,0.7,0.15)),
    
  ncol=1, rel_heights=c(0.45,0.55), labels=LETTERS[1:2]),
  
  rbind(Diff_alpha_reg, Diff_delta_reg) %>%
    ggplot(aes(x=slope, y=intercept, col=Variant)) +
    geom_hline(yintercept=0, linewidth=0.5, col='gray60') +
    geom_vline(xintercept=0, linewidth=0.5, col='gray60') +
    geom_smooth(method="lm", cex=0.75) +
    geom_segment(aes(x=slope_CI_lower, xend=slope_CI_upper, y=intercept, yend=intercept),
                 linewidth=0.2, alpha=0.75) +
    geom_segment(aes(x=slope, xend=slope, y=intercept_CI_lower, yend=intercept_CI_upper),
                 linewidth=0.2, alpha=0.75) +
    geom_point(cex=3, pch=18, alpha=0.75) +
    labs(x="\nSlope (estimated early decay rate, per day)",
         y="Intercept (estimated initial spatial log-differentition)\n")+
    scale_y_continuous(breaks=-2:4) +
    scale_color_manual(values=c('#d95f02', "#1b9e77")) +
    theme_minimal() +
    theme(axis.title.x=element_text(size=11), axis.text.x=element_text(size=10),
          axis.title.y=element_text(size=11), axis.text.y=element_text(size=10),
          legend.title=element_text(size=15), legend.text=element_text(size=12),
          legend.position=c(0.85,0.85)),
  
  ncol=2, rel_heights=c(0.1,0.9), labels=c(rep("",3), "C")) # 10 x 6

# ggsave("FigS8.png", width=10, heigh=6, dpi=600)

# Fig. S9. Change in visitors in the United Kingdom (Google mobility reports)

mobility <- read.csv("Data/changes-visitors-covid.csv", header=TRUE) %>% subset(Code=="GBR")

mobility <- mobility[,c("Day", "Retail", "Grocery", "Transit", "Workplaces") %>%
                       sapply(function(x){grep(x,colnames(mobility))})]

mobility$Day <- ymd(mobility$Day)

mobility$trend4 <- mobility[,-1] %>% apply(1,mean)

t_Alpha_min <- min(subset(Alpha, q>=0.1 & q<=0.99)$week)
t_Alpha_max <- max(subset(Alpha, q>=0.1 & q<=0.99)$week)+6

t_Delta_min <- min(subset(Delta2, q>=0.1 & q<=0.99)$date)
t_Delta_max <- max(subset(Delta2, q>=0.1 & q<=0.99)$date)

mobility %>% subset(Day<ymd("2021-09-01")) %>%
  gather(id, value, -Day) %>%
  ggplot(aes(x=Day, y=value)) +
  geom_rect(xmin= t_Alpha_min, xmax=t_Alpha_max, ymin=-Inf, ymax=20, fill="gray80") +
  geom_rect(xmin=t_Delta_min, xmax=t_Delta_max, ymin=-Inf, ymax=20, fill="gray80") +
  geom_hline(yintercept=0) +
  geom_line(aes(col=id)) +
  labs(y="Change in visitors\n") +
  scale_x_date(date_labels = "%Y %b.",
               breaks=c(ymd("2020-07-01"), ymd("2021-01-01"), ymd("2021-07-01"))) +
  scale_y_continuous(limits=c(-80,20), expand=c(0,0),breaks=seq(-80,20,20),
                     labels=label_percent(scale=1)) +
  scale_color_manual(values=c("forestgreen", "orange", "red2", "cornflowerblue", "black"),
                     labels=c("Grocery & pharmacy", "Retail & recreation",
                              "Transit stations", "Workplaces", "Mean")) +
  annotation_custom(
    grob=textGrob("Alpha sweep", gp=gpar(col="black", fontsize=12, fontface="italic")),
    xmin=t_Alpha_min + (t_Alpha_max-t_Alpha_min)/2,
    xmax=t_Alpha_min + (t_Alpha_max-t_Alpha_min)/2,
    ymin=25, ymax=22) +
  annotation_custom(
    grob=textGrob("Delta sweep", gp=gpar(col="black", fontsize=12, fontface="italic")),
    xmin=t_Delta_min + (t_Delta_max-t_Delta_min)/2, 
    xmax=t_Delta_min + (t_Delta_max-t_Delta_min)/2, 
    ymin=25, ymax=22) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=11),
        axis.title.y=element_text(size=12), axis.text.y=element_text(size=11),
        legend.title=element_blank(), legend.text=element_text(size=10), 
        plot.margin=margin(1,0.1,0.1,0.1,"cm")) # 8 x 5