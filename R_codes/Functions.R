# The interplay between migration and selection on the dynamics of pathogen variants
# R functions

SIR_metapop <- function(t, y, parms){ # ODEs for the full model in a 2-patch metapopulation
  
  # t, current time
  # y, current state variables
  # parms, vector of parameters
  
  # State variables
  
  ## Pop. A
  SA <- y[["SA"]] # Susceptible
  IAw <- y[["IAw"]] # Infected (wildtype)
  IAm <- y[["IAm"]] # Infected (variant)
  RA <- y[["RA"]] # Recovered
  IA <- IAw+IAm # Total infected density
  nA <- SA+IA+RA # Total population density
  
  ## Pop. B
  SB <- y[["SB"]] # Susceptible
  IBw <- y[["IBw"]] # Infected (wildtype)
  IBm <- y[["IBm"]] # Infected (variant)
  RB <- y[["RB"]] # Recovered
  IB <- IBw+IBm # Total infected density
  nB <- SB+IB+RB # Total population density
  
  # Parameters
  
  # Transmission rates
  bAw <- parms[["bAw"]]
  bBw <- parms[["bBw"]]
  bAm <- parms[["bAm"]]
  bBm <- parms[["bBm"]]
  
  # Recovery rates
  gAw <- parms[["gAw"]]
  gBw <- parms[["gBw"]]
  gAm <- parms[["gAm"]]
  gBm <- parms[["gBm"]]
  
  # Migration probabilities
  mA <- parms[["mA"]]
  mB <- parms[["mB"]]
  z <- ifelse("z" %in% names(parms),
              yes=1-parms[["z"]], no=1) # multiplier for the mobility of infected hosts
  
  # Waning immunity
  w <- ifelse("w" %in% names(parms), yes=parms[["w"]], no=0)
  
  # Forces of infection
  hAw <- bAw*((1-z*mA)*IAw+z*mB*IBw)/(nA-mA*(SA+z*IA+RA)+mB*(SB+z*IB+RB))
  hAm <- bAm*((1-z*mA)*IAm+z*mB*IBm)/(nA-mA*(SA+z*IA+RA)+mB*(SB+z*IB+RB))
  hA <- hAw+hAm
  hBw <- bBw*((1-z*mB)*IBw+z*mA*IAw)/(nB-mB*(SB+z*IB+RB)+mA*(SA+z*IA+RA))
  hBm <- bBm*((1-z*mB)*IBm+z*mA*IAm)/(nB-mB*(SB+z*IB+RB)+mA*(SA+z*IA+RA))
  hB <- hBw+hBm
  
  # Differentiation with respect to time
  
  ## Pop. A
  dSA <- -((1-mA)*hA+mA*hB)*SA + w*RA
  dIAw <- ((1-mA)*hAw+mA*hBw)*SA - gAw*IAw
  dIAm <- ((1-mA)*hAm+mA*hBm)*SA - gAm*IAm
  dRA <- gAw*IAw + gAm*IAm - w*RA
  
  ## Pop. B
  dSB <- -((1-mB)*hB+mB*hA)*SB + w*RB
  dIBw <- ((1-mB)*hBw+mB*hAw)*SB - gBw*IBw
  dIBm <- ((1-mB)*hBm+mB*hAm)*SB - gBm*IBm
  dRB <- gBw*IBw + gBm*IBm - w*RB
  
  return(list(c(dSA, dIAw, dIAm, dRA, dSB, dIBw, dIBm, dRB)))
}

# ODEs for the model in a 2-patch metapopulation assuming weak migration (WM) approx
SIR_metapop.WM <- function(t, y, parms){ 
  
  # t, current time
  # y, current state variables
  # parms, vector of parameters
  
  # State variables
  
  ## Pop. A
  SA <- y[["SA"]] # Susceptible
  IAw <- y[["IAw"]] # Infected (wildtype)
  IAm <- y[["IAm"]] # Infected (variant)
  RA <- y[["RA"]] # Recovered
  IA <- IAw+IAm # Total infected density
  nA <- SA+IA+RA # Total population density
  
  ## Pop. B
  SB <- y[["SB"]] # Susceptible
  IBw <- y[["IBw"]] # Infected (wildtype)
  IBm <- y[["IBm"]] # Infected (variant)
  RB <- y[["RB"]] # Recovered
  IB <- IBw+IBm # Total infected density
  nB <- SB+IB+RB # Total population density
  
  # Parameters
  
  # Transmission rates
  bAw <- parms[["bAw"]]
  bBw <- parms[["bBw"]]
  bAm <- parms[["bAm"]]
  bBm <- parms[["bBm"]]
  
  # Recovery rates
  gAw <- parms[["gAw"]]
  gBw <- parms[["gBw"]]
  gAm <- parms[["gAm"]]
  gBm <- parms[["gBm"]]
  
  # Migration probabilities
  mA <- parms[["mA"]]
  mB <- parms[["mB"]]
  z <- ifelse("z" %in% names(parms),
              yes=1-parms[["z"]], no=1) # multiplier for the mobility of infected hosts
  
  # Waning immunity
  w <- ifelse("w" %in% names(parms), yes=parms[["w"]], no=0)
  
  # New infections
  
  infection_Aw <- ((1-mA*(1-(1-z)*(1-IA/nA))-mB*(nB/nA-(1-z)*IB/nA))*bAw*IAw/nA+z*mB*bAw*IBw/nA+mA*bBw*IBw/nB)*SA
  infection_Am <- ((1-mA*(1-(1-z)*(1-IA/nA))-mB*(nB/nA-(1-z)*IB/nA))*bAm*IAm/nA+z*mB*bAm*IBm/nA+mA*bBm*IBm/nB)*SA
  
  infection_Bw <- ((1-mB*(1-(1-z)*(1-IB/nB))-mA*(nA/nB-(1-z)*IA/nB))*bBw*IBw/nB+z*mA*bBw*IAw/nB+mB*bAw*IAw/nA)*SB
  infection_Bm <- ((1-mB*(1-(1-z)*(1-IB/nB))-mA*(nA/nB-(1-z)*IA/nB))*bBm*IBm/nB+z*mA*bBm*IAm/nB+mB*bAm*IAm/nA)*SB
  
  # Differentiation with respect to time
  
  ## Pop. A
  dSA <- -infection_Aw-infection_Am + w*RA
  dIAw <- infection_Aw - gAw*IAw
  dIAm <- infection_Am - gAm*IAm
  dRA <- gAw*IAw + gAm*IAm - w*RA
  
  ## Pop. B
  dSB <- -infection_Bw-infection_Bm + w*RB
  dIBw <- infection_Bw - gBw*IBw
  dIBm <- infection_Bm - gBm*IBm
  dRB <- gBw*IBw + gBm*IBm - w*RB
  
  return(list(c(dSA, dIAw, dIAm, dRA, dSB, dIBw, dIBm, dRB)))
}

simulate <- function(init, times, parms, exact=TRUE, approx.WM=TRUE){ # Simulate models (see ODE systems above)
  
  # init = initial conditions for each state variable
  # parms = parameter values
  # times = vector of requested time points
  # exact = full model
  # approx.WM = weak migration approximation
  
  if(exact){
    simul <- lsoda(y=init, times=times, func=SIR_metapop, parms=parms) %>% as.data.frame
    
    IA <- simul$IAw+simul$IAm
    IB <- simul$IBw+simul$IBm
    
    simul$IA <- IA
    simul$IB <- IB
    simul$f <- IA/(IA+IB)
    simul$qA <- simul$IAm/IA
    simul$qB <- simul$IBm/IB
    simul$q <- (simul$IAm+simul$IBm)/(IA+IB)
    
    simul$logit_qA <- logit(simul$qA)
    simul$logit_qB <- logit(simul$qB)
    simul$logit_q <- logit(simul$q)
    
    simul$lnQ <- simul$logit_qA-simul$logit_qB
  }
  if(approx.WM){
    simul.WM <- lsoda(y=init, times=times, func=SIR_metapop.WM, parms=parms) %>% as.data.frame
    
    IA.WM <- simul.WM$IAw+simul.WM$IAm
    IB.WM <- simul.WM$IBw+simul.WM$IBm
    
    simul.WM$IA <- IA.WM
    simul.WM$IB <- IB.WM
    simul.WM$f <- IA.WM/(IA.WM+IB.WM)
    simul.WM$qA <- simul.WM$IAm/IA.WM
    simul.WM$qB <- simul.WM$IBm/IB.WM
    simul.WM$q <- (simul.WM$IAm+simul.WM$IBm)/(IA.WM+IB.WM)
    
    simul.WM$logit_qA <- logit(simul.WM$qA)
    simul.WM$logit_qB <- logit(simul.WM$qB)
    simul.WM$logit_q <- logit(simul.WM$q)
    
    simul.WM$lnQ <- simul.WM$logit_qA-simul.WM$logit_qB
  }
  if(exact & approx.WM){
    return(simul %>% data.frame("Approx_WM"=FALSE) %>% rbind(simul.WM %>% data.frame("Approx_WM"=TRUE)))
  }else if(exact){
    return(simul %>% data.frame("Approx_WM"=FALSE))
  }else if(approx.WM){
    return(simul.WM %>% data.frame("Approx_WM"=TRUE))
  }else{
    return(NULL)
  }
}

# Simulate models (see ODE systems above) with 2 phases (see below)
simulate2 <- function(init, times, parms, tM, qB0){ 
  
  # Same as simulate but simulate 2 phases :
  # (i) the variant is only present in pop. A and the two pop. are isolated (i.e., migration prob. mA=mB=0);
  # (ii) the variant is introduced in pop. B at t=tM and at frequency qB=qB0 (and pop. can be coupled through migration)
  
  if(!tM %in% times){
    times <- times %>% c(tM) %>% sort
  }
  
  times1 <- times[times<=tM]
  times2 <- times[times>=tM]
  
  parms2 <- parms
  parms[["mA"]] <- parms[["mB"]] <- 0
  
  simul <- simulate(init=init, times=times1, parms=parms)
  
  init2 <- simul[simul$time==tM & simul$Approx_WM==FALSE,2:9] %>% unlist
  init2[["IBm"]] <- qB0/(1-qB0)*init2[["IBw"]]
  
  init2.WM <- simul[simul$time==tM & simul$Approx_WM==TRUE,2:9] %>% unlist
  init2.WM[["IBm"]] <- qB0/(1-qB0)*init2.WM[["IBw"]]
  
  return(simul %>%
           rbind(simulate(init=init2, times=times2, parms=parms2, approx.WM=FALSE)) %>%
           rbind(simulate(init=init2, times=times2, parms=parms2, exact=FALSE)))
}

SIR_metapop.n <- function(t, y, parms){ # ODEs for the full model with n >= 2 pop.
  
  # t, current time
  # y, current state variables
  # parms, vector of parameters
  
  # State variables
  
  S <- y[grep("S", names(y))] # Susceptible
  Iw <- y[grep("Iw", names(y))] # Infected (wildtype)
  Im <- y[grep("Im", names(y))] # Infected (variant)
  R <- y[grep("R", names(y))] # Recovered
  N <- S+Iw+Im+R # Total population density
  
  # Parameters
  
  # Transmission rates
  bw <- parms[["bw"]]
  bm <- parms[["bm"]]
  
  # Recovery rates
  gw <- parms[["gw"]]
  gm <- parms[["gm"]]
  
  # Migration probabilities
  M <- parms[["M"]]
  
  # Waning immunity
  w <- ifelse("w" %in% names(parms), yes=parms[["w"]], no=0)
  
  # Forces of infection
  hw <- bw*M%*%Iw/(M%*%N)
  hm <- bm*M%*%Im/(M%*%N)
  h <- hw+hm
  
  # Differentiation with respect to time
  dS <- -t(M)%*%h*S + w*R
  dIw <- t(M)%*%hw*S - gw*Iw
  dIm <- t(M)%*%hm*S - gm*Im
  dR <- gw*Iw + gm*Im - w*R
  
  return(list(c(dS, dIw, dIm, dR)))
}

simulate.n <- function(init, times, parms){ # Simulate the ODE system for the full model with n >= 2 pop.
  
  n <- nrow(parms[["M"]])

  simul <- lsoda(y=init, times=times, func=SIR_metapop.n, parms=parms) %>% as.data.frame
  
  simul <- simul %>% cbind(simul[,grep("Iw", colnames(simul))]+simul[,grep("Im", colnames(simul))])
  colnames(simul)[(ncol(simul)-n+1):ncol(simul)] <- paste0("I_tot",1:n)
  
  simul <- simul %>% cbind(simul[,grep("Im", colnames(simul))]/simul[,grep("I_tot", colnames(simul))])
  colnames(simul)[(ncol(simul)-n+1):ncol(simul)] <- paste0("q",1:n)
  
  simul <- simul[,grep("q", colnames(simul))] %>% lapply(logit) %>% as.data.frame %>% cbind(simul,.)
  colnames(simul)[(ncol(simul)-n+1):ncol(simul)] <- paste0("logit_q",c(1:n))
  
  return(simul)
}

format_scientific <- function(x){
  
  x0 <- which(x==0)
  x1 <- which(x==1)
  x <- scales::scientific_format()(x) %>% gsub("e", " %*% 10^", .) %>% gsub("\\+", "", .)
  x[x0] <- "0"
  x[x1] <- "1"
  
  return(parse(text=x))
}

html_style <- function(x){
  return(paste0("<span style='color:",
                colors[sub("-.*", "", x) %>%
                         match(regions %>% str_replace_all(abbrev))],
                "'>", x, "</span>"))
}