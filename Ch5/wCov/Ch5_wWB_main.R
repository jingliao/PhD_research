###******** required packages/functions
library(MCMCpack)  # dirichlet distribution
library(tidyverse)
library(islandR)
library(forcats)
library(GGally)

source("new_functionsCov.R")

###******** initial settings
#' deviation of normal prior
sigma <- qsigma <- 1

#' number of simulating sampling distribution
simulation <- 100

##########********  WATER ATTRIBUTION  ********##########
## load the source file for island or dirichlet model
## island
source("Ch5_wWB_sourceI.R")

## Dirichlet
source("Ch5_wWB_sourceD.R")

##########********  MCMC ALGORITHM  ********##########
sampDistSTs <- split(sampDistST,
                   rep(1:simulation, each=nrow(freqtable)))

## STs from sources that are also observed from humans 
y <- match(humanST, sourceST)

###******** MCMC - Metropolis-Hasting Random walk 
Iteration <- 10000
Thinning <- 100
Burnin <- 3000

totalMH <- Burnin+Iteration*Thinning

for (s in 1:simulation){

    sampDistS <- as.matrix(sampDistSTs[[s]])
    
    #' water
    par_pW <- sampDistS %*% pW_j
    
    #' human attribution combined with water attribution
    comb_pH <- cbind(par_pW, sampDistS)
    
    par_pH <- as.numeric(rowSums(comb_pH[y,]*ph_capF_R[,c(4,1,2,3,5)])) 
    #' joint loglikelihood
    logLH_pW <- sum(tableWater$Water*log(par_pW))
    logLH_pH <- sum(log(par_pH))
    jointlogLH <- logLH_pW + logLH_pH

###******** MCMC - Metropolis-Hasting Random walk
    accept <- reject <- 0

    posttheta  <- mcmcalgmCovWB(qsigma=1, sigma=1,
                              Wtheta=thetaW,
                              Htheta=thetaH_R,
                              log_LH=jointlogLH,
                              post_pij=sampDistS,
                              design_mtrx=designR_mtrx,
                              y=y)

    ## choose file path for posterior samples given the sampling
    ## distribution is based one which model

    ## files for island model
    save(posttheta,
         file=paste("postthetaR_RI_WB_",s,".RData", sep=""))
    
    ## files for Dirichlet model
    save(posttheta,
         file=paste("postthetaR_RD_WB_",s,".RData", sep=""))

}


