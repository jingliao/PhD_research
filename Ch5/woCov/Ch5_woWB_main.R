###******** required packages/functions
library(MCMCpack)  # dirichlet distribution
library(tidyverse)
library(islandR)
library(forcats)
library(ggplot2)
library(GGally) # for ggpairs

source("new_functions_ch5woCov.R")

###******** load the data
#' the .R file must be located outside data-raw file
alltypes <- read.csv("human_types_age.csv")

###******** divide the data into human, source(P,R,O) and water
#' frequency table of ST vs. sources
freqtable <- as.data.frame.matrix(table(alltypes$ST,
                                          alltypes$Source))

SeqTypes <- rownames(freqtable)

#' human data, h_{i}
tableHum <- freqtable %>% select(Human)

#' source data, x_{ij}
tableSource <- freqtable %>% select(-Human, -Water)

#' water data
tableWater <- freqtable %>% select(Water)

###******** initial points of each parameter before and after delogit
#' initial settings
sigma <- qsigma <- 1

#' number of simulation for sampling distribution of genotypes
# simulation <- 1
simulation <- 100

##########********  WATER ATTRIBUTION  ********##########
## load the source file for island or dirichlet model
## island
source("Ch5_woWB_sourceI.R")

## Dirichlet
source("Ch5_woWB_sourceD.R")

###******** joint likelihoods of water and human via P(STi|sourcej)
sampDistSTs <- split(sampDistST,
                   rep(1:simulation, each=nrow(freqtable)))

postthetalist <- list()

Iteration <- 10000
Thinning <- 100
Burnin <- 3000

for (s in 1:simulation){

    sampDistS <- sampDistSTs[[s]][,-4]
    
    #' water
    par_pW <- as.numeric(as.matrix(sampDistS)%*%pW_j)

    #' human
    #' water part
    par_pH1 <- par_pW*pH_j[3]
    #' other sources part
    par_pH2 <- as.numeric(as.matrix(sampDistS)%*%pH_j[-3])
    
    par_pH <- par_pH1 + par_pH2

    #' joint loglikelihood
    logLH_pW <- sum(tableWater$Water*log(par_pW))
    logLH_pH <- sum(tableHum$Human*log(par_pH))
    jointlogLH <- logLH_pW + logLH_pH

###******** MCMC - Metropolis-Hasting Random walk
    totalMH <- Burnin+Iteration*Thinning

    accept <- reject <- 0

    postthetalist[[s]] <- mcmcalgm(qsigma=1, sigma=1,
                               thetaW=pw_source,
                               thetaH=ph_source,
                               log_LH=jointlogLH,
                               post_pij=sampDistS)

}

posttheta <- do.call(rbind, postthetalist)


#' columns for simulation and iteration
colsm <- data.frame(s=rep(1:simulation,
                          Iteration/simulation),
                    m=1:Iteration) %>%
    expand(s,m) ## for diagnostics section


##**** island 
load("newTheoryI(s100).RData")

head(posttheta)
dim(posttheta)

###******** posterior pW_j and pH_j
post_pw_source <- posttheta[,1:2]
post_ph_source <- posttheta[,3:5]

postWjI <- post_pw_source
colnames(postWjI) <- c("Poultry", "Other")

postHjI <- post_ph_source
colnames(postHjI) <- c("Poultry", "Other", "Water")

#' water attribution
post_pW_j <- delogit_mtrx(par_mtrx=post_pw_source)
colnames(post_pW_j) <- c("Poultry", "Other", "Ruminants")

smWI <- as.data.frame(post_pW_j) %>%
    mutate(s=colsm$s, m=colsm$m) %>% group_by(s, m) %>%
    gather(source, postW, Poultry:Ruminants) %>%
    mutate(model="Island")   ## for diagnostics section

postWI <- as.data.frame(post_pW_j) %>%
    gather(Source, pW, Poultry:Ruminants) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "Other")),
           model="Island") 

#' human attribution
post_pH_j <- delogit_mtrx(par_mtrx=post_ph_source)
colnames(post_pH_j) <- c("Poultry","Other","Water","Ruminants")

smHI <- as.data.frame(post_pH_j) %>%
    mutate(s=colsm$s, m=colsm$m) %>% group_by(s, m) %>%
    gather(source, postH, Poultry:Ruminants) %>%
    mutate(model="Island")   ## for diagnostics section

postHI <- as.data.frame(post_pH_j) %>%
    gather(Source, pH, Poultry:Ruminants) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other")),
           model="Island")

#' trace plots
## water
postpwI <- as.data.frame(post_pw_source) %>%
    rename("Poultry"=V1, "Other"=V2) %>%
    mutate(Simulation=as.factor(rep(seq(1, simulation),
                                    each=Iteration)),
           Iteration=rep(seq(1, Iteration), simulation)) %>%
    gather(Source, value, Poultry:Other) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Other")),
           Iteration=as.numeric(Iteration),
           model="Island")

## human
postphI <- as.data.frame(post_ph_source) %>%
    rename("Poultry"=V1, "Other"=V2, "Water"=V3) %>%
    mutate(Simulation=as.factor(rep(seq(1, simulation),
                                    each=Iteration)),
           Iteration=rep(seq(1, Iteration), simulation)) %>%
    gather(Source, value, Poultry:Water) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Water",
                                          "Other")),
           Iteration=as.numeric(Iteration),
           model="Island")

##**** Dirichlet
load("newTheoryD(s100).RData")

head(posttheta)
dim(posttheta)

###******** posterior pW_j and pH_j
post_pw_source <- posttheta[,1:2]
post_ph_source <- posttheta[,3:5]

postWjD <- post_pw_source
colnames(postWjD) <- c("Poultry", "Other")

postHjD <- post_ph_source
colnames(postHjD) <- c("Poultry", "Other", "Water")

#' water attribution
post_pW_j <- delogit_mtrx(par_mtrx=post_pw_source)
colnames(post_pW_j) <- c("Poultry", "Other", "Ruminants")

smWD <- as.data.frame(post_pW_j) %>%
    mutate(s=colsm$s, m=colsm$m) %>% group_by(s, m) %>%
    gather(source, postW, Poultry:Ruminants) %>%
    mutate(model="Dirichlet")   ## for diagnostics section

postWD <- as.data.frame(post_pW_j) %>%
    gather(Source, pW, Poultry:Ruminants) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "Other")),
           model="Dirichlet") 

#' human attribution
post_pH_j <- delogit_mtrx(par_mtrx=post_ph_source)
colnames(post_pH_j) <- c("Poultry","Other","Water","Ruminants")

smHD <- as.data.frame(post_pH_j) %>%
    mutate(s=colsm$s, m=colsm$m) %>% group_by(s, m) %>%
    gather(source, postH, Poultry:Ruminants) %>%
    mutate(model="Dirichlet")   ## for diagnostics section

postHD <- as.data.frame(post_pH_j) %>%
    gather(Source, pH, Poultry:Ruminants) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other")),
           model="Dirichlet")

#' trace plots
## water
postpwD <- as.data.frame(post_pw_source) %>%
    rename("Poultry"=V1, "Other"=V2) %>%
    mutate(Simulation=as.factor(rep(seq(1, simulation),
                                    each=Iteration)),
           Iteration=rep(seq(1, Iteration), simulation)) %>% 
    gather(Source, value, Poultry:Other) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Other")),
           Iteration=as.numeric(Iteration),
           model="Dirichlet") 
## human
postphD <- as.data.frame(post_ph_source) %>%
    rename("Poultry"=V1, "Other"=V2, "Water"=V3) %>%
    mutate(Simulation=as.factor(rep(seq(1, simulation),
                                   each=Iteration)),
           Iteration=rep(seq(1, Iteration), simulation)) %>%
    gather(Source, value, Poultry:Water) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Water",
                                          "Other")),
           Iteration=as.numeric(Iteration),
           model="Dirichlet") 

####################################################################
###  comparison plots for Water and Human attribution            ###
####################################################################

##*************** for result section
postW <- bind_rows(postWI, postWD) %>%
    mutate(pW=pW*100,
           model=factor(model, levels=c("Island", "Dirichlet")))

postH <- bind_rows(postHI, postHD) %>%
    mutate(model=factor(model, levels=c("Island", "Dirichlet")))

jpeg("fig_ch5_2_upper.jpg", width = 700, height =500, quality=100)

## comparison of water attribution between island and dirichlet model
ggplot(postW, aes(x=Source, y=pW)) +
    geom_violin() + stat_summary(fun.y=median, geom="point", size=1,
                                 color="red") +
    ylim(0,100) +
    ylab("percentage of water isolates (%)") + xlab("") +
    facet_wrap(~model) + theme_bw() +
    theme(text=element_text(size=12, family="serif"))

dev.off()

jpeg("fig_ch5_3_upper.jpg", width = 700, height =500, quality=100)
## comparison of human attribution between two models

ggplot(postH) +
    geom_boxplot(aes(x=Source, y=pH*100)) +
    ylim(0,100) +
    ylab("percentage of human cases (%)") + xlab("") +
    facet_wrap(~model) +
    theme_bw() +
    theme(text=element_text(size=12, family="serif"))

dev.off()

#' comparison of trace plot
jpeg("fig_ch5_5_upper.jpg", width = 700, height =500, quality=100)

## water
postpw <- bind_rows(postpwI, postpwD) %>%
    mutate(model=factor(model, levels=c("Island", "Dirichlet"))) %>%
    rename("simulation"=Simulation)

ggplot(postpw%>%filter(simulation==10|simulation==36|simulation==98),
       aes(x=Iteration, y=value, color=simulation)) +
    geom_line() +
    xlab("iteration") + ylab("posterior g") +
    facet_grid(model~Source) + 
    theme_bw() +
    theme(legend.position="bottom",
          legend.box.margin = margin(-10, 6, -8, 6),
          text=element_text(size=12, family="serif"))

dev.off()

jpeg("fig_ch5_6_upper.jpg", width = 700, height =500, quality=100)

postph <- bind_rows(postphI, postphD) %>%
    mutate(model=factor(model, levels=c("Island", "Dirichlet"))) %>%
    rename("simulation"=Simulation)

ggplot(postph%>%filter(simulation==10|simulation==36|simulation==98),
       aes(x=Iteration, y=value, color=simulation)) +
    geom_line() +
    xlab("iteration") + ylab("posterior f") +
    facet_grid(model~Source) + 
    theme_bw() +
    theme(legend.position="bottom",
          legend.box.margin = margin(-10, 6, -8, 6),
          text=element_text(size=12, family="serif"))

dev.off()

##*************** for diagnostics section
#' water attribution
## combine island and dirichlet models
smW_diag <- bind_rows(smWI, smWD)

group1 <- smW_diag %>% filter(s==1, m==1) %>% mutate(group=1)

group2 <- smW_diag %>% filter(s==1, m<=100) %>% mutate(group=2)

group3 <- smW_diag %>% filter(m<=100) %>% mutate(group=3)

postW_diag <- bind_rows(group1, group2, group3) %>%
    mutate(group=factor(group, levels=c("1", "2", "3"),
                        labels=c("s=1, m=1", "s=1, m=100",
                                 "s=100, m=100")),
           model=factor(model, levels=c("Island", "Dirichlet")),
           postW=postW*100,
           source=factor(source, levels=c("Poultry", "Ruminants",
                                          "Other")))

pdf(file="fig_ch5_7_upper.pdf", width=7, height=5)

ggplot(postW_diag, aes(x=source, y=postW)) + geom_violin() + 
    stat_summary(fun.y=median, geom="point", size=1, color="red") +
    ylab("percentage of water isolates (%)") + ylim(c(0,100)) +
    theme_bw() +
    theme(text=element_text(size=12, family="serif")) +
    facet_grid(model~group)

dev.off()

#' human attribution
## combine island and dirichlet models
smH_diag <- bind_rows(smHI, smHD)

group1 <- smH_diag %>% filter(s==1, m==1) %>% mutate(group=1)

group2 <- smH_diag %>% filter(s==1, m<=100) %>% mutate(group=2)

group3 <- smH_diag %>% filter(m<=100) %>% mutate(group=3)

postH_diag <- bind_rows(group1, group2, group3) %>%
    mutate(group=factor(group, levels=c("1", "2", "3"),
                        labels=c("s=1, m=1", "s=1, m=100",
                                 "s=100, m=100")),
           model=factor(model, levels=c("Island", "Dirichlet")),
           postH=postH*100,
           source=factor(source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other")))

pdf(file="fig_ch5_8_upper.pdf", width=7, height=5)

ggplot(postH_diag, aes(x=source, y=postH)) + geom_boxplot() + 
    ylab("percentage of human cases (%)") + ylim(c(0,100)) +
    theme_bw() +
    theme(text=element_text(size=12, family="serif"),
          axis.text.x=element_text(size=9)) +
    facet_grid(model~group)

dev.off()

