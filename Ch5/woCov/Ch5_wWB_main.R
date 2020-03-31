###******** required packages/functions
library(MCMCpack)  # dirichlet distribution
library(tidyverse)
library(islandR)
library(forcats)
library(GGally) ## for ggpairs

source("new_functions_ch5woCov.R")

##########********  WATER ATTRIBUTION  ********##########

###******** load the data
#' the .R file must be located outside data-raw file
alltypes <- read.csv("human_types_wbirds.csv")

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

#' initial settings
sigma <- qsigma <- 1

#' number of simulation for sampling distribution of genotypes
simulation <- 100

##########********  WATER ATTRIBUTION  ********##########
## load the source file for island or dirichlet model
## island
source("Ch5_wWB_sourceI.R")

## Dirichlet
source("Ch5_wWB_sourceD.R")


###******** joint likelihoods of water and human via P(STi|sourcej)
sampDistSTs <- split(sampDistST,
                   rep(1:simulation, each=nrow(freqtable)))

postthetalist <- list()

Iteration <- 10000
Thinning <- 100
Burnin <- 3000

#' MCMC - Metropolis-Hasting Random walk 

for (s in 1:simulation){

    sampDistS <- sampDistSTs[[s]][-5]

    #' water
    par_pW <- as.numeric(as.matrix(sampDistS)%*%pW_j)
    
    #' human
    #' water part
    par_pH1 <- par_pW*pH_j[4]
    #' other sources part
    par_pH2 <- as.numeric(as.matrix(sampDistS)%*%pH_j[-4])

    par_pH <- par_pH1 + par_pH2

    #' joint loglikelihood
    logLH_pW <- sum(tableWater$Water*log(par_pW))
    logLH_pH <- sum(tableHum$Human*log(par_pH))
    jointlogLH <- logLH_pW + logLH_pH

###******** MCMC - Metropolis-Hasting Random walk 
  
    totalMH <- Burnin+Iteration*Thinning

    accept <- reject <- 0

    postthetalist[[s]] <- mcmcalgmWB(qsigma=1, sigma=1,
                                     thetaW=pw_source,
                                     thetaH=ph_source,
                                     log_LH=jointlogLH,
                                     post_pij=sampDistS)

}

posttheta <- do.call(rbind, postthetalist)

############################ save/load samples, choose island or Dirichlet model
#' columns for simulation and iteration
colsm <- data.frame(s=rep(1:simulation,
                          Iteration/simulation),
                    m=1:Iteration) %>%
    expand(s,m) ## for diagnostics section

##**** island
load("waterBirdI(s100).RData")

head(posttheta)
dim(posttheta)

###******** posterior pW_j and pH_j
post_pw_source <- posttheta[,1:3]
post_ph_source <- posttheta[,4:7]

postWjI <- post_pw_source
colnames(postWjI) <- c("Poultry", "Other", "WBirds")

postHjI <- post_ph_source
colnames(postHjI) <- c("Poultry", "Other", "WBirds", "Water")

#' water attribution
post_pW_j <- delogit_mtrx(par_mtrx=post_pw_source)
colnames(post_pW_j) <- c("Poultry", "Other", "WBirds", "Ruminants")

smWI <- as.data.frame(post_pW_j) %>%
    mutate(s=colsm$s, m=colsm$m) %>% group_by(s, m) %>%
    gather(source, postW, Poultry:Ruminants) %>%
    mutate(model="Island")   ## for diagnostics section

postWI <- as.data.frame(post_pW_j) %>%
    gather(Source, pW, Poultry:Ruminants) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "Other", "WBirds")),
           model="Island") 


#' human attribution
post_pH_j <- delogit_mtrx(par_mtrx=post_ph_source)
colnames(post_pH_j) <- c("Poultry","Other","WBirds","Water",
                         "Ruminants")

smHI <- as.data.frame(post_pH_j) %>%
    mutate(s=colsm$s, m=colsm$m) %>% group_by(s, m) %>%
    gather(source, postH, Poultry:Ruminants) %>%
    mutate(model="Island")   ## for diagnostics section

postHI <- as.data.frame(post_pH_j) %>%
    gather(Source, pH, Poultry:Ruminants) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "WBirds","Water","Other")),
           model="Island")


#' trace plots
## water
postpwI <- as.data.frame(post_pw_source) %>%
    rename("Poultry"=V1, "Other"=V2, "WBirds"=V3) %>%
    mutate(Simulation=rep(seq(1, simulation), each=Iteration),
           Iteration=rep(seq(1, Iteration), simulation)) %>%
    gather(Source, value, Poultry:WBirds) %>%
    mutate(Source=factor(Source, levels=c("Poultry",
                                          "WBirds", "Other")),
           Iteration=as.numeric(Iteration),
           model="Island")

## human
postphI <- as.data.frame(post_ph_source) %>%
    rename("Poultry"=V1, "Other"=V2, "WBirds"=V3, "Water"=V4) %>%
    mutate(Simulation=rep(seq(1, simulation), each=Iteration),
           Iteration=rep(seq(1, Iteration), simulation)) %>%
    gather(Source, value, Poultry:Water) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "WBirds",
                                          "Water", "Other")),
           Iteration=as.numeric(Iteration),
           model="Island")

##**** Dirichlet
load("waterBirdD(s100).RData")

head(posttheta)
dim(posttheta)

###******** posterior pW_j and pH_j
post_pw_source <- posttheta[,1:3]
post_ph_source <- posttheta[,4:7]

postWjD <- post_pw_source
colnames(postWjD) <- c("Poultry", "Other", "WBirds")

postHjD <- post_ph_source
colnames(postHjD) <- c("Poultry", "Other", "WBirds", "Water")

#' water attribution
post_pW_j <- delogit_mtrx(par_mtrx=post_pw_source)
colnames(post_pW_j) <- c("Poultry", "Other", "WBirds", "Ruminants")

smWD <- as.data.frame(post_pW_j) %>%
    mutate(s=colsm$s, m=colsm$m) %>% group_by(s, m) %>%
    gather(source, postW, Poultry:Ruminants) %>%
    mutate(model="Dirichlet")   ## for diagnostics section

postWD <- as.data.frame(post_pW_j) %>%
    gather(Source, pW, Poultry:Ruminants) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "Other","WBirds")),
           model="Dirichlet") 


#' human attribution
post_pH_j <- delogit_mtrx(par_mtrx=post_ph_source)
colnames(post_pH_j) <- c("Poultry","Other","WBirds","Water",
                         "Ruminants")

smHD <- as.data.frame(post_pH_j) %>%
    mutate(s=colsm$s, m=colsm$m) %>% group_by(s, m) %>%
    gather(source, postH, Poultry:Ruminants) %>%
    mutate(model="Dirichlet")   ## for diagnostics section

postHD <- as.data.frame(post_pH_j) %>%
    gather(Source, pH, Poultry:Ruminants) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "WBirds","Water","Other")),
           model="Dirichlet")


#' trace plots
## water
postpwD <- as.data.frame(post_pw_source) %>%
    rename("Poultry"=V1, "Other"=V2, "WBirds"=V3) %>%
    mutate(Simulation=rep(seq(1, simulation), each=Iteration),
           Iteration=rep(seq(1, Iteration), simulation)) %>%
    gather(Source, value, Poultry:WBirds) %>%
    mutate(Source=factor(Source, levels=c("Poultry",
                                          "WBirds", "Other")),
           Iteration=as.numeric(Iteration),
           model="Dirichlet") 
## human
postphD <- as.data.frame(post_ph_source) %>%  
    rename("Poultry"=V1, "Other"=V2, "WBirds"=V3, "Water"=V4) %>%
    mutate(Simulation=rep(seq(1, simulation), each=Iteration),
           Iteration=rep(seq(1, Iteration), simulation)) %>%
    gather(Source, value, Poultry:Water) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "WBirds",
                                          "Water", "Other")),
           Iteration=as.numeric(Iteration),
           model="Dirichlet")

####################################################################
###  comparison plots for Water and Human attribution            ###
####################################################################
##*************** for simulation 100
postW <- bind_rows(postWI, postWD) %>%
    mutate(pW=pW*100,
           model=factor(model, levels=c("Island", "Dirichlet")))

postH <- bind_rows(postHI, postHD) %>%
    mutate(Source=factor(Source, labels=c("Poultry", "Ruminants",
                                          "WBirds", "Water",
                                          "Other")),
           model=factor(model, levels=c("Island", "Dirichlet")))

## comparison of water attribution between island and dirichlet model
jpeg("fig_ch5_2_bottom.jpg", width = 700, height =500, quality=100)

ggplot(postW, aes(x=Source, y=pW)) +
    geom_violin() + stat_summary(fun.y=median, geom="point", size=1,
                                 color="red") +
    ylab("percentage of water isolates (%)") + xlab("") +
    ylim(0,100) +
    facet_wrap(~model) + theme_bw() +
    theme(text=element_text(size=12, family="serif"))

dev.off()

## comparison of human attribution between two models
jpeg("fig_ch5_3_bottom.jpg", width = 700, height =500, quality=100)

ggplot(postH) +
    geom_boxplot(aes(x=Source, y=pH*100)) +
    ylim(0,100) + 
    ylab("percentage of human cases (%)") + xlab("") +
    facet_wrap(~model) +
    theme_bw() +
    theme(text=element_text(size=12, family="serif"))

dev.off()

## for simulation 100
#' comparison of trace plot
jpeg(file="fig_ch5_5_bottom.jpg", height=500, width=700, quality=100)
## water
postpw <- bind_rows(postpwI, postpwD) %>%
    rename("simulation"=Simulation) %>%
    mutate(model=factor(model, levels=c("Island", "Dirichlet")),
           simulation=as.factor(simulation))

## sample three simulation of pi
#  sample(x=seq(1,100,1), size=3)

ggplot(postpw%>%filter(simulation==10|simulation==36|simulation==98),
       aes(x=Iteration, y=value, color=simulation)) +
    geom_line() +
    xlab("iteration") + ylab("posterior g") +
    facet_grid(model~Source) + 
    theme_bw() + theme(legend.position="bottom",
                       legend.box.margin = margin(-10, 6, -8, 6),
                       text=element_text(size=12, family="serif"))

dev.off()

jpeg("fig_ch5_6_bottom.jpg", width = 700, height =500, quality=100)

postph <- bind_rows(postphI, postphD) %>%
    rename("simulation"=Simulation) %>%
    mutate(model=factor(model, levels=c("Island", "Dirichlet")),
           simulation=as.factor(simulation))

ggplot(postph%>%filter(simulation==10|simulation==36|simulation==98),
       aes(x=Iteration, y=value, color=simulation)) +
    geom_line() +
    xlab("iteration") + ylab("posterior f") +
    facet_grid(model~Source) + 
    theme_bw() + theme(legend.position="bottom",
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
                                          "Other", "WBirds")))

pdf(file="fig_ch5_7_bottom.pdf", width=7, height=5)

ggplot(postW_diag, aes(x=source, y=postW)) + geom_violin() + 
    stat_summary(fun.y=median, geom="point", size=1, color="red") +
    ylab("percentage of water isolates (%)") + ylim(c(0,100)) +
    theme_bw() +
    theme(text=element_text(size=12, family="serif"),
          axis.text.x=element_text(size=9)) +
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
                                          "WBirds", "Water",
                                          "Other")))

pdf(file="fig_ch5_8_bottom.pdf", width=7, height=5)

ggplot(postH_diag, aes(x=source, y=postH)) + geom_boxplot() + 
    ylab("percentage of human cases (%)") + ylim(c(0,100)) +
    theme_bw() +
    theme(text=element_text(size=12, family="serif"),
          axis.text.x=element_text(size=7)) +
    facet_grid(model~group)

dev.off()


