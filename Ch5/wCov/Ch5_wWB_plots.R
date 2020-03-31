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

###******** MCMC - Metropolis-Hasting Random walk 
Iteration <- 10000
Thinning <- 100
Burnin <- 3000

##**** load samples and combine .csv files as one data frame
base_path <- "samples_post/wWB"

rur_z <- sort(unique(na.omit(humanCov$Rurality)))

## island model
files <- list.files(base_path,
                    pattern = "post_sampleR_RI_WB_[0-9]+.csv")

# grab out the data for these chains
# for multiple simulations
allps <- list()

for (f in seq_along(files)) {

    ps <- read.csv(file.path(base_path, files[f]), header=T)
    
    allps[[f]] <- ps

}
    
post_sampleI <- do.call(rbind, allps) 

## Dirichlet model
files <- list.files(base_path,
                    pattern = "post_sampleR_RD_WB_[0-9]+.csv")


# grab out the data for these chains
# for multiple simulations
allps <- list()

for (f in seq_along(files)) {

    ps <- read.csv(file.path(base_path, files[f]), header=T)
    
    allps[[f]] <- ps

}
    
post_sampleD <- do.call(rbind, allps) 

#####******** RESULTS VISUALIZATION AND COMPARISON ********####
#' comparison of water parameter
## island
post_thetaWI <- post_sampleI[,c(1,5:7)]
names(post_thetaWI) <- c("simulation", "Other", "Poultry", "WBirds")
post_thetaWI$iteration <- rep(seq(1, Iteration), simulation)

post_thetaWlongI <- post_thetaWI %>%
    gather(source, value, Other:WBirds) %>%
    mutate(source=factor(source,
                         levels=c("Poultry", "Other", "WBirds")),
           model="Island")

## Dirichlet
post_thetaWD <- post_sampleD[,c(1,5:7)]
names(post_thetaWD) <- c("simulation", "Other", "Poultry", "WBirds")
post_thetaWD$iteration <- rep(seq(1, Iteration), simulation)

post_thetaWlongD <- post_thetaWD %>%
    gather(source, value, Other:WBirds) %>%
    mutate(source=factor(source,
                         levels=c("Poultry", "Other", "WBirds")),
           model="Dirichlet")

post_thetaWlong <- bind_rows(post_thetaWlongI, post_thetaWlongD) %>%
    mutate(model=factor(model, levels=c("Island", "Dirichlet")),
           simulation=as.factor(simulation))


#' comparison of human parameters
## island
postH_alpha <- post_sampleI[,c(1,8:11)]
postH_beta1 <- post_sampleI[,c(1,12:15)]

names(postH_alpha) <- names(postH_beta1) <- c("s","Other", "Poultry",
                                            "WBird", "Water")

postH_alpha$para <- "alpha"
postH_beta1$para <- "beta"
postH_alpha$iteration <- postH_beta1$iteration <- rep(seq(1,
                                                  Iteration), simulation)
post_thetaHI <- rbind(postH_alpha, postH_beta1)

post_thetaHlongI <- post_thetaHI %>%
    gather(source, value, Other:Water) %>%
    mutate(para_level=ifelse(para=="alpha", "intercept",
                             "slope of rurality"),
           model="Island")

## Dirichlet
postH_alpha <- post_sampleD[,c(1,8:11)]
postH_beta1 <- post_sampleD[,c(1,12:15)]

names(postH_alpha) <- names(postH_beta1) <- c("s","Other", "Poultry",
                                            "WBird", "Water")

postH_alpha$para <- "alpha"
postH_beta1$para <- "beta"
postH_alpha$iteration <- postH_beta1$iteration <- rep(seq(1,
                                                  Iteration), simulation)
post_thetaHD <- rbind(postH_alpha, postH_beta1)

post_thetaHlongD <- post_thetaHD %>%
    gather(source, value, Other:Water) %>%
    mutate(para_level=ifelse(para=="alpha", "intercept",
                             "slope of rurality"),
           model="Dirichlet")

post_thetaHlong <- bind_rows(post_thetaHlongI, post_thetaHlongD) %>%
    mutate(model=factor(model, levels=c("Island", "Dirichlet")),
           source=factor(source,
                         levels=c("WBird", "Water", "Other",
                                  "Poultry")),
           simulation=as.factor(s)) %>%
    select(simulation, iteration, source, value, model, para, para_level)

## model comparison
jpeg("fig_ch5_appxC15_bottom.jpg", width = 700, height =500, quality=100) 

#' trace plots for water parameter
ggplot(post_thetaWlong%>%
       filter(simulation==82|simulation==33|simulation==42),
       aes(x=iteration, y=value, color=simulation)) + 
    geom_line() + 
    ylab("posterior g") + 
    theme_bw() + facet_grid(model~source) +
    theme(legend.position="bottom",
          text=element_text(size=12, family="serif")) 

dev.off()

jpeg("fig_ch5_appxC16_bottom.jpg", width = 700, height =500, quality=100) 

#' trace plots for human parameters
ggplot(post_thetaHlong %>%
       filter(simulation==82|simulation==33|simulation==42 & para=="alpha"),
       aes(x=iteration, y=value, color=simulation)) + 
    geom_line() + 
    ylab("intercept") + 
    theme_bw() + facet_grid(model~source) +
    theme(legend.position="bottom",
          text=element_text(size=12, family="serif")) 

dev.off()

jpeg("fig_ch5_appxC17_bottom.jpg", width = 700, height =500, quality=100) 

ggplot(post_thetaHlong %>%
       filter(simulation==82|simulation==33|simulation==42 & para=="beta"),
       aes(x=iteration, y=value, color=simulation)) + 
    geom_line() + 
    ylab("slope of rurality") + 
    theme_bw() + facet_grid(model~source) +
    theme(legend.position="bottom",
          text=element_text(size=12, family="serif")) 

dev.off()
#################################### water final attribution
## island model
post_pW_jI <- data.frame(simulation=post_thetaWI$simulation,
                        iteration=post_thetaWI$iteration,
                        delogit_mtrx(par_mtrx=post_thetaWI %>%
                                         select(Other:WBirds))) %>%
    rename(Ruminants=X0)

post_pWjlongI <- post_pW_jI %>% 
    gather(source, postW, Other:Ruminants) %>%
    mutate(source=factor(source,
                         levels=c("Poultry", "Ruminants",
                                  "Other", "WBirds")),
           postW=postW*100,
           model="Island")


## Dirichlet model
post_pW_jD <- data.frame(simulation=post_thetaWD$simulation,
                        iteration=post_thetaWD$iteration,
                        delogit_mtrx(par_mtrx=post_thetaWD %>%
                                         select(Other:WBirds))) %>%
    rename(Ruminants=X0)

post_pWjlongD <- post_pW_jD %>%
    gather(source, postW, Other:Ruminants) %>%
    mutate(source=factor(source,
                         levels=c("Poultry", "Ruminants",
                                  "Other", "WBirds")),
           postW=postW*100,
           model="Dirichlet")

# combine water cases between Island and Dirichlet
post_pWjlong <- bind_rows(post_pWjlongI, post_pWjlongD) %>%
    mutate(model=factor(model, levels=c("Island", "Dirichlet")))

# for result section: water attribution
pdf(file="fig_ch5_appxC14_bottom.pdf", width=7, height=5)

ggplot(post_pWjlong, aes(x=source, y=postW)) +
    geom_violin() + stat_summary(fun.y=median, geom="point",
                                 size=1, color="red") +
    ylab("percentage of water isolates (%)") + ylim(c(0,100)) + 
    theme_bw() +
    theme(text=element_text(size=12,  family="serif")) +
    facet_grid(.~model)

dev.off()

# for diagnostics section: water attribution
group1 <- post_pWjlong %>% filter(simulation==2, iteration==1) %>%
    mutate(group=1)

group2 <- post_pWjlong %>% filter(simulation==2, iteration<=100) %>%
    mutate(group=2)

group3 <- post_pWjlong %>% filter(iteration<=100) %>% mutate(group=3)

postW_diag <- bind_rows(group1, group2, group3) %>%
    mutate(group=factor(group, levels=c("1", "2", "3"),
                        labels=c("s=1, m=1", "s=1, m=100",
                                 "s=100, m=100")))

pdf(file="fig_ch5_9_bottom.pdf", width=7, height=5)

ggplot(postW_diag, aes(x=source, y=postW)) + geom_violin() + 
    stat_summary(fun.y=median, geom="point", size=1, color="red") +
    ylab("percentage of water isolates (%)") + ylim(c(0,100)) +
    theme_bw() +
    theme(text=element_text(size=12, family="serif"),
          axis.text.x=element_text(size=9)) +
    facet_grid(model~group)

dev.off()

## load posterior samples from both models
load("WBH2_capF.RData")

#' for result section: human attribution
## island model
post_capFI <- postFIDlist[[1]] %>% select(-sim, -ite)

datanamecapF <- post_capFI %>% rename("V1"=Rurality,"V2"=Other,
                                      "V3"=Poultry,
                                     "V4"=WBirds, "V5"=Water,
                                     "V6"=Ruminants)

tidycapFI <- as.data.frame(tidydataWB(dataname=datanamecapF,
                     sourcename=names(post_capFI)[-1],
                     CIpercent=80)) %>%
    mutate(model="Island") %>% rename("Rurality"=V1)

## Dirichlet 
post_capFD <- postFIDlist[[2]] %>% select(-sim, -ite)

datanamecapF <- post_capFD %>% rename("V1"=Rurality,"V2"=Other,
                                      "V3"=Poultry,
                                     "V4"=WBirds, "V5"=Water,
                                     "V6"=Ruminants)

tidycapFD <- as.data.frame(tidydataWB(dataname=datanamecapF,
                     sourcename=names(post_capFD)[-1],
                     CIpercent=80)) %>%
    mutate(model="Dirichlet") %>% rename("Rurality"=V1)


## combine posterior F between island and Dirichlet models
tidycapF <- bind_rows(tidycapFI, tidycapFD) %>%
    mutate(model=factor(model, levels=c("Island", "Dirichlet")),
           Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other",
                                          "WBirds")))

pdf(file="fig_ch5_4_bottom.pdf", width=7, height=5)

ggplot(tidycapF) +
    geom_ribbon(aes(x=Rurality, ymin=lower*100, ymax=upper*100,
                    fill=Source), alpha=0.4) + 
    geom_line(aes(x=Rurality, y=mean*100, color=Source)) +
    scale_x_continuous(name=NULL, breaks=c(-3,3),
                       labels=c("Rural", "Urban")) +
    scale_y_continuous(name="percentage of human cases (%)",
                       limits=c(0,100)) +
    scale_color_manual(name="source",
                       values=c("red", "green", "blue",
                                "pink", "black")) +
    scale_fill_manual(name="source",
                      values = c("red", "green", "blue",
                                 "pink", "black")) +
    facet_grid(.~model) + 
    theme_bw() + xlab("rurality") +  
    theme(legend.position="bottom",
          text=element_text(size=12, family="serif"))

dev.off()


#' for diagnostics section: human attribution

newpostFIDlist <- list()
newpostFIDlist[[1]] <- postFIDlist[[1]] %>% filter(ite>500)
newpostFIDlist[[2]] <- postFIDlist[[2]] %>% filter(ite>500)

## island model
post_capFg1 <- newpostFIDlist[[1]] %>% filter(sim==1 & ite==501) %>%
    select(-sim, -ite) 

tidyFg1 <- post_capFg1 %>%
    gather(Source, mean, Other:Ruminants) %>%
    mutate(lower=mean, upper=mean, model="Island", group=1) %>%
    select(Rurality, mean, lower, upper, Source, model, group)
    
post_capFg2 <- newpostFIDlist[[1]] %>% filter(sim==1 & ite<=600) %>%
    select(-sim, -ite)

datanamecapFg2 <- post_capFg2 %>% rename(V1=Rurality,
                                     V2=Other,
                                     V3=Poultry,
                                     V4=WBirds,
                                     V5=Water,
                                     V6=Ruminants)

tidyFg2 <- as.data.frame(tidydataWB(dataname=datanamecapFg2,
                     sourcename=names(post_capFg2)[-1],
                     CIpercent=80)) %>%
    mutate(model="Island", group=2) %>% rename("Rurality"=V1)

post_capFg3 <- postFIDlist[[1]] %>% filter(ite<=600) %>%
    select(-sim, -ite)

datanamecapFg3 <- post_capFg3 %>% rename(V1=Rurality,
                                     V2=Other,
                                     V3=Poultry,
                                     V4=WBirds,
                                     V5=Water,
                                     V6=Ruminants)

tidyFg3 <- as.data.frame(tidydataWB(dataname=datanamecapFg3,
                     sourcename=names(post_capFg3)[-1],
                     CIpercent=80)) %>%
    mutate(model="Island", group=3) %>% rename("Rurality"=V1)

tidyFI <- bind_rows(tidyFg1, tidyFg2, tidyFg3)


## Dirichlet model
post_capFg1 <- newpostFIDlist[[2]] %>% filter(sim==1 & ite==501) %>%
    select(-sim, -ite) 

tidyFg1 <- post_capFg1 %>%
    gather(Source, mean, Other:Ruminants) %>%
    mutate(lower=mean, upper=mean, model="Dirichlet", group=1) %>%
    select(Rurality, mean, lower, upper, Source, model, group)
    
post_capFg2 <- newpostcapFDlist %>% filter(sim==1 & ite<=600) %>%
    select(-sim, -ite)

datanamecapFg2 <- post_capFg2 %>% rename(V1=Rurality,
                                     V2=Other,
                                     V3=Poultry,
                                     V4=WBirds,
                                     V5=Water,
                                     V6=Ruminants)

tidyFg2 <- as.data.frame(tidydataWB(dataname=datanamecapFg2,
                     sourcename=names(post_capFg2)[-1],
                     CIpercent=80)) %>%
    mutate(model="Dirichlet", group=2) %>% rename("Rurality"=V1)

post_capFg3 <- newpostcapFDlist %>% filter(ite<=600) %>%
    select(-sim, -ite)

datanamecapFg3 <- post_capFg3 %>% rename(V1=Rurality,
                                     V2=Other,
                                     V3=Poultry,
                                     V4=WBirds,
                                     V5=Water,
                                     V6=Ruminants)

tidyFg3 <- as.data.frame(tidydataWB(dataname=datanamecapFg3,
                     sourcename=names(post_capFg3)[-1],
                     CIpercent=80)) %>%
    mutate(model="Dirichlet", group=3) %>% rename("Rurality"=V1)

tidyFD <- bind_rows(tidyFg1, tidyFg2, tidyFg3)

## combine samples from island and Dirichlet models
tidyG <- bind_rows(tidyFI, tidyFD) %>%
    mutate(mean=mean*100, lower=lower*100, upper=upper*100,
           Source=factor(Source,
                         levels=c("Poultry", "Ruminants",
                                  "Water", "Other", "WBirds")),
           model=factor(model, levels=c("Island", "Dirichlet")),
           group=factor(group, levels=c(1,2,3),
                        labels=c("s=1, m=1", "s=1, m=100",
                                 "s=100, m=100")))

## human attribution with different simulation and iteration
pdf(file="fig_ch5_10_bottom.pdf", width=7, height=5)

ggplot(tidyG) +
    geom_ribbon(aes(x=Rurality, ymin=lower, ymax=upper,
                    fill=Source), alpha=0.4) + 
    geom_line(aes(x=Rurality, y=mean, color=Source)) +
    scale_x_continuous(name=NULL, breaks=c(-3,3),
                       labels=c("Rural", "Urban")) +
    scale_y_continuous(name="percentage of human cases (%)",
                       limits=c(0,100)) +
    scale_color_manual(name="source",
                       values=c("red", "green", "blue",
                                "pink", "black")) +
    scale_fill_manual(name="source",
                      values = c("red", "green", "blue",
                                 "pink", "black")) +
    facet_grid(model~group) + 
    theme_bw() +   
    theme(legend.position="bottom",
          text=element_text(size=12, family="serif"),
          axis.text.x=element_text(size=7))

dev.off()


