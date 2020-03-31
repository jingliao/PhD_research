###******** required packages/functions
library(GGally) # for ggpair scatter plot
library(MCMCpack) # before tidyverse library due to conflict 'select'
library(tidyverse)
library(islandR)

source("new_functionsCov.R")

###******** MCMC settings 
Iteration <- 10000

###******** initial settings
#' deviation of normal prior
sigma <- qsigma <- 1

#' number of simulating sampling distribution
# simulation <- 1
simulation <- 100

##########********  WATER ATTRIBUTION  ********##########
## load the source file for island or dirichlet model
## island
source("Ch5_woWB_sourceI.R")

## Dirichlet
source("Ch5_woWB_sourceD.R")

##**** load samples and combine .csv files as one data frame
base_path <- "samples_post/woWB"

rur_z <- sort(unique(na.omit(humanCov$Rurality)))

## island model
filesI <- list.files(base_path,
                    pattern = "post_sampleR_RI_s[0-9]+.csv") 

#' grab out the data for these chains
# for multiple simulation
allps <- list()

for (f in seq_along(filesI)) {

    ps <- read.csv(file.path(base_path, filesI[f]), header=T)
    
    allps[[f]] <- ps

}
    
post_sampleI <- do.call(rbind, allps)

## Dirichlet model
filesD <- list.files(base_path,
                    pattern = "post_sampleR_RD_s[0-9]+.csv") 

#' grab out the data for these chains
# for multiple simulations
allps <- list()

for (f in seq_along(filesD)) {

    ps <- read.csv(file.path(base_path, filesD[f]), header=T)
    
    allps[[f]] <- ps

}
    
post_sampleD <- do.call(rbind, allps)

#####******** RESULTS VISUALIZATION AND COMPARISON ********####
#' comparison of water parameter
## island
post_thetaWI <- post_sampleI[,c(1,5,6)]
names(post_thetaWI) <- c("simulation", "Poultry", "Other")
post_thetaWI$iteration <- rep(seq(1, Iteration), simulation)

post_thetaWlongI <- post_thetaWI %>%
    gather(source, value, Poultry:Other) %>%
    mutate(source=factor(source, levels=c("Poultry", "Other")),
           model="Island")


## Dirichlet
post_thetaWD <- post_sampleD[,c(1,5,6)]
names(post_thetaWD) <- c("simulation", "Poultry", "Other")
post_thetaWD$iteration <- rep(seq(1, Iteration), simulation)

post_thetaWlongD <- post_thetaWD %>%
    gather(source, value, Poultry:Other) %>%
    mutate(source=factor(source, levels=c("Poultry", "Other")),
           model="Dirichlet")

post_thetaWlong <- bind_rows(post_thetaWlongI, post_thetaWlongD) %>%
    mutate(model=factor(model, levels=c("Island", "Dirichlet")),
           simulation=as.factor(simulation))

#' comparison of human parameters
## island
postH_alpha <- post_sampleI[,c(1,7:9)]
postH_beta1 <- post_sampleI[,c(1,10:12)]

names(postH_alpha) <- names(postH_beta1) <- c("s","Poultry",
                                            "Other", "Water")

postH_alpha$para <- "alpha"
postH_beta1$para <- "beta"
postH_alpha$iteration <- postH_beta1$iteration <- rep(seq(1,
                                                  Iteration), simulation)
post_thetaHI <- rbind(postH_alpha, postH_beta1)

post_thetaHlongI <- post_thetaHI %>%
    gather(source, value, Poultry:Water) %>%
    mutate(source=factor(source,
                         levels=c("Water", "Other", "Poultry")),
           para_level=ifelse(para=="alpha", "intercept",
                             "slope of rurality"),
           model="Island")

## Dirichlet
postH_alpha <- post_sampleD[,c(1,7:9)]
postH_beta1 <- post_sampleD[,c(1,10:12)]

names(postH_alpha) <- names(postH_beta1) <- c("s","Poultry",
                                            "Other", "Water")

postH_alpha$para <- "alpha"
postH_beta1$para <- "beta"
postH_alpha$iteration <- postH_beta1$iteration <- rep(seq(1,
                                                  Iteration), simulation)
post_thetaHD <- rbind(postH_alpha, postH_beta1)

post_thetaHlongD <- post_thetaHD %>%
    gather(source, value, Poultry:Water) %>%
    mutate(source=factor(source,
                         levels=c("Water", "Other", "Poultry")),
           para_level=ifelse(para=="alpha", "intercept",
                             "slope of rurality"),
           model="Dirichlet")

post_thetaHlong <- bind_rows(post_thetaHlongI, post_thetaHlongD) %>%
    mutate(model=factor(model, levels=c("Island", "Dirichlet")),
           simulation=as.factor(s)) %>%
    select(simulation, iteration, source, value, model, para,
           para_level)

## model comparison
jpeg("fig_ch5_appxC15_upper.jpg", width = 700, height =500, quality=100) 

#' trace plots for water parameter

ggplot(post_thetaWlong %>%
       filter(simulation==45|simulation==62|simulation==1),
       aes(x=iteration, y=value, color=simulation)) + 
    geom_line() + 
    ylab("posterior g") + 
    theme_bw() + facet_grid(model~source) +
    theme(legend.position="bottom",
          text=element_text(size=12, family="serif")) 

dev.off()

#' trace plots for human parameters
jpeg("fig_ch5_appxC16_upper.jpg", width = 700, height =500, quality=100) 

ggplot(post_thetaHlong %>%
       filter(simulation==45|simulation==62|simulation==1 & para=="alpha"),
       aes(x=iteration, y=value, color=simulation)) + 
    geom_line() + 
    ylab("intercept") + 
    theme_bw() + facet_grid(model~source) +
    theme(legend.position="bottom",
          text=element_text(size=12, family="serif")) 

dev.off()

jpeg("fig_ch5_appxC17_upper.jpg", width = 700, height =500, quality=100) 

ggplot(post_thetaHlong %>%
       filter(simulation==45|simulation==62|simulation==1 & para=="beta"),
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
                                         select(Poultry, Other))) %>%
    rename(Ruminants=X0)

post_pWjlongI <- post_pW_jI %>%
    gather(source, postW, Poultry:Ruminants) %>%
    mutate(source=factor(source,
                         levels=c("Poultry", "Ruminants", "Other")),
           postW=postW*100,
           model="Island")

## Dirichlet model
post_pW_jD <- data.frame(simulation=post_thetaWD$simulation,
                        iteration=post_thetaWD$iteration,
                        delogit_mtrx(par_mtrx=post_thetaWD %>%
                                         select(Poultry, Other))) %>%
    rename(Ruminants=X0)

post_pWjlongD <- post_pW_jD %>%
    gather(source, postW, Poultry:Ruminants) %>%
    mutate(source=factor(source,
                         levels=c("Poultry", "Ruminants", "Other")),
           postW=postW*100,
           model="Dirichlet")

# combine water cases between Island and Dirichlet
post_pWjlong <- bind_rows(post_pWjlongI, post_pWjlongD) %>%
    mutate(model=factor(model, levels=c("Island", "Dirichlet")))

# for result section: water attribution
pdf(file="fig_ch5_appxC14_upper.pdf", width=7, height=5)

ggplot(post_pWjlong, aes(x=source, y=postW)) +
    geom_violin() + stat_summary(fun.y=median, geom="point",
                                 size=1, color="red") +
    ylab("percentage of water isolates (%)") + ylim(c(0,100)) + 
    theme_bw() +
    theme(text=element_text(size=12,  family="serif")) +
    facet_grid(.~model)

dev.off()

# for diagnostics section: water attribution
group1 <- post_pWjlong %>% filter(simulation==1, iteration==1) %>%
    mutate(group=1)

group2 <- post_pWjlong %>% filter(simulation==1, iteration<=100) %>%
    mutate(group=2)

group3 <- post_pWjlong %>% filter(iteration<=100) %>% mutate(group=3)

postW_diag <- bind_rows(group1, group2, group3) %>%
    mutate(group=factor(group, levels=c("1", "2", "3"),
                        labels=c("s=1, m=1", "s=1, m=100",
                                 "s=100, m=100")))

pdf(file="fig_ch5_9_upper.pdf", width=7, height=5)

ggplot(postW_diag, aes(x=source, y=postW)) + geom_violin() + 
    stat_summary(fun.y=median, geom="point", size=1, color="red") +
    ylab("percentage of water isolates (%)") + ylim(c(0,100)) +
    theme_bw() +
    theme(text=element_text(size=12, family="serif")) +
    facet_grid(model~group)

dev.off()

#################################### human final attribution
load("NTH2_capF.RData")

#' for result section: human attribution
## island model
post_capFI <- postFIDlist[[1]] %>% select(-sim, -ite)

datanamecapF <- post_capFI %>% rename(V1=Rurality,
                                     V2=Poultry,
                                     V3=Other,
                                     V4=Water,
                                     V5=Ruminants)

tidycapFI <- as.data.frame(tidydata(dataname=datanamecapF,
                     sourcename=names(post_capFI)[-1],
                     CIpercent=80)) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other")),
           model="Island") %>%
    rename("Rurality"=V1)

## Dirichlet model
post_capFD <- postFIDlist[[2]] %>% select(-sim, -ite)

datanamecapF <- post_capFD %>% rename(V1=Rurality,
                                     V2=Poultry,
                                     V3=Other,
                                     V4=Water,
                                     V5=Ruminants)

tidycapFD <- as.data.frame(tidydata(dataname=datanamecapF,
                     sourcename=names(post_capFD)[-1],
                     CIpercent=80)) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other")),
           model="Dirichlet") %>%
    rename("Rurality"=V1)

## combine posterior F between island and Dirichlet models
tidycapF <- bind_rows(tidycapFI, tidycapFD) %>%
    mutate(model=factor(model, levels=c("Island", "Dirichlet")))

## final human attribution
pdf(file="fig_ch5_4_upper.pdf", width=7, height=5)
ggplot(tidycapF) +
    geom_ribbon(aes(x=Rurality, ymin=lower*100, ymax=upper*100,
                    fill=Source), alpha=0.4) + 
    geom_line(aes(x=Rurality, y=mean*100, color=Source)) +
    scale_x_continuous(name=NULL, breaks=c(-3,3),
                       labels=c("Rural", "Urban")) +
    scale_y_continuous(name="percentage of human cases (%)",
                       limits=c(0,100)) +
    scale_color_manual(name="source",
                       values=c("red", "green", "blue","pink")) +
    scale_fill_manual(name="source",
                      values = c("red", "green", "blue", "pink")) +
    facet_grid(.~model) + 
    theme_bw() + xlab("rurality") +  
    theme(legend.position="none",
          text=element_text(size=12, family="serif"))
dev.off()

#' for diagnostics section: human attribution
## island model
post_capFg1 <- postFIDlist[[1]] %>% filter(sim==1 & ite==1) %>%
    select(-sim, -ite) 

tidyFg1 <- post_capFg1 %>%
    gather(Source, mean, Poultry:Ruminants) %>%
    mutate(lower=mean, upper=mean, model="Island", group=1) %>%
    select(Rurality, mean, lower, upper, Source, model, group)
    
post_capFg2 <- postFIDlist[[1]] %>% filter(sim==1 & ite<=100) %>%
    select(-sim, -ite)

datanamecapFg2 <- post_capFg2 %>% rename(V1=Rurality,
                                     V2=Poultry,
                                     V3=Other,
                                     V4=Water,
                                     V5=Ruminants)

tidyFg2 <- as.data.frame(tidydata(dataname=datanamecapFg2,
                     sourcename=names(post_capFg2)[-1],
                     CIpercent=80)) %>%
    mutate(model="Island", group=2) %>% rename("Rurality"=V1)

post_capFg3 <- postFIDlist[[1]] %>% filter(ite<=100) %>%
    select(-sim, -ite)

datanamecapFg3 <- post_capFg3 %>% rename(V1=Rurality,
                                     V2=Poultry,
                                     V3=Other,
                                     V4=Water,
                                     V5=Ruminants)

tidyFg3 <- as.data.frame(tidydata(dataname=datanamecapFg3,
                     sourcename=names(post_capFg3)[-1],
                     CIpercent=80)) %>%
    mutate(model="Island", group=3) %>% rename("Rurality"=V1)

tidyFI <- bind_rows(tidyFg1, tidyFg2, tidyFg3)

## Dirichlet model
post_capFg1 <- postFIDlist[[2]] %>% filter(sim==1 & ite==1) %>%
    select(-sim, -ite) 

tidyFg1 <- post_capFg1 %>%
    gather(Source, mean, Poultry:Ruminants) %>%
    mutate(lower=mean, upper=mean, model="Dirichlet", group=1) %>%
    select(Rurality, mean, lower, upper, Source, model, group)

post_capFg2 <- postFIDlist[[2]] %>% filter(sim==1 & ite<=100) %>%
    select(-sim, -ite)

datanamecapFg2 <- post_capFg2 %>% rename(V1=Rurality,
                                     V2=Poultry,
                                     V3=Other,
                                     V4=Water,
                                     V5=Ruminants)

tidyFg2 <- as.data.frame(tidydata(dataname=datanamecapFg2,
                     sourcename=names(post_capFg2)[-1],
                     CIpercent=80)) %>%
    mutate(model="Dirichlet", group=2) %>% rename("Rurality"=V1)

post_capFg3 <- postFIDlist[[1]] %>% filter(ite<=100) %>%
    select(-sim, -ite)

datanamecapFg3 <- post_capFg3 %>% rename(V1=Rurality,
                                     V2=Poultry,
                                     V3=Other,
                                     V4=Water,
                                     V5=Ruminants)

tidyFg3 <- as.data.frame(tidydata(dataname=datanamecapFg3,
                     sourcename=names(post_capFg3)[-1],
                     CIpercent=80)) %>%
    mutate(model="Dirichlet", group=3) %>% rename("Rurality"=V1)

tidyFD <- bind_rows(tidyFg1, tidyFg2, tidyFg3)

## combine samples from island and Dirichlet models
tidyG <- bind_rows(tidyFI, tidyFD) %>%
    mutate(mean=mean*100, lower=lower*100, upper=upper*100,
           Source=factor(Source,
                         levels=c("Poultry", "Ruminants",
                                  "Water", "Other")),
           model=factor(model, levels=c("Island", "Dirichlet")),
           group=factor(group, levels=c(1,2,3),
                        labels=c("s=1, m=1", "s=1, m=100",
                                 "s=100, m=100")))

## human attribution with different simulation and iteration
pdf(file="fig_ch5_10_upper.pdf", width=7, height=5)

ggplot(tidyG) +
    geom_ribbon(aes(x=Rurality, ymin=lower, ymax=upper,
                    fill=Source), alpha=0.4) + 
    geom_line(aes(x=Rurality, y=mean, color=Source)) +
    scale_x_continuous(name=NULL, breaks=c(-3,3),
                       labels=c("Rural", "Urban")) +
    scale_y_continuous(name="percentage of human cases (%)",
                       limits=c(0,100)) +
    scale_color_manual(name="source",
                       values=c("red", "green", "blue","pink")) +
    scale_fill_manual(name="source",
                      values = c("red", "green", "blue", "pink")) +
    facet_grid(model~group) + 
    theme_bw() +   
    theme(legend.position="none",
          text=element_text(size=12, family="serif"),
          axis.text.x=element_text(size=7))

dev.off()
