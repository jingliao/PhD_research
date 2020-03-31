##**** required libraries ----
library(tidyverse)
library(ggplot2)
library(coda)
library(cowplot)
library(gridExtra)
library(ggrepel) ## index STs with a line
library(tidyr)
library(ggjoy)
library(ggpubr) ## for command ggarrange

##**** functions required ----

getparaout <- function(data, parameter){

    return(output <- data %>% map(parameter) %>%
               map_dfr(.f=as.data.frame, .id="Iteration"))

}

####**** choose the data you want to load ----
## w/ covariate, c=0
load(file="wCov_c0H2.RData")
## chain 2
load(file="wCov_c0H2chain2.RData")
## chain 3
load(file="wCov_c0H2chain3.RData")

## w/ covariate, random c
load(file="wCov_randCH2SD0.05.RData")
## chain 2
load(file="wCov_randCSD0.05ch2.RData")
## chain 3
load(file="wCov_randCSD0.05ch3.RData")

####**** call parameters out of the posterior list ----
## w/ covariate, Rurality ----
#' new codes
## whole period
ru <- unique(H)[,2]

post_theta <- getparaout(data=posterior, parameter='theta') %>%
    mutate(Iteration=as.numeric(Iteration))

theta <- post_theta %>% select(-Iteration)

sname=colnames(theta)

smallf <- capitalF <- list()
oddno <- seq(1, dim(theta)[1], by=2)

for(i in order(oddno)){

    indt <- as.matrix(theta[oddno[i]:(oddno[i]+1),])
    smallf[[i]] <- data.frame(rurality=ru, unique(H)%*%indt)
    capitalF[[i]] <- data.frame(rurality=ru, computeF(theta=indt, H=unique(H)))

}

## for o=0
ch6F_o0 <- do.call(rbind, capitalF) %>% rename(Ruminants=V4) %>%
    group_by(rurality) %>% gather(source, Fprob, Other:Ruminants) %>%
    mutate(Source=factor(source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other")),
           UR2006_num=rurality) %>%
    group_by(UR2006_num, Source) %>%
    summarize(m=mean(Fprob), li=quantile(Fprob, 0.2),
              ui=quantile(Fprob, 0.8)) %>% ungroup

## for random o
ch6F_randO <- do.call(rbind, capitalF) %>% rename(Ruminants=V4) %>%
    group_by(rurality) %>% gather(source, Fprob, Other:Ruminants) %>%
    mutate(Source=factor(source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other")),
           UR2006_num=rurality) %>%
    group_by(UR2006_num, Source) %>%
    summarize(m=mean(Fprob), li=quantile(Fprob, 0.2),
              ui=quantile(Fprob, 0.8)) %>% ungroup

## if mixed o=0 and random o
combF_randO <- bind_rows(ch6F_o0%>%mutate(typeC="O=0"),
                         ch6F_randO%>%mutate(typeC="random O"))

## regression parameters beta
oddrowindex <- seq(from=1, to=dim(theta)[1], by=2)

alpha <- post_theta[oddrowindex,] %>%
    gather(source, alpha, Other:Water)

beta <- post_theta[oddrowindex+1,] %>%
    gather(source, beta, Other:Water)

## chain 1
thetalongch1 <- alpha %>% bind_cols(beta %>%
                                 select(-c(source, Iteration))) %>%
    mutate(chain="1")

## chain 2
thetalongch2 <- alpha %>% bind_cols(beta %>%
                                 select(-c(source, Iteration))) %>%
    mutate(chain="2")

## chain 3
thetalongch3 <- alpha %>% bind_cols(beta %>%
                                 select(-c(source, Iteration))) %>%
    mutate(chain="3")

## combine chains
thetalong <- bind_rows(thetalongch1, thetalongch2, thetalongch3)

############ comparison between ch3 and ch6 (o=0)
#' ch6, o=0
ch6capF <- ch6F_o0
#' ch3
ch3capF <- read.table(file="ch3_capFH2.txt", row.names=NULL) %>%
    rename(UR2006_num=V1, Other=V2, Poultry=V3, Water=V4,
           Ruminants=V5) %>% group_by(UR2006_num) %>%
    gather(Source, Fprob, Other:Ruminants) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other"))) %>%
    group_by(UR2006_num, Source) %>%
    summarize(m=mean(Fprob), li=quantile(Fprob, 0.2),
              ui=quantile(Fprob, 0.8)) %>% ungroup


combF_H2o0 <- bind_rows(ch3capF%>%mutate(model="chapter 3", H=2),
                      ch6capF%>%mutate(model="chapter 6", H=2))

## attribution capital F
pdf("fig_ch6_1.pdf", width=7, height=5)
ggplot(combF_H2o0) +
    geom_ribbon(aes(x=UR2006_num, ymin=li*100, ymax=ui*100,
                    fill=Source), alpha=0.3) +
    geom_line(aes(x=UR2006_num, y=m*100, col=Source), lwd=1) +
    scale_x_continuous(name=NULL, breaks=c(-3,3), labels=c("Rural", "Urban"), expand=c(0,0)) +
    scale_y_continuous(name="percentage of cases (%)", limits=c(0,100), expand=c(0,0)) +
    scale_color_manual(name="source", values = c(Poultry="red", Ruminants="green", Other="pink", Water="blue")) +
    scale_fill_manual(name="source", values = c(Poultry="red", Ruminants="green", Other="pink", Water="blue")) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(hjust=c(-0.5,1.1)),
          axis.text.y = element_text(vjust=c(-0.5,rep(0.5, 3), 1.1)),
          text=element_text(size=12, family="serif")) +
    facet_wrap(.~model)

dev.off()

############ comparison in ch6 between o=0 and random O
## chapter 6, o=0 vs. random Cauchy o
pdf("fig_ch6_2.pdf", width=7, height=5)

ggplot(combF_randO) +
    geom_ribbon(aes(x=UR2006_num, ymin=li*100, ymax=ui*100,
                    fill=Source), alpha=0.3) +
    geom_line(aes(x=UR2006_num, y=m*100, col=Source), lwd=1) +
    scale_x_continuous(name=NULL, breaks=c(-3,3), labels=c("Rural", "Urban"), expand=c(0,0)) +
    scale_y_continuous(name="percentage of cases (%)", limits=c(0,100), expand=c(0,0)) +
    scale_color_manual(name="source", values = c(Poultry="red", Ruminants="green", Other="pink", Water="blue")) +
    scale_fill_manual(name="source", values = c(Poultry="red", Ruminants="green", Other="pink", Water="blue")) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(hjust=c(-0.5,1.1)),
          axis.text.y = element_text(vjust=c(-0.5,rep(0.5, 3), 1.1)),
          text=element_text(size=12, family="serif")) +
    facet_wrap(.~typeC)

dev.off()


############################## plot for posterior pi_source, rand O
post_pis <- getparaout(data=posterior, parameter='pi_s') %>%
    mutate(Iteration=as.numeric(Iteration))

post_pis$ST <- rep(as.numeric(rownames(x)), n_iter)

postpis_long <- post_pis %>% group_by(Iteration, ST) %>%
    gather(source, pis, Other:Ruminants) %>%
    mutate(source=factor(source,
                         levels=c("Poultry","Ruminants","Water",
                                  "Other")))

############################## plot for posterior pi_human, rand O
post_pih <- getparaout(data=posterior, parameter='pi_h') %>%
    mutate(Iteration=as.numeric(Iteration))

post_pih$ST <- rep(as.numeric(rownames(x)), n_iter)

postpih_long <- post_pih %>% group_by(Iteration, ST) %>%
    gather(source, pih, Other:Ruminants) %>%
    mutate(source=factor(source,
                         levels=c("Poultry","Ruminants","Water",
                                  "Other")))

######################### comparison of pi between h and s, rand O
post_pisST <- posterior %>% map('pi_s') %>%
    map_dfr(.f = function(x) {
        as.data.frame(x) %>% tibble::rownames_to_column('ST') },
        .id = 'Iteration') %>%
    gather(Source, P, Other:Ruminants) %>%
    mutate(Iteration = as.numeric(Iteration)) %>%
    group_by(ST, Source) %>% summarize(mean=mean(P)) %>%
    mutate(group="s")

post_pihST <- posterior %>% map('pi_h') %>%
    map_dfr(.f = function(x) {
        as.data.frame(x) %>% tibble::rownames_to_column('ST') },
        .id = 'Iteration') %>%
    gather(Source, P, Other:Ruminants) %>%
    mutate(Iteration = as.numeric(Iteration)) %>%
    group_by(ST, Source) %>% summarize(mean=mean(P)) %>%
    mutate(group="h")

post_allpi <- rbind(post_pisST, post_pihST) %>%
    spread(group, mean) %>%
    mutate(Source=factor(Source,
                         levels=c("Poultry", "Ruminants",
                                  "Water", "Other")))

diff_pi <- post_allpi %>% group_by(ST, Source) %>%
    mutate(diff=(h-s)^2) 

top3diff <- diff_pi %>% group_by(Source) %>% 
    top_n(n=3, wt=diff) %>% arrange(Source, -diff) %>%
    mutate(Variable="most")
top3diff

bottom3diff <- diff_pi %>% group_by(Source) %>% 
    top_n(n=-3, wt=diff) %>% arrange(Source, diff) %>%
    mutate(Variable="least")

difftab <- rbind(top3diff, bottom3diff) %>% 
  select(Source, Variable, ST, diff, h, s) %>% 
  arrange(Source)

## STs with highest posteriore mean of pis for each source
topSTtab <- diff_pi %>% arrange(Source) %>%
    mutate(type=as.numeric(ST)) %>% group_by(Source, type) %>%
    filter(s>0.04 | h>0.04) %>% arrange(Source, desc(s)) %>% ungroup %>%
    select(Source, type, s, h, diff)
topSTtab

######################### boxplot of pi_s for top few ST
## boxplot for source
#' calculate posterior mean for sources and humans
postpis_mean <- postpis_long %>% group_by(source, ST) %>%
    summarise(mean=mean(pis))

postpih_mean <- postpih_long %>% group_by(source, ST) %>%
    summarise(mean=mean(pih))

#' with threshold of posterior mean of pi_s and pi_h, 0.02
pistop002 <- postpis_mean %>% group_by(source) %>%
    filter(mean>0.02) %>% select(source, ST)

pihtop002 <- postpih_mean %>% group_by(source) %>%
    filter(mean>0.02) %>% select(source, ST)

## to save time to select the top few ST, split it per source
#' split it per source 
stop002P <- pistop002 %>% filter(source=="Poultry")
stop002R <- pistop002 %>% filter(source=="Ruminants")
stop002W <- pistop002 %>% filter(source=="Water")
stop002O <- pistop002 %>% filter(source=="Other")

htop002P <- pihtop002 %>% filter(source=="Poultry")
htop002R <- pihtop002 %>% filter(source=="Ruminants")
htop002W <- pihtop002 %>% filter(source=="Water")
htop002O <- pihtop002 %>% filter(source=="Other")

#' select the top STs per source and combine them
pistopP <- post_pis %>% select(Iteration, Poultry, ST) %>%
    filter(ST %in% stop002P$ST) %>% gather(source, Pis, Poultry)

pistopR <- post_pis %>% select(Iteration, Ruminants, ST) %>%
    filter(ST %in% stop002R$ST) %>% gather(source, Pis, Ruminants)

pistopW <- post_pis %>% select(Iteration, Water, ST) %>%
    filter(ST %in% stop002W$ST) %>% gather(source, Pis, Water)

pistopO <- post_pis %>% select(Iteration, Other, ST) %>%
    filter(ST %in% stop002O$ST) %>% gather(source, Pis, Other)

pistopbox <- rbind(pistopP, pistopR, pistopW, pistopO) %>%
    mutate(source=factor(source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other")))

pihtopP <- post_pih %>% select(Iteration, Poultry, ST) %>%
    filter(ST %in% htop002P$ST) %>% gather(source, Pih, Poultry)
pihtopR <- post_pih %>% select(Iteration, Ruminants, ST) %>%
    filter(ST %in% htop002R$ST) %>% gather(source, Pih, Ruminants)
pihtopW <- post_pih %>% select(Iteration, Water, ST) %>%
    filter(ST %in% htop002W$ST) %>% gather(source, Pih, Water)
pihtopO <- post_pih %>% select(Iteration, Other, ST) %>%
    filter(ST %in% htop002O$ST) %>% gather(source, Pih, Other)
pihtopbox <- rbind(pihtopP, pihtopR, pihtopW, pihtopO) %>%
    mutate(source=factor(source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other")))

## plots
## given random c without seeding
pdf("fig_ch6_3.pdf", width=7, height=5)
## given random c with seeding
pdf("postPis_SD0.05seed.pdf", width=7, height=5)

ggplot(pistopbox, aes(x=as.factor(ST), y=Pis)) + geom_boxplot() +
    facet_wrap(~source, ncol=2, scales="free_x") + xlab("ST") +
    ylab(expression("posterior "*pi^s)) +
    scale_y_continuous(limits=c(0,1)) +
    theme_bw() +
    theme(legend.position="none",
          text=element_text(size=12, family="serif"))

ggplot(pihtopbox, aes(x=as.factor(ST), y=Pih)) + geom_boxplot() + 
    facet_wrap(~source, ncol=2, scales="free_x") + xlab("ST") +
    ylab(expression("posterior "*pi^h)) + theme_bw() +
    theme(legend.position="none",
          text=element_text(size=12, family="serif"))


dev.off()

## comparison of pi for each source
pdf("fig_ch6_4.pdf", width=7, height=5)
piP <- ggplot(diff_pi %>% filter(Source=="Poultry"),
       aes(x=s, y=h), label=ST) +  geom_point() +
    geom_abline(slope=1, intercept=0, linetype="dotted") +
    xlim(0,.4) + ylim(0,.4) +  xlab("source") + ylab("human") + 
    geom_text_repel(aes(label=ifelse(diff>=min(top3diff[top3diff$Source=="Poultry", "diff"]), as.character(ST),'')),
              hjust=1.5,vjust=-2.5) + 
    theme_bw() + 
    ggtitle(expression("Posterior "*mu[pi]*" for Poultry")) +
    theme(text=element_text(size=12, family="serif"))


piR <- ggplot(diff_pi %>% filter(Source=="Ruminants"),
       aes(x=s, y=h), label=ST) +  geom_point() +
    geom_abline(slope=1, intercept=0, linetype="dotted") +
    xlim(0,.4) + ylim(0,.4) +  xlab("source") + ylab("human") + 
    geom_text_repel(aes(label=ifelse(diff>=min(top3diff[top3diff$Source=="Ruminants", "diff"]), as.character(ST),'')),
              hjust=1.5,vjust=-2) + 
    theme_bw() + 
    ggtitle(expression("Posterior "*mu[pi]*" for Ruminants")) +
    theme(text=element_text(size=12, family="serif"))


piW <- ggplot(diff_pi %>% filter(Source=="Water"),
       aes(x=s, y=h), label=ST) +  geom_point() +
    geom_abline(slope=1, intercept=0, linetype="dotted") +
    xlim(0,0.1) + ylim(0,0.1) +  xlab("source") + ylab("human") + 
    geom_text_repel(aes(label=ifelse(diff>=min(top3diff[top3diff$Source=="Water", "diff"]), as.character(ST),'')),
              hjust=1,vjust=-2.5) + 
    theme_bw() + 
    ggtitle(expression("Posterior "*mu[pi]*" for Water")) +
    theme(text=element_text(size=12, family="serif"))


piO <- ggplot(diff_pi %>% filter(Source=="Other"),
       aes(x=s, y=h), label=ST) +  geom_point() +
    geom_abline(slope=1, intercept=0, linetype="dotted") +
    xlim(0,0.25) + ylim(0,0.25) + xlab("source") + ylab("human") + 
    geom_text_repel(aes(label=ifelse(diff>=min(top3diff[top3diff$Source=="Other", "diff"]), as.character(ST),'')),
              hjust=1.5,vjust=-2) + 
    theme_bw() + 
    ggtitle(expression("Posterior "*mu[pi]*" for Other")) +
    theme(text=element_text(size=12, family="serif"))

 
plot_grid(piP, piR, ncol=1)
plot_grid(piW, piO, ncol=1)

dev.off()

######################### Density and median with 80% HPD for O
post_o <- getparaout(data=posterior, parameter='c') %>%
    mutate(Iteration=as.numeric(Iteration))

post_o$ST <- rep(as.numeric(rownames(x)), n_iter)

posto_long <- post_o %>% group_by(Iteration, ST) %>%
    gather(source, posto, Other:Ruminants) %>%
    mutate(source=factor(source,
                         levels=c("Poultry", "Ruminants",
                                  "Water","Other")))

## posterior median of c vs. ST for all sources
## plots for posterior o ONLY for random c setting

summaryo <- posto_long %>% group_by(ST, source) %>%
    summarize(median=median(posto), sd=sd(posto),
              li=quantile(posto, 0.2), ui=quantile(posto, 0.8))

summaryo04 <- summaryo %>% filter(abs(median)>=0.4) %>%
    arrange(source, ST)

HPDmin <- floor(min(summaryo04$li))
HPDmax <- ceiling(max(summaryo04$ui))

postOmed04 <- summaryo %>% filter(abs(median)<0.4) %>%
    mutate(source=factor(source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other"),
                         labels=c("Poultry (|median|<0.4)",
                                          "Ruminants (|median|<0.4)",
                                          "Water",
                                          "Other")))


pdf("fig_ch6_5.pdf", width = 7, height =5)

ggplot(postOmed04) +
    geom_density(aes(x=median, fill=source), alpha=0.5) +
    theme_bw() +
    theme(text=element_text(size=12, family="serif"))


## poultry
oP <- ggplot(summaryo04%>%filter(source=="Poultry"),
             aes(x=as.factor(ST), y=median, label=median)) +
    geom_point(color="red", size=0.5) +
    geom_text(aes(label=round(median,3)), hjust=1.2,
              vjust=0.6, size=3) + 
    geom_linerange(aes(x=as.factor(ST), ymin=li, ymax=ui)) + 
    xlab("") + ylab("") + 
    scale_y_continuous(limits=c(HPDmin, HPDmax)) +
    ggtitle("Poultry") +
    theme_bw() +
    theme(text=element_text(size=12, family="serif"))

## Ruminants
oR <- ggplot(summaryo04%>%filter(source=="Ruminants"),
             aes(x=as.factor(ST), y=median, label=median)) +
    geom_point(color="red", size=1) +
    geom_text(aes(label=round(median, 3)), hjust=1.2,
              vjust=0.6, size=3) +
    geom_linerange(aes(x=as.factor(ST), ymin=li, ymax=ui)) + 
    xlab("") + ylab("") + 
    scale_y_continuous(limits=c(HPDmin, HPDmax)) +
    ggtitle("Ruminants") +
    theme_bw() +
    theme(text=element_text(size=12, family="serif"))

plotnoSeed <- ggarrange(oP, oR, ncol=1, nrow=2) # w/o seeding

annotate_figure(plotnoSeed,
                left = text_grob("posterior median of O", rot = 90,
                                 size=12, family="serif"),
                bottom = text_grob("ST", size=12, family="serif"))

dev.off()

## comparison of trace plots for regression parameters between O=0
## and random O, given three MCMC chains

## regression parameters beta
oddrowindex <- seq(from=1, to=dim(theta)[1], by=2)

alpha <- post_theta[oddrowindex,] %>%
    gather(source, alpha, Other:Water)

beta <- post_theta[oddrowindex+1,] %>%
    gather(source, beta, Other:Water)

## chain 1
thetalongch1 <- alpha %>% bind_cols(beta %>%
                                 select(-c(source, Iteration))) %>%
    mutate(chain="1")

## chain 2
thetalongch2 <- alpha %>% bind_cols(beta %>%
                                 select(-c(source, Iteration))) %>%
    mutate(chain="2")

## chain 3
thetalongch3 <- alpha %>% bind_cols(beta %>%
                                 select(-c(source, Iteration))) %>%
    mutate(chain="3")

## combine chains
thetalong <- bind_rows(thetalongch1, thetalongch2, thetalongch3)

## trace plot for alpha and beta
thetagather <- thetalong %>% group_by(chain, Iteration, source) %>%
    gather(theta, value, alpha:beta) %>%
    select(chain, Iteration, source, theta, value) %>%
    ungroup %>%
    mutate(theta=ifelse(theta=="alpha", "intercept",
                        "slope of rurality"),
           source=factor(source,
                         levels=c("Poultry", "Water", "Other")))

## O=0
pdf("fig_ch6_6_upper.pdf", width=7, height=5)
## O is random
pdf("fig_ch6_6_bottom.pdf", width=7, height=5)

ggplot(thetagather, aes(x=Iteration, y=value)) +
    geom_line(aes(color=chain), alpha=0.7) +
    scale_color_manual(name="chain",
                       values=c("grey", "blue", "red")) +
    xlab("iteration") + ylab("value") +
    theme_bw() + facet_grid(theta~source, scale="free") + 
    theme(legend.position="bottom",
          text=element_text(size=12, family="serif"))

dev.off()


## comparison for posterior o, given random o, w/ seeding
#** choose 1 data to load
load(file="wCov_randC_seed1SD0.05.RData")
load(file="wCov_randC_seed2SD0.05.RData")

### start reading out samples for posterior o
post_o <- getparaout(data=posterior, parameter='c') %>%
    mutate(Iteration=as.numeric(Iteration))

post_o$ST <- rep(as.numeric(rownames(x)), n_iter)


## seed 1 types
posto_longs1 <- post_o %>% group_by(Iteration, ST) %>%
    gather(source, posto, Other:Ruminants) %>%
    mutate(source=factor(source,
                         levels=c("Poultry", "Ruminants",
                                  "Water","Other")))

summaryos1 <- posto_longs1 %>% group_by(ST, source) %>%
    summarize(median=median(posto), sd=sd(posto),
              li=quantile(posto, 0.2), ui=quantile(posto, 0.8))

postOmed04s1 <- summaryos1 %>% filter(abs(median)<0.4) %>%
    mutate(source=factor(source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other"),
                         labels=c("Poultry (|median|<0.4)",
                                          "Ruminants (|median|<0.4)",
                                          "Water",
                                  "Other")),
           seed="seed 1")


typess1 <- posto_longs1 %>% group_by(ST, source) %>%
    summarize(median=median(posto)) %>% filter(abs(median)>=0.4) %>%
    select(source, ST) 

## seed 2 types
posto_longs2 <- post_o %>% group_by(Iteration, ST) %>%
    gather(source, posto, Other:Ruminants) %>%
    mutate(source=factor(source,
                         levels=c("Poultry", "Ruminants",
                                  "Water","Other")))

summaryos2 <- posto_longs2 %>% group_by(ST, source) %>%
    summarize(median=median(posto), sd=sd(posto),
              li=quantile(posto, 0.2), ui=quantile(posto, 0.8))

postOmed04s2 <- summaryos2 %>% filter(abs(median)<0.4) %>%
    mutate(source=factor(source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other"),
                         labels=c("Poultry (|median|<0.4)",
                                          "Ruminants (|median|<0.4)",
                                          "Water",
                                  "Other")),
           seed="seed 2")

typess2 <- posto_longs2 %>% group_by(ST, source) %>%
    summarize(median=median(posto)) %>% filter(abs(median)>=0.4) %>%
    select(source, ST) 

## combine types between seeds
types <- full_join(typess1, typess2) %>% arrange(source, ST)

## filter types from all posterior samples, seed 1
summaryo04s1 <- left_join(types, posto_longs1) %>%
    group_by(ST, source) %>% summarize(median=median(posto),
                                       sd=sd(posto),
                                       li=quantile(posto, 0.2),
                                       ui=quantile(posto, 0.8)) %>%
    mutate(seed="set.seed(123)") %>% ungroup

## filter types from all posterior samples, seed 2
summaryo04s2 <- left_join(types, posto_longs2) %>%
    group_by(ST, source) %>% summarize(median=median(posto),
                                       sd=sd(posto),
                                       li=quantile(posto, 0.2),
                                       ui=quantile(posto, 0.8)) %>%
    mutate(seed="set.seed(321)") %>% ungroup

## combine posteriors for these types 
comb_summaryo04 <- bind_rows(summaryo04s1, summaryo04s2)

pdf("fig_ch6_7.pdf", width=7, height=5)

ggplot(comb_summaryo04,
       aes(x=as.factor(ST), y=median, color=seed)) +
    geom_point(position=position_dodge(width=0.5), size=1) +
    geom_linerange(aes(x=as.factor(ST), ymin=li, ymax=ui),
                   position=position_dodge(width=0.5)) + 
    xlab("ST") + ylab("posterior median of O") +
    facet_grid(source~., scales="free") +
    theme_bw() +
    theme(legend.position="bottom",
          text=element_text(size=12, family="serif"))

dev.off()


