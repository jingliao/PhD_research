#### required packages
library(islandR)
library(dplyr)
library(tidyr)
library(coda)
library(tidyverse)

library(cowplot) # for plot_grid
library(ggplot2)
library(ggjoy)
library(gridExtra)

#### final attribution with prior SD=8
# filter out those where we don't have a location, and have a factor version
# of UR2006_num
manawatu <- read.csv("human_types_age.csv") %>% mutate(UR2006_num=UR_num)

dataset <- manawatu %>% filter(Source != "Human" | !is.na(UR2006_num)) %>%
  mutate(UR2006_fact = factor(UR2006_num))

num_samples = 100

# Fit the sequence type distribution using the island model
st_i = st_fit(formula = Source ~ ST,
         non_primary = "Human",
         data = dataset,
         method = "island",
         sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC,
         samples = num_samples)

# and using the Dirichlet
st_d = st_fit(formula = Source ~ ST,
              non_primary = "Human",
              data = dataset,
              method = "dirichlet",
              samples = num_samples)


# do the attribution
# sensitivity, only focus on linear fitted model
# change sigma from 1 to 4 where prec is variance/precision with default vlaue 1

humans <- dataset %>% filter(Source == "Human")

## calculate F
## default is SD=1, meaning: variance=1 = 1/prec, so prec=0.1
## SD=5, meaning variance=64, so prec=1/64=0.016
at_il <- attribution(ST ~ UR2006_num, st_i, humans,
                     priors=list(theta=list(mean=0, prec=0.016)))
at_dl <- attribution(ST ~ UR2006_num, st_d, humans,
                     priors=list(theta=list(mean=0, prec=0.016)))

df_il <- predict(at_il, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_num=([-0-9]+)", convert=TRUE) %>%
  mutate(GenotypeModel = 'Island', AttributionModel='Linear')

df_dl <- predict(at_dl, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_num=([-0-9]+)", convert=TRUE) %>%
  mutate(GenotypeModel = 'Dirichlet', AttributionModel='Linear')

df <- bind_rows(df_il, df_dl)

## calculate capital F
plot_df <- df %>%
  group_by(GenotypeModel, Source, UR2006_num) %>% summarize(
    m = mean(p),
    li = quantile(p, 0.1),
    ui = quantile(p, 0.8)
  ) %>% ungroup %>%
  mutate(Source = factor(Source, levels=rev(c("Other", "Water", "Ruminants", "Poultry"))))


## extract f
rurality <- as.numeric(unique(model.matrix(ST~UR2006_num,
                                           humans))[,"UR2006_num"])

source_name <- c("Other", "Poultry", "Ruminants")

smallf_il <-  at_il %>% pluck("posterior") %>% map("p") %>% map(.f=as.data.frame)

smallf_dl <-  at_dl %>% pluck("posterior") %>% map("p") %>%
    map(.f=as.data.frame)


flist <- function(data, GModel){

    a <- list()
    
    for (i in 1:length(data)){

        a[[i]]= data[[i]] %>% mutate(UR2006_num=rurality,
                                     GenotypeModel=GModel)
        
    }

    b <- a %>% map_dfr(.f=as.data.frame, .id="Iteration") %>%
        group_by(GenotypeModel, Iteration, UR2006_num) %>%
        gather(Source, f, V1:V3)

    return(b)
    
}

smallf <- bind_rows(flist(data=smallf_il, GModel="Island"),
                    flist(data=smallf_dl, GModel="Dirichlet")) %>%
    mutate(Source=ifelse(Source=="V1", source_name[1],
                  ifelse(Source=="V2", source_name[2],
                         source_name[3]))) 

## calculate f
plotf <- smallf %>% group_by(GenotypeModel, Source, UR2006_num) %>%
    summarize(m=mean(f), li=quantile(f, 0.2),
              ui=quantile(f, 0.8)) %>% ungroup %>%
    mutate(Source=factor(Source, levels=rev(c("Other", "Poultry", "Ruminants"))))


# plot
pdf("fig_ch4_9_bottom.pdf", width=7, height=5)

## f
plot1 <- ggplot(plotf)+
    geom_ribbon(aes(x=UR2006_num, ymin=li, ymax=ui, fill=Source),
                alpha=0.3) +
    geom_line(aes(x=UR2006_num, y=m, col=Source), lwd=1) +
    facet_wrap(GenotypeModel~.) +
    scale_x_continuous(name=NULL, breaks=c(-3,3),
                       labels=c("Rural", "Urban"), expand=c(0,0)) +
    scale_y_continuous(name="posterior mean of f",
                       limits=c(-20,20), expand=c(0,0)) +
    scale_color_manual(name="source", values = c(Poultry="red",
                                                 Ruminants="green",
                                                 Other="pink")) +
    scale_fill_manual(name="source", values = c(Poultry="red",
                                                Ruminants="green",
                                                Other="pink")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(hjust=c(-0.5,1.1)),
          axis.text.y = element_text(vjust=c(-0.5,rep(0.5, 3), 1.1)),
          text=element_text(size=12, family="serif"))


## capital F
plot2 <- ggplot(plot_df) +
    geom_ribbon(aes(x=UR2006_num, ymin=li*100, ymax=ui*100, fill=Source), alpha=0.3) +
    geom_line(aes(x=UR2006_num, y=m*100, col=Source), lwd=1) +
    facet_wrap(GenotypeModel~.) +
    scale_x_continuous(name=NULL, breaks=c(-3,3), labels=c("Rural", "Urban"), expand=c(0,0)) +
    scale_y_continuous(name="percentage of cases (%)", limits=c(0,100), expand=c(0,0)) +
    scale_color_manual(name="source", values = c(Poultry="red", Ruminants="green", Other="pink", Water="blue")) +
    scale_fill_manual(name="source", values = c(Poultry="red", Ruminants="green", Other="pink", Water="blue")) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "bottom",
          axis.text.x = element_text(hjust=c(-0.5,1.1)),
          axis.text.y = element_text(vjust=c(-0.5,rep(0.5, 3), 1.1)),
          text=element_text(size=12, family="serif"))

plot_grid(plot1, plot2, ncol=1, align="v")

dev.off()


################################## attribution based on SD=1
load('attribution_fits.Rdata')

rurality <- as.numeric(unique(model.matrix(ST~UR2006_num,
                                           humans))[,"UR2006_num"])

source_name <- c("Other", "Poultry", "Ruminants")

## calculate capital F
plot_dfo <- df %>%
  group_by(GenotypeModel, AttributionModel, Source, UR2006_num) %>% summarize(
    m = mean(p),
    li = quantile(p, 0.1),
    ui = quantile(p, 0.8)
  ) %>% ungroup %>%
  mutate(Source = factor(Source, levels=rev(c("Other", "Water", "Ruminants", "Poultry"))))


## extract f
smallf_ilo <-  at_il %>% pluck("posterior") %>% map("p") %>% map(.f=as.data.frame)

smallf_dlo <-  at_dl %>% pluck("posterior") %>% map("p") %>%
    map(.f=as.data.frame)

smallfo <- bind_rows(flist(data=smallf_ilo, GModel="Island"),
                    flist(data=smallf_dlo, GModel="Dirichlet")) %>%
    mutate(Source=ifelse(Source=="V1", source_name[1],
                  ifelse(Source=="V2", source_name[2],
                         source_name[3]))) 

## calculate f
plotfo <- smallfo %>% group_by(GenotypeModel, Source, UR2006_num) %>%
    summarize(m=mean(f), li=quantile(f, 0.2),
              ui=quantile(f, 0.8)) %>% ungroup %>%
    mutate(Source=factor(Source, levels=rev(c("Other", "Poultry", "Ruminants"))))


# plot
pdf("fig_ch4_9_upper.pdf", width=7, height=5)

## f
plot1 <- ggplot(plotfo)+
    geom_ribbon(aes(x=UR2006_num, ymin=li, ymax=ui, fill=Source),
                alpha=0.3) +
    geom_line(aes(x=UR2006_num, y=m, col=Source), lwd=1) +
    facet_wrap(GenotypeModel~.) +
    scale_x_continuous(name=NULL, breaks=c(-3,3),
                       labels=c("Rural", "Urban"), expand=c(0,0)) +
    scale_y_continuous(name="posterior mean of f",
                       limits=c(-20,20), expand=c(0,0)) +
    scale_color_manual(name="source", values = c(Poultry="red",
                                                 Ruminants="green",
                                                 Other="pink")) +
    scale_fill_manual(name="source", values = c(Poultry="red",
                                                Ruminants="green",
                                                Other="pink")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(hjust=c(-0.5,1.1)),
          axis.text.y = element_text(vjust=c(-0.5,rep(0.5, 3), 1.1)),
          text=element_text(size=12, family="serif"))


## capital F
plot2 <- ggplot(plot_dfo%>%filter(AttributionModel=="Linear")) +
    geom_ribbon(aes(x=UR2006_num, ymin=li*100, ymax=ui*100, fill=Source), alpha=0.3) +
    geom_line(aes(x=UR2006_num, y=m*100, col=Source), lwd=1) +
    facet_wrap(GenotypeModel~.) +
    scale_x_continuous(name=NULL, breaks=c(-3,3), labels=c("Rural", "Urban"), expand=c(0,0)) +
    scale_y_continuous(name="percentage of cases (%)", limits=c(0,100), expand=c(0,0)) +
    scale_color_manual(name="source", values = c(Poultry="red", Ruminants="green", Other="pink", Water="blue")) +
    scale_fill_manual(name="source", values = c(Poultry="red", Ruminants="green", Other="pink", Water="blue")) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "bottom",
          axis.text.x = element_text(hjust=c(-0.5,1.1)),
          axis.text.y = element_text(vjust=c(-0.5,rep(0.5, 3), 1.1)),
          text=element_text(size=12, family="serif"))

plot_grid(plot1, plot2, ncol=1, align="v")

dev.off()


################################ attribution with intervention
# filter out those where we don't have a location
dataset <- manawatu %>% filter(Source != "Human" | !is.na(UR2006_num)) %>%
  mutate(Intervention = ifelse(Year < 2008, "Before", "After"))

num_samples = 100

# Fit the sequence type distribution using the island model
st_i = st_fit(formula = Source ~ ST,
              non_primary = "Human",
              data = dataset,
              method = "island",
              sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC,
              samples = num_samples)

# and using the Dirichlet
st_d = st_fit(formula = Source ~ ST,
              non_primary = "Human",
              data = dataset,
              method = "dirichlet",
              samples = num_samples)

# do the attribution
humans <- dataset %>% filter(Source == "Human")
at_if <- attribution(ST ~ UR2006_num*Intervention, st_i, humans)
at_id <- attribution(ST ~ UR2006_num*Intervention, st_d, humans)

df_il <- predict(at_if, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_num=([-0-9]+)", convert=TRUE, remove=FALSE) %>%
  extract(X, into="Intervention", regex="Intervention([A-Za-z]+)=1") %>%
  replace_na(replace=list(Intervention="After")) %>%
  mutate(GenotypeModel = 'Island', AttributionModel='Linear')

df_dl <- predict(at_id, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_num=([-0-9]+)", convert=TRUE, remove=FALSE) %>%
  extract(X, into="Intervention", regex="Intervention([A-Za-z]+)=1") %>%
  replace_na(replace=list(Intervention="After")) %>%
  mutate(GenotypeModel = 'Dirichlet', AttributionModel='Linear')

df <- bind_rows(df_il, df_dl)

## capital F
plot_df <- df %>%
  group_by(GenotypeModel, AttributionModel, Source, UR2006_num, Intervention) %>% summarize(
    m = mean(p),
    li = quantile(p, 0.1),
    ui = quantile(p, 0.8)
  ) %>% ungroup %>%
  mutate(Source = factor(Source, levels=rev(c("Other", "Water", "Ruminants", "Poultry"))),
         Intervention = factor(Intervention, levels=c("Before", "After"), labels=c("2005-2007", "2008-2014")))

# plot
pdf("fig_ch4_8.pdf", width=7, height=5)
ggplot(plot_df) +
  geom_ribbon(aes(x=UR2006_num, ymin=li*100, ymax=ui*100, fill=Source), alpha=0.3) +
  geom_line(aes(x=UR2006_num, y=m*100, col=Source), lwd=1) +
  facet_grid(GenotypeModel~Intervention) +
  scale_x_continuous(name=NULL, breaks=c(-3,3), labels=c("Rural", "Urban"), expand=c(0,0)) +
  scale_y_continuous(name="percentage of cases (%)", limits=c(0,100), expand=c(0,0)) +
  scale_color_manual(name="source", values = c(Poultry="red", Ruminants="green", Other="pink", Water="blue")) +
  scale_fill_manual(name="source", values = c(Poultry="red", Ruminants="green", Other="pink", Water="blue")) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(hjust=c(-0.1,1.1)),
        axis.text.y = element_text(vjust=c(-0.1,rep(0.5, 3), 1.1)),
        text=element_text(size=12, family="serif"))
dev.off()


############################################# combine two chains

load('attribution_fits.Rdata')
ch1_il <-  at_il %>% pluck("posterior") %>% map("theta")
ch1_dl <-  at_dl %>% pluck("posterior") %>% map("theta")


load('attribution_fits_chain2.Rdata')
ch2_il <-  at_il %>% pluck("posterior") %>% map("theta")
ch2_dl <-  at_dl %>% pluck("posterior") %>% map("theta")


tidy_theta <- function(data){

    tlist <- list()

    for(i in 1:length(data)){

    tlist[[i]] <- as.numeric(t(data[[i]]))

    }

    output <- as.data.frame(do.call(rbind, tlist))
    colnames(output) <- c("intercept1", "intercept2", "intercept3",
                          "slope1", "slope2", "slope3")

    return(output)
}

tidych1_il <- tidy_theta(data=ch1_il)
tidych1_dl <- tidy_theta(data=ch1_dl)
tidych2_il <- tidy_theta(data=ch2_il)
tidych2_dl <- tidy_theta(data=ch2_dl)

## combine two chains for island and Dirichlet models
cmb_i <- mcmc.list(as.mcmc(tidych1_il), as.mcmc(tidych2_il))
cmb_d <- mcmc.list(as.mcmc(tidych1_dl), as.mcmc(tidych2_dl))

pdf("fig_ch4_4_ch4_6.pdf", width=7, height=5)

par(family='serif')
plot(cmb_i)
autocorr.plot(cmb_i)

dev.off()

pdf("fig_ch4_5_ch4_7.pdf", width=7, height=5)

par(family='serif')
plot(cmb_d)
autocorr.plot(cmb_d)

dev.off()

