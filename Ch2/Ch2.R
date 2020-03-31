################### for figure 2.2
library(meshblocknz)
library(islandR)
library(dplyr)
library(tidyr)
library(ggplot2)

## meshblock and population per rurality scale
mbMana <- mb2013 %>% filter(RC2013_name=="Manawatu-Wanganui Region",
                             DHB_name=="Mid Central") %>%
    select(MB2013, MB2006, Pop2001, Pop2006, Pop2013, UR2006_num,
           UR2013_num) %>%
    gather(Year, Population, Pop2001:Pop2013) %>%
    extract(Year, into="Year", regex="([0-9]+)", convert=TRUE)

pop13 <- split(mbMana, mbMana$MB2013)

predict_popn <- function(data){
    fit <- lm(Population~Year, data=data)
    data.frame(Year=2005:2016,
               Population=predict(fit, data.frame(Year=2005:2016)))
}

pop_extrap <- bind_rows(lapply(pop13, predict_popn),
                        .id="MB2013") %>%
    mutate(MB2013=as.numeric(MB2013)) %>%
    left_join(mb2013 %>% select(MB2006, MB2013, UR2006_num) %>%
    unique)
    

## load the latest manawatu data and cases per rurality scale
mana <- read.csv(file="data/human_types_age.csv")

mawacase <- mana %>% mutate(UR2006_num=UR_num) %>%
    group_by(UR2006_num, Year) %>% summarize(Cases=n()) 

## combine data
cases <- pop_extrap %>% group_by(UR2006_num, Year) %>% summarize(Population=sum(Population)) %>% left_join(mawacase)

caserates <- cases %>% ungroup %>% mutate(group=ifelse(UR2006_num>=0,
                                                       "urban",
                                                       "rural")) %>%
    group_by(group, Year) %>%
    summarize(pop=sum(Population, na.rm=TRUE),
              report=sum(Cases, na.rm=TRUE)) %>%
                                        # na.rm=T to ignore NA values
    mutate(rate=report/pop*100000) %>%
    mutate(rate=ifelse(Year==2005, rate/0.75, rate)) %>%
    ungroup() %>% mutate(Year=factor(Year))

pdf(file="fig_ch2_2.pdf", width=7, height=7)
## bar chart
ggplot(caserates) +
    geom_col(aes(x=Year, y=rate, fill=group), position="dodge2") +
    scale_fill_manual(values=c("grey30", "grey70")) +
    ylab("case rate per 100,000 population") +
    theme_bw() +
    theme(legend.position=c(0.06,0.9),
          legend.title=element_blank(),
          legend.background=element_rect(fill="transparent"))

dev.off()

################### for figure 2.3

# required package
library(tidyverse)
library(dplyr)   # tools working with data frames
library(tidyr)
library(ggpubr)

mana <- read.csv("data/human_types_age.csv")

## frequency table of ST per source
sts <- mana %>% group_by(Source, ST) %>%
    summarize(frequency=n()) %>% spread(Source, frequency) %>%
    replace_na(list(Human=0, Other=0, Poultry=0, Ruminants=0,
                    Water=0))

## set alpha color
alpha = function(col, alpha){
    rgb(t(col2rgb(col)/255), alpha=alpha) }

cols40 = c(rep(alpha("slateblue4", 0.9), 40), "grey70")
cols61 = c(rep(alpha("slateblue4", 0.9), 5), alpha("red4",0.8),
           rep(alpha("slateblue4",0.9),14), "grey70")

# find top 40 STs in human and sources
top40 = order(-sts$Human)[1:40]
top40p <- order(-sts$Poultry)[1:40]
top40r <- order(-sts$Ruminants)[1:40]
top40w <- order(-sts$Water)[1:40]
top40o <- order(-sts$Other)[1:40]

## calculate percentage
cal_percent <- function(data, category, toprange){

    top_rate <- data[[category]][toprange]/sum(data[[category]])
    bottom_rate <- sum(data[[category]][-toprange])/sum(data[[category]])

    return(c(top_rate, bottom_rate)*100)

}

tabrates <- data.frame(ST=c(sts$ST[top40],"Other"),
                       Human=cal_percent(sts, "Human", top40),
                       Poultry=cal_percent(sts, "Poultry", top40),
                       Ruminants=cal_percent(sts, "Ruminants",
                                             top40),
                       Water=cal_percent(sts, "Water", top40),
                       Other=cal_percent(sts, "Other", top40)) %>%
    gather(Source, rates, Human:Other)

## produce bar plot for each category
bar_plot <- function(data, category, toprange){

    barplot(cal_percent(data, category, toprange),
            names=c(sts$ST[toprange], "Other"), ylim=c(0,20),
            col=cols40, border=NA, yaxt="n", cex.names = 0.9,
            las=2, col.axis="grey10")
    axis(2, col="grey10", col.axis="grey10", las=1, cex.axis=0.8)
    title(main=category)
}


## ST distribution for each category 
pdf("fig_ch2_3.pdf", width=7, height=7)

layout(matrix(c(1,2,3,4), 4, 1, byrow=TRUE))
par(mai=c(0.6,0.6,0.2,0.1)) # margins settings

bar_plot(sts, "Human", top40)
bar_plot(sts, "Poultry", top40)
bar_plot(sts, "Ruminants", top40)
bar_plot(sts, "Other", top40)

mtext("sequence type", side=1, padj=-1, outer=TRUE)
mtext("percentage", side=2, padj=2, outer=TRUE)

dev.off()


