# required package
library(dplyr)  
library(ggplot2)

# load the data manually summarised from ESR reports

nzrates <- read.csv(file="NZrates.csv")

pdf("fig_ch1_1.pdf", width=7, height=5)
ggplot(data=nzrates, aes(x=year, y=rates)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = 2006, linetype="dashed") +
    scale_x_continuous(name="Reported year", breaks=seq(1990, 2016, 2)) +
    scale_y_continuous(name="Notification rates (100,000 population)", breaks=seq(0, 400, 50), limits=c(0,400)) +
    theme(text=element_text(size=12, family="serif"))
    
dev.off()
