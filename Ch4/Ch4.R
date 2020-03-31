#### required packages
library(islandR)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggjoy)

#### final attribution
manawatu <- read.csv("human_types_age.csv") %>% mutate(UR2006_num=UR_num)

dataset <- manawatu %>% filter(Source != "Human" | !is.na(UR2006_num)) %>%
  mutate(UR2006_fact = factor(UR2006_num))


# Fit the sequence type distribution using the island model
num_samples = 100

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
at_if <- attribution(ST ~ UR2006_fact, st_i, humans)
at_il <- attribution(ST ~ UR2006_num, st_i, humans)
at_df <- attribution(ST ~ UR2006_fact, st_d, humans)
at_dl <- attribution(ST ~ UR2006_num, st_d, humans)

df_if <- predict(at_if, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_fact([-0-9]+)=1", convert=TRUE) %>%
  replace_na(replace=list(UR2006_num=-3)) %>% mutate(GenotypeModel = 'Island', AttributionModel='Categorical')
df_il <- predict(at_il, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_num=([-0-9]+)", convert=TRUE) %>%
  mutate(GenotypeModel = 'Island', AttributionModel='Linear')
df_df <- predict(at_df, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_fact([-0-9]+)=1", convert=TRUE) %>%
  replace_na(replace=list(UR2006_num=-3)) %>% mutate(GenotypeModel = 'Dirichlet', AttributionModel='Categorical')
df_dl <- predict(at_dl, FUN=identity) %>% extract(X, into="UR2006_num", regex="UR2006_num=([-0-9]+)", convert=TRUE) %>%
  mutate(GenotypeModel = 'Dirichlet', AttributionModel='Linear')

df <- bind_rows(df_if, df_il, df_df, df_dl)

## calculate F
plot_df <- df %>%
  group_by(GenotypeModel, AttributionModel, Source, UR2006_num) %>% summarize(
    m = mean(p),
    li = quantile(p, 0.1),
    ui = quantile(p, 0.8)
  ) %>% ungroup %>%
  mutate(Source = factor(Source, levels=rev(c("Other", "Water", "Ruminants", "Poultry"))))

# plot
pdf("fig_ch4_2.pdf", width=7, height=5)
ggplot(plot_df) +
    geom_ribbon(aes(x=UR2006_num, ymin=li*100, ymax=ui*100, fill=Source), alpha=0.3) +
    geom_line(aes(x=UR2006_num, y=m*100, col=Source), lwd=1) +
    facet_grid(GenotypeModel~AttributionModel) +
    scale_x_continuous(name=NULL, breaks=c(-3,3), labels=c("Rural", "Urban"), expand=c(0,0)) +
    scale_y_continuous(name="percentage of cases (%)", limits=c(0,100), expand=c(0,0)) +
    scale_color_manual(name="source", values = c(Poultry="red", Ruminants="green", Other="pink", Water="blue")) +
    scale_fill_manual(name="source", values = c(Poultry="red", Ruminants="green", Other="pink", Water="blue")) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(hjust=c(-0.5,1.1)),
          axis.text.y = element_text(vjust=c(-0.5,rep(0.5, 3), 1.1)),
          text=element_text(size=12, family="serif"))

dev.off()

######################################### DIC calculation
# compute DIC's
dic_data <-
  rbind(data.frame(DIC = DIC(at_il), GenotypeModel='Island', AttributionModel='Linear'),
            b=data.frame(DIC = DIC(at_if), GenotypeModel='Island', AttributionModel='Categorical'),
            data.frame(DIC = DIC(at_dl), GenotypeModel='Dirichlet', AttributionModel='Linear'),
            data.frame(DIC = DIC(at_df), GenotypeModel='Dirichlet', AttributionModel='Categorical'))


# Table 4.1 : DIC values
dic_data %>% mutate(DIC=round(DIC, 1),
                    Model = factor(AttributionModel, levels=c("Linear", "Categorical"))) %>%
  select(Model, GenotypeModel, DIC) %>%
  spread(GenotypeModel, DIC) %>%
  knitr::kable('latex', booktabs=TRUE, linesep='', format.args=list(big.mark=','))

######################################### sampling distribution
# Fit the sequence type distribution using the island model
num_samples <- 10000

st_i = st_fit(formula = Source ~ ST,
              non_primary = "Human",
              data = manawatu,
              method = "island",
              sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC,
              samples=num_samples)

st_il <- lapply(1:num_samples, function(i) {
    data.frame(ST=as.numeric(st_i$types),
               st_i$sampling_distribution[,,i], Iteration=i) })
st_id <- do.call(rbind, st_il)

# and using the Dirichlet
###################################### Dir(alpha=1)
st_d = st_fit(formula = Source ~ ST,
              non_primary = "Human",
              data = manawatu,
              method = "dirichlet",
              samples = num_samples)

###################################### Dir(alpha=50)
st_d = st_fit(formula = Source ~ ST,
              non_primary = "Human",
              data = manawatu,
              method = "dirichlet",
              samples = num_samples,
              prior=50)

st_dl <- lapply(1:num_samples, function(i) {
    data.frame(ST=as.numeric(st_d$types),
               st_d$sampling_distribution[,,i], Iteration=i) })
st_dd <- do.call(rbind, st_dl)

# Join them up
st <- bind_rows(cbind(st_id, Model='Island'),
                cbind(st_dd, Model='Dirichlet')) %>%
    gather(Source, P, -ST, -Iteration, -Model)

final <- st %>% group_by(Iteration, Model, ST) %>%
  mutate(P = P/sum(P)) %>% spread(Model, P)

# Filter out the STs to plot
sts <- c(403, 2343, 2026, 474)
plot_dat <- final %>% filter(ST %in% sts) %>%
    gather(Model, P, Dirichlet, Island) %>% ungroup %>%
    mutate(ST = factor(ST, levels=sts, labels = paste0("ST-", sts)),
           Scale = ifelse(Model == "Dirichlet" & ST %in% c("ST-403","ST-2343"), 2, 1.2))

# and do the plot
###################################### Dir(alpha=1)
pdf("fig_ch4_3.pdf", width=7, height=5)

###################################### Dir(alpha=50)
pdf("fig_ch4_10.pdf", width=7, height=5)

ggplot(plot_dat) + geom_joy(aes(x=P, y=Source, fill=Model, scale=Scale), alpha=0.7, size=0.1, bandwidth=0.01) +
  facet_wrap(~ST, ncol=2) +
  ylab("") +
  scale_x_continuous(name = "p(source | ST)", limits=c(0,1), expand = c(0,0)) +
  scale_fill_manual(name="model", values=c("steelblue2", "brown")) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(hjust=c(0.2,rep(0.5, 3), 1.1)),
          text=element_text(size=12, family="serif"))
dev.off()

