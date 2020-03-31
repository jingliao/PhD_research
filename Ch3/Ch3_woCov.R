##########********  REQUIRED PACKAGES  ********##########
library(MCMCpack)  # dirichlet distribution
library(tidyverse)

##########********  Data processing  ********##########
###******** load the data
#' data path
alltypes <- read.csv("data/human_types_age.csv")

#' frequency table of ST vs. sources
freqtable <- as.data.frame.matrix(table(alltypes$ST,
                                          alltypes$Source))

freqtable <- freqtable %>% select(Human, Poultry, Ruminants, Water, Other)
SeqTypes <- rownames(freqtable)

###******** divide data into humans and sources
#' human data, y_{i}
y_i <- data.frame(Human=freqtable[,1])
rownames(y_i) <- SeqTypes

#' source data, x_{ij}
x_ij <- freqtable[,-1]

source_name <- colnames(x_ij)

#' number of sources
source_no <- dim(x_ij)[2]

#' number of unique STs found on humans and sources
ST_no <- dim(x_ij)[1]

#####********  Parameters in our model to update  ********#####
###******** priors
#' pi_ij and F_j ~Dir(alpha=1)
alpha_pij <- rep(1, ST_no)
pi <- t(rdirichlet(n=source_no, alpha=alpha_pij))

#' F
alpha_Fj <- rep(1,source_no)
F <- as.numeric(rdirichlet(n=1, alpha=alpha_Fj))

###******** gamma
numerator <- matrix(NA, nrow=ST_no, ncol=source_no)

for (i in 1:ST_no){

    numerator[i,] <- pi[i,]*F
}

gamma <- sweep(numerator, 1, FUN="/", rowSums(numerator))

###******** Z
Z <- matrix(NA, nrow=ST_no, ncol=source_no)

for(i in 1:ST_no){

    Z[i,] <- rmultinom(n=1, size=y_i[i,], prob=gamma[i,])

}

##########********  Gibbs sampling  ********##########
###******** updating order: pi_ij, z_ij, F_j
posterior <- list()

## settings
n_burnin <- 1e3
n_ite <- 1e4
n_thin <- 5
total_ite <- n_ite*n_thin+n_burnin

for (m in 1:total_ite){

    ## first update pi
    delta <- x_ij + alpha_pij + Z
    pi <- matrix(NA, ST_no, source_no)

    for (j in 1:source_no){

        pi[,j] <- rdirichlet(1, delta[,j])
    }

    ## update Z
    numerator <- matrix(NA, nrow=ST_no, ncol=source_no)

    for(i in 1:ST_no){

        numerator[i,] <- pi[i,]*F
    }

    gamma <- sweep(numerator, 1, FUN="/", rowSums(numerator))

    Z <- matrix(NA, ST_no, source_no)

    for (i in 1:ST_no){

        Z[i,] <- rmultinom(1, y_i[i,], gamma[i,])
    }

    ## update F

    F <- as.numeric(rdirichlet(n=1, alpha=colSums(Z)+1))
    
    if ((m > n_burnin) && (m%%n_thin==0)){

        posterior[[m]] <- list(Z=Z, pi=pi, F=F)

    }
        
}

#' load RData
load(file="OFbayes_woCov.RData")
load(file="OFbayes_woCov_rerun.RData")

post_F_untidy <- purrr::map(posterior, "F")
names(post_F_untidy) <- seq_along(post_F_untidy)
post_F_untidy[sapply(post_F_untidy, is.null)] <- NULL
post_F <- as.data.frame(do.call(rbind, post_F_untidy))
colnames(post_F) <- source_name

#### source attribution
post_F$Iteration <- seq(1,dim(post_F)[1],by=1)

#' source attribution
post_F_stack <- post_F %>% group_by(Iteration) %>% gather(Source, posteriorF, Poultry:Other) %>% mutate(Source=factor(Source, levels=c("Poultry", "Ruminants", "Water", "Other")))

#' boxplot and trace plot
pdf(file="fig_ch3_1_ch3_3_upper.pdf", height=5, width=7)
pdf(file="OldFbayes_woCov_rerun.pdf", height=7, width=7)

ggplot(post_F_stack) +
    geom_boxplot(aes(x=Source, y=posteriorF*100)) +
    xlab("source") +
    ylab("percentage of cases (%)") +
    theme_bw() + 
    theme(text=element_text(size=12,  family="serif"))

#' trace plot of posterior F_j
ggplot(post_F_stack, aes(x=Iteration, y=posteriorF,
                          group=Source)) +
    geom_line(aes(color=Source)) +
    scale_color_manual(name="source", values=c("red", "green", "blue", "pink")) +
    xlab("iteration") + ylab("posterior attribution F") +
    theme_bw() + theme(legend.position="bottom",
                       text=element_text(size=12, family="serif"))

dev.off()
