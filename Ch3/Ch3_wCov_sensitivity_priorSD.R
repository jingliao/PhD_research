##########********  REQUIRED PACKAGES  ********##########
library(MCMCpack)  # dirichlet distribution
library(dplyr)
library(tidyr)
library(tidyverse)
library(cowplot)
library(GGally) # for pairs plot in ggplot

source("functions_wCov.R")

##########********  Data processing  ********##########
###******** load the data
#' data path
base_path <- "data"

alltypes <- read.csv(file.path(base_path, "human_types_age.csv"))

#' frequency table of ST vs. sources
freqtable <- as.data.frame.matrix(table(alltypes$ST,
                                          alltypes$Source))

SeqTypes <- rownames(freqtable)

#' change order of sources to make the last one as the baseline
freqtable <- freqtable %>% select(Human, Poultry, Ruminants, Water, Other)

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

## human types data
humantypes <- alltypes %>% filter(Source=="Human")

############################## design matrix

## Setup the design matrix, and map from humans to ST
#' with rurality
humandata <- humantypes %>% dplyr::select(ST, UR_num) %>% na.omit
human_no <- nrow(humandata)
H <- model.matrix(ST ~ UR_num, data=humandata)

# Parameters in our model to update
Z <- matrix(NA, human_no, source_no)
pi <- matrix(NA, ST_no, source_no)

## fill from the priors
sample_theta <- function(para_rowno, norm_sd){ 
  
  matrix(rep(rnorm(source_no-1,0,norm_sd), para_rowno), nrow=para_rowno)
  
}

############################## choose overall SD for normal prior
#' smaller sigma than 1
sigma <- 0.5
theta <- sample_theta(para_rowno=ncol(H), norm_sd=sigma)

#' larger sigma than 1
sigma <- 2
theta <- sample_theta(para_rowno=ncol(H), norm_sd=sigma)

#' pi_ij and F_j ~Dir(alpha=1)
alpha_pij <- rep(1, ST_no)

pi <- t(rdirichlet(source_no, alpha_pij))

## To fill in Z, we need to sample using gamma, which involves
## computing F first
f <- H %*% theta
F <- delogit_mtrx(f)

## and then gamma
## Setup the map from each human row to each source type row 
st_index <- match(humandata$ST, SeqTypes)
num <- pi[st_index,] * F
gamma <- num / rowSums(num)

## FINALLY, we can sample Z:
## This is really the first iteration of our MCMC
for (i in 1:human_no) {
    Z[i,] <- rmultinom(1, 1, gamma[i,])
}

human_index <- matrix(FALSE, ST_no, human_no)

for (i in 1:ST_no) {
    human_index[i,] <- as.numeric(SeqTypes[i]) == humandata$ST
}

posterior <- list()

burnin <- 1000
thin <- 5
nite <- 10000

totalm <- nite*thin+burnin

accept <- reject <- 0

qsigma <- 1

## Now we iterate
for (m in 1:totalm) {
    
    ## First update pi (Gibbs)
    sumZ <- matrix(0, ST_no, source_no)
    
    for (i in 1:ST_no) {
        sumZ[i,] <- colSums(Z[human_index[i,],,drop=FALSE])
        }
    delta <- x_ij + alpha_pij + sumZ

    pi <- matrix(NA, ST_no, source_no)
    
    for (j in 1:source_no){

        pi[,j] <- rdirichlet(1, delta[,j])
    }
    
    ## Now update Z (you already have pi and F)

    num <- pi[st_index,]*F
    gamma <- num/rowSums(num)

    Z <- matrix(NA, human_no, source_no)
    
    for(i in 1:human_no){

        Z[i,] <- rmultinom(n=1, size=1, prob=gamma[i,])

    }

    ## Now update alpha
    order <- sample(length(theta), length(theta))

    for (id in order){

        theta_new <- theta
        ## don't bother qsigma, as long as it is the same onwards
        theta_new[id] <- rnorm(n=1, mean=theta[id], sd=qsigma)

        ## compute f, F
        f_new <- H %*% theta_new
        F_new <- delogit_mtrx(f_new)

        ## Compute LogLik
        diff_F <- matrix(log(F_new)-log(F), ncol=source_no)

        log_Fj_ratio <- sum(Z*diff_F)
        
        log_prior_ratio <- (theta[id]^2-theta_new[id]^2)/(2*sigma^2)

        log_alpha <- log_Fj_ratio+log_prior_ratio
        
        ## Decide if you accept
        if (log(runif(1)) < min(0, log_alpha)){

            theta <- theta_new
            F <- F_new
            accept <- accept + 1

        }

        else {

            reject <- reject + 1
            
        }

    }

    ## Thinning, storage etc.
    if( (m > burnin) && (m%%thin==0) ) {
        
        posterior[[m]] <- list(theta=as.numeric(t(theta)),
                               Z=Z,
                               pi=pi)
    }

}

############################## load file
#' smaller SD on normal prior
load(file="OFbayes_wCovH2_SD0.5.RData")

#' larger SD on normal prior
load(file="OFbayes_wCovH2_SD2.RData")


post_theta <- do.call(rbind, purrr::map(posterior, "theta"))

############################## choose way to calculate f and F
postf_list <- postcapF_list <- list()

for(i in 1:nite){
  
    ind_theta <- matrix(as.numeric(post_theta[i,]),
                        nrow=dim(theta)[1], byrow=TRUE)

    ### f
    ind_f <- unique(H) %*% ind_theta
   
    #' with covariate - rurality    
    postf_list[[i]] <- cbind(unique(H)[,2], ind_f)

    ### F
    ind_capF <- delogit_mtrx(ind_f)

    #' with covariate - rurality    
    postcapF_list[[i]] <- cbind(unique(H)[,2], ind_capF)
    
}

post_f <- as.data.frame(do.call(rbind, postf_list))
post_capF <- as.data.frame(as.data.frame(do.call(rbind, postcapF_list))) 


############################ source attribution
##' with covariate, rurality, f and F 
tidyf <- as.data.frame(tidydata(dataname=post_f,
                                sourcename=source_name[-4],
                                CIpercent=80,
                                H=H)) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "Water"))) %>%
    rename("Rurality"=V1)

tidycapF <- as.data.frame(tidydata(dataname=post_capF,
                     sourcename=source_name,
                     CIpercent=80,
                     H=H)) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other"))) %>%
    rename("Rurality"=V1)

######################### diagnosis with plots, with rurality
## f and F
saplotf <- ggplot(tidyf) +
    geom_ribbon(aes(x=Rurality, ymin=lower, ymax=upper,
                    fill=Source), alpha=0.4) + 
    geom_line(aes(x=Rurality, y=mean, color=Source)) + 
    scale_y_continuous(name="posterior mean of f", limits=c(-6,10)) +
    scale_x_continuous(name=NULL, breaks=c(-3,3),
                       labels=c("Rural", "Urban")) + 
    scale_color_manual(values=c("red", "green", "blue")) +
    scale_fill_manual(values = c("red", "green", "blue")) +
    theme_bw() + 
    theme(legend.position="none",
          text=element_text(size=12, family="serif")) 


saplotcapF <- ggplot(tidycapF) +
    geom_ribbon(aes(x=Rurality, ymin=lower*100, ymax=upper*100,
                    fill=Source), alpha=0.4) + 
    geom_line(aes(x=Rurality, y=mean*100, color=Source)) +
    scale_x_continuous(name=NULL,
                       breaks=c(-3,3), labels=c("Rural", "Urban")) +
    scale_y_continuous(name="percentage of cases (%)",
                       limits=c(0,100)) +
    scale_color_manual(name="source",
                       values=c("red", "green", "blue","pink")) +
    scale_fill_manual(name="source",
                      values = c("red", "green", "blue", "pink")) +
    theme_bw() + xlab("") +  
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position="bottom",
          legend.box.margin = margin(-20, 6, -8, 6),
          text=element_text(size=12, family="serif"))


## alpha and beta
post_theta <- as.data.frame(post_theta)
post_theta  <- post_theta %>% mutate(iteration=as.numeric(rownames(post_theta)))

post_beta0 <- post_theta %>% dplyr::select(iteration, V1:V3) %>%
    mutate(parameter="beta0") %>% gather(source, value, V1:V3)
post_beta1 <- post_theta %>% dplyr::select(iteration, V4:V6) %>%
        mutate(parameter="beta1") %>% gather(source, value, V4:V6)

postt_all <- rbind(post_beta0, post_beta1) %>% 
    mutate(source=factor(recode(source, "V1"="Poultry",
                                "V2"="Ruminants",
                                "V3"="Water", "V4"="Poultry",
                                "V5"="Ruminants", "V6"="Water"),
                         levels=c("Water", "Poultry", "Ruminants")),
           para_level=factor(recode(parameter, "beta0"="intercept",
                                   "beta1"="slope of rurality"),
                            levels=c("intercept", "slope of rurality")))

# scatter plot
newbeta0 <- postt_all %>% dplyr::select(-para_level) %>%
    filter(parameter=="beta0") %>% spread(source, value) %>%
    dplyr::select(Poultry, Ruminants, Water)

newbeta1 <- postt_all %>% dplyr::select(-para_level) %>%
    filter(parameter=="beta1") %>% spread(source, value) %>%
    dplyr::select(Poultry, Ruminants, Water)

#' smaller SD on normal prior
pdf(file="fig_ch3_appxC8_C9_C10_upper.pdf", height=7, width=7)

#' larger SD on normal prior
pdf(file="fig_ch3_appxC8_C9_C10_bottom.pdf", height=7, width=7)

# posterior mean and 80% CI
saplotf
saplotcapF

# traceplot -- alpha and beta
ggplot(postt_all,
       aes(x=iteration, y=value, color=source)) + 
    geom_line() + 
    scale_color_manual(name="source",
                       values=c("blue", "red", "green")) +
    ylab("values") + 
    theme_bw() + facet_grid(para_level~.) +
    theme(legend.position="bottom",
          text=element_text(size=12, family="serif"))

dev.off()

