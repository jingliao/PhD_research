##########********  REQUIRED PACKAGES  ********##########
library(dplyr)
library(MCMCpack)  # dirichlet distribution
library(tidyverse)
library(cowplot)
library(GGally)
library(ggpubr) ## combine ggplots with common legend

source("functions_wCov.R")

##########********  Data processing  ********##########
###******** load the data
base_path <- "data"

#' data path
alltypes <- read.csv(file.path(base_path, "human_types_age.csv"))

#' frequency table of ST vs. sources
freqtable <- as.data.frame.matrix(table(alltypes$ST,
                                          alltypes$Source))

freqtable <- freqtable %>% dplyr::select(Human, Poultry, Ruminants, Water, Other)
SeqTypes <- rownames(freqtable)

#' change order of sources to make the last one as the baseline
## freqtable <- freqtable %>% select(Human, Other, Ruminants, Water, Poultry)

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

############################## choose design matrix
#' without covariates
humandata <- humantypes %>% dplyr::select(ST)
human_no <- nrow(humandata)
H <- model.matrix(ST ~ 1, data=humandata)

#' with rurality
humandata <- humantypes %>% dplyr::select(ST, UR_num) %>% na.omit
human_no <- nrow(humandata)
H <- model.matrix(ST ~ UR_num, data=humandata)

#' with rurality and age
humantypes$age_group <- ifelse(humantypes$Age<5, 1, 0)
humandata <- humantypes %>% dplyr::select(ST, UR_num, age_group) %>% na.omit
human_no <- nrow(humandata)
H <- model.matrix(ST ~ UR_num*age_group, data=humandata)

############################## starting points
# Parameters in our model to update
Z <- matrix(NA, human_no, source_no)
pi <- matrix(NA, ST_no, source_no)

## fill from the priors
sample_theta <- function(para_rowno){ 
  
  matrix(rep(rnorm(source_no-1,0,1), para_rowno), nrow=para_rowno)
  
}

theta <- sample_theta(para_rowno=ncol(H))
                      
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

## FINALLY, sample Z: This is really the first iteration of our MCMC
for (i in 1:human_no) {
    Z[i,] <- rmultinom(1, 1, gamma[i,])
}

## Setup map from each source type row to the human rows that match it
human_index <- matrix(FALSE, ST_no, human_no)

for (i in 1:ST_no) {
    human_index[i,] <- as.numeric(SeqTypes[i]) == humandata$ST
}

## output
burnin <- 1000
thin <- 5
nite <- 10000

totalm <- nite*thin+burnin

accept <- reject <- 0

sigma <- qsigma <- 1

posterior <- list()

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
    if( (m > burnin) & (m%%thin==0) ){
        
        posterior[[m]] <- list(theta=as.numeric(t(theta)),
                               Z=Z,
                               pi=pi)
    }

}


############################## load file
#' without covariate
load(file="OFbayes_wCovH1.RData")

#' with covariate rurality
load(file="OFbayes_wCovH2.RData")

#' with covariates rurality and age
load(file="OFbayes_wCovH3.RData")

############################## drag posterior parameter out
post_theta <- do.call(rbind, purrr::map(posterior, "theta")) 

############################## choose way to calculate f and F
#' without or with covariates

postf_list <- postcapF_list <- list()

for(i in 1:nite){

    ind_theta <- matrix(as.numeric(post_theta[i,]),
                        nrow=dim(theta)[1], byrow=TRUE)

    ### f
    ind_f <- unique(H) %*% ind_theta

    #' without covariate
    postf_list[[i]] <- ind_f
    
    #' with covariate - rurality    
    postf_list[[i]] <- cbind(unique(H)[,2], ind_f)

    #' with covariates - rurality and age    
    postf_list[[i]] <- cbind(unique(H)[,c(2,3)], ind_f)

    ### F
    ind_capF <- delogit_mtrx(ind_f)

    #' without covariate
    postcapF_list[[i]] <- ind_capF
    
    #' with covariate - rurality    
    postcapF_list[[i]] <- cbind(unique(H)[,2], ind_capF)

    #' with covariates - rurality and age    
    postcapF_list[[i]] <- cbind(unique(H)[,c(2,3)], ind_capF)
    
}

post_f <- as.data.frame(do.call(rbind, postf_list))
post_capF <- as.data.frame(as.data.frame(do.call(rbind, postcapF_list))) 


############################ choose w/o or w/ cov, source attribution
##' without covariate, f and F 
names(post_f)  <-  source_name[-4]
names(post_capF) <- source_name

postf_long <- post_f %>% gather(source, postf, Poultry:Water) %>%
    mutate(source=factor(source, levels=c("Poultry", "Ruminants",
                                          "Water")))
postcapF_long <- post_capF %>% gather(source, postcapF,
                                      Poultry:Other) %>%
    mutate(source=factor(source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other")))

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


##' with covariates, rurality and age, f and F
names(post_f)[c(1,2)]  <- names(post_capF)[c(1,2)] <- c("V1", "V2") 

tidyf <- as.data.frame(tidydata(dataname=post_f,
                                sourcename=source_name[-4],
                                CIpercent=80,
                                H=H)) %>%
    rename("Rurality"=V1, "Age"=V2) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "Water")),
           age=factor(Age, labels=c("under 5", "5 or above"),
                      levels=c(1,0)))
           

tidycapF <- as.data.frame(tidydata(dataname=post_capF,
                     sourcename=source_name,
                     CIpercent=80,
                     H=H)) %>% rename("Rurality"=V1, "Age"=V2) %>%
    mutate(Source=factor(Source, levels=c("Poultry", "Ruminants",
                                          "Water", "Other")),
           age=factor(Age, labels=c("age < 5 years", "age >= 5 years"),
                      levels=c(1,0)))


######################### diagnosis with plots, without covariate
post_theta <- as.data.frame(post_theta)
colnames(post_theta) <- source_name[-4]
post_theta$iteration <- seq(1,nite)

postt_long <- as.data.frame(post_theta) %>%
    gather(source, alpha, Poultry:Water) %>%
    mutate(source=factor(source, levels=c("Poultry", "Ruminants",
                                          "Water")))

pdf(file="fig_ch3_1_ch3_3_bottom.pdf", height=5, width=7)
##' without covariate, boxplots for f and F and trace plot for alpha
ggplot(postcapF_long, aes(x=source, y=postcapF*100)) +
    geom_boxplot() +
    ylab("percentage of cases (%)") + 
    theme_bw() +
    theme(text=element_text(size=12,  family="serif"))


ggplot(postt_long, aes(x=iteration, y=alpha, color=source)) + 
  geom_line() + 
    scale_color_manual(name="source",
                       values=c("red", "green","blue")) +
    ylab("intercept") +
    theme_bw() +
    theme(legend.position="bottom",
          text=element_text(size=12,
                            family="serif"))

dev.off()

######################### diagnosis with plots, with rurality
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


pdf(file="fig_ch3_2upper_ch3_4_ch3_5_appxC2.pdf", height=5, width=7)

# posterior mean and 80% CI
ggplot(tidycapF) +
    geom_ribbon(aes(x=Rurality, ymin=lower*100, ymax=upper*100,
                    fill=Source), alpha=0.4) + 
    geom_line(aes(x=Rurality, y=mean*100, color=Source)) +
      scale_x_continuous(name=NULL, breaks=c(-3,3), labels=c("Rural", "Urban")) +
    scale_y_continuous(name="percentage of cases (%)", limits=c(0,100)) +
    scale_color_manual(name="source",
                       values=c("red", "green", "blue","pink")) +
    scale_fill_manual(name="source",
                      values = c("red", "green", "blue", "pink")) +
    theme_bw() + xlab("rurality") +  
    theme(legend.position="none", text=element_text(size=12,
                                                    family="serif"))

# traceplot -- alpha and beta
ggplot(postt_all,
       aes(x=iteration, y=value, color=source)) + 
    geom_line() + 
    scale_color_manual(name="source",
                       values=c("blue", "red", "green")) +
    ylab("values") + 
    theme_bw() + facet_grid(para_level~.) +
    theme(legend.position="bottom",
          text=element_text(size=12,
                            family="serif"))


# scatter plot
ggpairs(newbeta0, title="intercept") +
    theme(text=element_text(size=12, family="serif"))

ggpairs(newbeta1, title="slope of rurality") +
    theme(text=element_text(size=12, family="serif"))



postt_spread <- as.data.frame(postt_all %>%
                              dplyr::select(-para_level) %>%
    spread(parameter, value) %>%
    mutate(source=factor(source,
                         levels=c("Poultry",
                                  "Ruminants",
                                  "Water"))) %>%
    dplyr::select(-iteration) %>% arrange(source))

ggplot(postt_spread, aes(x=beta0, y=beta1, color=source)) +
    geom_point() +
    scale_color_manual(values=c("red", "green", "blue")) +
    theme_bw() +
    xlab("intercept") + ylab("slope of rurality") +
    theme(text=element_text(size=12, family="serif"),
          legend.position="bottom")


dev.off()


######################### diagnosis with plots, with rurality and age

## alpha and beta
post_theta <- as.data.frame(post_theta)
post_theta  <- post_theta %>% mutate(iteration=as.numeric(rownames(post_theta)))

post_beta0 <- post_theta %>% dplyr::select(iteration, V1:V3) %>%
    mutate(parameter="beta0") %>% gather(source, value, V1:V3)
post_beta1 <- post_theta %>% dplyr::select(iteration, V4:V6) %>%
        mutate(parameter="beta1") %>% gather(source, value, V4:V6)
post_beta2 <- post_theta %>% dplyr::select(iteration, V7:V9) %>%
        mutate(parameter="beta2") %>% gather(source, value, V7:V9)
post_beta12 <- post_theta %>% dplyr::select(iteration, V10:V12) %>%
        mutate(parameter="beta12") %>% gather(source, value, V10:V12)

postt_all <- rbind(post_beta0, post_beta1, post_beta2,
                             post_beta12) %>% 
    mutate(source=factor(recode(source, "V1"="Poultry",
                                "V2"="Ruminants",
                                "V3"="Water", "V4"="Poultry",
                                "V5"="Ruminants", "V6"="Water",
                                "V7"="Poultry", "V8"="Ruminants",
                                "V9"="Water", "V10"="Poultry",
                                "V11"="Ruminants", "V12"="Water"),
                         levels=c("Water", "Poultry", "Ruminants")),
           para_level=factor(recode(parameter, "beta0"="intercept",
                                   "beta1"="slope of rurality",
                                   "beta2"="slope of age",
                                   "beta12"="slope of interaction"),
                            levels=c("intercept", "slope of rurality", "slope of age", "slope of interaction")))


# scatter plot
newbeta0 <- postt_all %>% dplyr::select(-para_level) %>%
    filter(parameter=="beta0") %>% spread(source, value) %>%
    dplyr::select(Poultry, Ruminants, Water)

newbeta1 <- postt_all %>% dplyr::select(-para_level) %>%
    filter(parameter=="beta1") %>% spread(source, value) %>%
    dplyr::select(Poultry, Ruminants, Water)

newbeta2 <- postt_all %>% dplyr::select(-para_level) %>%
    filter(parameter=="beta2") %>% spread(source, value) %>%
    dplyr::select(Poultry, Ruminants, Water)

newbeta12 <- postt_all %>% dplyr::select(-para_level) %>%
    filter(parameter=="beta12") %>% spread(source, value) %>%
    dplyr::select(Poultry, Ruminants, Water)


pdf(file="fig_ch3_2bottom_ch3_6_appxC3_C4.pdf", height=5, width=7)

# posterior mean and 80% CI
ggplot(tidycapF) +
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
    theme_bw() + facet_wrap(.~age) +   
    theme(legend.position="bottom", text=element_text(size=12, family="serif"))


# traceplots
ggplot(postt_all,
       aes(x=iteration, y=value, color=source)) + 
    geom_line() + 
    scale_color_manual(name="source",
                       values=c("blue", "red", "green")) +
    ylab("values") + 
    theme_bw() + facet_wrap(.~para_level) +
    theme(legend.position="bottom",
          text=element_text(size=12, family="serif"))

# scatter plot
ggpairs(newbeta0, title="intercept") +
    theme(text=element_text(size=12, family="serif"))
ggpairs(newbeta1, title="slope of rurality") +
    theme(text=element_text(size=12, family="serif"))
ggpairs(newbeta2, title="slope of age") +
    theme(text=element_text(size=12, family="serif"))
ggpairs(newbeta12, title="slope of interaction") +
    theme(text=element_text(size=12, family="serif"))

dev.off()
